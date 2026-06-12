#include "TpcPolyTrackReco.h"

#include "FullTrack.h"
#include "FullTrackContainer.h"
#include "IdealPadMap.h"
#include "TpcPolyTrackContainerv1.h"
#include "TpcPolyTrackv1.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>

#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <phgarfield/PHGarfield.h>

#include <TPolyLine.h>
#include <TPolyLine3D.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

namespace
{
  double wrap_pi(double phi)
  {
    while (phi > M_PI) phi -= 2.0 * M_PI;
    while (phi <= -M_PI) phi += 2.0 * M_PI;
    return phi;
  }

  double unwrap_near(double phi, const double ref)
  {
    while (phi - ref > M_PI) phi -= 2.0 * M_PI;
    while (phi - ref < -M_PI) phi += 2.0 * M_PI;
    return phi;
  }

  bool solve_3x3(double A[3][3], double b[3], double x[3])
  {
    double M[3][4] = {
      {A[0][0], A[0][1], A[0][2], b[0]},
      {A[1][0], A[1][1], A[1][2], b[1]},
      {A[2][0], A[2][1], A[2][2], b[2]}
    };
    for (int col = 0; col < 3; ++col)
    {
      int pivot = col;
      for (int row = col + 1; row < 3; ++row)
        if (std::fabs(M[row][col]) > std::fabs(M[pivot][col])) pivot = row;
      if (std::fabs(M[pivot][col]) < 1.0e-20) return false;
      if (pivot != col)
        for (int k = col; k < 4; ++k) std::swap(M[col][k], M[pivot][k]);
      const double div = M[col][col];
      for (int k = col; k < 4; ++k) M[col][k] /= div;
      for (int row = 0; row < 3; ++row)
      {
        if (row == col) continue;
        const double factor = M[row][col];
        for (int k = col; k < 4; ++k) M[row][k] -= factor * M[col][k];
      }
    }
    x[0] = M[0][3];
    x[1] = M[1][3];
    x[2] = M[2][3];
    return true;
  }

  bool line_fit(const std::vector<double>& x, const std::vector<double>& y,
                double& slope, double& intercept, double& chi2, int& ndof)
  {
    if (x.size() < 2 || x.size() != y.size()) return false;
    double S = 0.0, Sx = 0.0, Sy = 0.0, Sxx = 0.0, Sxy = 0.0;
    for (unsigned int i = 0; i < x.size(); ++i)
    {
      S += 1.0;
      Sx += x[i];
      Sy += y[i];
      Sxx += x[i] * x[i];
      Sxy += x[i] * y[i];
    }
    const double den = S * Sxx - Sx * Sx;
    if (std::fabs(den) < 1.0e-20) return false;
    slope = (S * Sxy - Sx * Sy) / den;
    intercept = (Sy - slope * Sx) / S;
    chi2 = 0.0;
    for (unsigned int i = 0; i < x.size(); ++i)
    {
      const double r = y[i] - (slope * x[i] + intercept);
      chi2 += r * r;
    }
    ndof = static_cast<int>(x.size()) - 2;
    return true;
  }
}

TpcPolyTrackReco::TpcPolyTrackReco(const std::string& name)
  : SubsysReco(name)
  , m_inputNodeName("FULLTRACKS")
  , m_outputNodeName("TPCPOLYTRACKS")
{
}

TpcPolyTrackReco::~TpcPolyTrackReco()
{
  delete m_idealPadMap;
  m_idealPadMap = nullptr;
  delete m_garfield;
  m_garfield = nullptr;
}

int TpcPolyTrackReco::InitRun(PHCompositeNode* topNode)
{
  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK) return Fun4AllReturnCodes::ABORTRUN;
  if (createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK) return Fun4AllReturnCodes::ABORTRUN;

  delete m_idealPadMap;
  m_idealPadMap = new IdealPadMap();
  if (m_idealPadMap->load_from_cdb(Verbosity()) != 0 || !m_idealPadMap->is_loaded())
  {
    std::cerr << Name() << "::InitRun - failed to load IdealPadMap" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  delete m_garfield;
  m_garfield = new PHGarfield(Name() + "_PHGarfield");
  if (m_garfield->InitRun(topNode) != Fun4AllReturnCodes::EVENT_OK)
  {
    std::cerr << Name() << "::InitRun - PHGarfield InitRun failed" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_event = 0;
  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcPolyTrackReco::getNodes(PHCompositeNode* topNode)
{
  m_fullTracks = findNode::getClass<FullTrackContainer>(topNode, m_inputNodeName);
  if (!m_fullTracks)
  {
    std::cerr << Name() << "::getNodes - missing " << m_inputNodeName << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_hits = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!m_hits)
  {
    std::cerr << Name() << "::getNodes - missing TRKR_HITSET" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcPolyTrackReco::createNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    dstNode = new PHCompositeNode("DST");
    topNode->addNode(dstNode);
  }

  m_polyTracks = findNode::getClass<TpcPolyTrackContainer>(topNode, m_outputNodeName);
  if (!m_polyTracks)
  {
    m_polyTracks = new TpcPolyTrackContainerv1();
    PHIODataNode<PHObject>* node = new PHIODataNode<PHObject>(m_polyTracks, m_outputNodeName, "PHObject");
    dstNode->addNode(node);
    std::cout << Name() << "::createNodes - created " << m_outputNodeName << " node" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

bool TpcPolyTrackReco::make_xyz_point(TrkrDefs::hitsetkey hsk,
                                      TrkrDefs::hitkey hk,
                                      int side,
                                      Point& p) const
{
  if (!m_hits || !m_idealPadMap || !m_garfield) return false;
  TrkrHitSet* hitset = m_hits->findHitSet(hsk);
  if (!hitset) return false;
  TrkrHit* hit = hitset->getHit(hk);
  if (!hit) return false;

  const unsigned int layer = TrkrDefs::getLayer(hsk);
  const unsigned int pad = TpcDefs::getPad(hk);
  const unsigned int tbin = TpcDefs::getTBin(hk);
  if (layer < 7 || layer > 54) return false;

  const double radius = m_idealPadMap->get_radius(layer);
  const double phi = m_idealPadMap->get_phi(static_cast<unsigned int>(side), layer, pad);
  if (!std::isfinite(radius) || !std::isfinite(phi)) return false;

  const double corrected_tbin = static_cast<double>(tbin) - m_t0;
  const double target_time_ns = corrected_tbin * m_tpcAdcClock;
  if (target_time_ns <= 0.0 || !std::isfinite(target_time_ns)) return false;
  if (m_reverseDriftStepNs <= 0.0 || !std::isfinite(m_reverseDriftStepNs)) return false;

  const double x0 = radius * std::cos(phi);
  const double y0 = radius * std::sin(phi);
  const double z0 = (side == 0) ? m_startZSouth : m_startZNorth;

  TPolyLine3D* drift = m_garfield->ReverseDrift(x0, y0, z0, m_reverseDriftStepNs);
  if (!drift || drift->GetN() <= 0)
  {
    delete drift;
    return false;
  }

  const int npoints = drift->GetN();
  const Float_t* xyz = drift->GetP();
  if (!xyz || npoints <= 0)
  {
    delete drift;
    return false;
  }

  const double max_time_ns = static_cast<double>(npoints - 1) * m_reverseDriftStepNs;
  if (target_time_ns > max_time_ns)
  {
    delete drift;
    return false;
  }

  const double fbin = target_time_ns / m_reverseDriftStepNs;
  const int i0 = std::min(static_cast<int>(std::floor(fbin)), npoints - 1);
  const int i1 = std::min(i0 + 1, npoints - 1);
  const double frac = fbin - static_cast<double>(i0);

  const int idx0 = 3 * i0;
  const int idx1 = 3 * i1;
  const double x = xyz[idx0] + frac * (xyz[idx1] - xyz[idx0]);
  const double y = xyz[idx0 + 1] + frac * (xyz[idx1 + 1] - xyz[idx0 + 1]);
  const double z = xyz[idx0 + 2] + frac * (xyz[idx1 + 2] - xyz[idx0 + 2]);
  delete drift;

  if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) return false;
  p.hitsetkey = hsk;
  p.hitkey = hk;
  p.x = x;
  p.y = y;
  p.z = z;
  return true;
}

bool TpcPolyTrackReco::fit_points(const std::vector<Point>& points, FitResult& fit) const
{
  fit = FitResult();
  if (points.size() < 3) return false;

  double A[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  double b[3] = {0.0, 0.0, 0.0};
  for (const Point& p : points)
  {
    const double row[3] = {p.x, p.y, 1.0};
    const double rhs = -(p.x * p.x + p.y * p.y);
    for (int i = 0; i < 3; ++i)
    {
      b[i] += row[i] * rhs;
      for (int j = 0; j < 3; ++j) A[i][j] += row[i] * row[j];
    }
  }

  double sol[3] = {0.0, 0.0, 0.0};
  if (!solve_3x3(A, b, sol)) return false;

  const double xc = -0.5 * sol[0];
  const double yc = -0.5 * sol[1];
  const double r2 = xc * xc + yc * yc - sol[2];
  if (r2 <= 0.0) return false;
  const double R = std::sqrt(r2);
  const double dc = std::hypot(xc, yc);
  if (R <= 0.0 || dc <= 1.0e-12) return false;

  const Point& first = points.front();
  const Point& mid = points[points.size() / 2];
  const Point& last = points.back();
  const double cross = (mid.x - first.x) * (last.y - mid.y) - (mid.y - first.y) * (last.x - mid.x);
  const double sign = (cross >= 0.0) ? 1.0 : -1.0;

  fit.curvature = sign / R;
  fit.d0 = sign * (dc - R);

  const double px = xc * (1.0 - R / dc);
  const double py = yc * (1.0 - R / dc);
  const double rx = px - xc;
  const double ry = py - yc;
  double tx = -sign * ry / R;
  double ty =  sign * rx / R;
  const double dxfl = last.x - first.x;
  const double dyfl = last.y - first.y;
  if (tx * dxfl + ty * dyfl < 0.0)
  {
    tx = -tx;
    ty = -ty;
  }
  fit.phi0 = wrap_pi(std::atan2(ty, tx));

  const double phi_perigee = std::atan2(ry, rx);
  std::vector<double> svals;
  std::vector<double> zvals;
  svals.reserve(points.size());
  zvals.reserve(points.size());
  double prev_angle = phi_perigee;
  for (const Point& p : points)
  {
    double angle = std::atan2(p.y - yc, p.x - xc);
    angle = unwrap_near(angle, prev_angle);
    prev_angle = angle;
    double dangle = angle - phi_perigee;
    if (sign * dangle < 0.0) dangle += sign * 2.0 * M_PI;
    svals.push_back(R * dangle);
    zvals.push_back(p.z);

    const double resid = std::hypot(p.x - xc, p.y - yc) - R;
    fit.chi2_xy += resid * resid;
  }

  double dzds = 0.0;
  if (!line_fit(svals, zvals, dzds, fit.z0, fit.chi2_z, fit.ndof_z)) return false;
  fit.theta = std::atan2(1.0, dzds);
  fit.ndof_xy = static_cast<int>(points.size()) - 3;
  fit.ok = true;
  return true;
}

int TpcPolyTrackReco::process_event(PHCompositeNode*)
{
  if (!m_fullTracks || !m_polyTracks) return Fun4AllReturnCodes::EVENT_OK;
  m_polyTracks->Reset();

  const unsigned int nfull = m_fullTracks->size();
  for (unsigned int ifull = 0; ifull < nfull; ++ifull)
  {
    const FullTrack* full = m_fullTracks->get_track(ifull);
    if (!full) continue;

    std::vector<Point> points;
    points.reserve(full->size_hit_indices());
    for (unsigned int ih = 0; ih < full->size_hit_indices(); ++ih)
    {
      const FullTrack::HitIndex hi = full->get_hit_index(ih);
      Point p;
      if (make_xyz_point(hi.first, hi.second, full->get_side(), p)) points.push_back(p);
    }

    FitResult fit;
    const bool fit_ok = fit_points(points, fit);

    TpcPolyTrackv1* out = new TpcPolyTrackv1();
    out->set_event(m_event);
    out->set_track_id(m_polyTracks->size());
    out->set_source_full_track_id(full->get_track_id());
    out->set_side(full->get_side());
    out->set_fit_status(fit_ok ? 1 : 0);
    if (fit_ok)
    {
      out->set_d0(fit.d0);
      out->set_z0(fit.z0);
      out->set_phi0(fit.phi0);
      out->set_theta(fit.theta);
      out->set_curvature(fit.curvature);
      out->set_chi2_xy(fit.chi2_xy);
      out->set_chi2_z(fit.chi2_z);
      out->set_ndof_xy(fit.ndof_xy);
      out->set_ndof_z(fit.ndof_z);
    }
    for (const Point& p : points) out->add_hit(p.hitsetkey, p.hitkey, p.x, p.y, p.z);
    m_polyTracks->add_track(out);
  }

  if (Verbosity() > 0)
  {
    std::cout << Name() << "::process_event - event " << m_event
              << " full_tracks=" << nfull
              << " poly_tracks=" << m_polyTracks->size() << std::endl;
  }
  ++m_event;
  return Fun4AllReturnCodes::EVENT_OK;
}
