#include "FullTrackVertexer.h"

#include "Fitter.h"
#include "FullTrack.h"
#include "FullTrackContainer.h"
#include "FullTrackVertexv1.h"
#include "FullTrackVertexContainer.h"
#include "FullTrackVertexContainerv1.h"
#include "IdealPadMap.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>

#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <TMath.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

namespace
{
  double wrap_phi(double phi)
  {
    while (phi > TMath::Pi()) phi -= 2.0 * TMath::Pi();
    while (phi <= -TMath::Pi()) phi += 2.0 * TMath::Pi();
    return phi;
  }

  double unwrap_phi_near(double phi, const double reference)
  {
    while (phi - reference > TMath::Pi()) phi -= 2.0 * TMath::Pi();
    while (phi - reference < -TMath::Pi()) phi += 2.0 * TMath::Pi();
    return phi;
  }

  bool is_good_number(const double x)
  {
    return std::isfinite(x) && std::fabs(x) < 1.0e30;
  }

  double sagitta_model_derivative(const double xrot,
                                  const double x0,
                                  const double invR)
  {
    const double dx = xrot - x0;
    const double dx2 = dx * dx;
    const double invR2 = invR * invR;
    const double invR3 = invR2 * invR;
    const double invR5 = invR3 * invR2;
    return -invR * dx - 0.5 * invR3 * dx2 * dx - 0.375 * invR5 * dx2 * dx2 * dx;
  }

  double sagitta_phi_at_radius(const double radius,
                               const Fitter::SagittaFit& fit)
  {
    const double c = std::cos(fit.theta);
    const double s = std::sin(fit.theta);
    double yy = std::tan(fit.theta) * radius;

    for (unsigned int iter = 0; iter < 25; ++iter)
    {
      const double xrot = c * radius + s * yy;
      const double yrot = -s * radius + c * yy;
      const double f = Fitter::sagittaModel(xrot, fit.S, fit.x0, fit.invR);
      const double g = yrot - f;
      const double df = sagitta_model_derivative(xrot, fit.x0, fit.invR);
      const double dg = c - df * s;
      if (std::fabs(dg) < 1.0e-12) break;
      const double step = g / dg;
      yy -= step;
      if (std::fabs(step) < 1.0e-10) break;
    }

    return fit.b + yy;
  }

  struct RadiusSort
  {
    const std::vector<FullTrackVertexer::HitPoint>* pts;
    explicit RadiusSort(const std::vector<FullTrackVertexer::HitPoint>* p) : pts(p) {}
    bool operator()(unsigned int a, unsigned int b) const
    {
      const FullTrackVertexer::HitPoint& pa = (*pts)[a];
      const FullTrackVertexer::HitPoint& pb = (*pts)[b];
      if (pa.radius != pb.radius) return pa.radius < pb.radius;
      return pa.tbin < pb.tbin;
    }
  };

  struct VertexFit
  {
    bool ok {false};
    bool use_sagitta {false};
    double d0 {0.0};
    double timebin0 {0.0};
    double phi_slope {0.0};
    double phi_intercept {0.0};
    double tbin_slope {0.0};
    double tbin_intercept {0.0};
    Fitter::SagittaFit phi_sagitta;
    double rmin {0.0};
    double rmax {0.0};
  };

  struct PcaPoint
  {
    bool ok {false};
    double radius {0.0};
    double phi {0.0};
    double timebin {0.0};
  };

  struct CollisionTrack
  {
    VertexFit fit;
    unsigned int nlayers {0};
    int side {-1};
  };

  struct CollisionFit
  {
    bool ok {false};
    double radius {0.0};
    double phi {0.0};
    double timebin {0.0};
    double timebin_rms {0.0};
    unsigned int ntracks {0};
  };

  double median(std::vector<double> values)
  {
    if (values.empty()) return 0.0;
    std::sort(values.begin(), values.end());
    const std::size_t n = values.size();
    const std::size_t mid = n / 2;
    if (n % 2 == 1) return values[mid];
    return 0.5 * (values[mid - 1] + values[mid]);
  }

  unsigned int count_unique_layers(const std::vector<FullTrackVertexer::HitPoint>& pts)
  {
    std::vector<unsigned int> layers;
    layers.reserve(pts.size());
    for (const auto& p : pts)
    {
      if (p.ok) layers.push_back(p.layer);
    }
    std::sort(layers.begin(), layers.end());
    layers.erase(std::unique(layers.begin(), layers.end()), layers.end());
    return static_cast<unsigned int>(layers.size());
  }

  VertexFit fit_vertex_points(const std::vector<FullTrackVertexer::HitPoint>& pts,
                              const bool use_sagitta,
                              const double weight_power,
                              const double weight_floor_frac)
  {
    VertexFit result;
    if (pts.size() < 2) return result;

    std::vector<unsigned int> order;
    order.reserve(pts.size());
    for (unsigned int i = 0; i < pts.size(); ++i) order.push_back(i);
    std::sort(order.begin(), order.end(), RadiusSort(&pts));

    double max_adc = 0.0;
    for (const auto& p : pts) max_adc = std::max(max_adc, static_cast<double>(p.adc));
    if (max_adc <= 0.0) max_adc = 1.0;

    std::vector<Fitter::FitPoint> radius_phi_points;
    std::vector<Fitter::FitPoint> radius_tbin_points;
    radius_phi_points.reserve(pts.size());
    radius_tbin_points.reserve(pts.size());

    bool first = true;
    double phi_reference = 0.0;
    result.rmin = std::numeric_limits<double>::max();
    result.rmax = -std::numeric_limits<double>::max();
    for (unsigned int io = 0; io < order.size(); ++io)
    {
      const FullTrackVertexer::HitPoint& p = pts[order[io]];
      if (!is_good_number(p.radius) || !is_good_number(p.global_phi)) continue;

      double phi = p.global_phi;
      if (first)
      {
        phi_reference = phi;
        first = false;
      }
      else
      {
        phi = unwrap_phi_near(phi, phi_reference);
        phi_reference = phi;
      }

      const double w = Fitter::adcWeight(static_cast<double>(p.adc), max_adc,
                                         weight_power, weight_floor_frac);
      radius_phi_points.emplace_back(p.radius, phi, w);
      radius_tbin_points.emplace_back(p.radius, static_cast<double>(p.tbin), w);
      result.rmin = std::min(result.rmin, p.radius);
      result.rmax = std::max(result.rmax, p.radius);
    }

    if (radius_phi_points.size() < 2 || radius_tbin_points.size() < 2) return result;

    const Fitter::LineFit phi_line = Fitter::fitLine(radius_phi_points);
    const Fitter::LineFit tbin_line = Fitter::fitLine(radius_tbin_points);
    if (!phi_line.ok || !tbin_line.ok) return result;

    result.d0 = phi_line.intercept;
    result.timebin0 = tbin_line.intercept;
    result.phi_slope = phi_line.slope;
    result.phi_intercept = phi_line.intercept;
    result.tbin_slope = tbin_line.slope;
    result.tbin_intercept = tbin_line.intercept;

    if (use_sagitta && radius_phi_points.size() >= 3)
    {
      result.phi_sagitta = Fitter::fitSagitta(radius_phi_points);
      result.use_sagitta = result.phi_sagitta.ok;
      if (result.use_sagitta) result.d0 = sagitta_phi_at_radius(0.0, result.phi_sagitta);
    }

    if (!is_good_number(result.d0) || !is_good_number(result.timebin0)) return result;
    result.d0 = wrap_phi(result.d0);
    result.ok = true;
    return result;
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
      {
        if (std::fabs(M[row][col]) > std::fabs(M[pivot][col])) pivot = row;
      }
      if (std::fabs(M[pivot][col]) < 1.0e-20) return false;
      if (pivot != col)
      {
        for (int k = col; k < 4; ++k) std::swap(M[col][k], M[pivot][k]);
      }

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

  PcaPoint closest_approach_to_beamline(const std::vector<FullTrackVertexer::HitPoint>& pts,
                                        const VertexFit& fit)
  {
    PcaPoint pca;
    if (!fit.ok || pts.size() < 3) return pca;

    double A[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    double b[3] = {0.0, 0.0, 0.0};
    unsigned int nxy = 0;

    for (const auto& point : pts)
    {
      if (!point.ok || !is_good_number(point.radius) || !is_good_number(point.global_phi)) continue;
      const double x = point.radius * std::cos(point.global_phi);
      const double y = point.radius * std::sin(point.global_phi);
      const double row[3] = {x, y, 1.0};
      const double rhs = -(x * x + y * y);
      for (int i = 0; i < 3; ++i)
      {
        b[i] += row[i] * rhs;
        for (int j = 0; j < 3; ++j) A[i][j] += row[i] * row[j];
      }
      ++nxy;
    }

    if (nxy < 3) return pca;

    double sol[3] = {0.0, 0.0, 0.0};
    if (!solve_3x3(A, b, sol)) return pca;

    const double xc = -0.5 * sol[0];
    const double yc = -0.5 * sol[1];
    const double r2 = xc * xc + yc * yc - sol[2];
    if (r2 <= 0.0) return pca;

    const double R = std::sqrt(r2);
    const double dc = std::hypot(xc, yc);
    if (R <= 0.0 || dc <= 1.0e-12) return pca;

    const double px = xc * (1.0 - R / dc);
    const double py = yc * (1.0 - R / dc);
    pca.radius = std::hypot(px, py);
    pca.phi = wrap_phi(std::atan2(py, px));
    pca.timebin = fit.tbin_slope * pca.radius + fit.tbin_intercept;
    pca.ok = is_good_number(pca.radius) && is_good_number(pca.phi) && is_good_number(pca.timebin);
    return pca;
  }

  CollisionFit fit_collision_point(const std::vector<CollisionTrack>& tracks)
  {
    CollisionFit result;
    if (tracks.empty()) return result;

    std::vector<double> timebins;
    std::vector<double> phis;
    unsigned int nvalid = 0;

    for (const auto& trk : tracks)
    {
      if (!trk.fit.ok || !is_good_number(trk.fit.d0) || !is_good_number(trk.fit.timebin0)) continue;

      timebins.push_back(trk.fit.timebin0);
      phis.push_back(wrap_phi(trk.fit.d0));
      ++nvalid;
    }

    if (nvalid == 0) return result;

    double sum_sin_d0 = 0.0;
    double sum_cos_d0 = 0.0;
    for (const double phi : phis)
    {
      sum_sin_d0 += std::sin(phi);
      sum_cos_d0 += std::cos(phi);
    }

    const double phi_reference = (sum_sin_d0 == 0.0 && sum_cos_d0 == 0.0)
      ? phis.front()
      : std::atan2(sum_sin_d0, sum_cos_d0);
    for (double& phi : phis)
    {
      phi = unwrap_phi_near(phi, phi_reference);
    }

    result.radius = 0.0;
    result.phi = wrap_phi(median(phis));
    result.timebin = median(timebins);
    if (!is_good_number(result.phi) || !is_good_number(result.timebin)) return result;

    double sum_res2 = 0.0;
    unsigned int nrms = 0;
    for (const auto& trk : tracks)
    {
      if (!trk.fit.ok || !is_good_number(trk.fit.d0) || !is_good_number(trk.fit.timebin0)) continue;

      const double residual = trk.fit.timebin0 - result.timebin;
      sum_res2 += residual * residual;
      ++nrms;
    }

    result.timebin_rms = nrms > 0 ? std::sqrt(sum_res2 / static_cast<double>(nrms)) : 0.0;
    result.ntracks = nvalid;
    result.ok = true;
    return result;
  }

  std::vector<CollisionFit> fit_collision_points_per_side(std::vector<CollisionTrack> tracks)
  {
    std::vector<CollisionFit> collisions;
    if (tracks.empty()) return collisions;

    tracks.erase(std::remove_if(tracks.begin(), tracks.end(),
                                [](const CollisionTrack& trk)
                                {
                                  return trk.side < 0 || trk.side > 1 ||
                                         !trk.fit.ok ||
                                         !is_good_number(trk.fit.d0) ||
                                         !is_good_number(trk.fit.timebin0);
                                }),
                 tracks.end());
    if (tracks.empty()) return collisions;

    for (int side = 0; side < 2; ++side)
    {
      std::vector<CollisionTrack> side_tracks;
      side_tracks.reserve(tracks.size());
      for (const auto& trk : tracks)
      {
        if (trk.side == side) side_tracks.push_back(trk);
      }
      const CollisionFit collision = fit_collision_point(side_tracks);
      if (collision.ok) collisions.push_back(collision);
    }

    return collisions;
  }
}

FullTrackVertexer::FullTrackVertexer(const std::string& name)
  : SubsysReco(name)
  , m_inputNodeName("FULLTRACKS")
  , m_outputNodeName("FULLTRACKVERTICES")
  , m_fullTracks(nullptr)
  , m_vertices(nullptr)
  , m_hits(nullptr)
  , m_idealPadMap(nullptr)
  , m_fitWeightPower(1.0)
  , m_fitWeightFloorFrac(0.05)
  , m_useSagittaPhiFit(true)
  , m_collisionMinTrackLayers(10)
  , m_collisionTimebinSeparation(100.0)
{
}

FullTrackVertexer::~FullTrackVertexer()
{
  delete m_idealPadMap;
  m_idealPadMap = nullptr;
}

int FullTrackVertexer::InitRun(PHCompositeNode* topNode)
{
  if (!createNodes(topNode)) return Fun4AllReturnCodes::ABORTRUN;

  delete m_idealPadMap;
  m_idealPadMap = new IdealPadMap();
  if (m_idealPadMap->load_from_cdb(Verbosity()) != 0 || !m_idealPadMap->is_loaded())
  {
    std::cerr << Name() << "::InitRun - cannot load IdealPadMap from CDB" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

bool FullTrackVertexer::getNodes(PHCompositeNode* topNode)
{
  m_fullTracks = findNode::getClass<FullTrackContainer>(topNode, m_inputNodeName.c_str());
  if (!m_fullTracks)
  {
    const char* candidate_names[] = {
      "FULLTRACKS",
      "FullTracks",
      "FullTrackContainer",
      "FULLTRACKCONTAINER",
      "FULLTRACKS_CONTAINER"
    };

    for (unsigned int i = 0;
         i < sizeof(candidate_names) / sizeof(candidate_names[0]) && !m_fullTracks;
         ++i)
    {
      m_fullTracks = findNode::getClass<FullTrackContainer>(topNode, candidate_names[i]);
      if (m_fullTracks) m_inputNodeName = candidate_names[i];
    }
  }

  if (!m_fullTracks)
  {
    std::cerr << Name() << "::getNodes - missing FullTrackContainer node" << std::endl;
    return false;
  }

  m_hits = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!m_hits)
  {
    std::cerr << Name() << "::getNodes - missing TRKR_HITSET" << std::endl;
    return false;
  }

  return true;
}

bool FullTrackVertexer::createNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    dstNode = new PHCompositeNode("DST");
    topNode->addNode(dstNode);
  }

  m_vertices = findNode::getClass<FullTrackVertexContainer>(topNode, m_outputNodeName.c_str());
  if (!m_vertices)
  {
    m_vertices = new FullTrackVertexContainerv1();
    PHIODataNode<PHObject>* node = new PHIODataNode<PHObject>(m_vertices, m_outputNodeName, "PHObject");
    dstNode->addNode(node);
    std::cout << Name() << "::createNodes - created " << m_outputNodeName << " node" << std::endl;
  }

  return true;
}

int FullTrackVertexer::process_event(PHCompositeNode* topNode)
{
  if (!m_fullTracks || !m_hits || !m_vertices)
  {
    if (!getNodes(topNode) || !createNodes(topNode)) return Fun4AllReturnCodes::EVENT_OK;
  }

  if (!m_idealPadMap || !m_idealPadMap->is_loaded())
  {
    std::cerr << Name() << "::process_event - IdealPadMap is not loaded" << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  m_vertices->Reset();

  std::vector<CollisionTrack> collision_tracks;

  const unsigned int ntracks = m_fullTracks ? m_fullTracks->size() : 0;
  for (unsigned int itrk = 0; itrk < ntracks; ++itrk)
  {
    const FullTrack* trk = m_fullTracks->get_track(itrk);
    if (!trk) continue;

    std::vector<HitPoint> pts;
    pts.reserve(trk->size_hit_indices());
    for (unsigned int ih = 0; ih < trk->size_hit_indices(); ++ih)
    {
      const FullTrack::HitIndex idx = trk->get_hit_index(ih);
      const HitPoint p = make_hit_point(idx.first, idx.second);
      if (p.ok) pts.push_back(p);
    }

    const VertexFit fit = fit_vertex_points(pts, m_useSagittaPhiFit,
                                            m_fitWeightPower, m_fitWeightFloorFrac);
    if (!fit.ok) continue;

    FullTrackVertexv1* out = new FullTrackVertexv1();
    out->set_track_id(trk->get_track_id());
    out->set_d0(fit.d0);
    out->set_timebin0(fit.timebin0);
    m_vertices->add_vertex(out);
    const PcaPoint pca = closest_approach_to_beamline(pts, fit);
    if (pca.ok)
    {
      out->set_pca_valid(1);
      out->set_pca_radius(pca.radius);
      out->set_pca_phi(pca.phi);
      out->set_pca_timebin(pca.timebin);
    }

    const unsigned int nlayers = count_unique_layers(pts);
    if (nlayers >= m_collisionMinTrackLayers)
    {
      CollisionTrack accepted;
      accepted.fit = fit;
      accepted.nlayers = nlayers;
      accepted.side = trk->get_side();
      collision_tracks.push_back(accepted);
    }
  }

  const std::vector<CollisionFit> collisions = fit_collision_points_per_side(collision_tracks);
  m_vertices->set_collision_min_layers(m_collisionMinTrackLayers);
  m_vertices->clear_collision_vertices();
  for (const auto& collision : collisions)
  {
    m_vertices->add_collision_vertex(collision.radius, collision.phi, collision.timebin,
                                     collision.timebin_rms, collision.ntracks);
  }
  m_vertices->set_collision_vertex_valid(collisions.empty() ? 0 : 1);

  if (Verbosity() > 0)
  {
    std::cout << Name() << "::process_event - input tracks=" << ntracks
              << " vertices=" << m_vertices->size()
              << " collision_vertices=" << m_vertices->get_collision_vertex_count()
              << " first_collision_ntracks=" << m_vertices->get_collision_ntracks()
              << " first_collision_radius=" << m_vertices->get_collision_radius()
              << " first_collision_phi=" << m_vertices->get_collision_phi()
              << " first_collision_timebin=" << m_vertices->get_collision_timebin()
              << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

FullTrackVertexer::HitPoint
FullTrackVertexer::make_hit_point(const TrkrDefs::hitsetkey hsk,
                                  const TrkrDefs::hitkey hk) const
{
  HitPoint p;

  TrkrHitSet* hitset = m_hits ? m_hits->findHitSet(hsk) : nullptr;
  if (!hitset) return p;

  TrkrHit* hit = hitset->getHit(hk);
  if (!hit) return p;

  p.hitsetkey = hsk;
  p.hitkey = hk;
  p.layer = TrkrDefs::getLayer(hsk);
  p.pad = TpcDefs::getPad(hk);
  p.tbin = TpcDefs::getTBin(hk);
  p.adc = hit->getAdc();

  if (p.layer < 7 || p.layer > 54) return p;
  if (!m_idealPadMap) return p;

  p.radius = m_idealPadMap->get_radius(p.layer);
  p.global_phi = wrap_phi(m_idealPadMap->get_phi(static_cast<unsigned int>(TpcDefs::getSide(hsk)), p.layer, p.pad));

  if (!is_good_number(p.radius) || !is_good_number(p.global_phi)) return p;

  p.ok = true;
  return p;
}
