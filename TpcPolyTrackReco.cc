#include "TpcPolyTrackReco.h"

#include "FullTrack.h"
#include "FullTrackContainer.h"
#include "IdealPadMap.h"
#include "TpcPolyTrackContainerv1.h"
#include "TpcPolyTrackv1.h"
#include "TpcPolyHelixFitter.h"

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
  delete m_helixFitter;
  m_helixFitter = nullptr;
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

  delete m_helixFitter;
  m_helixFitter = new TpcPolyHelixFitter();
  if (!m_helixFitter->InitField(Verbosity()))
  {
    std::cerr << Name() << "::InitRun - failed to load FIELDMAP_TRACKING" << std::endl;
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
//  std::cout<<"    TpcPolyTrackReco::make_xyz_point"<<std::endl;
  if (!m_hits || !m_idealPadMap || !m_garfield)
  {
    if (!m_hits) std::cout<<"No m_hits "<<std::endl;
    if (!m_idealPadMap) std::cout<<"No m_idealPadMap"<<std::endl;
    if (!m_garfield)std::cout<<"No m_garfield"<<std::endl;
    return false;
  } 
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
//std::cout<<"Before Reverse Drift"<<std::endl;
  TPolyLine3D* drift = m_garfield->ReverseDrift(x0, y0, z0, m_reverseDriftStepNs);
  if (!drift || drift->GetN() <= 0)
  {
    delete drift;
    return false;
  }
//std::cout<<"After Reverse Drift"<<std::endl;
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


int TpcPolyTrackReco::process_event(PHCompositeNode*)
{
//  std::cout<<"TpcPolyTrackReco::process_event!!!!!!"<<std::endl;
  if (!m_fullTracks || !m_polyTracks || !m_helixFitter) return Fun4AllReturnCodes::EVENT_OK;
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

    std::vector<TpcPolyHelixFitter::Point> fit_points;
    fit_points.reserve(points.size());
    for (const Point& p : points)
    {
      TpcPolyHelixFitter::Point fp;
      fp.x = p.x;
      fp.y = p.y;
      fp.z = p.z;
      fit_points.push_back(fp);
    }

    TpcPolyHelixFitter::FitResult fit;
    const bool fit_ok = m_helixFitter->fit(fit_points, fit);

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
