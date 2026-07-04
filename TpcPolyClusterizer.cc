#include "TpcPolyClusterizer.h"

#include "FullTrack.h"
#include "FullTrackContainer.h"
#include "IdealPadMap.h"
#include "TpcPolyClusterTrackContainerv1.h"
#include "TpcPolyClusterTrackv1.h"

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

#include <TPolyLine3D.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <vector>

namespace
{
  double wrap_phi(double phi)
  {
    while (phi > M_PI) phi -= 2.0 * M_PI;
    while (phi <= -M_PI) phi += 2.0 * M_PI;
    return phi;
  }
}

TpcPolyClusterizer::TpcPolyClusterizer(const std::string& name)
  : SubsysReco(name)
  , m_inputNodeName("FULLTRACKS")
  , m_outputNodeName("TPCPOLYCLUSTERTRACKS")
{
}

TpcPolyClusterizer::~TpcPolyClusterizer()
{
  delete m_idealPadMap;
  m_idealPadMap = nullptr;
  delete m_garfield;
  m_garfield = nullptr;
}

int TpcPolyClusterizer::InitRun(PHCompositeNode* topNode)
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

int TpcPolyClusterizer::getNodes(PHCompositeNode* topNode)
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

int TpcPolyClusterizer::createNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    dstNode = new PHCompositeNode("DST");
    topNode->addNode(dstNode);
  }

  m_clusterTracks = findNode::getClass<TpcPolyClusterTrackContainer>(topNode, m_outputNodeName);
  if (!m_clusterTracks)
  {
    m_clusterTracks = new TpcPolyClusterTrackContainerv1();
    PHIODataNode<PHObject>* node = new PHIODataNode<PHObject>(m_clusterTracks, m_outputNodeName, "PHObject");
    dstNode->addNode(node);
    std::cout << Name() << "::createNodes - created " << m_outputNodeName << " node" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

bool TpcPolyClusterizer::make_xyz_point(TrkrDefs::hitsetkey hsk,
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
  const double adc = hit->getAdc();
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
  p.layer = layer;
  p.pad = pad;
  p.tbin = tbin;
  p.adc = adc;
  p.x = x;
  p.y = y;
  p.z = z;
  return true;
}

TpcPolyClusterizer::ClusterParameters
TpcPolyClusterizer::make_cluster_parameters(const std::vector<Point>& points,
                                            const Centroid& centroid,
                                            const int side) const
{
  ClusterParameters params;
  if (points.empty() || !centroid.ok || !m_idealPadMap) return params;

  std::set<unsigned int> pads;
  std::set<unsigned int> tbins;
  std::map<unsigned int, double> adc_by_pad;
  for (const Point& p : points)
  {
    params.adc += p.adc;
    pads.insert(p.pad);
    tbins.insert(p.tbin);
    adc_by_pad[p.pad] += p.adc;
  }

  params.phi_width = static_cast<unsigned int>(pads.size());
  params.time_width = static_cast<unsigned int>(tbins.size());

  unsigned int max_adc_pad = 0;
  double max_adc = -std::numeric_limits<double>::max();
  for (const auto& pad_adc : adc_by_pad)
  {
    if (pad_adc.second > max_adc)
    {
      max_adc = pad_adc.second;
      max_adc_pad = pad_adc.first;
    }
  }

  const unsigned int total_phibins = m_idealPadMap->get_total_phibins(centroid.layer);
  const double pad_phi_width = total_phibins > 0U ? 2.0 * M_PI / static_cast<double>(total_phibins) : 0.0;
  const double cluster_phi = std::atan2(centroid.y, centroid.x);
  const double max_adc_phi = m_idealPadMap->get_phi(static_cast<unsigned int>(side), centroid.layer, max_adc_pad);
  if (pad_phi_width > 0.0 && std::isfinite(cluster_phi) && std::isfinite(max_adc_phi))
  {
    params.phase = wrap_phi(cluster_phi - max_adc_phi) / pad_phi_width;
  }

  return params;
}

TpcPolyClusterizer::Centroid
TpcPolyClusterizer::make_centroid(const std::vector<Point>& points)
{
  Centroid c;
  if (points.empty()) return c;

  double sx = 0.0;
  double sy = 0.0;
  double sz = 0.0;
  for (const Point& p : points)
  {
    sx += p.x;
    sy += p.y;
    sz += p.z;
  }

  const double n = static_cast<double>(points.size());
  c.x = sx / n;
  c.y = sy / n;
  c.z = sz / n;

  double sxx = 0.0;
  double syy = 0.0;
  double szz = 0.0;
  for (const Point& p : points)
  {
    const double dx = p.x - c.x;
    const double dy = p.y - c.y;
    const double dz = p.z - c.z;
    sxx += dx * dx;
    syy += dy * dy;
    szz += dz * dz;
  }

  c.rms_x = std::sqrt(sxx / n);
  c.rms_y = std::sqrt(syy / n);
  c.rms_z = std::sqrt(szz / n);
  c.layer = points.front().layer;
  c.ok = std::isfinite(c.x) && std::isfinite(c.y) && std::isfinite(c.z);
  return c;
}

int TpcPolyClusterizer::process_event(PHCompositeNode*)
{
  if (!m_fullTracks || !m_clusterTracks) return Fun4AllReturnCodes::EVENT_OK;
  m_clusterTracks->Reset();

  const unsigned int nfull = m_fullTracks->size();
  unsigned int nclusters = 0;
  for (unsigned int ifull = 0; ifull < nfull; ++ifull)
  {
    const FullTrack* full = m_fullTracks->get_track(ifull);
    if (!full) continue;

    std::map<unsigned int, std::vector<Point>> points_by_layer;
    for (unsigned int ih = 0; ih < full->size_hit_indices(); ++ih)
    {
      const FullTrack::HitIndex hi = full->get_hit_index(ih);
      Point p;
      if (make_xyz_point(hi.first, hi.second, full->get_side(), p)) points_by_layer[p.layer].push_back(p);
    }
    if (points_by_layer.empty()) continue;

    TpcPolyClusterTrackv1* out = new TpcPolyClusterTrackv1();
    out->set_event(m_event);
    out->set_track_id(m_clusterTracks->size());
    out->set_source_full_track_id(full->get_track_id());
    out->set_side(full->get_side());

    for (const auto& layer_points : points_by_layer)
    {
      const unsigned int layer = layer_points.first;
      const std::vector<Point>& points = layer_points.second;
      const Centroid centroid = make_centroid(points);
      if (!centroid.ok) continue;

      const ClusterParameters params = make_cluster_parameters(points, centroid, full->get_side());
      out->add_cluster(layer, centroid.x, centroid.y, centroid.z,
                       centroid.rms_x, centroid.rms_y, centroid.rms_z,
                       params.adc, params.phi_width, params.time_width, params.phase);
      for (const Point& p : points) out->add_hit_to_last_cluster(p.hitsetkey, p.hitkey, p.x, p.y, p.z);
      ++nclusters;
    }

    if (out->size_clusters() == 0)
    {
      delete out;
      continue;
    }
    m_clusterTracks->add_track(out);
  }

  if (Verbosity() > 0)
  {
    std::cout << Name() << "::process_event - event " << m_event
              << " full_tracks=" << nfull
              << " cluster_tracks=" << m_clusterTracks->size()
              << " layer_clusters=" << nclusters << std::endl;
  }

  ++m_event;
  return Fun4AllReturnCodes::EVENT_OK;
}
