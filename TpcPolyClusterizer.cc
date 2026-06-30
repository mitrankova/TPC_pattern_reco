#include "TpcPolyClusterizer.h"

#include "FullTrack.h"
#include "FullTrackContainer.h"
#include "FinalTrackContainerv1.h"
#include "FinalTrackv1.h"
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
#include <phgenfit/Fitter.h>
#include <phgenfit/Measurement.h>
#include <phgenfit/SpacepointMeasurement.h>
#include <phgenfit/Track.h>
#include <phfield/PHFieldUtility.h>
#include <GenFit/MeasuredStateOnPlane.h>
#include <GenFit/RKTrackRep.h>

#include <TMatrixDSym.h>
#include <TGeoManager.h>
#include <TPolyLine3D.h>
#include <TVector3.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <vector>

TpcPolyClusterizer::TpcPolyClusterizer(const std::string& name)
  : SubsysReco(name)
  , m_inputNodeName("FULLTRACKS")
  , m_outputNodeName("TPCPOLYCLUSTERTRACKS")
  , m_finalTrackNodeName("FINALTRACKS")
{
}

TpcPolyClusterizer::~TpcPolyClusterizer()
{
  delete m_idealPadMap;
  m_idealPadMap = nullptr;
  delete m_garfield;
  m_garfield = nullptr;
  m_genfitFitter.reset();
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

  m_genfitFitter.reset();

  if (!gGeoManager)
  {
    std::cerr << Name() << "::InitRun - no active TGeoManager. "
              << "Load Tracking_Geometry before initializing GenFit." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  auto* fieldMap = PHFieldUtility::GetFieldMapNode(nullptr, topNode, Verbosity());
  if (!fieldMap)
  {
    std::cerr << Name() << "::InitRun - failed to get magnetic field map" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_genfitFitter.reset(PHGenFit::Fitter::getInstance(gGeoManager,
                                                     fieldMap,
                                                     "KalmanFitter",
                                                     "RKTrackRep",
                                                     false));
  if (m_genfitFitter) m_genfitFitter->set_verbosity(Verbosity());

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

  m_finalTracks = findNode::getClass<FinalTrackContainer>(topNode, m_finalTrackNodeName);
  if (!m_finalTracks)
  {
    m_finalTracks = new FinalTrackContainerv1();
    PHIODataNode<PHObject>* node = new PHIODataNode<PHObject>(m_finalTracks, m_finalTrackNodeName, "PHObject");
    dstNode->addNode(node);
    std::cout << Name() << "::createNodes - created " << m_finalTrackNodeName << " node" << std::endl;
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
  p.x = x;
  p.y = y;
  p.z = z;
  return true;
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
  c.ok = std::isfinite(c.x) && std::isfinite(c.y) && std::isfinite(c.z);
  return c;
}


FinalTrack*
TpcPolyClusterizer::fit_cluster_track(const FullTrack* full, const TpcPolyClusterTrack* cluster_track) const
{
  if (!full || !cluster_track || !m_genfitFitter) return nullptr;
  if (!full->get_seed_valid() || cluster_track->size_clusters() < 3) return nullptr;

  std::vector<TVector3> positions;
  positions.reserve(cluster_track->size_clusters());
  for (unsigned int icluster = 0; icluster < cluster_track->size_clusters(); ++icluster)
  {
    const double x = cluster_track->get_cluster_x(icluster);
    const double y = cluster_track->get_cluster_y(icluster);
    const double z = cluster_track->get_cluster_z(icluster);
    if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) continue;
    positions.emplace_back(x, y, z);
  }
  if (positions.size() < 3) return nullptr;

  TVector3 seed_pos = positions.front();
  TVector3 seed_mom(full->get_seed_px(), full->get_seed_py(), full->get_seed_pz());
  if (!std::isfinite(seed_mom.X()) || !std::isfinite(seed_mom.Y()) || !std::isfinite(seed_mom.Z()) || seed_mom.Mag() <= 0.0)
  {
    seed_mom = positions.back() - positions.front();
  }
  if (seed_mom.Mag() <= 0.0) return nullptr;
  seed_mom.SetMag(1.0);

  TMatrixDSym seed_cov(6);
  seed_cov.Zero();
  for (unsigned int i = 0; i < 6; ++i)
  {
    for (unsigned int j = 0; j < 6; ++j) seed_cov(i, j) = full->get_seed_cov(i, j);
  }
  bool have_cov = false;
  for (unsigned int i = 0; i < 6; ++i)
  {
    if (seed_cov(i, i) > 0.0 && std::isfinite(seed_cov(i, i))) have_cov = true;
  }
  if (!have_cov)
  {
    seed_cov(0, 0) = 25.0;
    seed_cov(1, 1) = 25.0;
    seed_cov(2, 2) = 25.0;
    seed_cov(3, 3) = 1.0;
    seed_cov(4, 4) = 1.0;
    seed_cov(5, 5) = 1.0;
  }

  std::unique_ptr<genfit::AbsTrackRep> rep(new genfit::RKTrackRep(m_genfitParticleId));
  std::unique_ptr<PHGenFit::Track> gf_track(new PHGenFit::Track(rep.release(), seed_pos, seed_mom, seed_cov));

  std::vector<PHGenFit::Measurement*> measurements;
  measurements.reserve(positions.size());
  for (const TVector3& pos : positions)
  {
    measurements.push_back(new PHGenFit::SpacepointMeasurement(pos, m_genfitMeasurementResolution));
  }
  gf_track->addMeasurements(measurements);

  int fit_status = -1;
  try
  {
    fit_status = m_genfitFitter->processTrack(gf_track.get(), false);
  }
  catch (...)
  {
    return nullptr;
  }
  if (fit_status != 0) return nullptr;

  TVector3 final_pos = seed_pos;
  TVector3 final_mom = gf_track->get_mom();
  TMatrixDSym final_cov(6);
  final_cov.Zero();

  try
  {
    std::unique_ptr<genfit::MeasuredStateOnPlane> state(gf_track->extrapolateToPoint(seed_pos));
    if (state)
    {
      state->getPosMomCov(final_pos, final_mom, final_cov);
    }
  }
  catch (...)
  {
    for (unsigned int i = 0; i < 6; ++i)
    {
      for (unsigned int j = 0; j < 6; ++j) final_cov(i, j) = seed_cov(i, j);
    }
  }

  if (!std::isfinite(final_pos.X()) || !std::isfinite(final_pos.Y()) || !std::isfinite(final_pos.Z()) ||
      !std::isfinite(final_mom.X()) || !std::isfinite(final_mom.Y()) || !std::isfinite(final_mom.Z()))
  {
    return nullptr;
  }

  FinalTrackv1* out = new FinalTrackv1();
  out->set_event(m_event);
  out->set_source_full_track_id(full->get_track_id());
  out->set_fit_status(1);
  out->set_nclusters(static_cast<unsigned int>(positions.size()));
  out->set_x(final_pos.X());
  out->set_y(final_pos.Y());
  out->set_z(final_pos.Z());
  out->set_px(final_mom.X());
  out->set_py(final_mom.Y());
  out->set_pz(final_mom.Z());
  out->set_charge(gf_track->get_charge());
  out->set_chi2(gf_track->get_chi2());
  out->set_ndf(gf_track->get_ndf());
  for (unsigned int i = 0; i < 6; ++i)
  {
    for (unsigned int j = 0; j < 6; ++j)
    {
      out->set_cov(i, j, final_cov(i, j));
    }
  }
  return out;
}

int TpcPolyClusterizer::process_event(PHCompositeNode*)
{
  if (!m_fullTracks || !m_clusterTracks || !m_finalTracks) return Fun4AllReturnCodes::EVENT_OK;
  m_clusterTracks->Reset();
  m_finalTracks->Reset();

  const unsigned int nfull = m_fullTracks->size();
  unsigned int nclusters = 0;
  unsigned int nfinal = 0;
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

      out->add_cluster(layer, centroid.x, centroid.y, centroid.z,
                       centroid.rms_x, centroid.rms_y, centroid.rms_z);
      for (const Point& p : points) out->add_hit_to_last_cluster(p.hitsetkey, p.hitkey, p.x, p.y, p.z);
      ++nclusters;
    }

    if (out->size_clusters() == 0)
    {
      delete out;
      continue;
    }

    FinalTrack* final_track = fit_cluster_track(full, out);
    if (final_track)
    {
      final_track->set_track_id(m_finalTracks->size());
      m_finalTracks->add_track(final_track);
      ++nfinal;
    }

    m_clusterTracks->add_track(out);
  }

  if (Verbosity() > 0)
  {
    std::cout << Name() << "::process_event - event " << m_event
              << " full_tracks=" << nfull
              << " cluster_tracks=" << m_clusterTracks->size()
              << " layer_clusters=" << nclusters
              << " final_tracks=" << nfinal << std::endl;
  }

  ++m_event;
  return Fun4AllReturnCodes::EVENT_OK;
}
