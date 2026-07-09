#include "TpcPolyClusterTrackReco.h"

#include "FinalTrackContainerv1.h"
#include "FinalTrackv1.h"
#include "TpcPolyClusterTrack.h"
#include "TpcPolyClusterTrackContainer.h"
#include "TpcPolyHelixFitter.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>

#include <cmath>
#include <iostream>
#include <vector>

TpcPolyClusterTrackReco::TpcPolyClusterTrackReco(const std::string& name)
  : SubsysReco(name)
  , m_inputNodeName("TPCPOLYCLUSTERTRACKS")
  , m_outputNodeName("FINALTRACKS")
{
}

int TpcPolyClusterTrackReco::InitRun(PHCompositeNode* topNode)
{
  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK) return Fun4AllReturnCodes::ABORTRUN;
  if (createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK) return Fun4AllReturnCodes::ABORTRUN;

  m_event = 0;
  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcPolyClusterTrackReco::getNodes(PHCompositeNode* topNode)
{
  m_clusterTracks = findNode::getClass<TpcPolyClusterTrackContainer>(topNode, m_inputNodeName);
  if (!m_clusterTracks)
  {
    std::cerr << Name() << "::getNodes - missing " << m_inputNodeName << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcPolyClusterTrackReco::createNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    dstNode = new PHCompositeNode("DST");
    topNode->addNode(dstNode);
  }

  m_finalTracks = findNode::getClass<FinalTrackContainer>(topNode, m_outputNodeName);
  if (!m_finalTracks)
  {
    m_finalTracks = new FinalTrackContainerv1();
    PHIODataNode<PHObject>* node = new PHIODataNode<PHObject>(m_finalTracks, m_outputNodeName, "PHObject");
    dstNode->addNode(node);
    std::cout << Name() << "::createNodes - created " << m_outputNodeName << " node" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void TpcPolyClusterTrackReco::fillFinalTrack(const TpcPolyClusterTrack* in,
                                             const TpcPolyHelixFitter::FitResult& fit,
                                             const bool fit_ok)
{
  FinalTrackv1* out = new FinalTrackv1();
  out->set_event(m_event);
  out->set_track_id(m_finalTracks->size());
  out->set_source_full_track_id(in->get_source_full_track_id());
  out->set_fit_status(fit_ok ? 1 : 0);
  out->set_nclusters(in->size_clusters());

  if (fit_ok)
  {
    if (fit.is_line)
    {
      out->set_x(fit.line_x);
      out->set_y(fit.line_y);
      out->set_z(fit.line_z);
      out->set_px(fit.line_dx);
      out->set_py(fit.line_dy);
      out->set_pz(fit.line_dz);
      out->set_charge(0.0);
    }
    else
    {
      const double sin_phi = std::sin(fit.phi0);
      const double cos_phi = std::cos(fit.phi0);
      out->set_x(-fit.d0 * sin_phi);
      out->set_y(fit.d0 * cos_phi);
      out->set_z(fit.z0);

      const double abs_curvature = std::fabs(fit.curvature);
      const double pt = abs_curvature > 0.0 ? 0.003 * std::fabs(m_magneticFieldTesla) / abs_curvature : 0.0;
      const double tan_theta = std::tan(fit.theta);
      out->set_px(pt * cos_phi);
      out->set_py(pt * sin_phi);
      const double pz = std::fabs(tan_theta) > 1.0e-12 ? (pt / tan_theta) : 0.0;
      out->set_pz(pz);
      out->set_charge(fit.curvature >= 0.0 ? -1.0 : 1.0);
    }
    out->set_chi2(fit.chi2_xy + fit.chi2_z);
    out->set_ndf(static_cast<double>(fit.ndof_xy + fit.ndof_z));
  }

  m_finalTracks->add_track(out);
}

int TpcPolyClusterTrackReco::process_event(PHCompositeNode* topNode)
{
  if (!m_clusterTracks || !m_finalTracks)
  {
    if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK ||
        createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    {
      return Fun4AllReturnCodes::EVENT_OK;
    }
  }

  m_finalTracks->Reset();

  const unsigned int ncluster_tracks = m_clusterTracks->size();
  for (unsigned int itrk = 0; itrk < ncluster_tracks; ++itrk)
  {
    const TpcPolyClusterTrack* cluster_track = m_clusterTracks->get_track(itrk);
    if (!cluster_track) continue;

    std::vector<TpcPolyHelixFitter::Point> fit_points;
    fit_points.reserve(cluster_track->size_clusters());
    for (unsigned int icluster = 0; icluster < cluster_track->size_clusters(); ++icluster)
    {
      TpcPolyHelixFitter::Point fp;
      fp.x = cluster_track->get_cluster_x(icluster);
      fp.y = cluster_track->get_cluster_y(icluster);
      fp.z = cluster_track->get_cluster_z(icluster);
      if (std::isfinite(fp.x) && std::isfinite(fp.y) && std::isfinite(fp.z)) fit_points.push_back(fp);
    }

    TpcPolyHelixFitter::FitResult fit;
    const bool fit_ok = (m_fitMode == FitMode::Line3D) ?
      TpcPolyHelixFitter::fitLine3D(fit_points, fit) :
      TpcPolyHelixFitter::fit(fit_points, fit);
    fillFinalTrack(cluster_track, fit, fit_ok);
  }

  if (Verbosity() > 0)
  {
    std::cout << Name() << "::process_event - event " << m_event
              << " cluster_tracks=" << ncluster_tracks
              << " final_tracks=" << m_finalTracks->size() << std::endl;
  }

  ++m_event;
  return Fun4AllReturnCodes::EVENT_OK;
}
