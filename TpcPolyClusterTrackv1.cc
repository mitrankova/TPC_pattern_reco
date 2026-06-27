#include "TpcPolyClusterTrackv1.h"

#include <cmath>
#include <limits>

ClassImp(TpcPolyClusterTrackv1)

TpcPolyClusterTrackv1::TpcPolyClusterTrackv1()
{
  Reset();
}

void TpcPolyClusterTrackv1::identify(std::ostream& os) const
{
  os << "TpcPolyClusterTrackv1:"
     << " event=" << m_event
     << " track_id=" << m_track_id
     << " source_full_track_id=" << m_source_full_track_id
     << " side=" << m_side
     << " clusters=" << m_cluster_layer.size()
     << " hits=" << m_hit_indices.size()
     << std::endl;
}

void TpcPolyClusterTrackv1::Reset()
{
  m_event = 0;
  m_track_id = 0;
  m_source_full_track_id = 0;
  m_side = 0;
  m_cluster_layer.clear();
  m_cluster_first_hit.clear();
  m_cluster_nhits.clear();
  m_cluster_x.clear();
  m_cluster_y.clear();
  m_cluster_z.clear();
  m_cluster_rms_x.clear();
  m_cluster_rms_y.clear();
  m_cluster_rms_z.clear();
  m_hit_indices.clear();
  m_hit_x.clear();
  m_hit_y.clear();
  m_hit_z.clear();
}

int TpcPolyClusterTrackv1::isValid() const
{
  return m_cluster_layer.empty() ? 0 : 1;
}

void TpcPolyClusterTrackv1::add_cluster(const unsigned int layer,
                                        const double x,
                                        const double y,
                                        const double z,
                                        const double rms_x,
                                        const double rms_y,
                                        const double rms_z)
{
  m_cluster_layer.push_back(layer);
  m_cluster_first_hit.push_back(static_cast<unsigned int>(m_hit_indices.size()));
  m_cluster_nhits.push_back(0);
  m_cluster_x.push_back(x);
  m_cluster_y.push_back(y);
  m_cluster_z.push_back(z);
  m_cluster_rms_x.push_back(rms_x);
  m_cluster_rms_y.push_back(rms_y);
  m_cluster_rms_z.push_back(rms_z);
}

void TpcPolyClusterTrackv1::add_hit_to_last_cluster(const TrkrDefs::hitsetkey hsk,
                                                    const TrkrDefs::hitkey hk,
                                                    const double x,
                                                    const double y,
                                                    const double z)
{
  if (m_cluster_nhits.empty()) return;
  m_hit_indices.emplace_back(hsk, hk);
  m_hit_x.push_back(x);
  m_hit_y.push_back(y);
  m_hit_z.push_back(z);
  ++m_cluster_nhits.back();
}

unsigned int TpcPolyClusterTrackv1::flat_hit_index(const unsigned int icluster,
                                                   const unsigned int ihit) const
{
  if (icluster >= m_cluster_first_hit.size() || icluster >= m_cluster_nhits.size()) return std::numeric_limits<unsigned int>::max();
  if (ihit >= m_cluster_nhits[icluster]) return std::numeric_limits<unsigned int>::max();
  return m_cluster_first_hit[icluster] + ihit;
}

TpcPolyClusterTrack::HitIndex
TpcPolyClusterTrackv1::get_cluster_hit_index(const unsigned int icluster,
                                             const unsigned int ihit) const
{
  const unsigned int idx = flat_hit_index(icluster, ihit);
  if (idx >= m_hit_indices.size()) return {0, 0};
  return m_hit_indices[idx];
}

double TpcPolyClusterTrackv1::get_cluster_hit_x(const unsigned int icluster,
                                                const unsigned int ihit) const
{
  const unsigned int idx = flat_hit_index(icluster, ihit);
  return idx < m_hit_x.size() ? m_hit_x[idx] : 0.0;
}

double TpcPolyClusterTrackv1::get_cluster_hit_y(const unsigned int icluster,
                                                const unsigned int ihit) const
{
  const unsigned int idx = flat_hit_index(icluster, ihit);
  return idx < m_hit_y.size() ? m_hit_y[idx] : 0.0;
}

double TpcPolyClusterTrackv1::get_cluster_hit_z(const unsigned int icluster,
                                                const unsigned int ihit) const
{
  const unsigned int idx = flat_hit_index(icluster, ihit);
  return idx < m_hit_z.size() ? m_hit_z[idx] : 0.0;
}
