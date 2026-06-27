#include "TpcPolyClusterv1.h"

#include <cmath>

ClassImp(TpcPolyClusterv1)

TpcPolyClusterv1::TpcPolyClusterv1()
{
  Reset();
}

void TpcPolyClusterv1::identify(std::ostream& os) const
{
  os << "TpcPolyClusterv1:"
     << " event=" << m_event
     << " cluster_id=" << m_cluster_id
     << " source_full_track_id=" << m_source_full_track_id
     << " side=" << m_side
     << " nhits=" << m_hit_indices.size()
     << " centroid_x=" << m_centroid_x
     << " centroid_y=" << m_centroid_y
     << " centroid_z=" << m_centroid_z
     << " rms_x=" << m_rms_x
     << " rms_y=" << m_rms_y
     << " rms_z=" << m_rms_z
     << std::endl;
}

void TpcPolyClusterv1::Reset()
{
  m_event = 0;
  m_cluster_id = 0;
  m_source_full_track_id = 0;
  m_side = 0;
  m_centroid_x = 0.0;
  m_centroid_y = 0.0;
  m_centroid_z = 0.0;
  m_rms_x = 0.0;
  m_rms_y = 0.0;
  m_rms_z = 0.0;
  m_hit_indices.clear();
  m_hit_x.clear();
  m_hit_y.clear();
  m_hit_z.clear();
}

int TpcPolyClusterv1::isValid() const
{
  return (!m_hit_indices.empty() &&
          std::isfinite(m_centroid_x) &&
          std::isfinite(m_centroid_y) &&
          std::isfinite(m_centroid_z)) ? 1 : 0;
}
