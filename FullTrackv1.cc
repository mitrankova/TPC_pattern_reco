#include "FullTrackv1.h"

ClassImp(FullTrackv1)

FullTrackv1::FullTrackv1()
{
  Reset();
}

void FullTrackv1::identify(std::ostream& os) const
{
  os << "FullTrackv1:"
     << " event=" << m_event
     << " track_id=" << m_track_id
     << " side=" << m_side
     << " nsegments=" << m_nsegments
     << " nrawhits=" << m_nrawhits
     << " first_layer=" << m_first_layer
     << " last_layer=" << m_last_layer
     << " first_region=" << m_first_region
     << " last_region=" << m_last_region
     << " first_sector=" << m_first_sector
     << " last_sector=" << m_last_sector
     << " phi_slope=" << m_phi_slope
     << " phi_intercept=" << m_phi_intercept
     << " tbin_slope=" << m_tbin_slope
     << " tbin_intercept=" << m_tbin_intercept
     << " chi2_phi=" << m_chi2_phi
     << " chi2_tbin=" << m_chi2_tbin
     << " ndof_phi=" << m_ndof_phi
     << " ndof_tbin=" << m_ndof_tbin
     << " source_tracks=" << m_source_track_ids.size()
     << " hit_indices=" << m_hit_indices.size()
     << std::endl;
}

void FullTrackv1::Reset()
{
  m_event = 0;
  m_track_id = 0;
  m_side = 0;

  m_nsegments = 0;
  m_nrawhits = 0;
  m_first_layer = 0;
  m_last_layer = 0;
  m_first_sector = 0;
  m_last_sector = 0;
  m_first_region = 0;
  m_last_region = 0;

  m_phi_slope = 0.0;
  m_phi_intercept = 0.0;
  m_tbin_slope = 0.0;
  m_tbin_intercept = 0.0;
  m_chi2_phi = 0.0;
  m_chi2_tbin = 0.0;
  m_ndof_phi = 0;
  m_ndof_tbin = 0;

  m_source_track_ids.clear();
  m_source_regions.clear();
  m_source_sectors.clear();
  m_hit_indices.clear();
}

int FullTrackv1::isValid() const
{
  return (m_nsegments > 0 && !m_hit_indices.empty()) ? 1 : 0;
}
