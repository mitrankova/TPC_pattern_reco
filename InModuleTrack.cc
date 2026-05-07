#include "InModuleTrack.h"

ClassImp(InModuleTrack)

InModuleTrack::InModuleTrack()
{
  Reset();
}

void InModuleTrack::identify(std::ostream& os) const
{
  os << "InModuleTrack: event=" << m_event
     << " track_id=" << m_track_id
     << " region=" << m_region
     << " sector=" << m_sector
     << " side=" << m_side
     << " nblobs=" << m_nblobs
     << " nrawhits=" << m_nrawhits
     << " first_layer=" << m_first_layer
     << " last_layer=" << m_last_layer
     << " hits=" << m_hits.size()
     << std::endl;
}

void InModuleTrack::Reset()
{
  m_event = 0;
  m_track_id = 0;

  m_region = 0;
  m_sector = 0;
  m_side = 0;

  m_nblobs = 0;
  m_nrawhits = 0;

  m_first_layer = 0;
  m_last_layer = 0;

  m_pad_slope = 0.0;
  m_pad_intercept = 0.0;

  m_tbin_slope = 0.0;
  m_tbin_intercept = 0.0;

  m_chi2_pad = 0.0;
  m_chi2_tbin = 0.0;

  m_ndof_pad = 0;
  m_ndof_tbin = 0;

  m_hits.clear();
}

int InModuleTrack::isValid() const
{
  return m_hits.empty() ? 0 : 1;
}