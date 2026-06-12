#include "TpcPolyTrackv1.h"

ClassImp(TpcPolyTrackv1)

TpcPolyTrackv1::TpcPolyTrackv1()
{
  Reset();
}

void TpcPolyTrackv1::identify(std::ostream& os) const
{
  os << "TpcPolyTrackv1:"
     << " event=" << m_event
     << " track_id=" << m_track_id
     << " source_full_track_id=" << m_source_full_track_id
     << " side=" << m_side
     << " fit_status=" << m_fit_status
     << " nhits=" << m_hit_indices.size()
     << " d0=" << m_d0
     << " z0=" << m_z0
     << " phi0=" << m_phi0
     << " theta=" << m_theta
     << " curvature=" << m_curvature
     << " chi2_xy=" << m_chi2_xy
     << " chi2_z=" << m_chi2_z
     << std::endl;
}

void TpcPolyTrackv1::Reset()
{
  m_event = 0;
  m_track_id = 0;
  m_source_full_track_id = 0;
  m_side = 0;
  m_fit_status = 0;
  m_d0 = 0.0;
  m_z0 = 0.0;
  m_phi0 = 0.0;
  m_theta = 0.0;
  m_curvature = 0.0;
  m_chi2_xy = 0.0;
  m_chi2_z = 0.0;
  m_ndof_xy = 0;
  m_ndof_z = 0;
  m_hit_indices.clear();
  m_hit_x.clear();
  m_hit_y.clear();
  m_hit_z.clear();
}

int TpcPolyTrackv1::isValid() const
{
  return m_hit_indices.size() >= 3 ? 1 : 0;
}
