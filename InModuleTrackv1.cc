#include "InModuleTrackv1.h"

ClassImp(InModuleTrackv1)

InModuleTrackv1::InModuleTrackv1()
{
  Reset();
}

void InModuleTrackv1::identify(std::ostream& os) const
{
  os << "InModuleTrackv1:"
     << " event=" << m_event
     << " track_id=" << m_track_id
     << " region=" << m_region
     << " sector=" << m_sector
     << " side=" << m_side
     << " nblobs=" << m_nblobs
     << " nrawhits=" << m_nrawhits
     << " first_layer=" << m_first_layer
     << " last_layer=" << m_last_layer
     << " nhit_indices=" << m_hit_indices.size()
     << " local_fit_mode=" << m_local_fit_mode
     << " radius_tbin_slope=" << m_radius_tbin_slope
     << " radius_tbin_intercept=" << m_radius_tbin_intercept
     << " chi2_radius_tbin=" << m_chi2_radius_tbin
     << " ndof_radius_tbin=" << m_ndof_radius_tbin
     << " radius_phi_slope=" << m_radius_phi_slope
     << " radius_phi_intercept=" << m_radius_phi_intercept
     << " chi2_radius_phi_line=" << m_chi2_radius_phi_line
     << " ndof_radius_phi_line=" << m_ndof_radius_phi_line
     << " circle_x0=" << m_circle_x0
     << " circle_y0=" << m_circle_y0
     << " circle_radius=" << m_circle_radius
     << " chi2_radius_phi_circle=" << m_chi2_radius_phi_circle
     << " ndof_radius_phi_circle=" << m_ndof_radius_phi_circle
     << " has_sagitta_fit=" << m_has_sagitta_fit
     << " sagitta_S=" << m_sagitta_S
     << " sagitta_x0=" << m_sagitta_x0
     << " sagitta_invR=" << m_sagitta_invR
     << " sagitta_theta=" << m_sagitta_theta
     << " sagitta_b=" << m_sagitta_b
     << " chi2_radius_phi_sagitta=" << m_chi2_radius_phi_sagitta
     << " ndof_radius_phi_sagitta=" << m_ndof_radius_phi_sagitta
     << std::endl;
}
     
    


void InModuleTrackv1::Reset()
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

  m_local_fit_mode = -1;

  m_radius_tbin_slope = 0.0;
  m_radius_tbin_intercept = 0.0;
  m_chi2_radius_tbin = 0.0;
  m_ndof_radius_tbin = 0;

  m_radius_phi_slope = 0.0;
  m_radius_phi_intercept = 0.0;
  m_chi2_radius_phi_line = 0.0;
  m_ndof_radius_phi_line = 0;

  m_circle_x0 = 0.0;
  m_circle_y0 = 0.0;
  m_circle_radius = 0.0;
  m_chi2_radius_phi_circle = 0.0;
  m_ndof_radius_phi_circle = 0;

  m_has_sagitta_fit = 0;
  m_sagitta_S = 0.0;
  m_sagitta_x0 = 0.0;
  m_sagitta_invR = 0.0;
  m_sagitta_theta = 0.0;
  m_sagitta_b = 0.0;
  m_chi2_radius_phi_sagitta = 0.0;
  m_ndof_radius_phi_sagitta = 0;

  m_hit_indices.clear();
}

int InModuleTrackv1::isValid() const
{
  return m_hit_indices.empty() ? 0 : 1;
}
