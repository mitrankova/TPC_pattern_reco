#include "InModuleTrackHit.h"

ClassImp(InModuleTrackHit)

InModuleTrackHit::InModuleTrackHit()
{
  Reset();
}

void InModuleTrackHit::Reset()
{
  m_event = 0;
  m_track_id = 0;

  m_region = 0;
  m_sector = 0;
  m_side = 0;

  m_layer = 0;

  m_hitsetkey = 0;
  m_hitkey = 0;

  m_pad = 0.0;
  m_tbin = 0.0;
  m_adc = 0.0;

  m_local_phi = 0.0;
  m_local_radius = 0.0;
}