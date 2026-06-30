#include "FinalTrackVertexv1.h"

#include <cmath>

ClassImp(FinalTrackVertexv1)

FinalTrackVertexv1::FinalTrackVertexv1()
{
  Reset();
}

void FinalTrackVertexv1::identify(std::ostream& os) const
{
  os << "FinalTrackVertexv1:"
     << " track_id=" << m_track_id
     << " source_full_track_id=" << m_source_full_track_id
     << " dca2d=" << m_dca2d
     << " z0=" << m_z0
     << " pca_valid=" << m_pca_valid
     << " pca=(" << m_pca_x << ", " << m_pca_y << ", " << m_pca_z << ")"
     << std::endl;
}

void FinalTrackVertexv1::Reset()
{
  m_track_id = 0;
  m_source_full_track_id = 0;
  m_dca2d = 0.0;
  m_z0 = 0.0;
  m_pca_valid = 0;
  m_pca_x = 0.0;
  m_pca_y = 0.0;
  m_pca_z = 0.0;
  m_pca_radius = 0.0;
  m_pca_phi = 0.0;
}

int FinalTrackVertexv1::isValid() const
{
  return m_pca_valid && std::isfinite(m_z0);
}
