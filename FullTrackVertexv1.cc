#include "FullTrackVertexv1.h"

#include <cmath>

ClassImp(FullTrackVertexv1)

FullTrackVertexv1::FullTrackVertexv1()
{
  Reset();
}

void FullTrackVertexv1::identify(std::ostream& os) const
{
  os << "FullTrackVertexv1:"
     << " track_id=" << m_track_id
     << " d0=" << m_d0
     << " timebin0=" << m_timebin0
     << std::endl;
}

void FullTrackVertexv1::Reset()
{
  m_track_id = 0;
  m_d0 = 0.0;
  m_timebin0 = 0.0;
}

int FullTrackVertexv1::isValid() const
{
  return (std::isfinite(m_d0) && std::isfinite(m_timebin0)) ? 1 : 0;
}
