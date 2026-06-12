#include "TpcPolyTrack.h"

const std::vector<TpcPolyTrack::HitIndex> TpcPolyTrack::s_empty_indices;

void TpcPolyTrack::identify(std::ostream& os) const
{
  os << "TpcPolyTrack base class" << std::endl;
}

void TpcPolyTrack::Reset()
{
}

int TpcPolyTrack::isValid() const
{
  return 0;
}

const std::vector<TpcPolyTrack::HitIndex>& TpcPolyTrack::get_hit_indices() const
{
  return s_empty_indices;
}
