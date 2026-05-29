#include "FullTrack.h"

const std::vector<FullTrack::HitIndex> FullTrack::s_empty_indices;

void FullTrack::identify(std::ostream& os) const
{
  os << "FullTrack base class" << std::endl;
}

void FullTrack::Reset()
{
}

int FullTrack::isValid() const
{
  return 0;
}

const std::vector<FullTrack::HitIndex>& FullTrack::get_hit_indices() const
{
  return s_empty_indices;
}
