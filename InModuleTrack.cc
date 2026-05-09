#include "InModuleTrack.h"

const std::vector<InModuleTrack::HitIndex> InModuleTrack::s_empty_indices;

void InModuleTrack::identify(std::ostream& os) const
{
  os << "InModuleTrack base class" << std::endl;
}

void InModuleTrack::Reset()
{
}

int InModuleTrack::isValid() const
{
  return 0;
}

const std::vector<InModuleTrack::HitIndex>& InModuleTrack::get_hit_indices() const
{
  return s_empty_indices;
}
