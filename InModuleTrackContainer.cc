#include "InModuleTrackContainer.h"

void InModuleTrackContainer::identify(std::ostream& os) const
{
  os << "InModuleTrackContainer base class" << std::endl;
}

void InModuleTrackContainer::Reset()
{
}

int InModuleTrackContainer::isValid() const
{
  return 0;
}
