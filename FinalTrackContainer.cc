#include "FinalTrackContainer.h"

void FinalTrackContainer::identify(std::ostream& os) const
{
  os << "FinalTrackContainer base class" << std::endl;
}

void FinalTrackContainer::Reset()
{
}

int FinalTrackContainer::isValid() const
{
  return 0;
}
