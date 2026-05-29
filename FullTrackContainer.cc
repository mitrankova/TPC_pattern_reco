#include "FullTrackContainer.h"

void FullTrackContainer::identify(std::ostream& os) const
{
  os << "FullTrackContainer base class" << std::endl;
}

void FullTrackContainer::Reset()
{
}

int FullTrackContainer::isValid() const
{
  return 0;
}
