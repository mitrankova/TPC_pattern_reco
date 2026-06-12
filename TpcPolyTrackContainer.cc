#include "TpcPolyTrackContainer.h"

void TpcPolyTrackContainer::identify(std::ostream& os) const
{
  os << "TpcPolyTrackContainer base class" << std::endl;
}

void TpcPolyTrackContainer::Reset()
{
}

int TpcPolyTrackContainer::isValid() const
{
  return 0;
}
