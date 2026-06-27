#include "TpcPolyClusterTrackContainer.h"

void TpcPolyClusterTrackContainer::identify(std::ostream& os) const
{
  os << "TpcPolyClusterTrackContainer base class" << std::endl;
}

void TpcPolyClusterTrackContainer::Reset()
{
}

int TpcPolyClusterTrackContainer::isValid() const
{
  return 0;
}
