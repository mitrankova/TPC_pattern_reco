#include "TpcPolyClusterTrack.h"

void TpcPolyClusterTrack::identify(std::ostream& os) const
{
  os << "TpcPolyClusterTrack base class" << std::endl;
}

void TpcPolyClusterTrack::Reset()
{
}

int TpcPolyClusterTrack::isValid() const
{
  return 0;
}
