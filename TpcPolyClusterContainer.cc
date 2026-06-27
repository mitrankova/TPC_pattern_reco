#include "TpcPolyClusterContainer.h"

void TpcPolyClusterContainer::identify(std::ostream& os) const
{
  os << "TpcPolyClusterContainer base class" << std::endl;
}

void TpcPolyClusterContainer::Reset()
{
}

int TpcPolyClusterContainer::isValid() const
{
  return 0;
}
