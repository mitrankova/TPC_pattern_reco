#include "FullTrackVertex.h"

void FullTrackVertex::identify(std::ostream& os) const
{
  os << "FullTrackVertex base class" << std::endl;
}

void FullTrackVertex::Reset()
{
}

int FullTrackVertex::isValid() const
{
  return 0;
}
