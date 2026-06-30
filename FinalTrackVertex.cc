#include "FinalTrackVertex.h"

void FinalTrackVertex::identify(std::ostream& os) const
{
  os << "FinalTrackVertex base class" << std::endl;
}

void FinalTrackVertex::Reset()
{
}

int FinalTrackVertex::isValid() const
{
  return 0;
}
