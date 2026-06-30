#include "FinalTrack.h"

void FinalTrack::identify(std::ostream& os) const
{
  os << "FinalTrack base class" << std::endl;
}

void FinalTrack::Reset()
{
}

int FinalTrack::isValid() const
{
  return 0;
}
