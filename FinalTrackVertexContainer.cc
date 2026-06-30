#include "FinalTrackVertexContainer.h"

void FinalTrackVertexContainer::identify(std::ostream& os) const
{
  os << "FinalTrackVertexContainer base class" << std::endl;
}

void FinalTrackVertexContainer::Reset()
{
}

int FinalTrackVertexContainer::isValid() const
{
  return 0;
}
