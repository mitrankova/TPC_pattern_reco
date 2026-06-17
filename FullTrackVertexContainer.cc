#include "FullTrackVertexContainer.h"

void FullTrackVertexContainer::identify(std::ostream& os) const
{
  os << "FullTrackVertexContainer base class" << std::endl;
}

void FullTrackVertexContainer::Reset()
{
}

int FullTrackVertexContainer::isValid() const
{
  return 0;
}
