#include "InModuleTrackContainer.h"

ClassImp(InModuleTrackContainer)

InModuleTrackContainer::InModuleTrackContainer()
{
  Reset();
}

void InModuleTrackContainer::identify(std::ostream& os) const
{
  os << "InModuleTrackContainer with "
     << m_tracks.size()
     << " in-module tracks"
     << std::endl;
}

void InModuleTrackContainer::Reset()
{
  m_tracks.clear();
}

int InModuleTrackContainer::isValid() const
{
  return m_tracks.empty() ? 0 : 1;
}