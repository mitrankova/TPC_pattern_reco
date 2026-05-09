#include "InModuleTrackContainerv1.h"
#include "InModuleTrack.h"

ClassImp(InModuleTrackContainerv1)

InModuleTrackContainerv1::InModuleTrackContainerv1()
{
  Reset();
}

InModuleTrackContainerv1::~InModuleTrackContainerv1()
{
  Reset();
}

void InModuleTrackContainerv1::identify(std::ostream& os) const
{
  os << "InModuleTrackContainerv1 with "
     << m_tracks.size()
     << " in-module tracks"
     << std::endl;
}

void InModuleTrackContainerv1::Reset()
{
  for (auto* trk : m_tracks)
  {
    delete trk;
  }
  m_tracks.clear();
}

int InModuleTrackContainerv1::isValid() const
{
  return m_tracks.empty() ? 0 : 1;
}

PHObject* InModuleTrackContainerv1::CloneMe() const
{
  // Deep copy: clone each owned track.
  auto* copy = new InModuleTrackContainerv1();
  for (const auto* trk : m_tracks)
  {
    copy->m_tracks.push_back(static_cast<InModuleTrack*>(trk->CloneMe()));
  }
  return copy;
}
