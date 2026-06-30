#include "FinalTrackContainerv1.h"
#include "FinalTrack.h"

ClassImp(FinalTrackContainerv1)

FinalTrackContainerv1::FinalTrackContainerv1()
{
  Reset();
}

FinalTrackContainerv1::~FinalTrackContainerv1()
{
  Reset();
}

void FinalTrackContainerv1::identify(std::ostream& os) const
{
  os << "FinalTrackContainerv1 with " << m_tracks.size() << " tracks" << std::endl;
}

void FinalTrackContainerv1::Reset()
{
  for (unsigned int i = 0; i < m_tracks.size(); ++i) delete m_tracks[i];
  m_tracks.clear();
}

int FinalTrackContainerv1::isValid() const
{
  return m_tracks.empty() ? 0 : 1;
}

PHObject* FinalTrackContainerv1::CloneMe() const
{
  FinalTrackContainerv1* copy = new FinalTrackContainerv1();
  for (unsigned int i = 0; i < m_tracks.size(); ++i)
  {
    if (m_tracks[i]) copy->m_tracks.push_back(static_cast<FinalTrack*>(m_tracks[i]->CloneMe()));
  }
  return copy;
}
