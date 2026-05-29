#include "FullTrackContainerv1.h"
#include "FullTrack.h"

ClassImp(FullTrackContainerv1)

FullTrackContainerv1::FullTrackContainerv1()
{
  Reset();
}

FullTrackContainerv1::~FullTrackContainerv1()
{
  Reset();
}

void FullTrackContainerv1::identify(std::ostream& os) const
{
  os << "FullTrackContainerv1 with "
     << m_tracks.size()
     << " full tracks"
     << std::endl;
}

void FullTrackContainerv1::Reset()
{
  for (unsigned int i = 0; i < m_tracks.size(); ++i)
  {
    delete m_tracks[i];
  }
  m_tracks.clear();
}

int FullTrackContainerv1::isValid() const
{
  return m_tracks.empty() ? 0 : 1;
}

PHObject* FullTrackContainerv1::CloneMe() const
{
  FullTrackContainerv1* copy = new FullTrackContainerv1();
  for (unsigned int i = 0; i < m_tracks.size(); ++i)
  {
    copy->m_tracks.push_back(static_cast<FullTrack*>(m_tracks[i]->CloneMe()));
  }
  return copy;
}
