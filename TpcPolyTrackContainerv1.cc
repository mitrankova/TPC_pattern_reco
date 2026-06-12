#include "TpcPolyTrackContainerv1.h"
#include "TpcPolyTrack.h"

ClassImp(TpcPolyTrackContainerv1)

TpcPolyTrackContainerv1::TpcPolyTrackContainerv1()
{
  Reset();
}

TpcPolyTrackContainerv1::~TpcPolyTrackContainerv1()
{
  Reset();
}

void TpcPolyTrackContainerv1::identify(std::ostream& os) const
{
  os << "TpcPolyTrackContainerv1 with " << m_tracks.size() << " tracks" << std::endl;
}

void TpcPolyTrackContainerv1::Reset()
{
  for (unsigned int i = 0; i < m_tracks.size(); ++i) delete m_tracks[i];
  m_tracks.clear();
}

int TpcPolyTrackContainerv1::isValid() const
{
  return m_tracks.empty() ? 0 : 1;
}

PHObject* TpcPolyTrackContainerv1::CloneMe() const
{
  TpcPolyTrackContainerv1* copy = new TpcPolyTrackContainerv1();
  for (unsigned int i = 0; i < m_tracks.size(); ++i)
  {
    copy->m_tracks.push_back(static_cast<TpcPolyTrack*>(m_tracks[i]->CloneMe()));
  }
  return copy;
}
