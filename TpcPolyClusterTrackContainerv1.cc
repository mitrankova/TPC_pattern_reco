#include "TpcPolyClusterTrackContainerv1.h"
#include "TpcPolyClusterTrack.h"

ClassImp(TpcPolyClusterTrackContainerv1)

TpcPolyClusterTrackContainerv1::TpcPolyClusterTrackContainerv1()
{
  Reset();
}

TpcPolyClusterTrackContainerv1::~TpcPolyClusterTrackContainerv1()
{
  Reset();
}

void TpcPolyClusterTrackContainerv1::identify(std::ostream& os) const
{
  os << "TpcPolyClusterTrackContainerv1 with " << m_tracks.size() << " cluster tracks" << std::endl;
}

void TpcPolyClusterTrackContainerv1::Reset()
{
  for (unsigned int i = 0; i < m_tracks.size(); ++i) delete m_tracks[i];
  m_tracks.clear();
}

int TpcPolyClusterTrackContainerv1::isValid() const
{
  return m_tracks.empty() ? 0 : 1;
}

PHObject* TpcPolyClusterTrackContainerv1::CloneMe() const
{
  TpcPolyClusterTrackContainerv1* copy = new TpcPolyClusterTrackContainerv1();
  for (unsigned int i = 0; i < m_tracks.size(); ++i)
  {
    if (m_tracks[i]) copy->m_tracks.push_back(static_cast<TpcPolyClusterTrack*>(m_tracks[i]->CloneMe()));
  }
  return copy;
}
