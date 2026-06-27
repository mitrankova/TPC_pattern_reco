#pragma once

#include "TpcPolyClusterTrackContainer.h"

#include <iostream>
#include <vector>

class TpcPolyClusterTrack;

class TpcPolyClusterTrackContainerv1 : public TpcPolyClusterTrackContainer
{
 public:
  TpcPolyClusterTrackContainerv1();
  ~TpcPolyClusterTrackContainerv1() override;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;
  PHObject* CloneMe() const override;

  unsigned int size() const override { return static_cast<unsigned int>(m_tracks.size()); }
  void add_track(TpcPolyClusterTrack* trk) override { m_tracks.push_back(trk); }
  const TpcPolyClusterTrack* get_track(unsigned int i) const override
  {
    if (i >= m_tracks.size()) return nullptr;
    return m_tracks[i];
  }
  TpcPolyClusterTrack* get_track(unsigned int i) override
  {
    if (i >= m_tracks.size()) return nullptr;
    return m_tracks[i];
  }

 private:
  std::vector<TpcPolyClusterTrack*> m_tracks;

  ClassDefOverride(TpcPolyClusterTrackContainerv1, 1)
};
