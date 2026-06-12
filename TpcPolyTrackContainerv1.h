#pragma once

#include "TpcPolyTrackContainer.h"

#include <iostream>
#include <vector>

class TpcPolyTrack;

class TpcPolyTrackContainerv1 : public TpcPolyTrackContainer
{
 public:
  TpcPolyTrackContainerv1();
  ~TpcPolyTrackContainerv1() override;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;
  PHObject* CloneMe() const override;

  unsigned int size() const override { return static_cast<unsigned int>(m_tracks.size()); }
  void add_track(TpcPolyTrack* trk) override { m_tracks.push_back(trk); }
  const TpcPolyTrack* get_track(unsigned int i) const override
  {
    if (i >= m_tracks.size()) return nullptr;
    return m_tracks[i];
  }
  TpcPolyTrack* get_track(unsigned int i) override
  {
    if (i >= m_tracks.size()) return nullptr;
    return m_tracks[i];
  }

 private:
  std::vector<TpcPolyTrack*> m_tracks;

  ClassDefOverride(TpcPolyTrackContainerv1, 1)
};
