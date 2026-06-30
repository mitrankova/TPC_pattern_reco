#pragma once

#include "FinalTrackContainer.h"

#include <iostream>
#include <vector>

class FinalTrack;

class FinalTrackContainerv1 : public FinalTrackContainer
{
 public:
  FinalTrackContainerv1();
  ~FinalTrackContainerv1() override;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;
  PHObject* CloneMe() const override;

  unsigned int size() const override { return static_cast<unsigned int>(m_tracks.size()); }
  void add_track(FinalTrack* trk) override { m_tracks.push_back(trk); }
  const FinalTrack* get_track(unsigned int i) const override
  {
    if (i >= m_tracks.size()) return nullptr;
    return m_tracks[i];
  }
  FinalTrack* get_track(unsigned int i) override
  {
    if (i >= m_tracks.size()) return nullptr;
    return m_tracks[i];
  }

 private:
  std::vector<FinalTrack*> m_tracks;

  ClassDefOverride(FinalTrackContainerv1, 1)
};
