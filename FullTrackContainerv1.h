#pragma once

#include "FullTrackContainer.h"

#include <iostream>
#include <vector>

class FullTrack;

class FullTrackContainerv1 : public FullTrackContainer
{
 public:
  FullTrackContainerv1();
  ~FullTrackContainerv1() override;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;
  PHObject* CloneMe() const override;

  unsigned int size() const override
  {
    return static_cast<unsigned int>(m_tracks.size());
  }

  void add_track(FullTrack* trk) override
  {
    m_tracks.push_back(trk);
  }

  const FullTrack* get_track(unsigned int i) const override
  {
    if (i >= m_tracks.size()) return nullptr;
    return m_tracks[i];
  }

  FullTrack* get_track(unsigned int i) override
  {
    if (i >= m_tracks.size()) return nullptr;
    return m_tracks[i];
  }

 private:
  std::vector<FullTrack*> m_tracks;

  ClassDefOverride(FullTrackContainerv1, 1)
};
