#pragma once

#include "InModuleTrackContainer.h"

#include <iostream>
#include <vector>

class InModuleTrack;

class InModuleTrackContainerv1 : public InModuleTrackContainer
{
 public:
  InModuleTrackContainerv1();
  ~InModuleTrackContainerv1() override;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;
  PHObject* CloneMe() const override;

  unsigned int size() const override
  {
    return static_cast<unsigned int>(m_tracks.size());
  }

  // Takes ownership of the track pointer.
  void add_track(InModuleTrack* trk) override
  {
    m_tracks.push_back(trk);
  }

  const InModuleTrack* get_track(unsigned int i) const override
  {
    if (i >= m_tracks.size()) return nullptr;
    return m_tracks[i];
  }

  InModuleTrack* get_track(unsigned int i) override
  {
    if (i >= m_tracks.size()) return nullptr;
    return m_tracks[i];
  }

 private:
  std::vector<InModuleTrack*> m_tracks;

  ClassDefOverride(InModuleTrackContainerv1, 1)
};
