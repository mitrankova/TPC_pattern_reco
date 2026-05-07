#pragma once

#include "InModuleTrack.h"

#include <phool/PHObject.h>

#include <vector>
#include <iostream>

class InModuleTrackContainer : public PHObject
{
 public:
  InModuleTrackContainer();
  virtual ~InModuleTrackContainer() {}

  void identify(std::ostream& os = std::cout) const;
  void Reset();
  int isValid() const;
  PHObject* CloneMe() const { return new InModuleTrackContainer(*this); }

  unsigned int size() const
  {
    return static_cast<unsigned int>(m_tracks.size());
  }

  void add_track(const InModuleTrack& trk)
  {
    m_tracks.push_back(trk);
  }

  const InModuleTrack* get_track(unsigned int i) const
  {
    if (i >= m_tracks.size()) return 0;
    return &m_tracks[i];
  }

  InModuleTrack* get_track(unsigned int i)
  {
    if (i >= m_tracks.size()) return 0;
    return &m_tracks[i];
  }

  const std::vector<InModuleTrack>& get_tracks() const
  {
    return m_tracks;
  }

 private:
  std::vector<InModuleTrack> m_tracks;

  ClassDef(InModuleTrackContainer, 1)
};