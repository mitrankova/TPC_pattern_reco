#pragma once

#include <phool/PHObject.h>

#include <iostream>

class InModuleTrack;

// Abstract container for in-module tracks.
// Concrete storage is in versioned subclasses (InModuleTrackContainerv1, ...).
class InModuleTrackContainer : public PHObject
{
 public:
  InModuleTrackContainer() = default;
  ~InModuleTrackContainer() override = default;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;

  virtual unsigned int size() const { return 0; }

  virtual void add_track(InModuleTrack* /*trk*/) {}

  virtual const InModuleTrack* get_track(unsigned int /*i*/) const { return nullptr; }
  virtual InModuleTrack* get_track(unsigned int /*i*/) { return nullptr; }

 private:
  ClassDefOverride(InModuleTrackContainer, 0)
};
