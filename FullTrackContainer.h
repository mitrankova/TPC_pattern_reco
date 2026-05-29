#pragma once

#include <phool/PHObject.h>

#include <iostream>

class FullTrack;

class FullTrackContainer : public PHObject
{
 public:
  FullTrackContainer() = default;
  ~FullTrackContainer() override = default;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;

  virtual unsigned int size() const { return 0; }
  virtual void add_track(FullTrack*) {}
  virtual const FullTrack* get_track(unsigned int) const { return nullptr; }
  virtual FullTrack* get_track(unsigned int) { return nullptr; }

 private:
  ClassDefOverride(FullTrackContainer, 0)
};
