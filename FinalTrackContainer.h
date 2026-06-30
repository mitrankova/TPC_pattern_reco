#pragma once

#include <phool/PHObject.h>

#include <iostream>

class FinalTrack;

class FinalTrackContainer : public PHObject
{
 public:
  FinalTrackContainer() = default;
  ~FinalTrackContainer() override = default;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;

  virtual unsigned int size() const { return 0; }
  virtual void add_track(FinalTrack*) {}
  virtual const FinalTrack* get_track(unsigned int) const { return nullptr; }
  virtual FinalTrack* get_track(unsigned int) { return nullptr; }

 private:
  ClassDefOverride(FinalTrackContainer, 0)
};
