#pragma once

#include <phool/PHObject.h>

#include <iostream>

class TpcPolyTrack;

class TpcPolyTrackContainer : public PHObject
{
 public:
  TpcPolyTrackContainer() = default;
  ~TpcPolyTrackContainer() override = default;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;

  virtual unsigned int size() const { return 0; }
  virtual void add_track(TpcPolyTrack*) {}
  virtual const TpcPolyTrack* get_track(unsigned int) const { return nullptr; }
  virtual TpcPolyTrack* get_track(unsigned int) { return nullptr; }

 private:
  ClassDefOverride(TpcPolyTrackContainer, 0)
};
