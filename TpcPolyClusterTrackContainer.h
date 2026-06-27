#pragma once

#include <phool/PHObject.h>

#include <iostream>

class TpcPolyClusterTrack;

class TpcPolyClusterTrackContainer : public PHObject
{
 public:
  TpcPolyClusterTrackContainer() = default;
  ~TpcPolyClusterTrackContainer() override = default;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;

  virtual unsigned int size() const { return 0; }
  virtual void add_track(TpcPolyClusterTrack*) {}
  virtual const TpcPolyClusterTrack* get_track(unsigned int) const { return nullptr; }
  virtual TpcPolyClusterTrack* get_track(unsigned int) { return nullptr; }

 private:
  ClassDefOverride(TpcPolyClusterTrackContainer, 0)
};
