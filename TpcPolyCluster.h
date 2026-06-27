#pragma once

#include <phool/PHObject.h>
#include <trackbase/TrkrDefs.h>

#include <iostream>
#include <utility>
#include <vector>

class TpcPolyCluster : public PHObject
{
 public:
  TpcPolyCluster() = default;
  ~TpcPolyCluster() override = default;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;

  using HitIndex = std::pair<TrkrDefs::hitsetkey, TrkrDefs::hitkey>;

  virtual unsigned int get_event() const { return 0; }
  virtual unsigned int get_cluster_id() const { return 0; }
  virtual unsigned int get_source_full_track_id() const { return 0; }
  virtual int get_side() const { return 0; }
  virtual unsigned int get_nhits() const { return 0; }

  virtual double get_centroid_x() const { return 0.0; }
  virtual double get_centroid_y() const { return 0.0; }
  virtual double get_centroid_z() const { return 0.0; }
  virtual double get_rms_x() const { return 0.0; }
  virtual double get_rms_y() const { return 0.0; }
  virtual double get_rms_z() const { return 0.0; }

  virtual void set_event(unsigned int) {}
  virtual void set_cluster_id(unsigned int) {}
  virtual void set_source_full_track_id(unsigned int) {}
  virtual void set_side(int) {}
  virtual void set_centroid_x(double) {}
  virtual void set_centroid_y(double) {}
  virtual void set_centroid_z(double) {}
  virtual void set_rms_x(double) {}
  virtual void set_rms_y(double) {}
  virtual void set_rms_z(double) {}

  virtual void add_hit(TrkrDefs::hitsetkey, TrkrDefs::hitkey,
                       double /*x*/, double /*y*/, double /*z*/) {}
  virtual unsigned int size_hits() const { return 0; }
  virtual HitIndex get_hit_index(unsigned int) const { return {0, 0}; }
  virtual double get_hit_x(unsigned int) const { return 0.0; }
  virtual double get_hit_y(unsigned int) const { return 0.0; }
  virtual double get_hit_z(unsigned int) const { return 0.0; }
  virtual const std::vector<HitIndex>& get_hit_indices() const;

 private:
  static const std::vector<HitIndex> s_empty_indices;

  ClassDefOverride(TpcPolyCluster, 0)
};
