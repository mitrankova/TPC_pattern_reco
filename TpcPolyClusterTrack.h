#pragma once

#include <phool/PHObject.h>
#include <trackbase/TrkrDefs.h>

#include <iostream>
#include <utility>
#include <vector>

class TpcPolyClusterTrack : public PHObject
{
 public:
  using HitIndex = std::pair<TrkrDefs::hitsetkey, TrkrDefs::hitkey>;

  TpcPolyClusterTrack() = default;
  ~TpcPolyClusterTrack() override = default;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;

  virtual unsigned int get_event() const { return 0; }
  virtual unsigned int get_track_id() const { return 0; }
  virtual unsigned int get_source_full_track_id() const { return 0; }
  virtual int get_side() const { return 0; }
  virtual unsigned int size_clusters() const { return 0; }
  virtual unsigned int size_hits() const { return 0; }

  virtual unsigned int get_cluster_layer(unsigned int) const { return 0; }
  virtual unsigned int get_cluster_nhits(unsigned int) const { return 0; }
  virtual double get_cluster_x(unsigned int) const { return 0.0; }
  virtual double get_cluster_y(unsigned int) const { return 0.0; }
  virtual double get_cluster_z(unsigned int) const { return 0.0; }
  virtual double get_cluster_rms_x(unsigned int) const { return 0.0; }
  virtual double get_cluster_rms_y(unsigned int) const { return 0.0; }
  virtual double get_cluster_rms_z(unsigned int) const { return 0.0; }
  virtual double get_cluster_adc(unsigned int) const { return 0.0; }
  virtual unsigned int get_cluster_phi_width(unsigned int) const { return 0; }
  virtual unsigned int get_cluster_time_width(unsigned int) const { return 0; }
  virtual double get_cluster_phase(unsigned int) const { return 0.0; }
  virtual HitIndex get_cluster_hit_index(unsigned int, unsigned int) const { return {0, 0}; }
  virtual double get_cluster_hit_x(unsigned int, unsigned int) const { return 0.0; }
  virtual double get_cluster_hit_y(unsigned int, unsigned int) const { return 0.0; }
  virtual double get_cluster_hit_z(unsigned int, unsigned int) const { return 0.0; }

  virtual void set_event(unsigned int) {}
  virtual void set_track_id(unsigned int) {}
  virtual void set_source_full_track_id(unsigned int) {}
  virtual void set_side(int) {}
  virtual void add_cluster(unsigned int /*layer*/, double /*x*/, double /*y*/, double /*z*/,
                           double /*rms_x*/, double /*rms_y*/, double /*rms_z*/,
                           double /*adc*/ = 0.0, unsigned int /*phi_width*/ = 0,
                           unsigned int /*time_width*/ = 0, double /*phase*/ = 0.0) {}
  virtual void add_hit_to_last_cluster(TrkrDefs::hitsetkey, TrkrDefs::hitkey,
                                       double /*x*/, double /*y*/, double /*z*/) {}

 private:
  ClassDefOverride(TpcPolyClusterTrack, 0)
};
