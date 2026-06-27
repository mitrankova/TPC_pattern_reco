#pragma once

#include "TpcPolyClusterTrack.h"

#include <iostream>
#include <vector>

class TpcPolyClusterTrackv1 : public TpcPolyClusterTrack
{
 public:
  TpcPolyClusterTrackv1();
  ~TpcPolyClusterTrackv1() override = default;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;
  PHObject* CloneMe() const override { return new TpcPolyClusterTrackv1(*this); }

  unsigned int get_event() const override { return m_event; }
  unsigned int get_track_id() const override { return m_track_id; }
  unsigned int get_source_full_track_id() const override { return m_source_full_track_id; }
  int get_side() const override { return m_side; }
  unsigned int size_clusters() const override { return static_cast<unsigned int>(m_cluster_layer.size()); }
  unsigned int size_hits() const override { return static_cast<unsigned int>(m_hit_indices.size()); }

  unsigned int get_cluster_layer(unsigned int i) const override { return i < m_cluster_layer.size() ? m_cluster_layer[i] : 0; }
  unsigned int get_cluster_nhits(unsigned int i) const override { return i < m_cluster_nhits.size() ? m_cluster_nhits[i] : 0; }
  double get_cluster_x(unsigned int i) const override { return i < m_cluster_x.size() ? m_cluster_x[i] : 0.0; }
  double get_cluster_y(unsigned int i) const override { return i < m_cluster_y.size() ? m_cluster_y[i] : 0.0; }
  double get_cluster_z(unsigned int i) const override { return i < m_cluster_z.size() ? m_cluster_z[i] : 0.0; }
  double get_cluster_rms_x(unsigned int i) const override { return i < m_cluster_rms_x.size() ? m_cluster_rms_x[i] : 0.0; }
  double get_cluster_rms_y(unsigned int i) const override { return i < m_cluster_rms_y.size() ? m_cluster_rms_y[i] : 0.0; }
  double get_cluster_rms_z(unsigned int i) const override { return i < m_cluster_rms_z.size() ? m_cluster_rms_z[i] : 0.0; }
  HitIndex get_cluster_hit_index(unsigned int icluster, unsigned int ihit) const override;
  double get_cluster_hit_x(unsigned int icluster, unsigned int ihit) const override;
  double get_cluster_hit_y(unsigned int icluster, unsigned int ihit) const override;
  double get_cluster_hit_z(unsigned int icluster, unsigned int ihit) const override;

  void set_event(unsigned int v) override { m_event = v; }
  void set_track_id(unsigned int v) override { m_track_id = v; }
  void set_source_full_track_id(unsigned int v) override { m_source_full_track_id = v; }
  void set_side(int v) override { m_side = v; }
  void add_cluster(unsigned int layer, double x, double y, double z,
                   double rms_x, double rms_y, double rms_z) override;
  void add_hit_to_last_cluster(TrkrDefs::hitsetkey hsk, TrkrDefs::hitkey hk,
                               double x, double y, double z) override;

 private:
  unsigned int flat_hit_index(unsigned int icluster, unsigned int ihit) const;

  unsigned int m_event {0};
  unsigned int m_track_id {0};
  unsigned int m_source_full_track_id {0};
  int m_side {0};

  std::vector<unsigned int> m_cluster_layer;
  std::vector<unsigned int> m_cluster_first_hit;
  std::vector<unsigned int> m_cluster_nhits;
  std::vector<double> m_cluster_x;
  std::vector<double> m_cluster_y;
  std::vector<double> m_cluster_z;
  std::vector<double> m_cluster_rms_x;
  std::vector<double> m_cluster_rms_y;
  std::vector<double> m_cluster_rms_z;

  std::vector<HitIndex> m_hit_indices;
  std::vector<double> m_hit_x;
  std::vector<double> m_hit_y;
  std::vector<double> m_hit_z;

  ClassDefOverride(TpcPolyClusterTrackv1, 1)
};
