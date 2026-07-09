#pragma once

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

#include <string>
#include <vector>

class FullTrackContainer;
class IdealPadMap;
class PHCompositeNode;
class PHGarfield;
class TpcPolyClusterTrackContainer;
class TrkrHitSetContainer;

class TpcPolyClusterizer : public SubsysReco
{
 public:
  explicit TpcPolyClusterizer(const std::string& name = "TpcPolyClusterizer");
  ~TpcPolyClusterizer() override;

  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;

  void setInputNodeName(const std::string& n) { m_inputNodeName = n; }
  void setOutputNodeName(const std::string& n) { m_outputNodeName = n; }
  void setT0(double v) { m_t0 = v; }
  void setTpcAdcClock(double v) { m_tpcAdcClock = v; }
  void setReverseDriftStepNs(double v) { m_reverseDriftStepNs = v; }
  void setStartZ(double south_z, double north_z)
  {
    m_startZSouth = south_z;
    m_startZNorth = north_z;
  }

 private:
  struct Point
  {
    TrkrDefs::hitsetkey hitsetkey {0};
    TrkrDefs::hitkey hitkey {0};
    unsigned int layer {0};
    unsigned int side {0};
    unsigned int pad {0};
    unsigned int tbin {0};
    double adc {0.0};
    double x {0.0};
    double y {0.0};
    double z {0.0};
  };

  struct Centroid
  {
    bool ok {false};
    unsigned int layer {0};
    double x {0.0};
    double y {0.0};
    double z {0.0};
    double rms_x {0.0};
    double rms_y {0.0};
    double rms_z {0.0};
  };

  struct ClusterParameters
  {
    double adc {0.0};
    unsigned int phi_width {0};
    unsigned int time_width {0};
    double phase {0.0};
  };

  struct ClusterizedTrack
  {
    unsigned int source_full_track_id {0};
    int side {0};
    std::vector<unsigned int> layers;
    std::vector<Centroid> centroids;
    std::vector<ClusterParameters> parameters;
    std::vector<std::vector<Point>> cluster_points;
  };

  int getNodes(PHCompositeNode*);
  int createNodes(PHCompositeNode*);
  bool make_xyz_point(TrkrDefs::hitsetkey hsk, TrkrDefs::hitkey hk, PHGarfield* garfield, Point& p) const;
  ClusterParameters make_cluster_parameters(const std::vector<Point>& points, const Centroid& centroid, int side) const;
  static Centroid make_centroid(const std::vector<Point>& points);
  void cluster_sector_side(unsigned int sector,
                           int side,
                           PHGarfield* garfield,
                           std::vector<ClusterizedTrack>& output,
                           unsigned int& nclusters) const;

  std::string m_inputNodeName;
  std::string m_outputNodeName;
  FullTrackContainer* m_fullTracks {nullptr};
  TpcPolyClusterTrackContainer* m_clusterTracks {nullptr};
  TrkrHitSetContainer* m_hits {nullptr};
  IdealPadMap* m_idealPadMap {nullptr};
  PHGarfield* m_garfield {nullptr};
  std::vector<PHGarfield*> m_workerGarfields;
  unsigned int m_event {0};
  double m_t0 {6};
  double m_tpcAdcClock {56.881262};
  double m_reverseDriftStepNs {56.881262};
  double m_startZSouth {-102.325};
  double m_startZNorth {102.325};
};
