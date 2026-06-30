#pragma once

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

#include <memory>
#include <string>
#include <vector>

class FinalTrack;
class FullTrack;
class FullTrackContainer;
class IdealPadMap;
class PHCompositeNode;
class PHGarfield;
class FinalTrackContainer;
namespace PHGenFit { class Fitter; }
class TpcPolyClusterTrack;
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
  void setFinalTrackNodeName(const std::string& n) { m_finalTrackNodeName = n; }
  void setGenFitParticleId(int pid) { m_genfitParticleId = pid; }
  void setGenFitMeasurementResolution(double v) { m_genfitMeasurementResolution = v; }
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

  int getNodes(PHCompositeNode*);
  int createNodes(PHCompositeNode*);
  bool make_xyz_point(TrkrDefs::hitsetkey hsk, TrkrDefs::hitkey hk, int side, Point& p) const;
  static Centroid make_centroid(const std::vector<Point>& points);
  FinalTrack* fit_cluster_track(const FullTrack* full, const TpcPolyClusterTrack* cluster_track) const;

  std::string m_inputNodeName;
  std::string m_outputNodeName;
  std::string m_finalTrackNodeName;
  FullTrackContainer* m_fullTracks {nullptr};
  TpcPolyClusterTrackContainer* m_clusterTracks {nullptr};
  FinalTrackContainer* m_finalTracks {nullptr};
  TrkrHitSetContainer* m_hits {nullptr};
  IdealPadMap* m_idealPadMap {nullptr};
  PHGarfield* m_garfield {nullptr};
  std::shared_ptr<PHGenFit::Fitter> m_genfitFitter;
  unsigned int m_event {0};
  double m_t0 {6};
  double m_tpcAdcClock {56.881262};
  double m_reverseDriftStepNs {56.881262};
  double m_startZSouth {-102.325};
  double m_startZNorth {102.325};
  int m_genfitParticleId {-13};
  double m_genfitMeasurementResolution {0.2};
};
