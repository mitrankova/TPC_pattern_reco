#pragma once

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

#include <string>
#include <vector>

class FullTrackContainer;
class IdealPadMap;
class PHCompositeNode;
class PHGarfield;
class TpcPolyTrackContainer;
class TrkrHitSetContainer;

class TpcPolyTrackReco : public SubsysReco
{
 public:
  explicit TpcPolyTrackReco(const std::string& name = "TpcPolyTrackReco");
  ~TpcPolyTrackReco() override;

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
    double x {0.0};
    double y {0.0};
    double z {0.0};
  };

  int getNodes(PHCompositeNode*);
  int createNodes(PHCompositeNode*);
  bool make_xyz_point(TrkrDefs::hitsetkey hsk, TrkrDefs::hitkey hk, int side, Point& p) const;

  std::string m_inputNodeName;
  std::string m_outputNodeName;
  FullTrackContainer* m_fullTracks {nullptr};
  TpcPolyTrackContainer* m_polyTracks {nullptr};
  TrkrHitSetContainer* m_hits {nullptr};
  IdealPadMap* m_idealPadMap {nullptr};
  PHGarfield* m_garfield {nullptr};
  unsigned int m_event {0};
  double m_t0 {6};//329.0
  double m_tpcAdcClock {56.881262};//50.037280//53.326184//56.881262
  double m_reverseDriftStepNs {56.881262};
  double m_startZSouth {-102.325};
  double m_startZNorth {102.325};
};
