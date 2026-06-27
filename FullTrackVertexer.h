#pragma once

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

#include <string>

class FullTrackContainer;
class FullTrackVertexContainer;
class IdealPadMap;
class PHCompositeNode;
class TrkrHitSetContainer;

class FullTrackVertexer : public SubsysReco
{
 public:
  explicit FullTrackVertexer(const std::string& name = "FullTrackVertexer");
  ~FullTrackVertexer() override;

  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;

  void setInputNodeName(const std::string& n) { m_inputNodeName = n; }
  void setOutputNodeName(const std::string& n) { m_outputNodeName = n; }
  void setFitWeightPower(double v) { m_fitWeightPower = v; }
  void setFitWeightFloorFrac(double v) { m_fitWeightFloorFrac = v; }
  void setUseSagittaPhiFit(bool v) { m_useSagittaPhiFit = v; }
  void setCollisionMinTrackLayers(unsigned int v) { m_collisionMinTrackLayers = v; }
  void setCollisionTimebinSeparation(double v) { m_collisionTimebinSeparation = v; }

 public:
  struct HitPoint
  {
    bool ok {false};
    TrkrDefs::hitsetkey hitsetkey {0};
    TrkrDefs::hitkey hitkey {0};
    unsigned int layer {0};
    unsigned int pad {0};
    unsigned int tbin {0};
    unsigned short adc {0};
    double radius {0.0};
    double global_phi {0.0};
  };

 private:
  bool getNodes(PHCompositeNode* topNode);
  bool createNodes(PHCompositeNode* topNode);
  HitPoint make_hit_point(TrkrDefs::hitsetkey hsk, TrkrDefs::hitkey hk) const;

  std::string m_inputNodeName;
  std::string m_outputNodeName;
  FullTrackContainer* m_fullTracks;
  FullTrackVertexContainer* m_vertices;
  TrkrHitSetContainer* m_hits;
  IdealPadMap* m_idealPadMap;
  double m_fitWeightPower;
  double m_fitWeightFloorFrac;
  bool m_useSagittaPhiFit;
  unsigned int m_collisionMinTrackLayers;
  double m_collisionTimebinSeparation;
};
