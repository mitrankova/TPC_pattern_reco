#pragma once

#include "TpcPolyHelixFitter.h"

#include <fun4all/SubsysReco.h>

#include <string>

class FinalTrackContainer;
class PHCompositeNode;
class TpcPolyClusterTrack;
class TpcPolyClusterTrackContainer;

class TpcPolyClusterTrackReco : public SubsysReco
{
 public:
  explicit TpcPolyClusterTrackReco(const std::string& name = "TpcPolyClusterTrackReco");
  ~TpcPolyClusterTrackReco() override = default;

  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;

  void setInputNodeName(const std::string& n) { m_inputNodeName = n; }
  void setOutputNodeName(const std::string& n) { m_outputNodeName = n; }
  void setMagneticFieldTesla(double v) { m_magneticFieldTesla = v; }

 private:
  int getNodes(PHCompositeNode*);
  int createNodes(PHCompositeNode*);
  void fillFinalTrack(const TpcPolyClusterTrack* in,
                      const TpcPolyHelixFitter::FitResult& fit,
                      bool fit_ok);

  std::string m_inputNodeName;
  std::string m_outputNodeName;
  TpcPolyClusterTrackContainer* m_clusterTracks {nullptr};
  FinalTrackContainer* m_finalTracks {nullptr};
  unsigned int m_event {0};
  double m_magneticFieldTesla {1.4};
};
