#pragma once

#include "TpcPolyHelixFitter.h"

#include <fun4all/SubsysReco.h>

#include <string>

class FinalTrackContainer;
class IdealPadMap;
class PHCompositeNode;
class TpcPolyClusterTrack;
class TpcPolyClusterTrackContainer;

class TpcPolyClusterTrackReco : public SubsysReco
{
 public:
  enum class FitMode
  {
    Helix,
    Line3D
  };

  explicit TpcPolyClusterTrackReco(const std::string& name = "TpcPolyClusterTrackReco");
  ~TpcPolyClusterTrackReco() override;

  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;

  void setInputNodeName(const std::string& n) { m_inputNodeName = n; }
  void setOutputNodeName(const std::string& n) { m_outputNodeName = n; }
  void setMagneticFieldTesla(double v) { m_magneticFieldTesla = v; }
  void setFitMode(FitMode mode) { m_fitMode = mode; }
  void setUseLine3DFit(bool v) { m_fitMode = v ? FitMode::Line3D : FitMode::Helix; }

 private:
  int getNodes(PHCompositeNode*);
  int createNodes(PHCompositeNode*);
  double calc_dedx(const TpcPolyClusterTrack* track,
                   const TpcPolyHelixFitter::FitResult& fit,
                   bool fit_ok) const;
  void fillFinalTrack(const TpcPolyClusterTrack* in,
                      const TpcPolyHelixFitter::FitResult& fit,
                      bool fit_ok);

  std::string m_inputNodeName;
  std::string m_outputNodeName;
  TpcPolyClusterTrackContainer* m_clusterTracks {nullptr};
  FinalTrackContainer* m_finalTracks {nullptr};
  IdealPadMap* m_idealPadMap {nullptr};
  unsigned int m_event {0};
  double m_magneticFieldTesla {1.4};
  FitMode m_fitMode {FitMode::Helix};
};
