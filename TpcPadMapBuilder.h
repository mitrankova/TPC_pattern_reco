#ifndef TPCPADMAPBUILDER_H
#define TPCPADMAPBUILDER_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;

class TpcPadMapBuilder : public SubsysReco
{
 public:
  TpcPadMapBuilder(const std::string& name = "TpcPadMapBuilder",
                   const std::string& nodeName = "TPC_PADMAP");
  ~TpcPadMapBuilder() override {}

  int InitRun(PHCompositeNode* topNode) override;

 private:
  std::string m_nodeName;
};

#endif
