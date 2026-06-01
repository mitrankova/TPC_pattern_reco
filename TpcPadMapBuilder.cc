#include "TpcPadMapBuilder.h"

#include "TpcPadMap.h"
#include "TpcPadMapv1.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>

#include <iostream>

typedef PHIODataNode<PHObject> PHObjectNode_t;

TpcPadMapBuilder::TpcPadMapBuilder(const std::string& name,
                                   const std::string& nodeName)
  : SubsysReco(name)
  , m_nodeName(nodeName)
{
}

int TpcPadMapBuilder::InitRun(PHCompositeNode* topNode)
{
  if (!topNode)
  {
    std::cerr << "TpcPadMapBuilder::InitRun - null topNode" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  TpcPadMap* existing = findNode::getClass<TpcPadMap>(topNode, m_nodeName.c_str());
  if (existing)
  {
    if (existing->isValid())
    {
      std::cout << "TpcPadMapBuilder::InitRun - node " << m_nodeName
                << " already exists and is valid" << std::endl;
      return Fun4AllReturnCodes::EVENT_OK;
    }

    std::cerr << "TpcPadMapBuilder::InitRun - node " << m_nodeName
              << " already exists but is invalid" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  PHNodeIterator iter(topNode);
  PHCompositeNode* runNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "RUN"));
  if (!runNode)
  {
    runNode = new PHCompositeNode("RUN");
    topNode->addNode(runNode);
  }

  TpcPadMapv1* padMap = new TpcPadMapv1();
  if (padMap->load_from_cdb(1) != 0)
  {
    std::cerr << "TpcPadMapBuilder::InitRun - failed to load TPC pad map" << std::endl;
    delete padMap;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  PHObjectNode_t* padMapNode = new PHObjectNode_t(padMap, m_nodeName.c_str(), "PHObject");
  runNode->addNode(padMapNode);

  std::cout << "TpcPadMapBuilder::InitRun - added " << m_nodeName
            << " to RUN node" << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}
