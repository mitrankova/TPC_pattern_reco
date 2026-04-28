#include "InModuleTracks.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <TFile.h>
#include <TTree.h>

#include <pthread.h>

#include <algorithm>
#include <iostream>
#include <stdint.h>
#include <vector>

// ====================================================================
// InModuleThreadData constructors
// ====================================================================

InModuleThreadData::InModuleThreadData()
  : layer(0)
  , region(0)
  , sector(0)
  , side(0)
  , hitsetkey(0)
  , module_key(0)
  , hitset(0)
  , tGeometry(0)
  , pedestal(0.0)
  , verbosity(0)
{
}

InModuleThreadData::HitInfo::HitInfo()
  : layer(0)
  , region(0)
  , sector(0)
  , side(0)
  , pad(0)
  , tbin(0)
  , adc(0)
{
}

// ====================================================================
// pthread worker
// ====================================================================

namespace
{
  void* ProcessModule(void* arg)
  {
    InModuleThreadData* d = static_cast<InModuleThreadData*>(arg);

    if (!d)
    {
      return 0;
    }

    if (!d->hitset)
    {
      return 0;
    }

    TrkrHitSet::ConstRange range = d->hitset->getHits();

    for (TrkrHitSet::ConstIterator hitr = range.first;
         hitr != range.second;
         ++hitr)
    {
      const TrkrDefs::hitkey hitkey = hitr->first;
      const TrkrHit* hit = hitr->second;

      if (!hit)
      {
        continue;
      }

      const unsigned short pad = TpcDefs::getPad(hitkey);
      const unsigned short tbin = TpcDefs::getTBin(hitkey);
      const unsigned short rawAdc = hit->getAdc();

      const double fadc = static_cast<double>(rawAdc) - d->pedestal;

      if (fadc <= 0.0)
      {
        continue;
      }

      InModuleThreadData::HitInfo hi;

      hi.region = d->region;
      hi.sector = d->sector;
      hi.side = d->side;
      hi.pad = pad;
      hi.tbin = tbin;
      hi.adc = static_cast<unsigned short>(fadc);

      d->hits.push_back(hi);
    }

    return 0;
  }
}

// ====================================================================
// InModuleTracks implementation
// ====================================================================

InModuleTracks::InModuleTracks(const std::string& name,
                               const std::string& filename)
  : SubsysReco(name)
  , m_outputFileName(filename)
  , m_outputFile(0)
  , m_tree(0)
  , m_hits(0)
  , m_tGeometry(0)
  , m_event(0)
  , m_maxThreads(72)
  , m_tree_event(0)
{
}

InModuleTracks::~InModuleTracks()
{
  if (m_outputFile)
  {
    m_outputFile->Close();
    delete m_outputFile;
    m_outputFile = 0;
  }
}

void InModuleTracks::setMaxThreads(unsigned int n)
{
  if (n == 0)
  {
    m_maxThreads = 1;
  }
  else
  {
    m_maxThreads = n;
  }
}

int InModuleTracks::Init(PHCompositeNode* /*topNode*/)
{
  m_outputFile = new TFile(m_outputFileName.c_str(), "RECREATE");

  if (!m_outputFile || m_outputFile->IsZombie())
  {
    std::cerr << Name()
              << "::Init - cannot create output file "
              << m_outputFileName
              << std::endl;

    return Fun4AllReturnCodes::ABORTRUN;
  }
  else{
    std::cout << Name()
              << "!!!!!!!!!!::Init - output file "
              << m_outputFileName
              << " created"
              << std::endl;
  }

  m_tree = new TTree("InModuleTracks", "TPC in-module hit information");

  m_tree->Branch("event", &m_tree_event, "event/I");

  m_tree->Branch("region", &m_tree_region);
  m_tree->Branch("sector", &m_tree_sector);
  m_tree->Branch("side", &m_tree_side);

  m_tree->Branch("pad", &m_tree_pad);
  m_tree->Branch("tbin", &m_tree_tbin);
  m_tree->Branch("adc", &m_tree_adc);

  return Fun4AllReturnCodes::EVENT_OK;
}

int InModuleTracks::InitRun(PHCompositeNode* topNode)
{
  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_event = 0;

  return Fun4AllReturnCodes::EVENT_OK;
}

int InModuleTracks::End(PHCompositeNode* /*topNode*/)
{
  if (m_outputFile)
  {
    m_outputFile->cd();

    if (m_tree)
    {
      m_tree->Write();
    }

    m_outputFile->Close();
    delete m_outputFile;
    m_outputFile = 0;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int InModuleTracks::getNodes(PHCompositeNode* topNode)
{
  m_hits = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");

  if (!m_hits)
  {
    std::cerr << Name()
              << "::getNodes - missing TRKR_HITSET"
              << std::endl;

    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");

  if (!m_tGeometry)
  {
    std::cerr << Name()
              << "::getNodes - missing ActsGeometry"
              << std::endl;

    return Fun4AllReturnCodes::ABORTRUN;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void InModuleTracks::reset_tree_vars()
{
  m_tree_event = m_event;

  m_tree_region.clear();
  m_tree_sector.clear();
  m_tree_side.clear();

  m_tree_pad.clear();
  m_tree_tbin.clear();
  m_tree_adc.clear();
}

// ====================================================================
// process_event
// ====================================================================

int InModuleTracks::process_event(PHCompositeNode* /*topNode*/)
{
  std::cout<< Name()
           << "!!!!!!!!!!!::process_event - processing event " << m_event
           << std::endl;
  reset_tree_vars();

  std::vector<InModuleThreadData> tdata;
  tdata.reserve(72);

  // ------------------------------------------------------------
  // Build one work package per module key:
  //
  //   3 regions/modules x 12 sectors x 2 sides = 72 jobs
  //
  // The key is generated exactly as:
  //
  //   TpcDefs::genModuleHitSetKey(region, sector, side)
  //
  // ------------------------------------------------------------

  for (unsigned int side = 0; side < 2; ++side)
  {
    for (unsigned int sector = 0; sector < 12; ++sector)
    {
      for (unsigned int region = 0; region < 3; ++region)
      {
        for (unsigned int l = 0; l < 16; ++l)
        {

            int layer = region * 16 + l + 7;
            const TrkrDefs::hitsetkey hitset_key = TpcDefs::genHitSetKey(layer, sector, side);
            const TrkrDefs::hitsetkey key =
              TpcDefs::genModuleHitSetKey(static_cast<uint8_t>(region),
                                          static_cast<uint8_t>(sector),
                                          static_cast<uint8_t>(side));

            TrkrHitSet* hitset = m_hits->findHitSet(hitset_key);

            if (!hitset)
            {
              if (Verbosity() > 2)
              {
                std::cout << Name()
                          << "::process_event - no hitset for"
                          << " layer=" << layer
                          << " region=" << region
                          << " sector=" << sector
                          << " side=" << side
                          << " key=" << hitset_key
                          << " (module key " << key << ")"
                          << std::endl;
              }

              continue;
            }

            InModuleThreadData td;

            td.layer = layer;
            td.region = region;
            td.sector = sector;
            td.side = static_cast<int>(side);
            td.hitsetkey = hitset_key;
            td.module_key = key;

            td.hitset = hitset;
            td.tGeometry = m_tGeometry;

            // Same default pedestal value used in many TPC raw-hit contexts.
            // Replace by your calibrated pedestal if available.
            td.pedestal = 74.4;

            td.verbosity = Verbosity();

            tdata.push_back(td);
        }
      }
    }
  }

 // if (Verbosity() > 0)
  //{
    std::cout << Name()
              << "::process_event - event " << m_event
              << " has " << tdata.size()
              << " non-empty  hitsets"
              << std::endl;
  //}

  if (tdata.empty())
  {
    if (m_tree)
    {
      m_tree->Fill();
    }

    ++m_event;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  // ------------------------------------------------------------
  // Launch pthreads in batches.
  //
  // This is important:
  // if m_maxThreads < tdata.size(), we must not process only the
  // first m_maxThreads modules and silently skip the rest.
  // ------------------------------------------------------------

  const unsigned int maxLive =
    std::max(1u, std::min(m_maxThreads,
                          static_cast<unsigned int>(tdata.size())));

  for (unsigned int start = 0;
       start < static_cast<unsigned int>(tdata.size());
       start += maxLive)
  {
    const unsigned int end =
      std::min(start + maxLive,
               static_cast<unsigned int>(tdata.size()));

    const unsigned int nLive = end - start;

    std::vector<pthread_t> threads;
    threads.resize(nLive);

    std::vector<int> thread_ok;
    thread_ok.resize(nLive, 0);

    for (unsigned int i = 0; i < nLive; ++i)
    {
      const unsigned int idx = start + i;

      const int rc =
        pthread_create(&threads[i],
                       0,
                       ProcessModule,
                       static_cast<void*>(&tdata[idx]));

      if (rc != 0)
      {
        std::cerr << Name()
                  << "::process_event - pthread_create failed for"
                  << " region=" << tdata[idx].region
                  << " sector=" << tdata[idx].sector
                  << " side=" << tdata[idx].side
                  << std::endl;

        thread_ok[i] = 0;
      }
      else
      {
        thread_ok[i] = 1;
      }
    }

    for (unsigned int i = 0; i < nLive; ++i)
    {
      if (thread_ok[i])
      {
        pthread_join(threads[i], 0);
      }
    }
  }

  // ------------------------------------------------------------
  // Serial merge.
  //
  // ROOT output is filled only here, after all pthreads are joined.
  // This avoids thread-unsafe ROOT writes.
  // ------------------------------------------------------------

  for (std::vector<InModuleThreadData>::const_iterator td = tdata.begin();
       td != tdata.end();
       ++td)
  {
    for (std::vector<InModuleThreadData::HitInfo>::const_iterator hi = td->hits.begin();
         hi != td->hits.end();
         ++hi)
    {
      m_tree_region.push_back(hi->region);
      m_tree_sector.push_back(hi->sector);
      m_tree_side.push_back(hi->side);

      m_tree_pad.push_back(hi->pad);
      m_tree_tbin.push_back(hi->tbin);
      m_tree_adc.push_back(hi->adc);

      if (Verbosity() > 2)
      {
        std::cout << "  event=" << m_event
                  << " region=" << hi->region
                  << " sector=" << hi->sector
                  << " side=" << hi->side
                  << " pad=" << hi->pad
                  << " tbin=" << hi->tbin
                  << " adc=" << hi->adc
                  << std::endl;
      }
    }
  }

  if (m_tree)
  {
    m_tree->Fill();
  }

  ++m_event;

  return Fun4AllReturnCodes::EVENT_OK;
}