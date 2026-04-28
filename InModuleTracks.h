#pragma once

#include <fun4all/SubsysReco.h>

#include <trackbase/TrkrDefs.h>

#include <string>
#include <vector>

// forward declarations
class PHCompositeNode;
class TrkrHitSetContainer;
class TrkrHitSet;
class ActsGeometry;
class TFile;
class TTree;

// ── per-module work package ─────────────────────────────────────────
//
// One object corresponds to one TPC module key:
//
//   region/module = 0,1,2
//   sector        = 0..11
//   side          = 0,1
//
// All input fields are filled before pthread_create.
// Output vector is filled only by its own thread.
// Main thread reads output only after pthread_join.
//
struct InModuleThreadData
{
  InModuleThreadData();

  // module identity
  unsigned int layer;
  unsigned int region;
  unsigned int sector;
  int side;
  TrkrDefs::hitsetkey hitsetkey;
  TrkrDefs::hitsetkey module_key;

  // input
  TrkrHitSet* hitset;
  ActsGeometry* tGeometry;

  double pedestal;
  int verbosity;

  // thread-private output
  struct HitInfo
  {
    HitInfo();
    unsigned int layer;
    unsigned int region;
    unsigned int sector;
    int side;

    unsigned short pad;
    unsigned short tbin;
    unsigned short adc;
  };

  std::vector<HitInfo> hits;
};

// ── module ──────────────────────────────────────────────────────────

class InModuleTracks : public SubsysReco
{
 public:
  explicit InModuleTracks(const std::string& name = "InModuleTracks",
                          const std::string& filename = "InModuleTracks.root");

  virtual ~InModuleTracks();

  int Init(PHCompositeNode*);
  int InitRun(PHCompositeNode*);
  int process_event(PHCompositeNode*);
  int End(PHCompositeNode*);

  // cap number of simultaneously live pthreads
  void setMaxThreads(unsigned int n);

 private:
  int getNodes(PHCompositeNode*);
  void reset_tree_vars();

  // output
  std::string m_outputFileName;
  TFile* m_outputFile;
  TTree* m_tree;

  // node-tree pointers
  TrkrHitSetContainer* m_hits;
  ActsGeometry* m_tGeometry;

  // event bookkeeping
  int m_event;
  unsigned int m_maxThreads;

  // tree variables, filled once per event
  int m_tree_event;

  std::vector<unsigned int> m_tree_region;
  std::vector<unsigned int> m_tree_sector;
  std::vector<int> m_tree_side;

  std::vector<unsigned short> m_tree_pad;
  std::vector<unsigned short> m_tree_tbin;
  std::vector<unsigned short> m_tree_adc;
};