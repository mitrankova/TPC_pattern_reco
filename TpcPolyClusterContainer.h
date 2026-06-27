#pragma once

#include <phool/PHObject.h>

#include <iostream>

class TpcPolyCluster;

class TpcPolyClusterContainer : public PHObject
{
 public:
  TpcPolyClusterContainer() = default;
  ~TpcPolyClusterContainer() override = default;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;

  virtual unsigned int size() const { return 0; }
  virtual void add_cluster(TpcPolyCluster*) {}
  virtual const TpcPolyCluster* get_cluster(unsigned int) const { return nullptr; }
  virtual TpcPolyCluster* get_cluster(unsigned int) { return nullptr; }

 private:
  ClassDefOverride(TpcPolyClusterContainer, 0)
};
