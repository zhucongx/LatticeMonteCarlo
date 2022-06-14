#ifndef LKMC_LKMC_ANSYS_INCLUDE_ITERATOR_H_
#define LKMC_LKMC_ANSYS_INCLUDE_ITERATOR_H_
#include "EnergyEstimator.h"
#include "ClustersFinder.h"
namespace ansys {

class Iterator {
  public:
    virtual void IterateToRun() = 0;
    Iterator(unsigned long long int initial_number,
                unsigned long long int increment_number,
                unsigned long long int finial_number);
    virtual ~Iterator();

  protected:
    const unsigned long long initial_number_;
    const unsigned long long increment_number_;
    const unsigned long long finial_number_;

};

} // namespace ansys

#endif //LKMC_LKMC_ANSYS_INCLUDE_ITERATOR_H_
