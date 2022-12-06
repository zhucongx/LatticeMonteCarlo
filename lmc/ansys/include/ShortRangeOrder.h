#ifndef LMC_LMC_ANSYS_INCLUDE_SHORTRANGEORDER_H_
#define LMC_LMC_ANSYS_INCLUDE_SHORTRANGEORDER_H_

#include "Config.h"

namespace ansys {

class ShortRangeOrder {
  public:
    explicit ShortRangeOrder(cfg::Config config, std::set<Element>  element_set);
    [[nodiscard]] std::map<std::string, double> FindWarrenCowley(size_t shell_number) const;
  protected:
    cfg::Config config_;
    const std::set<Element> element_set_;
};

} // ansys

#endif //LMC_LMC_ANSYS_INCLUDE_SHORTRANGEORDER_H_
