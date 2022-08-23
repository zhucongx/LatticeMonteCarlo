#include "Iterator.h"

#include <utility>
namespace ansys {

Iterator::Iterator(unsigned long long int initial_steps,
                   unsigned long long int increment_steps,
                   Element solvent_element,
                   std::set<Element> element_set,
                   size_t smallest_cluster_criteria,
                   // size_t solvent_bond_criteria,
                   const std::string &predictor_filename)
    : initial_steps_(initial_steps),
      increment_steps_(increment_steps),
      final_number_(increment_steps),
      solvent_element_(solvent_element),
      smallest_cluster_criteria_(smallest_cluster_criteria),
      // solvent_bond_criteria_(solvent_bond_criteria),
      energy_estimator_(predictor_filename, std::move(element_set)) {
  std::ifstream ifs("kmc_log_modified.txt", std::ifstream::in);
  if (!ifs.is_open()) {
    std::cerr << "Cannot open kmc_log_modified.txt\n";
    return;
  }
  unsigned long long filename;
  double time;
  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  while (true) {
    if (ifs.eof() || ifs.bad()) {
      break;
    }
    if (ifs.peek() == '#') {
      ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      continue;
    }

    ifs >> filename >> time;
    ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    if (filename >= initial_steps_ && (filename - initial_steps_) % increment_steps == 0) {
      filename_time_hashset_[filename] = time;
      final_number_ = filename;
    }
  }
  final_number_ -= increment_steps;
}
Iterator::~Iterator() = default;
void Iterator::RunCluster() const {
  // start
  std::ofstream ofs("clusters_info.json", std::ofstream::out);
  ofs << "[ \n";

  for (unsigned long long i = 0; i <= final_number_; i += increment_steps_) {
    std::cout << i << " / " << final_number_ << std::endl;

    ClustersFinder cluster_finder(cfg::Config::ReadMap(
                                      "lattice.txt", "element.txt",
                                      "map" + std::to_string(i) + ".txt"),
                                  solvent_element_,
                                  smallest_cluster_criteria_,
                                  // solvent_bond_criteria_,
                                  energy_estimator_);
    auto num_different_element =
        cluster_finder.FindClustersAndOutput("cluster", std::to_string(i) + "_cluster.cfg");

    ofs << "{ \n"
        << "\"index\" : "
        << "\"" << std::to_string(i) << "\",\n"
        << "\"time\" : " << filename_time_hashset_.at(i) << ",\n"
        << "\"clusters\" : [ \n";
    for (auto it = num_different_element.cbegin(); it < num_different_element.cend(); ++it) {
      ofs << "[ ";
      const auto &cluster = *it;
      std::for_each(cluster.first.cbegin(), cluster.first.cend(), [&ofs](auto ii) {
        ofs << ii.second << ", ";
      });
      ofs << std::setprecision(16) << cluster.second;
      if (it == num_different_element.cend() - 1) {
        ofs << "] \n";
      } else {
        ofs << "], \n";
      }
    }
    ofs << "]\n";
    if (i != final_number_) {
      ofs << "}, \n";
    } else {
      ofs << "} \n";
    }
  }
  ofs << " ]" << std::endl;
}
void Iterator::RunReformat() const {
  for (unsigned long long i = 0; i <= final_number_; i += increment_steps_) {
    std::cout << i << " / " << final_number_ << std::endl;
    auto config = cfg::Config::ReadCfg(std::to_string(i) + ".cfg");
    config.ReassignLatticeVector();
    if (i == 0) {
      config.WriteLattice("lattice.txt");
      config.WriteElement("element.txt");
    }
    config.WriteMap("map" + std::to_string(i) + ".txt");
  }
}
} // namespace ansys