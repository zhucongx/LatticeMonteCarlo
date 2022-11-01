#ifndef LKMC_LKMC_PRED_INCLUDE_VACANCYMIGRATIONPREDICTORQUARTICLRU_H_
#define LKMC_LKMC_PRED_INCLUDE_VACANCYMIGRATIONPREDICTORQUARTICLRU_H_
#include "VacancyMigrationPredictorQuartic.h"
#include <thread>
#include <mutex>
#include <shared_mutex>
namespace pred {

template<class T>
class LruCache {
  public:
    LruCache(size_t capacity) : capacity_(capacity) {}
    bool Exist(size_t key) {
      std::scoped_lock<std::shared_mutex> lock(mu_);
      return (cache_hashmap_.find(key) != cache_hashmap_.end());
    }
    void Add(size_t key, const T &value) {
      std::scoped_lock<std::shared_mutex> lock(mu_);
      auto it = cache_hashmap_.find(key);
      if (it != cache_hashmap_.end()) {
        cache_list_.erase(it->second);
      }
      cache_list_.push_front(std::make_pair(key, std::move(value)));
      cache_hashmap_[key] = cache_list_.begin();
      if (cache_hashmap_.size() > capacity_) {
        auto last = cache_list_.rbegin()->first;
        cache_list_.pop_back();
        cache_hashmap_.erase(last);
      }
    }
    T Get(size_t key) {
      std::scoped_lock<std::shared_mutex> lock(mu_);
      auto it = cache_hashmap_.find(key);
      if (it == cache_hashmap_.end()) {
        throw std::runtime_error("LRUCache::Get() : key not found");
      }
      cache_list_.splice(cache_list_.begin(), cache_list_, it->second);
      return it->second->second;
    }
  private:
    size_t capacity_;
    // cache list of key and value
    std::list<std::pair<size_t, T>> cache_list_{};
    // hashmap of key and iterator of cache list
    std::unordered_map<size_t, decltype(cache_list_.begin()) > cache_hashmap_{};
    std::shared_mutex mu_;
};

class VacancyMigrationPredictorQuarticLru : public VacancyMigrationPredictorQuartic {
  public:
    VacancyMigrationPredictorQuarticLru(const std::string &predictor_filename,
                                        const cfg::Config &reference_config,
                                        const std::set<Element> &element_set,
                                        size_t cache_size);
    ~VacancyMigrationPredictorQuarticLru() override;
  private:
    [[nodiscard]] std::pair<double, double> GetBarrierAndDiffFromLatticeIdPair(
        const cfg::Config &config,
        const std::pair<size_t, size_t> &lattice_id_jump_pair) const override;
    [[nodiscard]] size_t GetHashFromConfigAndLatticeIdPair(
        const cfg::Config &config,
        const std::pair<size_t, size_t> &lattice_id_jump_pair) const;
    mutable LruCache<std::pair<double, double>> lru_cache_;
};
} // namespace pred

#endif //LKMC_LKMC_PRED_INCLUDE_VACANCYMIGRATIONPREDICTORQUARTICLRU_H_
