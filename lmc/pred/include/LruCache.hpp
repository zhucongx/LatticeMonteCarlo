#ifndef LMC_LMC_PRED_INCLUDE_LRUCACHE_HPP
#define LMC_LMC_PRED_INCLUDE_LRUCACHE_HPP
#include <thread>
#include <mutex>
#include <shared_mutex>
#include <list>
#include <unordered_map>
template<typename K, class V>
class LruCache {
public:
  explicit LruCache(size_t capacity) : capacity_(capacity) {}
  void Add(const K key, const V &value) {
    std::unique_lock<std::shared_mutex> lock(mu_);
    auto it = cache_hashmap_.find(key);
    if (it != cache_hashmap_.end()) {
      cache_list_.erase(it->second);
    }
    cache_list_.push_front(std::make_pair(key, value));
    cache_hashmap_[key] = cache_list_.begin();
    if (cache_hashmap_.size() > capacity_) {
      auto last = cache_list_.rbegin()->first;
      cache_list_.pop_back();
      cache_hashmap_.erase(last);
    }
  }
  // Use shared locking for reads to reduce contention under OpenMP.
  // Option A: shared_lock + no reordering on hit (or defer reordering).
  // Option B: thread-local front cache and merge misses to the global LRU outside
  // parallel regions to avoid lock contention entirely. Flat hash maps can also
  // lower constant factors if K is small and dense.
  // Note: Profiling shows LRU tweaks (shared lock / FIFO) only gave ~2.5% gain,
  // so we keep the simple write-lock path unless cache contention grows.
  [[nodiscard]] bool Get(const K key, V &value) {
    std::unique_lock<std::shared_mutex> lock(mu_);
    auto it = cache_hashmap_.find(key);
    if (it == cache_hashmap_.end()) {
      return false;
    }
    cache_list_.splice(cache_list_.begin(), cache_list_, it->second);
    value = it->second->second;
    return true;
  }
private:
  size_t capacity_;
  // cache list of key and value
  std::list<std::pair<K, V> > cache_list_{};
  // hashmap of key and iterator of cache list
  std::unordered_map<K, decltype(cache_list_.begin())> cache_hashmap_{};
  std::shared_mutex mu_{};
};
#endif //LMC_LMC_PRED_INCLUDE_LRUCACHE_HPP
