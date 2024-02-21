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
