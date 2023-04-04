#ifndef INTERN_CACHE_H
#define INTERN_CACHE_H

#include <functional>
#include <shared_mutex>

CCL_NAMESPACE_BEGIN
/** @brief Intern_Cache.
 * @details
 */

template<typename T, typename Hash = std::hash<T>> class Intern_Cache {
 public:
  /// @name Initialization
  ///@{
  Intern_Cache() : _hash_table(256)
  {
  }

  ///@}

  /// @name Access
  ///@{
  template<typename F> const T *lookup(const T &item, F create)
  {
    size_t offset = Hash()(item) % _hash_table.size();
    int step = 1;
    _mutex.lock_shared();
    while (true) {
      // Check _hashTable[offset]_ for provided item
      if (!_hash_table[offset]) {
        // Insert item into open hash table entry
        _mutex.unlock_shared();
        _mutex.lock();
        // Double check that another thread hasn't inserted _item_
        size_t offset = Hash()(item) % _hash_table.size();
        int step = 1;
        while (true) {
          if (!_hash_table[offset])
            // fine--it's definitely not there
            break;
          else if (*_hash_table[offset] == item) {
            // Another thread inserted it
            const T *ret = _hash_table[offset];
            _mutex.unlock();
            return ret;
          }
          else {
            // collision
            offset += step;
            ++step;
            offset %= _hash_table.size();
          }
        }

        // Grow the hash table if needed
        if (4 * _num_entries > _hash_table.size()) {
          std::vector<const T *> newHash(2 * _hash_table.size());
          for (const T *ptr : _hash_table)
            if (ptr)
              Insert(ptr, &newHash);

          _hash_table.swap(newHash);
        }

        // Allocate new hash table entry and add it to the hash table
        ++_num_entries;
        T *newPtr = create(item);
        Insert(newPtr, &_hash_table);
        _mutex.unlock();
        return newPtr;
      }
      else if (*_hash_table[offset] == item) {
        // Return pointer for found _item_ in hash table
        const T *ret = _hash_table[offset];
        _mutex.unlock_shared();
        return ret;
      }
      else {
        // Advance _offset_ after hash table collision
        offset += step;
        ++step;
        offset %= _hash_table.size();
      }
    }
  }

  const T *lookup(const T &item)
  {
    return lookup(item, [](const T &item) { return new T(item); });
  }

  ///@}
  /// @name Measurement
  ///@{
  size_t size() const
  {
    return _num_entries;
  }
  size_t capacity() const
  {
    return _hash_table.size();
  }

  ///@}

 private:
  void Insert(const T *ptr, std::vector<const T *> *table)
  {
    size_t offset = Hash()(*ptr) % table->size();
    int step = 1;
    // Advance _offset_ to next free entry in hash table
    while ((*table)[offset]) {
      offset += step;
      ++step;
      offset %= table->size();
    }

    (*table)[offset] = ptr;
  }

  // InternCache Private Members
  size_t _num_entries = 0;
  std::vector<const T *> _hash_table;
  std::shared_mutex _mutex;

};  // end of class Intern_Cache

CCL_NAMESPACE_END

#endif  // INTERN_CACHE_H
