#pragma once

#include <boost/optional/optional.hpp>
#include <condition_variable>
#include <cstdint>
#include <functional>
#include <glog/logging.h>
#include <mutex>
#include <shared_mutex>
#include <sstream>
#include <string>
#include <thread>
#include <type_traits>
#include <vector>

CCL_NAMESPACE_BEGIN

// Parallel Function Declarations
void parallel_init(int nThreads = -1);
void parallel_cleanup();

int available_cores();
int running_threads();

template<typename T> class Thread_Local {
 public:
  Thread_Local() : _hash_table(4 * running_threads()), _create([]() { return T(); }) {}
  Thread_Local(std::function<T(void)> &&c) : _hash_table(4 * running_threads()), _create(c) {}

  T &get();

  template<typename F> void for_all(F &&func);

 private:
  // Thread_Local Private Members
  struct Entry {
    std::thread::id tid;
    T value;
  };
  std::shared_mutex _mutex;
  std::vector<boost::optional<Entry>> _hash_table;
  std::function<T(void)> _create;
};

// Thread_Local Inline Methods
template<typename T> inline T &Thread_Local<T>::get()
{
  std::thread::id tid = std::this_thread::get_id();
  uint32_t hash = std::hash<std::thread::id>()(tid);
  hash %= _hash_table.size();
  int step = 1;
  int tries = 0;

  _mutex.lock_shared();
  while (true) {
    CHECK_LT(++tries, _hash_table.size());  // full hash table

    if (_hash_table[hash] && _hash_table[hash]->tid == tid) {
      // Found it
      T &threadLocal = _hash_table[hash]->value;
      _mutex.unlock_shared();
      return threadLocal;
    }
    if (!_hash_table[hash]) {
      _mutex.unlock_shared();

      // Get reader-writer lock before calling the callback so that the user
      // doesn't have to worry about writing a thread-safe callback.
      _mutex.lock();
      T newItem = _create();

      if (_hash_table[hash]) {
        // someone else got there first--keep looking, but now
        // with a writer lock.
        while (true) {
          hash += step;
          ++step;
          if (hash >= _hash_table.size()) {
            hash %= _hash_table.size();
          }

          if (!_hash_table[hash]) {
            break;
          }
        }
      }

      _hash_table[hash] = Entry{tid, std::move(newItem)};
      T &threadLocal = _hash_table[hash]->value;
      _mutex.unlock();
      return threadLocal;
    }

    hash += step;
    ++step;
    if (hash >= _hash_table.size()) {
      hash %= _hash_table.size();
    }
  }
}

template<typename T> template<typename F> inline void Thread_Local<T>::for_all(F &&func)
{
  _mutex.lock();
  for (auto &entry : _hash_table) {
    if (entry) {
      func(entry->value);
    }
  }
  _mutex.unlock();
}

// Barrier Definition
class Barrier {
 public:
  explicit Barrier(int n) : _num_to_block(n), _num_to_exit(n) {}

  Barrier(const Barrier &) = delete;
  Barrier &operator=(const Barrier &) = delete;

  // All block. Returns true to only one thread (which should delete the
  // barrier).
  bool block();

 private:
  std::mutex _mutex;
  std::condition_variable _cv;
  int _num_to_block, _num_to_exit;
};

void Parallel_For(int64_t start, int64_t end, std::function<void(int64_t, int64_t)> func);

// Parallel Inline Functions
inline void Parallel_For(int64_t start, int64_t end, std::function<void(int64_t)> func)
{
  Parallel_For(start, end, [&func](int64_t start, int64_t end) {
    for (int64_t i = start; i < end; ++i) {
      func(i);
    }
  });
}

class Thread_Pool;

class Parallel_Job {
 public:
  // Parallel_Job Public Methods
  virtual ~Parallel_Job()
  {
    DCHECK(removed);
  }

  virtual bool have_work() const = 0;
  virtual void run_step(std::unique_lock<std::mutex> *lock) = 0;

  bool finished() const
  {
    return !have_work() && active_workers == 0;
  }

  virtual std::string to_string() const = 0;

  // Parallel_Job Public Members
  static Thread_Pool *thread_pool;

 protected:
  std::string base_to_string() const
  {
    std::stringstream ss;
    ss << "active_workers: " << active_workers;
    ss << " removed: " << removed;
    return ss.str();
  }

 private:
  // Parallel_Job Private Members
  friend class Thread_Pool;
  int active_workers = 0;
  Parallel_Job *prev = nullptr, *next = nullptr;
  bool removed = false;
};

// Thread_Pool Definition
class Thread_Pool {
 public:
  // Thread_Pool Public Methods
  explicit Thread_Pool(int nThreads);

  ~Thread_Pool();

  size_t size() const
  {
    return threads.size();
  }

  std::unique_lock<std::mutex> add_to_job_list(Parallel_Job *job);
  void remove_from_job_list(Parallel_Job *job);

  void work_or_wait(std::unique_lock<std::mutex> *lock, bool is_enqueuing_thread);
  bool work_or_return();

  void disable();
  void reenable();

  void for_each_thread(std::function<void(void)> func);

  std::string to_string() const;

 private:
  // Thread_Pool Private Methods
  void worker();

  // Thread_Pool Private Members
  std::vector<std::thread> threads;
  mutable std::mutex mutex;
  bool shutdownThreads = false;
  bool disabled = false;
  Parallel_Job *jobList = nullptr;
  std::condition_variable job_list_condition;
};

bool do_parallel_work();

template<typename T> class Async_Job : public Parallel_Job {
 public:
  Async_Job(std::function<T(void)> w) : func(std::move(w)) {}

  bool have_work() const override
  {
    return !started;
  }

  void run_step(std::unique_lock<std::mutex> *lock) override
  {
    thread_pool->remove_from_job_list(this);
    started = true;
    lock->unlock();

    // Execute asynchronous work and notify waiting threads of its completion
    T r = func();
    std::unique_lock<std::mutex> ul(_mutex);
    _result = r;
    _cv.notify_all();
  }

  bool is_ready() const
  {
    std::lock_guard<std::mutex> lock(_mutex);
    return _result.has_value();
  }

  T get_result()
  {
    wait();
    std::lock_guard<std::mutex> lock(_mutex);
    return *_result;
  }

  boost::optional<T> try_get_result(std::mutex *extMutex)
  {
    {
      std::lock_guard<std::mutex> lock(_mutex);
      if (_result) {
        return _result;
      }
    }

    extMutex->unlock();
    do_parallel_work();
    extMutex->lock();
    return {};
  }

  void wait()
  {
    while (!is_ready() && do_parallel_work()) {
    }
    std::unique_lock<std::mutex> lock(_mutex);
    if (!_result.has_value()) {
      _cv.wait(lock, [this]() { return _result.has_value(); });
    }
  }

  void do_work()
  {
    T r = func();
    std::unique_lock<std::mutex> l(_mutex);
    CHECK(!_result.has_value());
    _result = r;
    _cv.notify_all();
  }

  std::string to_string() const override
  {
    return string_printf("[ Async_Job started: %d ]", started);
  }

 private:
  std::function<T(void)> func;
  bool started = false;
  boost::optional<T> _result;
  mutable std::mutex _mutex;
  std::condition_variable _cv;
};

void for_each_thread(std::function<void(void)> func);

void disable_thread_pool();
void reenable_thread_pool();

template<typename F, typename... Args> inline auto run_async(F func, Args &&...args)
{
  // Create _AsyncJob_ for _func_ and _args_
  auto fvoid = std::bind(func, std::forward<Args>(args)...);
  using R = typename std::invoke_result_t<F, Args...>;
  Async_Job<R> *job = new Async_Job<R>(std::move(fvoid));

  // Enqueue _job_ or run it immediately
  std::unique_lock<std::mutex> lock;
  if (running_threads() == 1) {
    job->do_work();
  }
  else {
    lock = Parallel_Job::thread_pool->add_to_job_list(job);
  }
  return job;
}

CCL_NAMESPACE_END
