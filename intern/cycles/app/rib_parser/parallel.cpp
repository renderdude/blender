#include <iterator>
#include <list>
#include <sstream>
#include <thread>
#include <vector>

#include "parallel.h"
#include "util/log.h"

CCL_NAMESPACE_BEGIN

// Barrier Method Definitions
bool Barrier::block()
{
  std::unique_lock<std::mutex> lock(_mutex);

  --_num_to_block;
  CHECK_GE(_num_to_block, 0);

  if (_num_to_block > 0) {
    _cv.wait(lock, [this]() { return _num_to_block == 0; });
  }
  else
    _cv.notify_all();

  return --_num_to_exit == 0;
}

Thread_Pool *Parallel_Job::thread_pool;

// Thread_Pool Method Definitions
Thread_Pool::Thread_Pool(int nThreads)
{
  for (int i = 0; i < nThreads - 1; ++i)
    threads.push_back(std::thread(&Thread_Pool::worker, this));
}

void Thread_Pool::worker()
{
  VLOG(1) << "Started execution in worker thread";

  std::unique_lock<std::mutex> lock(mutex);
  while (!shutdownThreads)
    work_or_wait(&lock, false);

  VLOG(1) << "Exiting worker thread";
}

std::unique_lock<std::mutex> Thread_Pool::add_to_job_list(Parallel_Job *job)
{
  std::unique_lock<std::mutex> lock(mutex);
  // Add _job_ to head of _jobList_
  if (jobList)
    jobList->prev = job;
  job->next = jobList;
  jobList = job;

  job_list_condition.notify_all();
  return lock;
}

void Thread_Pool::work_or_wait(std::unique_lock<std::mutex> *lock, bool is_enqueuing_thread)
{
  DCHECK(lock->owns_lock());
  // Return if this is a worker thread and the thread pool is disabled
  if (!is_enqueuing_thread && disabled) {
    job_list_condition.wait(*lock);
    return;
  }

  Parallel_Job *job = jobList;
  while (job && !job->have_work())
    job = job->next;
  if (job) {
    // Execute work for _job_
    job->active_workers++;
    job->run_step(lock);
    // Handle post-job-execution details
    DCHECK(!lock->owns_lock());
    lock->lock();
    job->active_workers--;
    if (job->finished())
      job_list_condition.notify_all();
  }
  else
    // Wait for new work to arrive or the job to finish
    job_list_condition.wait(*lock);
}

bool Thread_Pool::work_or_return()
{
  std::unique_lock<std::mutex> lock(mutex);

  Parallel_Job *job = jobList;
  while (job && !job->have_work())
    job = job->next;
  if (!job)
    return false;

  // Execute work for _job_
  job->active_workers++;
  job->run_step(&lock);
  DCHECK(!lock.owns_lock());
  lock.lock();
  job->active_workers--;
  if (job->finished())
    job_list_condition.notify_all();

  return true;
}

void Thread_Pool::remove_from_job_list(Parallel_Job *job)
{
  DCHECK(!job->removed);

  if (job->prev)
    job->prev->next = job->next;
  else {
    DCHECK(jobList == job);
    jobList = job->next;
  }
  if (job->next)
    job->next->prev = job->prev;

  job->removed = true;
}

void Thread_Pool::for_each_thread(std::function<void(void)> func)
{
  Barrier *barrier = new Barrier(threads.size() + 1);

  parallel_for(0, threads.size() + 1, [barrier, &func](int64_t) {
    func();
    if (barrier->block())
      delete barrier;
  });
}

void Thread_Pool::disable()
{
  CHECK(!disabled);
  disabled = true;
  CHECK(jobList == nullptr);  // Nothing should be running when disable() is called.
}

void Thread_Pool::reenable()
{
  CHECK(disabled);
  disabled = false;
}

Thread_Pool::~Thread_Pool()
{
  if (threads.empty())
    return;

  {
    std::lock_guard<std::mutex> lock(mutex);
    shutdownThreads = true;
    job_list_condition.notify_all();
  }

  for (std::thread &thread : threads)
    thread.join();
}

std::string Thread_Pool::to_string() const
{
  std::stringstream ss;
  ss << "[ Thread_Pool threads.size(): " << threads.size();
  ss << " shutdownThreads: " << shutdownThreads;
  std::string s = ss.str();
  if (mutex.try_lock()) {
    s += "jobList: [ ";
    Parallel_Job *job = jobList;
    while (job) {
      s += job->to_string() + " ";
      job = job->next;
    }
    s += "] ";
    mutex.unlock();
  }
  else
    s += "(job list mutex locked) ";
  return s + "]";
}

bool do_parallel_work()
{
  CHECK(Parallel_Job::thread_pool);
  // lock should be held when this is called...
  return Parallel_Job::thread_pool->work_or_return();
}

// Parallel_For_Loop_1D Definition
class Parallel_For_Loop_1D : public Parallel_Job {
 public:
  // Parallel_For_Loop_1D Public Methods
  Parallel_For_Loop_1D(int64_t start_index,
                       int64_t end_index,
                       int chunk_size,
                       std::function<void(int64_t, int64_t)> func)
      : _func(std::move(func)),
        _next_index(start_index),
        _end_index(end_index),
        _chunk_size(chunk_size)
  {
  }

  bool have_work() const
  {
    return _next_index < _end_index;
  }

  void run_step(std::unique_lock<std::mutex> *lock);

  std::string to_string() const
  {
    std::stringstream ss;
    ss << "[ Parallel_For_Loop_1D next_index: " << _next_index;
    ss << " end_index: " << _end_index << " chunk_size: " << _chunk_size << " ]";

    return ss.str();
  }

 private:
  // Parallel_For_Loop_1D Private Members
  std::function<void(int64_t, int64_t)> _func;
  int64_t _next_index, _end_index;
  int _chunk_size;
};

// Parallel_For_Loop_1D Method Definitions
void Parallel_For_Loop_1D::run_step(std::unique_lock<std::mutex> *lock)
{
  // Determine the range of loop iterations to run in this step
  int64_t index_start = _next_index;
  int64_t indexEnd = std::min(index_start + _chunk_size, _end_index);
  _next_index = indexEnd;

  // Remove job from list if all work has been started
  if (!have_work())
    thread_pool->remove_from_job_list(this);

  // Release lock and execute loop iterations in _[index_start, indexEnd)_
  lock->unlock();
  _func(index_start, indexEnd);
}

// Parallel Function Definitions
void parallel_for(int64_t start, int64_t end, std::function<void(int64_t, int64_t)> func)
{
  CHECK(Parallel_Job::thread_pool);
  if (start == end)
    return;

  // Compute chunk size for parallel loop
  int64_t chunk_size = std::max<int64_t>(1, (end - start) / (8 * running_threads()));

  // Create and enqueue _Parallel_For_Loop_1D_ for this loop
  Parallel_For_Loop_1D loop(start, end, chunk_size, std::move(func));
  std::unique_lock<std::mutex> lock = Parallel_Job::thread_pool->add_to_job_list(&loop);

  // Help out with parallel loop iterations in the current thread
  while (!loop.finished())
    Parallel_Job::thread_pool->work_or_wait(&lock, true);
}

///////////////////////////////////////////////////////////////////////////

int available_cores()
{
  return std::max<int>(1, std::thread::hardware_concurrency());
}

int running_threads()
{
  return Parallel_Job::thread_pool ? (1 + Parallel_Job::thread_pool->size()) : 1;
}

void parallel_init(int nThreads)
{
  CHECK(!Parallel_Job::thread_pool);
  if (nThreads <= 0)
    nThreads = available_cores();
  Parallel_Job::thread_pool = new Thread_Pool(nThreads);
}

void parallel_cleanup()
{
  delete Parallel_Job::thread_pool;
  Parallel_Job::thread_pool = nullptr;
}

void for_each_thread(std::function<void(void)> func)
{
  if (Parallel_Job::thread_pool)
    Parallel_Job::thread_pool->for_each_thread(std::move(func));
}

void disable_thread_pool()
{
  CHECK(Parallel_Job::thread_pool);
  Parallel_Job::thread_pool->disable();
}

void reenable_thread_pool()
{
  CHECK(Parallel_Job::thread_pool);
  Parallel_Job::thread_pool->reenable();
}

CCL_NAMESPACE_END
