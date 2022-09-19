#include <aocommon/barrier.h>
#include <aocommon/parallelfor.h>

#include <iostream>

#include <cmath>
#include <mutex>

#include <unistd.h>  // for sleep

#include <boost/test/unit_test.hpp>

using namespace aocommon;

BOOST_AUTO_TEST_SUITE(parallelfor)

void func() {}

BOOST_AUTO_TEST_CASE(barrier) {
  BOOST_CHECK_NO_THROW(Barrier b(1, func); b.wait();

                       Barrier c(2, func); std::thread t([&]() { c.wait(); });
                       c.wait(); t.join(););
}

BOOST_AUTO_TEST_CASE(unused) { BOOST_CHECK_NO_THROW(ParallelFor<size_t>(4);); }

BOOST_AUTO_TEST_CASE(single_thread) {
  ParallelFor<size_t> loop(1);
  std::vector<size_t> counts(10, 0);
  loop.Run(0, 10, [&](size_t iter) { counts[iter]++; });

  const std::vector<size_t> ref(10, 1);
  BOOST_CHECK_EQUAL_COLLECTIONS(counts.begin(), counts.end(), ref.begin(),
                                ref.end());
}

BOOST_AUTO_TEST_CASE(parallel_with_thread_id) {
  ParallelFor<size_t> loop(4);
  std::mutex mutex;
  std::vector<size_t> counts(10, 0);
  loop.Run(0, 10, [&](size_t iter, size_t thread_id) {
    std::unique_lock<std::mutex> lock(mutex);
    BOOST_CHECK_LT(thread_id, 4);
    counts[iter]++;
  });

  const std::vector<size_t> ref(10, 1);
  BOOST_CHECK_EQUAL_COLLECTIONS(counts.begin(), counts.end(), ref.begin(),
                                ref.end());
}

BOOST_AUTO_TEST_CASE(parallel_without_thread_id) {
  ParallelFor<size_t> loop(4);
  std::mutex mutex;
  std::vector<size_t> counts(10, 0);
  loop.Run(0, 10, [&](size_t iter) {
    std::unique_lock<std::mutex> lock(mutex);
    counts[iter]++;
  });

  const std::vector<size_t> ref(10, 1);
  BOOST_CHECK_EQUAL_COLLECTIONS(counts.begin(), counts.end(), ref.begin(),
                                ref.end());
}

BOOST_AUTO_TEST_CASE(reuse_with_thread_id) {
  std::vector<size_t> counts(20, 0);
  std::mutex mutex;
  ParallelFor<size_t> loop(4);
  loop.Run(0, 10, [&](size_t iter, size_t thread_id) {
    std::unique_lock<std::mutex> lock(mutex);
    BOOST_CHECK_LT(thread_id, 4);
    counts[iter]++;
  });
  std::vector<size_t> ref(20, 0);
  std::fill(ref.begin(), ref.begin() + 10, 1);
  BOOST_CHECK_EQUAL_COLLECTIONS(counts.begin(), counts.end(), ref.begin(),
                                ref.end());

  loop.Run(10, 20, [&](size_t iter, size_t thread_id) {
    std::unique_lock<std::mutex> lock(mutex);
    BOOST_CHECK_LT(thread_id, 4);
    counts[iter]++;
  });
  ref = std::vector<size_t>(20, 1);
  BOOST_CHECK_EQUAL_COLLECTIONS(counts.begin(), counts.end(), ref.begin(),
                                ref.end());
}

BOOST_AUTO_TEST_CASE(reuse_without_thread_id) {
  std::vector<size_t> counts(20, 0);
  std::mutex mutex;
  ParallelFor<size_t> loop(4);
  loop.Run(0, 10, [&](size_t iter) {
    std::unique_lock<std::mutex> lock(mutex);
    counts[iter]++;
  });
  std::vector<size_t> ref(20, 0);
  std::fill(ref.begin(), ref.begin() + 10, 1);
  BOOST_CHECK_EQUAL_COLLECTIONS(counts.begin(), counts.end(), ref.begin(),
                                ref.end());

  loop.Run(10, 20, [&](size_t iter) {
    std::unique_lock<std::mutex> lock(mutex);
    counts[iter]++;
  });
  ref = std::vector<size_t>(20, 1);
  BOOST_CHECK_EQUAL_COLLECTIONS(counts.begin(), counts.end(), ref.begin(),
                                ref.end());
}

BOOST_AUTO_TEST_CASE(single_threaded_with_thread_id) {
  ParallelFor<size_t> pFor(1);
  std::vector<size_t> counts(10, 0);
  pFor.Run(0, 10, [&](size_t iter, size_t) { counts[iter]++; });

  std::vector<size_t> ref(10, 1);
  BOOST_CHECK_EQUAL_COLLECTIONS(counts.begin(), counts.end(), ref.begin(),
                                ref.end());
}

BOOST_AUTO_TEST_SUITE_END()
