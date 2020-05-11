#include "diploid.hpp"

#include "ram/minimizer_engine.hpp"
#include "biosoup/overlap.hpp"

#include <iostream>
#include <future>

namespace raven {

DiploidSequences Partition(
  ::std::vector<::std::unique_ptr<::biosoup::Sequence>> const& sequences,
  ::std::shared_ptr<::thread_pool::ThreadPool> thread_pool) {
  
  auto const k = 15;
  auto const w = 5;
  auto const f = 0.001;

  ::ram::MinimizerEngine minimizer_engine{k, w, thread_pool};
  
  minimizer_engine.Minimize(sequences.begin(), sequences.end());
  minimizer_engine.Filter(f);
  
  using Overlaps = ::std::vector<::biosoup::Overlap>;
  
  ::std::vector<::std::future<Overlaps>> overlap_futures;
  
  for (::std::uint_fast32_t i = 0; i < sequences.size(); ++i)
    overlap_futures.emplace_back(thread_pool->Submit(
      [&, i] {
        return minimizer_engine.Map(sequences[i], true, true, false);
      }
    ));
  
  for (auto&& f : overlap_futures)
    for (auto&& overlap : f.get())
      ::std::cout << overlap.lhs_id << " " << overlap.rhs_id << "\n";
  
  ::std::cerr << "[ram::diploid::Partition] unflattened overlaps size "
              << overlap_futures.size()
              << "\n";
  
  return {};
}

}  // namespace raven