#ifndef DIPLOID_HPP_
#define DIPLOID_HPP_

#include <vector>
#include <memory>
#include <utility>

#include "biosoup/sequence.hpp"

namespace thread_pool {
class ThreadPool;
}  // namespace thread_pool

namespace raven {
  
using DiploidSequences =
  ::std::pair<
    ::std::vector<::std::unique_ptr<::biosoup::Sequence>>,
    ::std::vector<::std::unique_ptr<::biosoup::Sequence>>
  >;

DiploidSequences Partition(
  ::std::vector<::std::unique_ptr<::biosoup::Sequence>> const& sequences,
  ::std::shared_ptr<::thread_pool::ThreadPool> thread_pool);

}  // namespace raven

#endif  // DIPLOID_HPP_