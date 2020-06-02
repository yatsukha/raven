#ifndef DIPLOID_HPP_
#define DIPLOID_HPP_

#include <memory>
#include <utility>
#include <vector>

#include "biosoup/sequence.hpp"

namespace thread_pool {
class ThreadPool;
}  // namespace thread_pool

namespace raven {

namespace diploid {

using DiploidSequences =
    ::std::pair<::std::vector<::std::unique_ptr<::biosoup::Sequence>>,
                ::std::vector<::std::unique_ptr<::biosoup::Sequence>>>;

DiploidSequences Partition(
    ::std::vector<::std::unique_ptr<::biosoup::Sequence>> const& sequences,
    ::std::shared_ptr<::thread_pool::ThreadPool> thread_pool,
    ::std::int8_t const m, ::std::int8_t const n, ::std::int8_t const g);

}  // namespace diploid

}  // namespace raven

#endif  // DIPLOID_HPP_