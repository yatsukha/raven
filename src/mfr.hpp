#pragma once

#include "biosoup/progress_bar.hpp"

#include <iostream>
#include <algorithm>
#include <future>
#include <string>
#include <tuple>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <memory>
#include <random>

namespace raven {

using Node       = ::std::uint_fast32_t;
using Neighbours = ::std::unordered_set<Node>;
using Removed    = ::std::unordered_set<Node>;
using Graph      = ::std::unordered_map<Node, Neighbours>;
using Cycle      = ::std::vector<Node>;

class NullOpt {};

template<typename T>
class Optional {
  ::std::shared_ptr<T> data_;
  
 public:
  Optional() = default;
  Optional(NullOpt&&) {};
  
  Optional(T const& data): data_(new (new char[sizeof(T)]) T{data}) {}
  Optional(T&& data): data_(new (new char[sizeof(T)]) T{::std::forward<T>(data)}) {}
  
  Optional(Optional const&) = default;
  Optional& operator=(Optional const&) = default;
  
  Optional(Optional&&) = default;
  Optional& operator=(Optional&&) = default;
  
  operator bool() const noexcept {
    return data_.operator bool();
  }
  
  T const& operator*() const noexcept {
    return *data_;  
  }
};

using VisitedDephts =
  ::std::unordered_map<
    Node,
    ::std::pair<Node, Optional<Node>>>;
using OptionalCycle = Optional<Cycle>;

OptionalCycle OddCycle(Graph const& g, Removed const& r = {});

// memoization

using ZobristArray = ::std::vector<::std::uint_fast32_t>;
using HashType     = ZobristArray::value_type;

inline ZobristArray GenZobrist(::std::size_t const sz) {
  ::std::mt19937 mt{::std::random_device{}()};
  ::std::uniform_int_distribution<HashType> dist(
    0, ::std::numeric_limits<HashType>::max());
  
  ZobristArray za(sz);
  for (auto&& n : za) {
    n = dist(mt);
  }
  return za;
}

inline ZobristArray::value_type CalcHash(ZobristArray const&, Removed const& r)
    noexcept {
  return ::std::accumulate(
    r.begin(), r.end(), 0,
    [](HashType const total, HashType const current) {
      return total ^ current;
    }
  );
}

using RemovedMem = ::std::unordered_map<
  ZobristArray::value_type,
  ::std::vector<::std::pair<Removed, ::std::size_t>>>;

::std::size_t Optima(Graph const& g, Removed& r, RemovedMem& mem,
                     ZobristArray const& z, HashType hash,
                     ::std::size_t b =
                       ::std::numeric_limits<::std::size_t>::max());

Graph& FragmentIntersection(Graph& g);


}