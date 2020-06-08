
#pragma once

#include <algorithm>
#include <functional>
#include <future>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace raven {

namespace diploid {

class NullOpt {};

template <typename T>
class Optional {
  ::std::shared_ptr<T> data_;

 public:
  Optional() = default;
  Optional(NullOpt&&){};

  Optional(T const& data) : data_(new T{data}) {}
  Optional(T&& data) : data_(new T{::std::forward<T>(data)}) {}

  Optional(Optional const&) = default;
  Optional& operator=(Optional const&) = default;

  Optional(Optional&&) = default;
  Optional& operator=(Optional&&) = default;

  operator bool() const noexcept { return data_.operator bool(); }

  T const& operator*() const noexcept { return *data_; }
};

class ConflictGraph {
 public:
  using Node = ::std::uint_fast32_t;
  using Neighbours = ::std::unordered_set<Node>;
  using Graph = ::std::unordered_map<Node, Neighbours>;

  ::std::size_t const max_size;

  ConflictGraph(::std::size_t const max_size) : max_size(max_size) {}
  ConflictGraph(Graph&& g, ::std::size_t const max_size)
      : max_size(max_size), g(::std::forward<Graph>(g)) {}

  Graph const& graph() const noexcept { return g; }

  Graph& graph() noexcept { return g; }

  void connect(Node const u, Node const v) {
    g[u].insert(v);
    g[v].insert(u);
  }

  template <typename Iterable>
  void connect(Node const u, Iterable const& c) {
    ::std::for_each(
        ::std::begin(c), ::std::end(c),
        ::std::bind(
            static_cast<void (ConflictGraph::*)(Node const, Node const)>(
                &ConflictGraph::connect),
            this, u, ::std::placeholders::_1));
  }

  void remove(Node const u) {
    for (auto&& v : g[u]) {
      g[v].erase(g[v].find(u));
    }
    g.erase(g.find(u));
  }

  using Removed = ::std::unordered_set<Node>;
  using Cycle = ::std::vector<Node>;

  using OptionalCycle = Optional<Cycle>;

  OptionalCycle OddCycle(Removed const& r = Removed{}) const;

 private:
  using VisitedDephts =
      ::std::unordered_map<Node, ::std::pair<Node, Optional<Node>>>;

  OptionalCycle OddCycleImpl(VisitedDephts& v, Removed const& r,
                             Optional<Node> p = NullOpt{}, Node c = 0u,
                             ::std::uint_fast32_t d = 0u) const;

  Graph g;
};

ConflictGraph& FragmentIntersection(ConflictGraph& g);

}  // namespace diploid

}  // namespace raven