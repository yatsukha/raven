#include "mfr.hpp"

#include <cassert>

#include "biosoup/progress_bar.hpp"

namespace raven {

namespace diploid {

namespace detail {

using ZobristArray = ::std::vector<::std::uint_fast32_t>;
using HashType = ZobristArray::value_type;

ZobristArray GenZobrist(::std::size_t const sz) {
  ::std::mt19937 mt{::std::random_device{}()};
  ::std::uniform_int_distribution<HashType> dist(
      0, ::std::numeric_limits<HashType>::max());

  ZobristArray za(sz);
  for (auto&& n : za) {
    n = dist(mt);
  }
  return za;
}

ZobristArray::value_type CalcHash(ZobristArray const& z,
                                  ConflictGraph::Removed const& r) noexcept {
  return ::std::accumulate(r.begin(), r.end(), 0,
                           [&z](HashType const total, HashType const current) {
                             return total ^ z[current];
                           });
}

using RemovedMem = ::std::unordered_map<
    ZobristArray::value_type,
    ::std::vector<::std::pair<ConflictGraph::Removed, ::std::size_t>>>;

ConflictGraph::Cycle merge(ConflictGraph::Cycle const& u,
                           ConflictGraph::Cycle const& v) {
  return u.size() < v.size()
             ? ConflictGraph::Cycle{v.rbegin() + u.size(), v.rend()}
             : ConflictGraph::Cycle{u.rbegin() + v.size(), u.rend()};
}

}  // namespace detail

ConflictGraph::OptionalCycle ConflictGraph::OddCycleImpl(
    VisitedDephts& v, Removed const& r, Optional<Node> p, Node c,
    ::std::uint_fast32_t d) const {
  if (r.count(c)) {
    return {};
  }

  v[c] = {d++, p};

  for (auto&& nb : g.find(c)->second) {
    auto iter = v.find(nb);
    if (iter != v.end()) {
      if ((d - iter->second.first) % 2) {
        auto cycle_from = [&v](Optional<Node> const& opt_node,
                               decltype(d) const& depth) {
          Cycle cc;
          cc.reserve(depth);
          auto jumper = opt_node;
          while (jumper) {
            cc.push_back(*jumper);
            jumper = v[*jumper].second;
          }
          return cc;
        };

        auto a = cycle_from(decltype(p){c}, d);
        auto b = cycle_from(iter->second.second, iter->second.first);
        return detail::merge(a, b);
      }
    } else {
      auto ret = OddCycleImpl(v, r, c, nb, d);
      if (ret) {
        return ret;
      }
    }
  }

  return {};
}

ConflictGraph::OptionalCycle ConflictGraph::OddCycle(Removed const& r) const {
  VisitedDephts v;

  ::std::unordered_map<
      detail::HashType,
      ::std::vector<::std::pair<Removed, OptionalCycle>>> static mem;
  detail::ZobristArray static z = detail::GenZobrist(g.size());

  auto hash = detail::CalcHash(z, r);

  if (mem.count(hash)) {
    for (auto&& m : mem[hash]) {
      if (r == m.first) {
        return m.second;
      }
    }
  }

  for (auto&& pair : g) {
    if (!r.count(pair.first) && !v.count(pair.first)) {
      auto ret = OddCycleImpl(v, r, {}, pair.first);
      if (ret) {
        mem[hash].emplace_back(r, ret);
        return ret;
      }
    }
  }

  mem[hash].emplace_back(r, OptionalCycle{});
  return {};
}

::std::size_t Optima(
    ConflictGraph const& g, ConflictGraph::Removed& r, detail::RemovedMem& mem,
    detail::ZobristArray const& z, detail::HashType hash,
    ::std::size_t b = ::std::numeric_limits<::std::size_t>::max()) {
  if (r.size() >= b) {
    return ::std::numeric_limits<decltype(b)>::max();
  } else {
    auto cycle_opt = g.OddCycle(r);
    if (cycle_opt) {
      if (mem.count(hash)) {
        for (auto&& e : mem[hash]) {
          if (e.first == r) {
            return e.second;
          }
        }
      }
      auto n = b;
      for (auto&& node : *cycle_opt) {
        auto loc = r.insert(node).second;
        n = ::std::min(n, Optima(g, r, mem, z, hash ^ z[node], b));
        r.erase(loc);
      }
      mem[hash].emplace_back(r, n);
      return n;
    }
    return r.size();
  }
}

ConflictGraph& FragmentIntersection(ConflictGraph& cg) {
  auto& g = cg.graph();
  ConflictGraph::Removed d, tmp;

  detail::ZobristArray z = detail::GenZobrist(g.size());
  detail::RemovedMem mem;

  auto s = Optima(cg, tmp, mem, z, 0);

  ::std::cerr << "[raven::diploid::FragmentIntersection] calculated an optimal "
                 "solution of "
              << s << " removed vertices"
              << "\n";

  ::std::cerr << "[raven::diploid::FragmentIntersection] second run optimal "
              << Optima(cg, tmp, mem, z, 0) << "\n";

  ::std::cerr << "[raven::diploid::FragmentIntersection] candidate optimals: ";

  for (auto&& pair : g) {
    auto iter = tmp.insert(pair.first).first;
    auto ss = Optima(cg, tmp, mem, z, z[pair.first]);
    if (ss == s) {
      d.insert(pair.first);
    }
    ::std::cerr << ss << ", ";
    tmp.erase(iter);
  }

  ::std::cerr << "\n";
  ::std::cerr << "[raven::diploid::FragmentIntersection] deletion set size: "
              << d.size() << "\n";

  ::std::for_each(
      d.begin(), d.end(),
      ::std::bind(&ConflictGraph::remove, &cg, ::std::placeholders::_1));

  return cg;
}

}  // namespace diploid

}  // namespace raven