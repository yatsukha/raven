#include "mfr.hpp"

namespace raven {

namespace detail {

Cycle merge(Cycle const& u, Cycle const& v) {
  return u.size() < v.size() ? Cycle{v.rbegin() + u.size(), v.rend()}
                             : Cycle{u.rbegin() + v.size(), u.rend()};
}

OptionalCycle OddCycle(Graph const& g, VisitedDephts& v, Removed const& r,
                           Optional<Node> p = NullOpt{}, Node c = 0u,
                           ::std::uint_fast32_t d = 0u) {
  if (r.count(c)) {
    return {};
  }
  
  v[c] = {d++, p};
  
  for (auto&& nb : g.find(c)->second) {
    auto iter = v.find(nb);
    if (iter != v.end()) {
      if ((d - iter->second.first) % 2) {
        auto cycle_from = [&v](Optional<Node> const& opt_node, decltype(d) const& depth) {
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
        return merge(a, b);
      }
    } else {
      auto ret = OddCycle(g, v, r, c, nb, d);
      if (ret) {
        return ret;
      }
    }
  }
  
  return {};
}
  
}

OptionalCycle OddCycle(Graph const& g, Removed const& r) {
  VisitedDephts v;
  for (auto&& pair : g) {
    if (!r.count(pair.first) && !v.count(pair.first)) {
      auto ret = detail::OddCycle(g, v, r, {}, pair.first);
      if (ret) {
        return ret;
      }
    }
  }
  return {};
}

::std::size_t Optima(Graph const& g, Removed& r, RemovedMem& mem,
                     ZobristArray const& z, HashType hash, ::std::size_t b) {
  if (r.size() >= b) {
    return ::std::numeric_limits<decltype(b)>::max();
  } else {
    auto cycle_opt = OddCycle(g, r);
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
        auto loc = r.insert(node).first;
        n = ::std::min(n, Optima(g, r, mem, z, hash ^ z[node], b));
        r.erase(loc);
      }
      mem[hash].emplace_back(r, n);
      return n;  
    }
    return r.size();
  }
}

Graph& FragmentIntersection(Graph& g) {
    Removed d, tmp;
    
    ZobristArray z = GenZobrist(g.size());
    RemovedMem mem;
    
    ::std::cerr << "aaa" << "\n";
    
    auto s = Optima(g, tmp, mem, z, 0);
    
    ::std::cerr << "aaa" << ::std::endl;
    
    ::biosoup::ProgressBar bar{
      static_cast<::std::uint32_t>(g.size()), 80ul};
    ::std::cerr << "\n[raven::diploid::FragmentIntersection] finding maximum fragment set "
                << "[" << bar << "]";
    
    auto idx = 0;
    for (auto&& pair : g) {
      auto iter = tmp.insert(pair.first).first;
      if (Optima(g, tmp, mem, z, z[pair.first]) == s) {
        d.insert(pair.first);
      }
      tmp.erase(iter);
      ::std::cerr << ++idx << ::std::endl;
      if (++bar) {
        ::std::cerr << "\r[raven::diploid::FragmentIntersection] finding maximum fragment set " << "[" << bar << "]";
      }
    }
    
    for (auto&& node : d) {
      for (auto&& nb : g[node]) {
        g[nb].erase(g[nb].find(node));
      }
      g.erase(g.find(node));
    }
      
    return g;
}

}