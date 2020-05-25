#include "diploid.hpp"

#include "spoa/spoa.hpp"
#include "biosoup/progress_bar.hpp"

#include <iostream>
#include <future>
#include <string>
#include <tuple>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <memory>

namespace raven {
  
namespace detail {
  
struct PercentageInfo {
  char primary;
  char secondary;
  
  float primary_p;
  float secondary_p;
  
  float err_p;
  
  ::std::uint_fast32_t depth;
};
  
PercentageInfo ComputePercentages(::std::vector<::std::uint_fast32_t> const& v) {
  using CountBase = ::std::pair<::std::uint_fast32_t, char>;
  
  auto primary = ::std::max<CountBase>(
    ::std::max<CountBase>({v['A'], 'A'}, {v['T'], 'T'}),
    ::std::max<CountBase>({v['G'], 'G'}, {v['C'], 'C'})
  );
  
  auto base = primary.second;
  
  auto secondary =
    base == 'A'
      ? ::std::max<CountBase>(::std::max<CountBase>({v['G'], 'G'}, {v['C'], 'C'}), {v['T'], 'T'})
    : base == 'T'
      ? ::std::max<CountBase>(::std::max<CountBase>({v['G'], 'G'}, {v['C'], 'C'}), {v['A'], 'A'})
    : base == 'G'
      ? ::std::max<CountBase>(::std::max<CountBase>({v['A'], 'A'}, {v['C'], 'C'}), {v['T'], 'T'})
      : ::std::max<CountBase>(::std::max<CountBase>({v['G'], 'G'}, {v['A'], 'A'}), {v['T'], 'T'});
      
  PercentageInfo pi;
  pi.depth = v['A'] + v['T'] + v['G'] + v['C'];
  
  pi.primary = base;
  pi.secondary = secondary.second;
  
  pi.primary_p = static_cast<float>(primary.first) / pi.depth;
  pi.secondary_p = static_cast<float>(secondary.first) / pi.depth;
  
  pi.err_p = 1.0f - pi.primary_p - pi.secondary_p;
  
  return pi;
}
  
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

Cycle merge(Cycle const& u, Cycle const& v) {
  return u.size() < v.size() ? Cycle{v.rbegin() + u.size(), v.rend()}
                             : Cycle{u.rbegin() + v.size(), u.rend()};
}

OptionalCycle OddCycleImpl(Graph const& g, VisitedDephts& v, Removed const& r,
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
      auto ret = OddCycleImpl(g, v, r, c, nb, d);
      if (ret) {
        return ret;
      }
    }
  }
  
  return {};
}

OptionalCycle OddCycle(Graph const& g, Removed const& r = {}) {
  VisitedDephts v;
  for (auto&& pair : g) {
    if (!r.count(pair.first) && !v.count(pair.first)) {
      auto ret = OddCycleImpl(g, v, r, {}, pair.first);
      if (ret) {
        return ret;
      }
    }
  }
  return {};
}

::std::size_t Optima(Graph const& g, Removed const& r,
                     ::std::size_t b =
                       ::std::numeric_limits<::std::size_t>::max()) {
  if (r.size() >= b) {
    return ::std::numeric_limits<decltype(b)>::max();
  } else {
    auto cycle_opt = OddCycle(g, r);
    if (cycle_opt) {
      auto n = b;
      auto r_copy = r;
      for (auto&& node : *cycle_opt) {
        auto loc = r_copy.insert(node).first;
        n = ::std::min(n, Optima(g, r_copy, b));
        r_copy.erase(loc);
      } 
      return n;  
    }
    return r.size();
  }
}

Graph& FragmentIntersection(Graph& g) {
    ::biosoup::ProgressBar bar{
      static_cast<::std::uint32_t>(g.size()), 80ul};
    ::std::cerr << "[raven::diploid::FragmentIntersection] finding maximum fragment set"
                << "\t" << "[" << bar << "]";
    Removed d;
    auto s = Optima(g, {});
    for (auto&& pair : g) {
      if (Optima(g, {pair.first}) == s) {
        d.insert(pair.first);
      }
      if (++bar) {
        ::std::cerr << "\r\t" << "[" << bar << "]";
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

}  // namespace detail

DiploidSequences Partition(
  ::std::vector<::std::unique_ptr<::biosoup::Sequence>> const& sequences,
  ::std::shared_ptr<::thread_pool::ThreadPool>) {
  
  
  
  auto alignment_engine = ::spoa::createAlignmentEngine(
    static_cast<::spoa::AlignmentType>(0), 5, -4, -8);
  
  auto graph = ::spoa::createGraph();
  ::biosoup::ProgressBar bar{static_cast<::std::uint32_t>(sequences.size()), 80ul};
  
  ::std::cerr << "[raven::diploid::Partition] aligning sequences "
              << "[" << bar << "]";
  
  for (auto&& s : sequences) {
    if (++bar) {
      ::std::cerr << "\r"
                  << "[raven::diploid::Partition] aligning sequences "
                  << "[" << bar << "]";
    }
    
    graph->add_alignment(alignment_engine->align(s->data, graph), s->data);
  }
  
  ::std::vector<::std::string> msa;
  graph->generate_multiple_sequence_alignment(msa);
  
  ::std::cerr << "\n"
              << "[raven::diploid::Partition] "
              << "generated multiple sequence alignment vector"
              << "\n";
  
  ::std::vector<::std::uint_fast32_t> offsets(sequences.size());
  ::std::vector<::std::uint_fast32_t> base_count('T' + 1, 0);
  
  auto clear_base_count = [&base_count]{
    base_count['A'] = base_count['T'] = base_count['G'] = base_count['C'] = 0;
  };
  
  auto const min_depth     = 5;
  auto const max_err       = 0.1f;
  auto const min_secondary = 0.3f;
  
  // TODO: grouping columns into batches and then using thread pool
  
  ::std::vector<::std::vector<::std::int_fast8_t>> snp_m{msa.size()};
  
  for (::std::uint_fast32_t i = 0; i < msa.front().size(); ++i) {
    for (decltype(i) j = 0; j < msa.size(); ++j) {
      ++base_count[msa[j][i]];
      if (msa[j][i] != '-')
        ++offsets[j];
    }
    
    auto info = detail::ComputePercentages(base_count);
    clear_base_count();
    
    if (info.depth < min_depth || info.err_p > max_err
          || info.secondary_p < min_secondary) {
      continue;
    }
    
    for (decltype(i) j = 0; j < msa.size(); ++j) {
      snp_m[j].push_back(msa[j][i] == info.primary
                           ? 1
                           : (msa[j][i] == info.secondary ? -1 : 0));
      
      /*if (snp_m[j].back() != 2) {
        ::std::cerr << "[raven::Diploid::Partition] snip on " << j << " at "
                    << i << " -> " << snp_m[j].back()
                    << "\n";
      }*/
    }
  }
  
  ::std::cerr << "[raven::Diploid::Partition] built SNP matrix" << "\n";
  
  detail::Graph conflict_graph;
  
  ::biosoup::ProgressBar graph_bar{
    static_cast<::std::uint32_t>(snp_m.size()), 80ul};
  
  ::std::cerr << "[raven::diploid::Partition] building fragment conflict graph "
              << "[" << graph_bar << "]";
  
  for (::std::uint_fast32_t i = 0; i < snp_m.size(); ++i) {
    if (++graph_bar) {
      ::std::cerr << "\r"
                  << "[raven::diploid::Partition] building fragment conflict graph "
                  << "[" << graph_bar << "]";
    }
    
    for (decltype(i) j = 0; j < snp_m.size(); ++j) {
      if (i != j) {
        for (decltype(i) k = 0; k < snp_m[i].size(); ++k) {
          if (snp_m[i][k] + snp_m[j][k] == 0) {
            conflict_graph[i].insert(j);
            conflict_graph[j].insert(i);
          }
        }
      }
    }
  }
  
  msa.clear();
  msa.shrink_to_fit();
  
  snp_m.clear();
  snp_m.shrink_to_fit();
  
  auto old     = conflict_graph.size();
  auto reduced = detail::FragmentIntersection(conflict_graph).size();
  
  ::std::cerr << "\n" << "Regular graph size: " << old
              << ", reduced: " << reduced
              << "\n";
  
  return {};
}

}  // namespace raven