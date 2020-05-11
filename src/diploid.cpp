#include "diploid.hpp"

#include "spoa/spoa.hpp"

#include <iostream>
#include <future>
#include <string>
#include <tuple>

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
  
}  // namespace detail

DiploidSequences Partition(
  ::std::vector<::std::unique_ptr<::biosoup::Sequence>> const& sequences,
  ::std::shared_ptr<::thread_pool::ThreadPool>) {
  
  auto alignment_engine = ::spoa::createAlignmentEngine(
    static_cast<::spoa::AlignmentType>(0), 5, -4, -8);
  
  auto graph = ::spoa::createGraph();
  
  for (auto&& s : sequences) {
    graph->add_alignment(alignment_engine->align(s->data, graph), s->data);
  }
  
  ::std::vector<::std::string> msa;
  graph->generate_multiple_sequence_alignment(msa);
  
  ::std::vector<::std::uint_fast32_t> offsets;
  ::std::vector<::std::uint_fast32_t> base_count('T' + 1, 0);
  
  auto clear_base_count = [&base_count]{
    base_count['A'] = base_count['T'] = base_count['G'] = base_count['C'] = 0;
  };
  
  auto const min_depth     = 0;
  auto const max_err       = 0.1f;
  auto const min_secondary = 0.3f;
  
  // TODO: grouping columns into batches and then using thread pool
  
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
      // TODO: put actual SNP marking here
      if (msa[j][i] == info.primary || msa[j][i] == info.secondary)
        ::std::cerr << "[raven::diploid::Partition] "
                    << "found SNP in \"" << sequences[j]->name << "\" "
                    << "at offset " << offsets[j] - 1
                    << "\n";
    }
  }
  
  return {};
}

}  // namespace raven