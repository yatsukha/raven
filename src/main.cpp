// Copyright (c) 2020 Robert Vaser

#include <getopt.h>

#include <iostream>

#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "biosoup/timer.hpp"
#include "diploid.hpp"
#include "graph.hpp"

std::atomic<std::uint32_t> biosoup::Sequence::num_objects{0};

namespace {

const char* raven_version = RAVEN_VERSION;

static struct option options[] = {
    {"diploid", no_argument, nullptr, 'd'},
    {"polishing-rounds", required_argument, nullptr, 'p'},
    {"match", required_argument, nullptr, 'm'},
    {"mismatch", required_argument, nullptr, 'n'},
    {"gap", required_argument, nullptr, 'g'},
#ifdef CUDA_ENABLED
    {"cuda-poa-batches", optional_argument, nullptr, 'c'},
    {"cuda-banded-alignment", no_argument, nullptr, 'b'},
    {"cuda-alignment-batches", required_argument, nullptr, 'a'},
#endif
    {"graphical-fragment-assembly", required_argument, nullptr, 'f'},
    {"resume", no_argument, nullptr, 'r'},
    {"threads", required_argument, nullptr, 't'},
    {"version", no_argument, nullptr, 'v'},
    {"help", no_argument, nullptr, 'h'},
    {nullptr, 0, nullptr, 0}};

std::unique_ptr<bioparser::Parser<biosoup::Sequence>> CreateParser(
    const std::string& path) {
  auto is_suffix = [](const std::string& s, const std::string& suff) {
    return s.size() < suff.size()
               ? false
               : s.compare(s.size() - suff.size(), suff.size(), suff) == 0;
  };

  if (is_suffix(path, ".fasta") || is_suffix(path, ".fa") ||
      is_suffix(path, ".fasta.gz") || is_suffix(path, ".fa.gz")) {
    try {
      return bioparser::Parser<biosoup::Sequence>::Create<
          bioparser::FastaParser>(path);  // NOLINT
    } catch (const std::invalid_argument& exception) {
      std::cerr << exception.what() << std::endl;
      return nullptr;
    }
  }
  if (is_suffix(path, ".fastq") || is_suffix(path, ".fq") ||
      is_suffix(path, ".fastq.gz") || is_suffix(path, ".fq.gz")) {
    try {
      return bioparser::Parser<biosoup::Sequence>::Create<
          bioparser::FastqParser>(path);  // NOLINT
    } catch (const std::invalid_argument& exception) {
      std::cerr << exception.what() << std::endl;
      return nullptr;
    }
  }

  std::cerr << "[raven::CreateParser] error: file " << path
            << " has unsupported format extension (valid extensions: .fasta, "
            << ".fasta.gz, .fa, .fa.gz, .fastq, .fastq.gz, .fq, .fq.gz)"
            << std::endl;
  return nullptr;
}

void Help() {
  std::cout
      << "usage: raven [options ...] <sequences>\n"
         "\n"
         "  # default output is stdout\n"
         "  <sequences>\n"
         "    input file in FASTA/FASTQ format (can be compressed with gzip)\n"
         "\n"
         "  options:\n"
         "    -p, --polishing-rounds <int>\n"
         "      default: 2\n"
         "      number of times racon is invoked\n"
         "    -m, --match <int>\n"
         "      default: 3\n"
         "      score for matching bases\n"
         "    -n, --mismatch <int>\n"
         "      default: -5\n"
         "      score for mismatching bases\n"
         "    -g, --gap <int>\n"
         "      default: -4\n"
         "      gap penalty (must be negative)\n"
#ifdef CUDA_ENABLED
         "    -c, --cuda-poa-batches <int>\n"
         "       default: 0\n"
         "       number of batches for CUDA accelerated polishing\n"
         "    -b, --cuda-banded-alignment\n"
         "       use banding approximation for polishing on GPU\n"
         "       (only applicable when -c is used)\n"
         "    -a, --cuda-alignment-batches <int>\n"
         "       default: 0\n"
         "       number of batches for CUDA accelerated alignment\n"
#endif
         "    --graphical-fragment-assembly <string>\n"
         "      prints the assemblg graph in GFA format\n"
         "    -d, --diploid\n"
         "      partitions fragments into two haplotype sets,\n"
         "      and assembles each set separately\n"
         "    --resume\n"
         "      resume previous run from last checkpoint\n"
         "    -t, --threads <int>\n"
         "      default: 1\n"
         "      number of threads\n"
         "    --version\n"
         "      prints the version number\n"
         "    -h, --help\n"
         "      prints the usage\n";
}

}  // namespace

int main(int argc, char** argv) {
  std::int32_t num_polishing_rounds = 2;
  std::int8_t m = 3;
  std::int8_t n = -5;
  std::int8_t g = -4;

  std::string gfa_path = "";
  bool resume = false;

  std::uint32_t num_threads = 1;

  std::uint32_t cuda_poa_batches = 0;
  std::uint32_t cuda_alignment_batches = 0;
  bool cuda_banded_alignment = false;

  bool diploid = false;

  std::string optstr = "dp:m:n:g:t:h";
#ifdef CUDA_ENABLED
  optstr += "c:b:a:";
#endif
  int arg;
  while ((arg = getopt_long(argc, argv, optstr.c_str(), options, nullptr)) !=
         -1) {  // NOLINT
    switch (arg) {
      case 'p':
        num_polishing_rounds = atoi(optarg);
        break;
      case 'm':
        m = atoi(optarg);
        break;
      case 'n':
        n = atoi(optarg);
        break;
      case 'g':
        g = atoi(optarg);
        break;
#ifdef CUDA_ENABLED
      case 'c':
        cuda_poa_batches = 1;
        // next text entry is not an option, assuming it's the arg for option c
        if (optarg == NULL && argv[optind] != NULL && argv[optind][0] != '-') {
          cuda_poa_batches = atoi(argv[optind++]);
        }
        // optional argument provided in the ususal way
        if (optarg != NULL) {
          cuda_poa_batches = atoi(optarg);
        }
        break;
      case 'b':
        cuda_banded_alignment = true;
        break;
      case 'a':
        cuda_alignment_batches = atoi(optarg);
        break;
#endif
      case 'f':
        gfa_path = optarg;
        break;
      case 'd':
        diploid = true;
        break;
      case 'r':
        resume = true;
        break;
      case 't':
        num_threads = atoi(optarg);
        break;
      case 'v':
        std::cout << raven_version << std::endl;
        return 0;
      case 'h':
        Help();
        return 0;
      default:
        return 1;
    }
  }

  if (argc == 1) {
    Help();
    return 0;
  }

  if (optind >= argc) {
    std::cerr << "[raven::] error: missing input file!" << std::endl;
    return 1;
  }

  auto sparser = CreateParser(argv[optind]);
  if (sparser == nullptr) {
    return 1;
  }

  biosoup::Timer timer{};
  timer.Start();

  auto thread_pool = std::make_shared<thread_pool::ThreadPool>(num_threads);

  raven::Graph graph{thread_pool};
  if (resume) {
    try {
      graph.Load();
    } catch (std::exception& exception) {
      std::cerr << exception.what() << std::endl;
      return 1;
    }

    std::cerr << "[raven::] loaded previous run " << std::fixed << timer.Stop()
              << "s" << std::endl;

    timer.Start();
  }

  std::vector<std::unique_ptr<biosoup::Sequence>> sequences;
  if ((graph.stage() < -3) ||
      (graph.stage() >= 0 && graph.stage() < num_polishing_rounds)) {
    try {
      sequences = sparser->Parse(-1);
    } catch (const std::invalid_argument& exception) {
      std::cerr << exception.what() << std::endl;
      return 1;
    }

    if (sequences.empty()) {
      std::cerr << "[raven::] error: empty sequences set" << std::endl;
      return 1;
    }

    std::cerr << "[raven::] loaded " << sequences.size() << " sequences "
              << std::fixed << timer.Stop() << "s" << std::endl;

    timer.Start();
  }

  if (diploid) {
    ::std::cerr << "[raven::] starting diploid partitioning"
                << "\n";
    ::raven::diploid::Partition(sequences, thread_pool);
    return 0;
  }

  // TODO: separate this step into two separate ones
  // TODO: possibly means separating gfa output file
  //       and adding two output file parameters for regular output
  graph.Construct(sequences);
  graph.Assemble();
  graph.Polish(sequences, m, n, g, cuda_poa_batches, cuda_banded_alignment,
               cuda_alignment_batches, num_polishing_rounds);
  graph.PrintGFA(gfa_path);

  for (const auto& it : graph.GetUnitigs(num_polishing_rounds > 0)) {
    std::cout << ">" << it->name << std::endl;
    std::cout << it->data << std::endl;
  }

  timer.Stop();
  std::cerr << "[raven::] " << std::fixed << timer.elapsed_time() << "s"
            << std::endl;

  return 0;
}
