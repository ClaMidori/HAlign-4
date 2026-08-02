#ifndef CONFIG_HPP
#define CONFIG_HPP

// ------------------------------------------------------------------
// config.hpp
// Description (detailed English comments)
//
// This file centrally defines global configuration constants, logger initialization functions, and auxiliary types and tools related to command-line parsing (CLI11) (such as custom formatters and validators).
//
// Purpose:
// - Provide a single-entry "configuration library" for the project, making it easy to reference consistent symbols throughout the code (e.g., working directory structure, default command templates, log filenames, etc.);
// - Provide convenient logger initialization functions `setupLogger` / `setupLoggerWithFile`, making it easy to uniformly configure log output to console and file in main;
// - Provide CLI beautification and input trimming (trim_whitespace), improving command-line experience and fault tolerance.
//
// Notes:
// - The default MSA command template (DEFAULT_MSA_CMD) here is for placeholder and demonstration purposes only; in production, please replace with actual available MSA tools and parameters.
// - All path constants are string literals (relative paths), and when used, they are usually concatenated with `workdir` to form absolute or work-relative paths.
// ------------------------------------------------------------------

// ------------------------------------------------------------------
// Include headers: functional modules include thread pool, command-line parsing, logging system, serialization library, etc.
// ------------------------------------------------------------------
#include <CLI/CLI.hpp>                       // CLI11 command-line parsing library (local version)
#include "spdlog/spdlog.h"                       // spdlog main header file
#include "spdlog/sinks/stdout_color_sinks.h"     // Console color output sink
#include "spdlog/sinks/basic_file_sink.h"        // File output sink
#include "spdlog/async.h"                        // Asynchronous logging support

#include <filesystem>
#include <sstream>
#include <cinttypes>
#include <random>
#include <chrono>
#include <iomanip>

// ------------------------------------------------------------------
// General configuration constants
// ------------------------------------------------------------------
#define VERSION "2.0.0"                   // Version number, can be printed at program startup for tracking
#define LOGGER_NAME "logger"              // Default logger name (for spdlog registration)
#define LOGGER_FILE "halign4.log"         // Default log filename (relative to working directory)
#define CONFIG_FILE "config.json"         // Default config file path (if external config is supported in the future)

const std::string MINIPOA_CMD = "minipoa {input} -S -t {thread} -r1 > {output}"; // Minipoa multi-sequence alignment command template example
const std::string MAFFT_MSA_CMD = "mafft --thread {thread} --auto {input} > {output}"; // MAFFT multi-sequence alignment command template example
const std::string CLUSTALO_MSA_CMD = "clustalo -i {input} -o {output} --threads {thread}"; // Clustal Omega multi-sequence alignment command template example

const std::string DEFAULT_MSA_CMD = MINIPOA_CMD; // Default multi-sequence alignment command template

// ------------------------------------------------------------------
// resolveMsaCmdTemplate: Parse the content entered by the user in -p/--msa-cmd into the "final command template".
//
// Requirements:
// - When the user enters minipoa / mafft / clustalo, automatically use the corresponding built-in template command;
// - When the user enters a custom template (containing placeholders like {input}/{output}), keep it as is;
// - When the user does not enter -p, use DEFAULT_MSA_CMD, without changing existing default behavior.

// Design notes (correctness/usability):
// - Cannot treat -p as a "file path" for validation (ExistingFile / requireRegularFile), because these tool names usually depend on PATH.
// - Here, only map cases that are "exactly equal to keywords", to avoid accidentally affecting user custom commands (e.g., "mafft --auto ...").
// ------------------------------------------------------------------
inline std::string resolveMsaCmdTemplate(const std::string& user_value) {
    // trim: remove leading/trailing whitespace to avoid user accidentally inputting spaces causing keyword match failure
    const auto start = user_value.find_first_not_of(" \t\n\r");
    if (start == std::string::npos) {
        return DEFAULT_MSA_CMD;
    }
    const auto end = user_value.find_last_not_of(" \t\n\r");
    std::string v = user_value.substr(start, end - start + 1);

    // tolower: keywords are case insensitive
    for (char& c : v) {
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }

    if (v == "minipoa") return MINIPOA_CMD;
    if (v == "mafft") return MAFFT_MSA_CMD;
    if (v == "clustalo") return CLUSTALO_MSA_CMD;

    // Other cases: considered user custom template, return as is
    return user_value;
}

// Working directory system
const std::string WORKDIR_DATA = "data";         // Raw data directory
const std::string WORKDIR_TMP = "temp";         // Temporary directory (short lifecycle files)
const std::string RESULTS_DIR = "result";      // Final results directory

const std::string DATA_RAW = "raw_data";        // Raw data subdirectory (original input downloaded/copied)
const std::string DATA_CLEAN = "clean_data";    // Cleaned data subdirectory (sequences written after preprocessing)

// Consensus related filenames (relative to DATA_CLEAN)
const std::string CLEAN_CONS_UNALIGNED = "consensus_unaligned.fasta"; // Consensus sequence filename (unaligned)
const std::string CLEAN_CONS_ALIGNED = "consensus_aligned.fasta";     // Consensus sequence filename (aligned)

const std::string CLEAN_CONS_FASTA = "consensus.fasta"; // Final consensus sequence FASTA filename
const std::string CLEAN_CONS_JSON = "consensus.json";   // Consensus statistics/count output (JSON)

// ------------------------------------------------------------------
// Alignment output related filenames (relative to RESULTS_DIR)
// ------------------------------------------------------------------
// Notes: These filenames are used for multi-sequence alignment (MSA) and insertion sequence processing intermediate files and final output files
// Purpose: Centralize filename constants management, avoid hardcoding string literals in code, easy to uniformly modify and maintain

// Final alignment result filename (MSA output after all sequences aligned)
#define FINAL_ALIGNED_FASTA "final_aligned.fasta"

// Insertion sequence related filenames
#define ALL_INSERTION_FASTA "all_insertion.fasta"           // Merged all insertion sequences (unaligned)
#define ALIGNED_INSERTION_FASTA "aligned_insertion.fasta"   // Aligned insertion sequences (MSA result)

// Thread level output filename templates (for parallel writing)
// Notes: Each thread writes SAM files independently to avoid thread competition; finally merged by main thread
#define THREAD_SAM_PREFIX "thread"                          // Thread SAM file prefix ("thread" + tid + ".sam")
#define THREAD_SAM_SUFFIX ".sam"                            // Thread SAM file suffix
#define THREAD_INSERTION_SAM_SUFFIX "_insertion.sam"        // Thread insertion sequence SAM file suffix ("thread" + tid + "_insertion.sam")

// ------------------------------------------------------------------
// Debug and integer precision configuration
// ------------------------------------------------------------------
#ifndef DEBUG
#define DEBUG 0
#endif

#ifndef M64
#define M64 0
#endif

// Switch between 32-bit or 64-bit integers based on macros
// Purpose: When processing very large data counts, using 64-bit can avoid overflow; for development/lightweight runs, 32-bit saves memory.
#if M64
typedef int64_t	int_t;
typedef uint64_t uint_t;
#define PRIdN	PRId64
#define U_MAX	UINT64_MAX
#define I_MAX	INT64_MAX
#define I_MIN	INT64_MIN
#else
typedef int32_t int_t;
typedef uint32_t uint_t;
#define PRIdN	PRId32
#define U_MAX	UINT32_MAX
#define I_MAX	INT32_MAX
#define I_MIN	INT32_MIN
#endif

// Get hardware concurrency thread count (fallback to 1).
static int get_default_threads() {
    unsigned int hc = std::thread::hardware_concurrency();
    return static_cast<int>(hc ? hc : 1u);
}

// ------------------------------------------------------------------
// Default workdir generator (key logic added)
//
// Requirements: When user does not pass -w/--workdir, automatically use "./tmp-random number".
// Design points:
// 1) What is returned here is a relative path string ("./tmp-..."), consistent with CLI input;
// 2) Use "timestamp + random number" concatenation to reduce collision probability in concurrent/repeated runs;
// 3) Do not create directories here: Directory creation/cleanup strategy is still handled uniformly by checkOption()->file_io::prepareEmptydir,
//    to ensure existing processes and error handling logic remain unchanged.
// ------------------------------------------------------------------
static std::string makeDefaultWorkdir() {
    using Clock = std::chrono::high_resolution_clock;
    const auto now = Clock::now().time_since_epoch();
    const auto ns = std::chrono::duration_cast<std::chrono::nanoseconds>(now).count();

    // Use random_device as seed, generate small random perturbation, further reduce collision probability when starting in the same nanosecond.
    std::random_device rd;
    std::mt19937_64 gen(static_cast<uint64_t>(rd()) ^ static_cast<uint64_t>(ns));
    std::uniform_int_distribution<uint32_t> dist(0u, 0xFFFFFFFFu);
    const uint32_t r = dist(gen);

    std::ostringstream oss;
    oss << "./tmp-" << ns << "-" << std::hex << std::setw(8) << std::setfill('0') << r;
    return oss.str();
}

struct Options {
	// Input/Output and Working Directory
	std::string input;          // -i: input sequence file (path or copy file)
	std::string output;         // -o: final output file (write location)
	std::string workdir;        // -w: working directory, all intermediate files (data/raw, data/clean, etc.) are placed in this directory

	// Optional parameters: center sequence, MSA command template
	std::string center_path;    // -c: optional, specify center sequence file path, if specified, bypass automatic selection
	std::string msa_cmd;        // -p: command template for MSA on consensus sequence (can contain {input} {output} {thread} placeholders)

	// Parallel and algorithm parameters
	int threads = get_default_threads(); // -t: number of threads, default is CPU core count
	int kmer_size = 15;         // --kmer-size: k-mer size for classification/clustering (used in subsequent steps)
	int kmer_window = 10;       // --kmer-window: minimizer window size w (in number of k-mers)
	int cons_n = 1000;          // --cons-n: number of sequences selected for consensus calculation (Top-K by length)
	int sketch_size = 2000;     // --sketch-size: size for sketch (default 2000)

	// keep length related switch:
	// - keep_length: keep the length of the "first/center sequence" unchanged (other sequences can change/fill according to alignment results), suitable for scenarios where only the output consensus/center sequence length is concerned.
	bool keep_length = false; // --keep-length
	// workdir cleanup switch:
	// - save_workdir: if true, keep the working directory after alignment/output completion; if false (default), delete the working directory after successful completion.
	bool save_workdir = false;      // --save-workdir
};


// setupCli: Define CLI parameters and bind to Options
// Notes:
// - Use CLI11 library for parameter parsing, support short/long parameters and basic validation (e.g., ExistingFile)
// - For parameters like -p/--msa-cmd that might be executable names (not full paths), ExistingFile will reject just command names;
//   If you want to allow command names (resolved in PATH), remove check(CLI::ExistingFile) or change to more lenient user-level judgment.
static void setupCli(CLI::App& app, Options& opt) {
    app.description("HAlign 4: A New Strategy for Rapidly Aligning Millions of Sequences");

    // Set version flag: --version prints version and exits, -h/--help will also display version information
    app.set_version_flag("-v,--version", std::string("halign4 version ") + VERSION);

   // Required parameters (support both short and long parameters)
    // -i/--input: input sequences (FASTA).
    // Note:
    // - Must be a local file path (current implementation checks existence during parameter validation);
    // - File format recommended .fasta/.fa (internal reading uses kseq).
    // Example: -i test/data/mt1x/mt1x.fasta
    app.add_option("-i,--input", opt.input,
                   "Input sequences in FASTA format (local file path).")
        ->required()
        ->check(CLI::ExistingFile);

    // -o/--output: final output aligned sequences (FASTA).
    // Note:
    // - Output path does not need to exist in advance; parent directory is recommended to exist (if not, writing may fail later).
    // - Output content is multi-sequence aligned FASTA.
    // Example: -o out/aligned.fasta
    app.add_option("-o,--output", opt.output,
                   "Output aligned sequences (FASTA file path).")
        ->required();

    // workdir: changed to optional.
    // If user does not provide -w, generate default value with makeDefaultWorkdir() after CLI parsing in main().
    // Reason for doing this: CLI11's default_val will display the default value directly in -h; but here the default value has random numbers,
    // displaying it would make the help information "different every time", and mislead users into thinking they must specify a fixed path.
    // -w/--workdir: working directory (stores intermediate files/logs/temporary results).
    // Note:
    // - If not provided, program will automatically generate ./tmp-<random>;
    // - Directory will create subdirectories like data/raw_data, data/clean_data, temp, result;
    // - In Release mode, requires workdir to be an empty directory (to avoid overwriting old results); in Debug mode, allows reuse (for iteration).
    app.add_option("-w,--workdir", opt.workdir,
                   "Working directory for intermediate files (default: ./tmp-<random>). ")
        ->capture_default_str();

    // Optional parameters (add long parameter forms)
    // -c/--center-path: specify center/reference sequence (FASTA).
    // Note:
    // - Not provided: program will automatically select and generate consensus/center sequence in preprocessing stage;
    // - Provided: will use this sequence as reference (and manage uniformly in workdir).
    // Typical use: in COVID data sets, can use covid-ref's first (Wuhan reference) as center.
    app.add_option("-c,--center-path", opt.center_path,
                   "Center/reference sequence in FASTA (optional). If not set, a consensus/center is generated.")
        ->check(CLI::ExistingFile);

    // If -p is "executable file path", ExistingFile usually works;
    // If you want to allow only command names (resolved in PATH), don't check
    // msa-cmd: supports keyword or custom command template.
    // - Keywords: minipoa / mafft / clustalo
    // - Custom template: e.g., "mafft --auto {input} > {output}"
    // Note: cannot use ExistingFile validation here, otherwise keywords/command names will be incorrectly rejected.
    // -p/--msa-cmd: high-quality MSA tool (for high-quality alignment of consensus/insertion sequences).
    // Supports two forms:
    // 1) Keywords: minipoa / mafft / clustalo
    //    - After entering keywords, program will automatically expand to built-in templates (see MINIPOA_CMD/MAFFT_MSA_CMD/CLUSTALO_MSA_CMD);
    // 2) Custom "command template string": must contain at least {input} and {output}; optionally {thread}.
    //    - E.g., "mafft --thread {thread} --auto {input} > {output}"
    // Note:
    // - Default uses minipoa (equivalent to -p minipoa if not passed);
    // - This command will do a tiny.fasta smoke test during parameter validation; if the environment lacks this tool, it will directly error.
    app.add_option("-p,--msa-cmd", opt.msa_cmd,
                   "High-quality MSA method: keyword {minipoa|mafft|clustalo} or a custom command template containing {input} and {output} (optional {thread}).");

    // -t/--thread: number of threads.
    // Note:
    // - Default value is hardware concurrency (std::thread::hardware_concurrency);
    // - Affects preprocessing, alignment, and {thread} replacement in external MSA commands.
    app.add_option("-t,--thread", opt.threads, "Number of threads.")
        ->default_val(get_default_threads())
        ->check(CLI::Range(1, 100000));

    // --kmer-size: k-mer size.
    // Note:
    // - Parameter for minimizer/hash related processes; generally no need to change.
    // - Legal range [4,31] (consistent with some bit operations/encoding implementation constraints).
    app.add_option("--kmer-size", opt.kmer_size, "K-mer size used in sketch/minimizer.")
        ->default_val(15)
        ->check(CLI::Range(4, 31));

    // --kmer-window: minimizer window size w (unit: number of k-mers).
    // Note: w larger, minimizer sparser; w smaller, seeds denser but possibly slower.
    app.add_option("--kmer-window", opt.kmer_window,
                   "Minimizer window size w (in number of k-mers).")
        ->default_val(10)
        ->check(CLI::Range(1, 1000000));

    // --cons-n: Top-N for generating consensus/center (selected by length).
    // Note:
    // - If input sequence count <= cons_n, program will directly call external MSA to do one alignment on all sequences (fast path);
    // - If input sequence count >> cons_n, first generate consensus with Top-N, then do batch/reference alignment.
    app.add_option("--cons-n", opt.cons_n,
                   "Number of sequences used to build the consensus/center (Top-N by length).")
        ->default_val(1000)
        ->check(CLI::Range(1, 1000000));

    // --sketch-size: sketch (minhash) size.
    // Note: larger, more robust but slower/more memory; default 2000 is usually sufficient.
    app.add_option("--sketch-size", opt.sketch_size, "Sketch size (minhash count).")
        ->default_val(2000)
        ->check(CLI::Range(1, 10000000));


    app.add_flag("--keep-length", opt.keep_length,
        "Keep all reference sequences lengths unchanged. ");

    // workdir management: whether to keep working directory after completion
    // --save-workdir: keep working directory (default will delete).
    app.add_flag("--save-workdir", opt.save_workdir,
        "Keep the working directory after completion (default: remove). Useful for debugging.");
}

// logParsedOptions: Output parsed parameters in a nice table format to log
// Notes: The output here is for user help and debugging (print truncated long strings, boolean friendly display, etc.),
// does not affect program behavior. If the program runs in headless environment (service/container), logs are also convenient for auditing and reproducing run parameters.
static void logParsedOptions(const Options& opt) {
    // Helper to convert values and truncate long strings for tidy display
    auto toString = [](const std::string& s, size_t maxLen) -> std::string {
        if (s.empty()) return "(empty)";
        if (s.size() <= maxLen) return s;
        return s.substr(0, maxLen - 3) + "...";
    };

    auto boolToStr = [](bool b) { return b ? "true" : "false"; };

    const size_t keyW = 14;
    const size_t valW = 60;
    const size_t innerW = keyW + 3 + valW; // "key : value"

    std::vector<std::pair<std::string, std::string>> rows = {
        {"input", toString(opt.input, valW)},
        {"output", toString(opt.output, valW)},
        {"workdir", toString(opt.workdir, valW)},
        {"center-path", toString(opt.center_path, valW)},
        {"msa_cmd", toString(opt.msa_cmd, valW)},
        {"threads", std::to_string(opt.threads)},
        {"kmer-size", std::to_string(opt.kmer_size)},
        {"kmer-window", std::to_string(opt.kmer_window)},
        {"cons_n", std::to_string(opt.cons_n)},
        {"sketch_size", std::to_string(opt.sketch_size)},
        {"keep-length", boolToStr(opt.keep_length)},
        {"save-workdir", boolToStr(opt.save_workdir)}
    };

    std::ostringstream oss;

    // top border
    oss << "+" << std::string(innerW, '-') << "+\n";

    // title centered
    const std::string title = " Parsed options ";
    size_t paddingLeft = 0;
    if (innerW > title.size()) paddingLeft = (innerW - title.size()) / 2;
    oss << "|" << std::string(paddingLeft, ' ') << title
        << std::string(innerW - paddingLeft - title.size(), ' ') << "|\n";

    // separator
    oss << "+" << std::string(innerW, '-') << "+\n";

    // rows
    for (auto &kv : rows) {
        oss << "| " << std::left << std::setw(keyW) << kv.first << " : "
            << std::setw(valW) << kv.second << "|\n";
    }

    // bottom border
    oss << "+" << std::string(innerW, '-') << "+";

    spdlog::info("\n{}", oss.str());
}

// ------------------------------------------------------------------
// CLI11 custom formatter (beautify option output)
// Notes:
// - Custom `make_option_opts` can display parameter type and default value in help, easy for users to understand;
// - Custom `make_usage` provides friendlier usage examples and descriptions, convenient for beginners to get started quickly.
// ------------------------------------------------------------------
class CustomFormatter : public CLI::Formatter {
public:
	CustomFormatter() : Formatter() {}

	// Customize parameter display style (with default values)
	std::string make_option_opts(const CLI::Option* opt) const override {
		if (opt->get_type_size() == 0) return "";
		std::ostringstream out;
		out << " " << opt->get_type_name();
		if (!opt->get_default_str().empty())
			out << " (default: " << opt->get_default_str() << ")";
		return out.str();
	}

	// A usage example is provided; this example may be updated based on the actual executable name of your project.
	std::string make_usage(const CLI::App* app, std::string name) const override {
		std::ostringstream out;
		out << "Usage:\n"
			<< "  ./halign4 -i <ref.fa> -o <output.fa> -w </path/to/workdir> [options]\n\n"
			<< "Example:\n"
			<< "  ./halign4 -i ref.fa -o output.fa -w ./tmp -t 8\n\n";
		return out.str();
	}
};

// ------------------------------------------------------------------
// CLI11 custom validator: automatically trim whitespace from parameter sides
// Notes: Sometimes users accidentally add spaces or copy-paste with newlines in command line, trim_whitespace can improve robustness
// ------------------------------------------------------------------
inline CLI::Validator trim_whitespace = CLI::Validator(
	[](std::string& s) {
		auto start = s.find_first_not_of(" \t\n\r");
		auto end = s.find_last_not_of(" \t\n\r");
		s = (start == std::string::npos) ? "" : s.substr(start, end - start + 1);
		return std::string();  // An empty string indicates that the verification passed.
	}, ""
);



// ------------------------------------------------------------------
// Logger system initialization
// Notes:
// - `setupLoggerWithFile(path)` will create an asynchronous spdlog logger, output to console and log file in specified directory;
// - Asynchronous logging (async_logger) can reduce blocking and I/O latency in high concurrency logging scenarios, but needs to configure thread pool at program startup;
// - `setupLogger()` only outputs to console, suitable for interactive debugging or short-term running;
// - Both functions will set default log level (Debug/Info) and periodic flush (flush_every), which helps retain logs in case of crash.
// ------------------------------------------------------------------
inline void setupLoggerWithFile(std::filesystem::path log_dir) {
	auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
	console_sink->set_pattern("%^[%Y-%m-%d %H:%M:%S] [%l] %v%$");

	std::filesystem::path log_file = log_dir / LOGGER_FILE;
	auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(log_file.string(), true);
	file_sink->set_pattern("[%Y-%m-%d %H:%M:%S] [%l] %v");

	spdlog::sinks_init_list sinks = { console_sink, file_sink };
	auto logger = std::make_shared<spdlog::async_logger>(
		LOGGER_NAME, sinks.begin(), sinks.end(), spdlog::thread_pool(), spdlog::async_overflow_policy::block);

	spdlog::set_default_logger(logger);
#ifdef _DEBUG
	spdlog::set_level(spdlog::level::debug);
#else
	spdlog::set_level(spdlog::level::info);
#endif
	spdlog::flush_every(std::chrono::seconds(3));
}

// Console logging (for development/debugging)
inline void setupLogger() {
	auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
	console_sink->set_pattern("%^[%Y-%m-%d %H:%M:%S.%e] [%l] %v%$");

	spdlog::sinks_init_list sinks = { console_sink };
	auto logger = std::make_shared<spdlog::async_logger>(
		LOGGER_NAME, sinks.begin(), sinks.end(), spdlog::thread_pool(), spdlog::async_overflow_policy::block);

	spdlog::set_default_logger(logger);
#ifdef _DEBUG
	spdlog::set_level(spdlog::level::debug);
#else
	spdlog::set_level(spdlog::level::info);
#endif

	spdlog::flush_every(std::chrono::seconds(3));
}

// Get full command line string (for logging and reproduction)
inline std::string getCommandLine(int argc, char** argv) {
	std::ostringstream cmd;
	for (int i = 0; i < argc; ++i) {
		cmd << argv[i];
		if (i != argc - 1) cmd << " ";
	}
	return cmd.str();
}

#endif // CONFIG_HPP
