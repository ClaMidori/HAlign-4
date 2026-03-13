#ifndef HALIGN4_PREPROCESS_H
#define HALIGN4_PREPROCESS_H

#include <cstddef>
#include "config.hpp"
#include "utils.h"
#include "consensus.h"

// ==============================================================
// Preprocessing module (preprocess) header file description (detailed English comments)
//
// This module is responsible for preparing the user's raw FASTA input (which can be a local path or remote URL) into standardized data for subsequent analysis,
// and preparing necessary intermediate files for consensus calculation/alignment. Main responsibilities include:
//  1) Copy or download the input file to `data/raw` under the working directory;
//  2) Read input sequences one by one and perform "cleaning/normalization" (e.g., uppercase, replace non-A/C/G/T/U with N, remove illegal characters, etc.);
//  3) Write the cleaned sequences to `data/clean`;
//  4) Maintain a Top-K selector (length priority) to select candidate sequence sets for building consensus (write to `consensus_unaligned.fasta`);
//  5) Return the total number of processed sequences for the upper layer to decide if further merging/more processing is needed.
//
// Important semantics and conventions:
// - `workdir`: Working directory path provided by the caller; this module will create necessary subdirectories (data/raw, data/clean, etc.) under this directory,
//   if the directory does not exist, it will try to create it; if required to be empty (passed and checked by upper layer), data structures will be created in the empty directory;
// - I/O behavior: If `input_path` is a remote URL (e.g., http(s):// or starting with //), this module will download to local;
//   otherwise, copy the local file to the working directory. Download/copy failure will throw an exception (std::runtime_error).
// - Exception and error handling: The function throws an exception (std::runtime_error) when encountering serious I/O or parsing errors; the upper layer should catch and log.
// - Return value: The function returns the number of processed records (uint_t), if the number exceeds the project configuration limit (U_MAX defined in config.hpp),
//   the value will be truncated to U_MAX and a warning will be logged (caller should note).
//
// Performance and concurrency considerations:
// - This function is I/O intensive: For large files (GB level), pay attention to disk bandwidth and buffering (can be tuned via utils::seq_io's io buffering);
// - In high concurrency environments, do not call this function in parallel to write to the same `workdir` to avoid race conditions; if parallel is needed, use different working directories or external coordination.
//
// ==============================================================

// Preprocess input FASTA and return the number of processed sequences (total records processed)
//
// Parameters:
//  - input_path: Input FASTA file path, supports local path or remote URL (string).
//  - workdir: Working directory (string), will create data/raw and data/clean subdirectories under this directory and write intermediate files.
//  - cons_n: Number of sequences to select for subsequent consensus (Top-K, selected by sequence length), default 1000.
//
// Return value:
//  - Returns the actual number of processed sequences (uint_t); if too many entries exceed U_MAX, will truncate to U_MAX and log warning.
//
// Exceptions:
//  - Throws std::runtime_error when unable to create working directory, unable to download/copy input, or failed to read FASTA.
//
// Output (side effects):
//  - Saves original input copy (or downloaded file) in workdir/data/raw;
//  - Writes cleaned FASTA file and selected consensus candidates in workdir/data/clean (filenames see config.hpp);
//
// Conventions:
//  - The function will try to modify and write data in place to reduce memory peaks; Top-K selector will keep K complete records in memory.
uint_t preprocessInputFasta(const std::string input_path, const std::string workdir, const int cons_n = 1000);


// ==============================================================
// alignConsensusSequence
//
// Description: Externally exposed utility function for performing multiple sequence alignment (MSA) on unaligned consensus sequences in `input_file`,
// writing the alignment results to `output_file`. This function is usually called after `preprocessInputFasta`, processing flow is:
//  1) Use `msa_cmd` template to construct command (template can contain {input} {output} {thread} placeholders);
//  2) Execute the command under the specified working directory `workdir` (via shell or cmd module), and wait for completion;
//  3) The function will log the running time and perform basic checks on the result file (existence and size).
//
// Parameters:
//  - input_file: Unaligned FASTA to align (FilePath)
//  - output_file: File path to write alignment results (FilePath)
//  - msa_cmd: Multiple sequence alignment command template string (e.g., "mafft --auto {input} > {output}")
//  - workdir: Current working directory for running the command (relative paths in command are based on this)
//  - threads: Number of threads allocated to MSA command (passed to {thread} in template), whether it takes effect depends on the MSA tool used
//
// Performance tips:
//  - MSA is usually CPU-intensive and memory-sensitive step, please adjust `threads` and MSA tool parameters according to the target machine;
//  - If MSA tool supports streaming interface, consider changing to streaming pipeline in future to reduce disk I/O.
//
// ==============================================================
void alignConsensusSequence(const FilePath& input_file, const FilePath& output_file,
                            const std::string& msa_cmd, int threads);


#endif //HALIGN4_PREPROCESS_H
