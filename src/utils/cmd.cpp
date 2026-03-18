#include "utils.h"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include <sys/wait.h> // WIFEXITED/WEXITSTATUS/WIFSIGNALED/WTERMSIG

#include "config.hpp"

// Do not use namespace detail: helpers visible only within this cpp (internal linkage)
// replaceAll: replace all occurrences of a substring in a string with another string.
// Note: if 'from' is empty, does nothing; this implementation modifies the original string and supports changing target lengths.
static void replaceAll(std::string& s, const std::string& from, const std::string& to)
{
    if (from.empty()) return;
    std::size_t pos = 0;
    while ((pos = s.find(from, pos)) != std::string::npos) {
        s.replace(pos, from.size(), to);
        pos += to.size();
    }
}

// containsToken: check if a string contains a specified substring, mainly for quickly determining if a template contains placeholders
// Lightweight tool, does not perform full template parsing, only existence check.
static bool containsToken(const std::string& s, const std::string& tok)
{
    return s.find(tok) != std::string::npos;
}

namespace cmd
{
    // buildCommand:
    // - cmd_template: user-provided command template, must contain at least {input} and {output} placeholders.
    // - input_path/output_path: will replace {input}/{output} in the template.
    // - thread: if template contains {thread} and thread >= 0, replaces with thread count; pass -1 to not replace {thread}.
    // - opt: BuildOptions controls quiet mode and whether to close stdin, etc.
    // Return: final command line string for system() (note: it's a shell command line, may include redirections).
    // Security note: this function does simple string replacement, no shell escaping of paths; if template comes from untrusted source,
    // use proper escaping or avoid shell execution (e.g., use execv family to pass args directly).
    std::string buildCommand(std::string cmd_template,
                             const std::string& input_path,
                             const std::string& output_path,
                             int thread,
                             const BuildOptions& opt)
    {
        // Must contain {input} and {output}
        if (!containsToken(cmd_template, "{input}")) {
            throw std::runtime_error("cmd template missing {input}");
        }
        if (!containsToken(cmd_template, "{output}")) {
            throw std::runtime_error("cmd template missing {output}");
        }

        // Placeholder replacement
        replaceAll(cmd_template, "{input}", input_path);
        replaceAll(cmd_template, "{output}", output_path);

        if (containsToken(cmd_template, "{thread}") && thread >= 0) {
            replaceAll(cmd_template, "{thread}", std::to_string(thread));
        }

        // If user requests quiet or close stdin, append redirection symbols after the command.
        // note: detect_stdout_redirect just checks for '>' char, simple heuristic, not fully reliable.
        if (!opt.quiet && !opt.close_stdin) {
            return cmd_template;
        }

        const bool has_stdout_redirect =
            opt.detect_stdout_redirect ? (cmd_template.find('>') != std::string::npos) : false;

        // Quiet strategy (Linux):
        // - If user template already redirects stdout: discard stderr only (2>/dev/null)
        // - Otherwise: discard stdout + stderr to stdout (> /dev/null 2>&1)
        if (opt.quiet) {
            if (has_stdout_redirect) {
                cmd_template.append(" 2>/dev/null");
            } else {
                cmd_template.append(" > /dev/null 2>&1");
            }
        }

        // Close stdin to prevent external command from waiting for interactive input
        if (opt.close_stdin) {
            cmd_template.append(" < /dev/null");
        }

        return cmd_template;
    }

    // runCommand:
    // - Directly uses std::system to execute shell command and returns command's exit code (not boolean)
    // - Return semantics:
    //     * -1: system() call itself failed (e.g., fork failed);
    //     * If WIFEXITED(status) is true, returns child's exit code (WEXITSTATUS(status));
    //     * Otherwise returns raw status (e.g., encoding for signal termination).
    // Note: std::system calls /bin/sh -c "command", has shell injection risk, command args should be properly escaped by caller.
    int runCommand(const std::string& command)
    {
        const int status = std::system(command.c_str());

        if (status == -1) {
            return status;
        }
        if (WIFEXITED(status)) {
            return WEXITSTATUS(status);
        }
        // If terminated by signal or other unrecognized return, treat as failure
        return status;
    }

    // testCommandTemplate:
    // - Creates a temporary directory WORKDIR_TMP under the given workdir, writes a small FASTA (tiny.fasta) as input,
    //   then constructs and runs the user-provided command template, checks if command succeeds (exit code 0) and output file is generated.
    // - Parameters:
    //     * cmd_template: user command template (will be replaced by buildCommand for {input}/{output})
    //     * workdir: base working directory (function creates temp files under workdir/WORKDIR_TMP)
    //     * thread: thread replacement value passed to buildCommand
    // - Return: bool, true means command executed successfully and produced output file; false means any step failed.
    // - Cleanup: function will try to delete temp directory at the end (even if command failed, will attempt cleanup), deletion failure only logs warning.
    // Safety and robustness notes:
    // - This function is for "self-check" if template can run normally, command execution has side effects (writes to local filesystem), run in trusted environment.
    // - For long-running or interactive commands, suggest adding timeout in template or using safer execution methods.
    bool testCommandTemplate(const std::string& cmd_template, const FilePath& workdir, int thread)
    {
        const FilePath temp_dir = workdir / WORKDIR_TMP;
        const FilePath in_path  = temp_dir / "tiny.fasta";
        const FilePath out_path = temp_dir / "aligned.fasta";

        bool ok = true;

        do {
            // Use file_io to ensure temp directory exists (will create or throw on error)
            try {
                file_io::ensureDirectoryExists(temp_dir, WORKDIR_TMP);
            } catch (const std::exception& e) {
                spdlog::error("Failed to create {}: {}", temp_dir.string(), e.what());
                ok = false; break;
            }

            // Write a very small FASTA
            {
                // Ensure parent directory of file to write exists
                file_io::ensureParentDirExists(in_path);
                std::ofstream ofs(in_path);
                if (!ofs) {
                    spdlog::error("Cannot open {} for writing.", in_path.string());
                    ok = false; break;
                }
                ofs <<
                    ">seq1\nACGTACGTGA\n"
                    ">seq2\nACGTTGCA\n"
                    ">seq3\nACGTACGA\n";
             }

             std::string cmd_line;
            try {
                // Ensure parent directory of output file exists (in case user template writes to subdirectory)
                file_io::ensureParentDirExists(out_path);
                cmd_line = buildCommand(cmd_template, in_path.string(), out_path.string(), thread);
            } catch (const std::exception& e) {
                spdlog::error("buildCommand failed: {}", e.what());
                ok = false; break;
            }

            spdlog::info("Running: {}", cmd_line);
            const int rc = runCommand(cmd_line);

            if (rc != 0) {
                spdlog::error("cmd failed");
                ok = false;
                break;
            }

            // Use file_io::requireExists to check if output was produced (will throw on failure)
            try {
                file_io::requireExists(out_path, "command output");
            } catch (const std::exception& e) {
                spdlog::error("Output file not found: {} ({})", out_path.string(), e.what());
                ok = false;
                break;
            }

            spdlog::info("cmd finished successfully, output: {}", out_path.string());

        } while (false);

        // Delete temp directory, using file_io::removeAll
        try {
            file_io::removeAll(temp_dir);
        } catch (const std::exception& e) {
            spdlog::warn("Failed to remove {}: {}", temp_dir.string(), e.what());
        }

        return ok;
    }

} // namespace cmd
