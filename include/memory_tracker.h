#ifndef HALIGN4_MEMORY_TRACKER_H
#define HALIGN4_MEMORY_TRACKER_H

#include <cstdint>
#include <string>
#include <fstream>
#include <unistd.h>
#include <spdlog/spdlog.h>

// Lightweight MemoryTracker: reads /proc/self/statm and reports RSS in bytes/MB.
// Note: Linux-specific (works on WSL). Returns 0 on failure.
class MemoryTracker {
public:
    // Return resident set size in bytes
    static std::uint64_t rss_bytes() {
        std::ifstream f("/proc/self/statm");
        if (!f.is_open()) return 0;
        long total_pages = 0;
        long rss_pages = 0;
        // statm: size resident shared text lib data dt
        if (!(f >> total_pages >> rss_pages)) return 0;
        const long page_size = sysconf(_SC_PAGESIZE);
        if (page_size <= 0) return 0;
        return static_cast<std::uint64_t>(rss_pages) * static_cast<std::uint64_t>(page_size);
    }

    // Return RSS in megabytes (floating, for nicer logging)
    static double rss_mb() {
        const std::uint64_t bytes = rss_bytes();
        return (bytes == 0) ? 0.0 : static_cast<double>(bytes) / (1024.0 * 1024.0);
    }

    // Helper: log current RSS with a message
    static void log_rss(const std::string& tag) {
        const double mb = rss_mb();
        if (mb <= 0.0) {
            spdlog::warn("{}: MemoryTracker: failed to read /proc/self/statm", tag);
        } else {
            spdlog::info("{}: RSS = {:.2f} MB", tag, mb);
        }
    }
};

#endif // HALIGN4_MEMORY_TRACKER_H
