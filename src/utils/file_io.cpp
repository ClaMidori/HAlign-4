#include "../include/utils.h"

#include <system_error>
#include <iostream>

namespace file_io {

    // Format filesystem error information into a readable string (with path and error message)
    std::string formatFsError(std::string_view msg,
                                     const FilePath& p,
                                     const std::error_code& ec) {
        std::string s;
        s.reserve(msg.size() + p.string().size() + 64);
        s.append(msg);
        s.append(": ");
        s.append(p.string());
        if (ec) {
            s.append(" (");
            s.append(ec.message());
            s.append(")");
        }
        return s;
    }

    // Ensure path exists, otherwise throw exception (for input validation)
    void requireExists(const FilePath& p, std::string_view what) {
        std::error_code ec;
        const bool ok = fs::exists(p, ec);
        if (ec || !ok) {
            throw std::runtime_error(formatFsError(std::string(what) + " does not exist", p, ec));
        }
    }

    // Ensure path is a regular file, otherwise throw exception
    void requireRegularFile(const FilePath& p, std::string_view what) {
        std::error_code ec;
        const bool ok = fs::is_regular_file(p, ec);
        if (ec || !ok) {
            throw std::runtime_error(formatFsError(std::string(what) + " is not a regular file", p, ec));
        }
    }

    // Ensure path is a directory, otherwise throw exception
    void requireDirectory(const FilePath& p, std::string_view what) {
        std::error_code ec;
        const bool ok = fs::is_directory(p, ec);
        if (ec || !ok) {
            throw std::runtime_error(formatFsError(std::string(what) + " is not a directory", p, ec));
        }
    }

    // Ensure directory exists; if not, create it
    void ensureDirectoryExists(const FilePath& p, std::string_view what) {
        std::error_code ec;
        if (fs::exists(p, ec)) {
            if (ec) {
                throw std::runtime_error(formatFsError(std::string(what) + " status check failed", p, ec));
            }
            requireDirectory(p, what);
            return;
        }
        if (ec) {
            throw std::runtime_error(formatFsError(std::string(what) + " status check failed", p, ec));
        }

        fs::create_directories(p, ec);
        if (ec) {
            throw std::runtime_error(formatFsError(std::string("failed to create ") + std::string(what), p, ec));
        }
    }

    // Check if a path is empty (file has zero bytes)
    bool isEmpty(const FilePath& p) {
        std::error_code ec;
        const bool empty = fs::is_empty(p, ec);
        if (ec) {
            throw std::runtime_error(formatFsError("failed to check emptiness", p, ec));
        }
        return empty;
    }

    // Prepare an empty directory: create it if missing; throw if must_be_empty and directory is not empty
    void prepareEmptydir(const FilePath& workdir, bool must_be_empty) {
        if (workdir.empty()) {
            throw std::runtime_error("workdir is empty");
        }

        ensureDirectoryExists(workdir, "workdir");

        if (must_be_empty && !isEmpty(workdir)) {
            throw std::runtime_error("workdir must be empty: " + workdir.string());
        }
    }

    // Ensure the output file's parent directory exists (call before writing)
    void ensureParentDirExists(const FilePath& out_file) {
        if (out_file.empty()) return;
        auto parent = out_file.parent_path();
        if (parent.empty()) return;
        ensureDirectoryExists(parent, "output parent dir");
    }

    // Determine whether a given path is a URL (simple check for scheme:// or // prefix)
    bool isUrl(const FilePath& p) {
        const std::string s = p.string();
        if (s.size() >= 2 && s[0] == '/' && s[1] == '/') return true;
        static const std::regex scheme_re(R"(^[a-zA-Z][a-zA-Z0-9+.\-]*://)");
        return std::regex_search(s, scheme_re);
    }

    // ------------------------------------------------------------------
    // Function: copyFile
    // Purpose: Copy a file, with support for cross-filesystem copy
    //
    // Implementation notes:
    // 1. Prefer std::filesystem::copy_file (efficient, preserves metadata)
    // 2. Fall back to stream copy on cross-device errors (better compatibility)
    // 3. Special case: if source and destination refer to the same file, return early
    //
    // Parameters:
    //   - src: source file path (must be a regular file)
    //   - dst: destination file path (will overwrite existing file)
    //
    // Exceptions:
    //   - throws if source doesn't exist or isn't a regular file
    //   - throws on copy failures (with detailed error info)
    // ------------------------------------------------------------------
    void copyFile(const FilePath& src, const FilePath& dst) {
        requireRegularFile(src, "source file");
        ensureParentDirExists(dst);

        // Edge case: source and destination are the same file; no need to copy
        // Note: use std::filesystem::equivalent to check if two paths refer to the same file
        //       This is more reliable than string comparison and handles symlinks/relative paths.
        std::error_code equiv_ec;
        if (fs::equivalent(src, dst, equiv_ec)) {
            // If both paths refer to the same file, return early (no copy needed)
            return;
        }
        // Note: if equivalent fails (e.g. destination doesn't exist), equiv_ec will be set,
        //       but we ignore it because copy_file will handle that case.

        // Try std::filesystem::copy_file first (fast)
        std::error_code ec;
        if (fs::copy_file(src, dst, fs::copy_options::overwrite_existing, ec)) {
            // copy succeeded
            return;
        }

        // copy_file returned false; check why
        if (!ec) {
            // No error code but returned false, this is unexpected.
            // Possible reasons:
            // 1. Source and destination are the same (but equivalent check should have handled this)
            // 2. Destination file already exists and is identical (some implementations return false)
            // 3. Other unknown reasons
            //
            // In this case, we check if the destination file exists and is readable; if so, we treat it as success.
            std::error_code exists_ec;
            if (fs::exists(dst, exists_ec) && !exists_ec && fs::is_regular_file(dst, exists_ec) && !exists_ec) {
                // Destination exists and is a regular file; treat as success (likely already present and identical)
                return;
            }
            // Otherwise throw (unknown reason)
            throw std::runtime_error(formatFsError("failed to copy file (unknown reason)", dst, exists_ec));
        }

        // There is an error code: check if it is a cross-device error and fall back to stream copy
        if (ec == std::make_error_code(std::errc::cross_device_link)) {
            // Cross-filesystem copy: use stream copy
            std::ifstream in(src, std::ios::binary);
            if (!in) {
                throw std::runtime_error(formatFsError("failed to open source for reading", src, std::make_error_code(std::errc::io_error)));
            }
            std::ofstream out(dst, std::ios::binary | std::ios::trunc);
            if (!out) {
                throw std::runtime_error(formatFsError("failed to open destination for writing", dst, std::make_error_code(std::errc::io_error)));
            }
            out << in.rdbuf();
            if (!out) {
                throw std::runtime_error(formatFsError("failed to write file", dst, std::make_error_code(std::errc::io_error)));
            }
            return;
        }

        // Other errors: throw
        throw std::runtime_error(formatFsError("failed to copy file", dst, ec));
    }

    // Downloading remote files to local machine: use libcurl first (if available), otherwise fall back to calling curl/wget command line.
    void downloadFile(const std::string& url, const FilePath& dst) {
        if (url.empty()) {
            throw std::runtime_error("download url is empty");
        }
        ensureParentDirExists(dst);

    #if defined(HALIGN4_HAVE_LIBCURL)
        auto write_cb = [](char* ptr, size_t size, size_t nmemb, void* userdata) -> size_t {
            auto out = static_cast<std::ofstream*>(userdata);
            out->write(ptr, static_cast<std::streamsize>(size * nmemb));
            return out->good() ? size * nmemb : 0;
        };

        std::ofstream out(dst, std::ios::binary | std::ios::trunc);
        if (!out) {
            throw std::runtime_error("failed to open destination for writing: " + dst.string());
        }

        CURL* curl = curl_easy_init();
        if (!curl) {
            throw std::runtime_error("failed to initialize libcurl");
        }
        CURLcode res;
        curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
        curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
        curl_easy_setopt(curl, CURLOPT_FAILONERROR, 1L);
        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, +write_cb);
        curl_easy_setopt(curl, CURLOPT_WRITEDATA, &out);
        res = curl_easy_perform(curl);
        curl_easy_cleanup(curl);
        if (res != CURLE_OK) {
            out.close();
            std::error_code rem_ec;
            fs::remove(dst, rem_ec);
            throw std::runtime_error(std::string("libcurl error: ") + curl_easy_strerror(res));
        }
        out.close();
    #else
        std::string cmd = "curl -L -sSf -o \"" + dst.string() + "\" \"" + url + "\"";
        int rc = std::system(cmd.c_str());
        if (rc == 0) return;
        cmd = "wget -q -O \"" + dst.string() + "\" \"" + url + "\"";
        rc = std::system(cmd.c_str());
        if (rc == 0) return;

        std::error_code rem_ec;
        fs::remove(dst, rem_ec);
        throw std::runtime_error("failed to download file using curl/wget (rc=" + std::to_string(rc) + ")");
    #endif
    }

    // Download if srcOrUrl is a URL, otherwise copy a local file to destination
    void fetchFile(const FilePath& srcOrUrl, const FilePath& dst) {
        if (isUrl(srcOrUrl)) {
            downloadFile(srcOrUrl.string(), dst);
        } else {
            copyFile(srcOrUrl, dst);
        }
    }

    // Recursively remove a path (file or directory); throw on failure
    void removeAll(const FilePath& p) {
        std::error_code ec;
        fs::remove_all(p, ec);
        if (ec) {
            throw std::runtime_error(formatFsError("remove_all failed", p, ec));
        }
    }

    // Read an entire file into a std::string (for small files)
    std::string readFileToString(const FilePath& p) {
        requireRegularFile(p, "input file");

        std::ifstream in(p, std::ios::binary);
        if (!in) {
            throw std::runtime_error("failed to open file for reading: " + p.string());
        }

        std::string content;
        in.seekg(0, std::ios::end);
        content.resize(static_cast<std::size_t>(in.tellg()));
        in.seekg(0, std::ios::beg);
        in.read(&content[0], static_cast<std::streamsize>(content.size()));
        if (!in) {
            throw std::runtime_error("failed to read file: " + p.string());
        }
        return content;
    }

} // namespace file_io

