#include "align.h"

#include <stdexcept>
#include <cctype>
#include <limits>
#include <string>

// CIGAR encoding/decoding/conversion and sequence projection
// - CIGAR compressed format: unit=(len<<4)|op, encoding 0=M,1=I,2=D,3=N,4=S,5=H,6=P,7==,8=X
// - Key operations: string conversion, encode/decode, insertion/deletion detection, sequence projection alignment
// - Performance optimization: avoid O(N^2), pre-allocation, backward filling

namespace cigar
{
    // Operator and encoding conversion
    static inline uint32_t opCharToCode(const char op)
    {
        switch (op) {
        case 'M': return 0;
        case 'I': return 1;
        case 'D': return 2;
        case 'N': return 3;
        case 'S': return 4;
        case 'H': return 5;
        case 'P': return 6;
        case '=': return 7;
        case 'X': return 8;
        default:
            throw std::runtime_error(std::string("Unknown CIGAR op char: ") + op);
        }
    }

    static inline char opCodeToChar(const uint32_t code)
    {
        switch (code) {
        case 0: return 'M';
        case 1: return 'I';
        case 2: return 'D';
        case 3: return 'N';
        case 4: return 'S';
        case 5: return 'H';
        case 6: return 'P';
        case 7: return '=';
        case 8: return 'X';
        default:
            throw std::runtime_error("Unknown CIGAR op code: " + std::to_string(code));
        }
    }

    // Encode (operation,len) -> CigarUnit: (len<<4)|op
    CigarUnit cigarToInt(char operation, uint32_t len)
    {
        constexpr uint32_t kMaxLen = (1u << 28) - 1u;
        if (len == 0 || len > kMaxLen) {
            throw std::runtime_error("cigarToInt: invalid length=" + std::to_string(len));
        }
        const uint32_t op = opCharToCode(operation);
        return (len << 4) | (op & 0x0Fu);
    }

    // Decode CigarUnit -> (operation,len)
    void intToCigar(CigarUnit cigar, char& operation, uint32_t& len)
    {
        const uint32_t op = (cigar & 0x0Fu);
        len = (cigar >> 4);
        operation = opCodeToChar(op);
    }

    // Check if insertion operation 'I' exists (O(n) with short-circuit optimization)
    bool hasInsertion(const Cigar_t& cigar)
    {
        for (const CigarUnit cu : cigar) {
            char op_char;
            uint32_t len;
            intToCigar(cu, op_char, len);
            if (op_char == 'I') return true;
        }
        return false;
    }

    // Compress Cigar_t -> SAM CIGAR string (e.g. "10M5I3D")
    std::string cigarToString(const Cigar_t& cigar)
    {
        if (cigar.empty()) return "";

        std::string out;
        out.reserve(cigar.size() * 5);  // pre-allocate to avoid expansion

        for (const CigarUnit cu : cigar) {
            char op_char;
            uint32_t len = 0;
            intToCigar(cu, op_char, len);
            out.append(std::to_string(len));
            out.push_back(op_char);
        }
        return out;
    }

    // SAM CIGAR string -> compressed Cigar_t (e.g. "10M5I3D" parse to compressed vector)
    Cigar_t stringToCigar(const std::string& cigar_str)
    {
        Cigar_t result;
        if (cigar_str.empty() || cigar_str == "*") return result;

        result.reserve(cigar_str.size() / 2 + 1);

        uint64_t len_acc = 0;
        bool has_number = false;

        for (std::size_t i = 0; i < cigar_str.size(); ++i) {
            const unsigned char c = static_cast<unsigned char>(cigar_str[i]);

            if (std::isspace(c)) {
                continue;  // skip whitespace
            }

            if (std::isdigit(c)) {
                has_number = true;
                len_acc = len_acc * 10 + (c - '0');

                // overflow protection
                if (len_acc > ((1ull << 28) - 1ull)) {
                    throw std::runtime_error("stringToCigar: op length overflow in '" + cigar_str + "'");
                }
                continue;
            }

            // here: c is not digit nor whitespace, treat as op character
            if (!has_number || len_acc == 0) {
                throw std::runtime_error("stringToCigar: missing/invalid length before op in '" + cigar_str + "'");
            }

            const uint32_t len = static_cast<uint32_t>(len_acc);
            result.push_back(cigarToInt(static_cast<char>(c), len));

            // Reset state for next op
            len_acc = 0;
            has_number = false;
        }

        // string ends with number (missing op) is format error
        if (has_number) {
            throw std::runtime_error("stringToCigar: trailing number without op in '" + cigar_str + "'");
        }

        return result;
    }

    // Project query to ref coordinate system: insert gap '-' (backward fill to avoid O(N^2))
    void padQueryToRefByCigar(std::string& query, const Cigar_t& cigar)
    {
        if (cigar.empty()) {
            return;
        }

        // Count result length and consumed query length
        std::size_t out_len = 0;
        std::size_t consume_query = 0;

        for (const CigarUnit cu : cigar) {
            char op_char;
            uint32_t len = 0;
            intToCigar(cu, op_char, len);

            if (op_char == 'D' || op_char == 'N') {
                out_len += len;
            } else if (op_char == 'H' || op_char == 'P') {
                continue;
            } else {
                out_len += len;
                consume_query += len;
            }
        }

        assert(consume_query == query.size());

        // Backup original query, allocate output space
        std::string old = std::move(query);
        query.assign(out_len, '-');

        // Fill backward
        std::size_t w = out_len;
        std::size_t r = old.size();

        for (auto it = cigar.rbegin(); it != cigar.rend(); ++it) {
            char op_char;
            uint32_t len = 0;
            intToCigar(*it, op_char, len);

            if (op_char == 'D' || op_char == 'N') {
                for (uint32_t i = 0; i < len; ++i) {
                    query[--w] = '-';
                }
                continue;
            }

            if (op_char == 'H' || op_char == 'P') {
                continue;
            }

            for (uint32_t i = 0; i < len; ++i) {
                assert(r > 0);
                query[--w] = old[--r];
            }
        }

        assert(w == 0);
        assert(r == 0);
    }

    // Adjust query by CIGAR: delete I operation bases, add gap for D operation
    void delQueryToRefByCigar(std::string& query, const Cigar_t& cigar)
    {
        if (cigar.empty()) {
            return;
        }

        // Count result length and query consumed length
        std::size_t out_len = 0;
        std::size_t consume_query = 0;

        for (const CigarUnit cu : cigar) {
            char op_char;
            uint32_t len = 0;
            intToCigar(cu, op_char, len);

            if (op_char == 'I') {
                consume_query += len;
            } else if (op_char == 'H' || op_char == 'P') {
                continue;
            } else if (op_char == 'D' || op_char == 'N') {
                out_len += len;
            } else {
                out_len += len;
                consume_query += len;
            }
        }

        assert(consume_query == query.size());

        // Fast path: if no I/D/N operations, return directly
        if (out_len == query.size()) {
            bool has_indel = false;
            for (const CigarUnit cu : cigar) {
                char op_char;
                uint32_t len = 0;
                intToCigar(cu, op_char, len);
                if (op_char == 'I' || op_char == 'D' || op_char == 'N') {
                    has_indel = true;
                    break;
                }
            }
            if (!has_indel) return;
        }

        // Backup original query, allocate output space
        std::string old = std::move(query);
        query.assign(out_len, '-');

        // Fill backward
        std::size_t w = out_len;
        std::size_t r = old.size();

        for (auto it = cigar.rbegin(); it != cigar.rend(); ++it) {
            char op_char;
            uint32_t len = 0;
            intToCigar(*it, op_char, len);

            if (op_char == 'I') {
                // skip I operation characters
                for (uint32_t i = 0; i < len; ++i) {
                    assert(r > 0);
                    --r;
                }
                continue;
            }

            if (op_char == 'H' || op_char == 'P') {
                continue;
            }

            if (op_char == 'D' || op_char == 'N') {
                // D/N corresponds to gap (already initialized to '-')
                assert(w >= len);
                w -= len;
                continue;
            }

            // M/S/=/X: Copy characters
            for (uint32_t i = 0; i < len; ++i) {
                assert(r > 0);
                assert(w > 0);
                query[--w] = old[--r];
            }
        }

        assert(w == 0);
        assert(r == 0);
    }

    // Append CIGAR smartly and merge adjacent same-type operations
    void appendCigar(Cigar_t& result, const Cigar_t& cigar_to_add)
    {
        for (const CigarUnit cu : cigar_to_add) {
            char op_char;
            uint32_t len = 0;
            intToCigar(cu, op_char, len);

            if (len == 0) continue;

            if (result.empty()) {
                result.push_back(cu);
            } else {
                char last_op_char;
                uint32_t last_len = 0;
                intToCigar(result.back(), last_op_char, last_len);

                if (last_op_char == op_char) {
                    constexpr uint32_t kMaxLen = (1u << 28) - 1u;
                    if (static_cast<uint64_t>(last_len) + len > kMaxLen) {
                        throw std::runtime_error("appendCigar: merged length overflow");
                    }
                    result.back() = cigarToInt(op_char, last_len + len);
                } else {
                    result.push_back(cu);
                }
            }
        }
    }

    // Calculate reference sequence length consumed by CIGAR
    std::size_t getRefLength(const Cigar_t& cigar)
    {
        std::size_t total = 0;
        for (const CigarUnit cu : cigar) {
            char op_char;
            uint32_t len = 0;
            intToCigar(cu, op_char, len);

            // M/D/N/=/X consume ref
            if (op_char == 'M' || op_char == 'D' || op_char == 'N' ||
                op_char == '=' || op_char == 'X') {
                total += len;
            }
        }
        return total;
    }

    // Calculate query sequence length consumed by CIGAR
    std::size_t getQueryLength(const Cigar_t& cigar)
    {
        std::size_t total = 0;
        for (const CigarUnit cu : cigar) {
            char op_char;
            uint32_t len = 0;
            intToCigar(cu, op_char, len);

            // M/I/S/=/X consume query
            if (op_char == 'M' || op_char == 'I' || op_char == 'S' ||
                op_char == '=' || op_char == 'X') {
                total += len;
            }
        }
        return total;
    }

} // namespace cigar

