#ifndef HALIGN4_SEED_H
#define HALIGN4_SEED_H

#include <cstdint>
#include <type_traits>
#include <string_view>
#include <stdexcept>
#include <utility>
#include <vector>
#include <string>
#include <cstddef>
#include <algorithm>
#include <limits>
#include <cmath>
#include "hash.h"
#include "anchor.h"

// ================================================================
// Abstract seed interface (not bound to specific minimizer/syncmer/strobemer implementation)
//
// Goals:
// - You may support multiple seeds in the future (minimizer and various variants, such as syncmer/strobemer)
// - At the same time, want to choose different storage at different stages: store only hash or hash+position
// - Here only provide the "most abstract" interface and tools (traits/comparators), no specific implementation.
// ================================================================
namespace seed
{
    // ------------------------------------------------------------
    // SeedKind: Recommended to be shared at the "container/batch" level (save memory)
    // ------------------------------------------------------------
    enum class SeedKind : std::uint8_t
    {
        minimizer = 0,
        syncmer   = 1,
        strobemer = 2,
    };

    inline constexpr const char* seedKindToString(SeedKind k) noexcept
    {
        switch (k) {
        case SeedKind::minimizer: return "minimizer";
        case SeedKind::syncmer:   return "syncmer";
        case SeedKind::strobemer: return "strobemer";
        default:                  return "unknown";
        }
    }


    // ------------------------------------------------------------
    // CRTP base class: hash+position seed (hit)
    // ------------------------------------------------------------
    // Derived class needs to provide:
    //   hash_t hash() const noexcept;     // hash perspective
    //   uint32_t pos() const noexcept;      // position (0-based)
    //   uint32_t rid() const noexcept;      // sequence id (multi-sequence scenario)
    //   bool strand() const noexcept;       // direction
    //   uint32_t span() const noexcept;     // coverage range (minimizer generally equals k; strobemer can be larger)
    template <typename Derived>
    struct SeedHitBase
    {
        constexpr hash_t hash() const noexcept { return static_cast<const Derived&>(*this).hash(); }
        constexpr std::uint32_t pos() const noexcept { return static_cast<const Derived&>(*this).pos(); }
        constexpr std::uint32_t rid() const noexcept { return static_cast<const Derived&>(*this).rid(); }
        constexpr bool strand() const noexcept { return static_cast<const Derived&>(*this).strand(); }
        constexpr std::uint32_t span() const noexcept { return static_cast<const Derived&>(*this).span(); }

        // A general sorting: hash first, then rid/pos/strand
        friend constexpr bool operator<(const SeedHitBase& a, const SeedHitBase& b) noexcept
        {
            if (a.hash() != b.hash()) return a.hash() < b.hash();
            if (a.rid() != b.rid()) return a.rid() < b.rid();
            if (a.pos() != b.pos()) return a.pos() < b.pos();
            return a.strand() < b.strand();
        }

        friend constexpr bool operator==(const SeedHitBase& a, const SeedHitBase& b) noexcept
        {
            return a.hash() == b.hash() && a.rid() == b.rid() && a.pos() == b.pos() && a.strand() == b.strand() && a.span() == b.span();
        }
    };

    // ------------------------------------------------------------
    // Unified access interface (free-functions)
    // ------------------------------------------------------------

    // Default trait: no position information (you can still manually specialize to override)
    template <typename T>
    struct has_position : std::false_type {};

    // Auto deduction: if T inherits from SeedHitBase<T>, then it is considered to "have position".
    // This way, every time you add a Hit type later, as long as it inherits SeedHitBase, it automatically takes effect, no need to write specialization again.
    template <typename T>
    struct has_position_auto : std::is_base_of<SeedHitBase<T>, T> {};

    template <typename T>
    inline constexpr bool has_position_v = has_position<T>::value || has_position_auto<T>::value;

    // hash-only: requires type to provide hash()
    template <typename SeedT>
    inline constexpr hash_t hash_value(const SeedT& s) noexcept
    {
        return s.hash();
    }

    // hit: requires type to provide hash()/pos()/rid()/strand()/span()
    template <typename HitT>
    inline constexpr std::uint32_t get_pos(const HitT& h) noexcept { return h.pos(); }
    template <typename HitT>
    inline constexpr std::uint32_t get_rid(const HitT& h) noexcept { return h.rid(); }
    template <typename HitT>
    inline constexpr bool get_strand(const HitT& h) noexcept { return h.strand(); }
    template <typename HitT>
    inline constexpr std::uint32_t get_span(const HitT& h) noexcept { return h.span(); }

    // Comparator that only looks at hash (applicable to sort+unique / Jaccard / containment etc.)
    struct HashOnlyLess
    {
        template <typename M>
        constexpr bool operator()(const M& a, const M& b) const noexcept
        {
            return hash_value(a) < hash_value(b);
        }
    };

    struct HashOnlyEqual
    {
        template <typename M>
        constexpr bool operator()(const M& a, const M& b) const noexcept
        {
            return hash_value(a) == hash_value(b);
        }
    };

} // namespace seed

// =============================================================
namespace minimizer
{

    // =============================================================
    // nt4_table
    // -------------------------------------------------------------
    // Map input characters to 2-bit encoding (0..3), other characters are invalid(4).
    // - 'A'/'a' -> 0
    // - 'C'/'c' -> 1
    // - 'G'/'g' -> 2
    // - 'T'/'t' -> 3
    // - 'U'/'u' -> 3   (RNA U treated as T)
    // - Others (including 'N'/'n', '-', etc.) -> 4
    //
    // Design purpose:
    // - Allow high-performance implementation (rolling k-mer/minimap2 style) to O(1) lookup table, avoid switch branches.
    // - Table placed in header file, reused by multiple .cpp, no need to repeat 256 item initialization.
    //
    // ASCII values: A=65, C=67, G=71, T=84, U=85, a=97, c=99, g=103, t=116, u=117
    // =============================================================
    inline constexpr std::uint8_t nt4_table[256] = {
        // 0-15
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        // 16-31
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        // 32-47
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        // 48-63
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        // 64-79: A=65->0, C=67->1, G=71->2
        4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,
        // 80-95: T=84->3, U=85->3
        4,4,4,4,3,3,4,4,4,4,4,4,4,4,4,4,
        // 96-111: a=97->0, c=99->1, g=103->2
        4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,
        // 112-127: t=116->3, u=117->3
        4,4,4,4,3,3,4,4,4,4,4,4,4,4,4,4,
        // 128-143
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        // 144-159
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        // 160-175
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        // 176-191
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        // 192-207
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        // 208-223
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        // 224-239
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        // 240-255
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4
    };
    // ------------------------- minimizer hit (with position) -------------------------
    // Used for chaining/positioning: need to know which sequence (rid), what position (pos), direction (strand), and coverage length (span).
    //
    // Memory optimization: Use minimap2's mm128_t packing method, store all information in only 16 bytes:
    //   x = (hash << 8) | span
    //     - bit[0..7]   : span (8 bit)
    //     - bit[8..63]  : hash (56 bit)
    //   y = (rid_with_strand << 32) | pos
    //     - bit[0..31]  : pos (32 bit)
    //     - bit[32..62] : rid (31 bit)
    //     - bit[63]     : strand (1 bit)
    //
    // Notes:
    // - 56bit hash is sufficient for most scenarios; if you want to keep 64bit hash, need to change packing layout.
    // - span uses 8bit, suitable for k<=255; if larger k is needed later, layout can be adjusted.
    struct MinimizerHit : public seed::SeedHitBase<MinimizerHit>
    {
        hash_t x{0};
        hash_t y{0};

        constexpr MinimizerHit() = default;
        constexpr MinimizerHit(hash_t x_, hash_t y_) : x(x_), y(y_) {}

        // pack/unpack x
        static constexpr hash_t pack_x(hash_t hash56, std::uint8_t span) noexcept
        {
            return (hash56 << 8) | static_cast<hash_t>(span);
        }
        static constexpr std::uint8_t span_from_x(hash_t x) noexcept
        {
            return static_cast<std::uint8_t>(x & 0xffULL);
        }
        static constexpr hash_t hash_from_x(hash_t x) noexcept
        {
            return (x >> 8);
        }

        // pack/unpack y
        static constexpr hash_t pack_y(std::uint32_t pos, std::uint32_t rid, bool strand) noexcept
        {
            const std::uint32_t rid_with_strand = (rid & 0x7fffffffU) | (strand ? 0x80000000U : 0U);
            return (static_cast<hash_t>(rid_with_strand) << 32) | static_cast<hash_t>(pos);
        }
        static constexpr std::uint32_t pos_from_y(hash_t y) noexcept
        {
            return static_cast<std::uint32_t>(y & 0xffffffffULL);
        }
        static constexpr std::uint32_t rid_with_strand_from_y(hash_t y) noexcept
        {
            return static_cast<std::uint32_t>((y >> 32) & 0xffffffffULL);
        }
        static constexpr std::uint32_t rid_from_y(hash_t y) noexcept
        {
            return rid_with_strand_from_y(y) & 0x7fffffffU;
        }
        static constexpr bool strand_from_y(hash_t y) noexcept
        {
            return (rid_with_strand_from_y(y) & 0x80000000U) != 0U;
        }

        // Convenient constructor: input "complete semantic fields", internal automatic packing
        constexpr MinimizerHit(hash_t hash56, std::uint32_t pos, std::uint32_t rid, bool strand, std::uint8_t span) noexcept
            : x(pack_x(hash56, span)), y(pack_y(pos, rid, strand))
        {
        }

        // SeedHitBase required API (still provide hash/pos/rid/strand/span)
        constexpr hash_t hash() const noexcept { return hash_from_x(x); }
        constexpr std::uint32_t pos() const noexcept { return pos_from_y(y); }
        constexpr std::uint32_t rid() const noexcept { return rid_from_y(y); }
        constexpr bool strand() const noexcept { return strand_from_y(y); }
        constexpr std::uint32_t span() const noexcept { return span_from_x(x); }
    };

    static_assert(sizeof(MinimizerHit) == 16, "MinimizerHit should be 16 bytes when packed as (x,y)");
    using MinimizerHits = std::vector<MinimizerHit>;

    // =====================================================================
    // extractMinimizer
    // ---------------------------------------------------------------------
    // Extract minimizer hit list (hash+position) from an input sequence.
    //
    // Parameters:
    //   seq       : input sequence
    //   k         : k-mer size
    //   w         : window size (in k-mer units)
    //
    // Returns:
    //   minimizer hit list (in scan order).
    //
    // Note: Only declaration provided here; implementation should be placed in corresponding .cpp file.
    // =====================================================================
    MinimizerHits extractMinimizer(const std::string& seq,
                                   std::size_t k,
                                   std::size_t w,
                                   bool non_canonical);

    // =====================================================================
    // collect_anchors - collect anchors (refer to minimap2 implementation)
    // ---------------------------------------------------------------------
    // Returns:
    //   anchor::Anchors - anchor list (unsorted)
    // =====================================================================
    anchor::Anchors collect_anchors(const MinimizerHits& ref_hits, const MinimizerHits& qry_hits, anchor::SeedFilterParams params = anchor::default_mm2_params());




} // namespace minimizer


#endif //HALIGN4_SEED_H
