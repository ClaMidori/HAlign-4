#include "hash.h"


// The code2bit is fixed to "big-endian" input (consistent across platforms).
// Perform a bswap on the little-endian machine to make the memory byte order equivalent to the original big-endian bytes.
static inline std::uint64_t to_be64(std::uint64_t x) noexcept {
#if defined(__BYTE_ORDER__) && defined(__ORDER_LITTLE_ENDIAN__) && (__BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__)
    return __builtin_bswap64(x);
#else
    return x;
#endif
}

hash_t getHash(const char* seq, int length, std::uint32_t seed)
{
    // The seed type of XXH3_64bits_withSeed is XXH64_hash_t (64-bit), which can be directly extended.
    return static_cast<hash_t>(
        XXH3_64bits_withSeed(seq, static_cast<size_t>(length), static_cast<XXH64_hash_t>(seed))
    );
}

hash_t getHash2bit(std::uint64_t code2bit, std::uint32_t seed)
{
    const std::uint64_t be = to_be64(code2bit);
    return static_cast<hash_t>(
        XXH3_64bits_withSeed(&be, sizeof(be), static_cast<XXH64_hash_t>(seed))
    );
}
