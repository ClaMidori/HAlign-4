#include "consensus.h"

// ---------------- TopKLongestSelector Implementation Notes ----------------
// TopKLongestSelector responsibility:in one linear scan (streaming),select top K longest sequences,
// maintain stability at equal lengths (preserve earlier first).
//
// Design goals and motivation:
// - Streaming: no load all to memory; suitable N large (millions+).
// - Space-limited: O(K) extra; K << N (e.g. hundreds/thousands).
// - Time-efficient: O(log K) per sequence; O(N log K) total.
// - Stability: at equal length preserve earlier (aids reproducibility/debug).
//
// Data structures and comparison:
// - Min-heap: top stores worst of K (shortest or latest at same len),
//   new better candidates compare with top; if better replace and sift, keep size K.
// - Item: len (length), order (order counter for stability), rec (SeqRecord).
// - Comparison: lex (len, -order) or two functions:
//     * worseThan(a,b): if a worse than b (worse elem at top first), true;
//       def: shorter -> worse; if equal len, larger order (later) -> worse.
//     * betterThan(a,b): opposite, if cand better than top.
//
// Stability:
// - order_counter_ increments to record read order.
// - At equal len, earlier rank better (smaller order). So at same len, earlier reserved first,
//   ensuring stability (same order gives same result).
//
// Memory and performance:
// - Memory: K SeqRecords use O(K * avg_len). If avg_len large and K large,
//   high memory; if sensitive, indices/metadata or external storage.
// - Performance: O(log K) heap ops. Small K -> fast. Large K -> other strategies (partial sort/external).
//
// Concurrency/thread-safe:
// - Not thread-safe; TopKLongestSelector assumes single-thread context (main reading file).
// - For multi-thread (concurrent reads/calls), use locking or thread-local TopK+merge (recommended):
//     * each thread owns TopK (thread-local), after batch merge to global (cost O(T*K log K)).
//
// Alternatives/extensions:
// - For unstable TopK only, use std::priority_queue + simple comparator.
// - If K large+low memory, external merge or two-pass (sample threshold, then filter).
// - For other metrics (quality/composite), adjust Item comparator.
//
// Boundaries and exceptions:
// - If k_==0: selector returns immediately; takeSortedDesc returns empty.
// - If ultra-long seq causes overflow (unrealistic), handled as size_t.
// - Assumes seq_io::SeqRecord move semantics work.


TopKLongestSelector::TopKLongestSelector(std::size_t k)
        : k_(k)
    {
        heap_.reserve(k_);
    }

    // Reset to new K, clear state
    void TopKLongestSelector::reset(std::size_t k)
    {
        k_ = k;
        order_counter_ = 0; // order_counter_ records order for stability
        heap_.clear();
        heap_.reserve(k_);
    }

    // Current heap size
    std::size_t TopKLongestSelector::size() const
    {
        return heap_.size();
    }

    // Capacity K
    std::size_t TopKLongestSelector::capacity() const
    {
        return k_;
    }

    // Is empty
    bool TopKLongestSelector::empty() const
    {
        return heap_.empty();
    }

    // Compare: is a worse than b
    // Min-heap: worse elem at top (replace/pop).
    bool TopKLongestSelector::worseThan(const Item& a, const Item& b)
    {
        if (a.len != b.len) return a.len < b.len;      // shorter is worse
        return a.order > b.order;                      // same len: later (larger order) worse
    }

    // Is cand better than worst (replace top?)
    bool TopKLongestSelector::betterThan(const Item& cand, const Item& worst)
    {
        if (cand.len != worst.len) return cand.len > worst.len; // longer is better
        // same len prefer earlier (smaller order)
        return cand.order < worst.order;
    }

    // Sift up: swap with parent if worse
    // Min-heap sift up: maintain worst at top
    void TopKLongestSelector::siftUp(std::size_t idx)
    {
        while (idx > 0) {
            const std::size_t parent = (idx - 1) / 2;
            // if current worse, sift up
            if (worseThan(heap_[idx], heap_[parent])) {
                std::swap(heap_[idx], heap_[parent]);
                idx = parent;
            } else {
                break;
            }
        }
    }

    // Sift down: sink top to position
    // find worse child, swap if needed
    void TopKLongestSelector::siftDown(std::size_t idx)
    {
        const std::size_t n = heap_.size();
        while (true) {
            const std::size_t left = idx * 2 + 1;
            if (left >= n) break; // no children, stop

            const std::size_t right = left + 1;

            // select worse child
            std::size_t worst_child = left;
            if (right < n && worseThan(heap_[right], heap_[left])) {
                worst_child = right;
            }

            // if child worse, sift
            if (worseThan(heap_[worst_child], heap_[idx])) {
                std::swap(heap_[idx], heap_[worst_child]);
                idx = worst_child;
            } else {
                break;
            }
        }
    }

    // Consider new record
    // Logic:
    // - if k_==0, return (no save)
    // - if not full, push+sift
    // - else compare with top; if better replace+sift
    //
    // Key perf opt:
    // - check len/order first, then move
    // - avoid moving 99.99% records
    // - 1M records K=100: avoid ~999k moves
    void TopKLongestSelector::consider(seq_io::SeqRecord rec)
    {
        if (k_ == 0) return;

        const std::size_t cand_len = rec.seq.size();
        const std::uint64_t cand_order = order_counter_++;

        if ((double)rec.n_num / (double)rec.seq.size() > 0.01) {
            return;
        }

        if ((double)rec.n_num / (double)rec.seq.size() > 0.01) {
            return;
        }

        // not full: accept
        if (heap_.size() < k_) {
            Item item;
            item.len = cand_len;
            item.order = cand_order;
            item.rec = std::move(rec);
            heap_.push_back(std::move(item));
            siftUp(heap_.size() - 1);
            return;
        }

        // full: check if better
        // avoid wrong moves
        const Item& worst = heap_[0];

        // quick: better?
        bool is_better;
        if (cand_len != worst.len) {
            is_better = (cand_len > worst.len);  // longer is better
        } else {
            is_better = (cand_order < worst.order);  // equal len: earlier better
        }

        // only if better, move+replace
        if (is_better) {
            heap_[0].len = cand_len;
            heap_[0].order = cand_order;
            heap_[0].rec = std::move(rec);
            siftDown(0);
        }
        // else: destroyed, avoid waste
    }

    // Output by descended len (stable)
    // Return: vector of selected records
    std::vector<seq_io::SeqRecord> TopKLongestSelector::takeSortedDesc()
    {
        // avoid copy, move to sort
        std::vector<Item> items = std::move(heap_);
        heap_.clear();
        heap_.shrink_to_fit();
        heap_.reserve(k_);

        // custom: len desc; equal len asc order
        std::sort(items.begin(), items.end(),
                  [](const Item& a, const Item& b) {
                      if (a.len != b.len) return a.len > b.len;
                      return a.order < b.order;
                  });

        // move items to result
        std::vector<seq_io::SeqRecord> out;
        out.reserve(items.size());
        for (auto& it : items) {
            out.push_back(std::move(it.rec));
        }
        return out;
    }
