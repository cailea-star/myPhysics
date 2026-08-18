/**
 * @file    tools_hashtable.hpp
 * @author  cailea
 * @date    2026-05-22
 * @brief   Sparse hash table for bounded integer-key tuples.
 */

#pragma once

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <unordered_map>
#include <vector>

/**
 * @brief Store sparse values indexed by bounded integer-key tuples.
 * @note  Mixed-radix packing is collision-free inside the configured closed intervals.
 */
template <typename Value_T = double>
class MyHashTable {
private:
    std::vector<int> keymin_I1D_i;
    std::vector<int> keymax_I1D_i;
    std::vector<std::uint64_t> keystride_I1D_i;
    std::unordered_map<std::uint64_t, Value_T> value_T1D_key;
    Value_T defaultvalue_T;

public:
    /**
     * @brief  Construct a mixed-radix hash table over closed key intervals.
     * @math   stride_i = ∏_{j<i}(keymax_j - keymin_j + 1)
     * @output Empty hash table with configured bounds and default value.
     * @note   Requires valid nonempty bounds and a packed-key range within uint64_t.
     */
    MyHashTable(const std::vector<int>& keymin_I1D_i, const std::vector<int>& keymax_I1D_i, const Value_T& defaultvalue_T = Value_T{})
    : keymin_I1D_i(keymin_I1D_i), keymax_I1D_i(keymax_I1D_i), keystride_I1D_i(keymin_I1D_i.size(), 1), defaultvalue_T(defaultvalue_T) {
        // dim(keymin) = dim(keymax) > 0.
        assert(keymin_I1D_i.size() == keymax_I1D_i.size());
        assert(!keymin_I1D_i.empty());

        // stride_i = ∏_{j<i}(keymax_j - keymin_j + 1).
        std::uint64_t capacity_I = 1;
        for (std::size_t i_I = 0; i_I < keymin_I1D_i.size(); ++i_I) {
            // keymin_i ≤ keymax_i.
            assert(keymin_I1D_i[i_I] <= keymax_I1D_i[i_I]);

            // width_i = keymax_i - keymin_i + 1.
            const std::uint64_t width_I = static_cast<std::uint64_t>(static_cast<long long>(keymax_I1D_i[i_I]) - static_cast<long long>(keymin_I1D_i[i_I]) + 1LL);

            // ∏_{j≤i} w_j ≤ uint64_t^{max}.
            assert(capacity_I <= std::numeric_limits<std::uint64_t>::max() / width_I);

            // stride_i ← ∏_{j<i} width_j; capacity ← ∏_{j≤i} width_j.
            keystride_I1D_i[i_I] = capacity_I;
            capacity_I *= width_I;
        }
    }

    /**
     * @brief  Test whether a bounded integer-key tuple is stored.
     * @math   contains(key) = [packedkey(key) ∈ keys]
     * @output True if the integer-key tuple has a stored value.
     * @note   Requires a key inside the configured intervals.
     */
    template <typename Key_T>
    bool contains(const Key_T& key_I1D_i) const {
        // packedkey(key) → lookup.
        return value_T1D_key.find(to_packedkey(key_I1D_i)) != value_T1D_key.end();
    }

    /**
     * @brief  Access a stored value by bounded integer-key tuple.
     * @math   value = M[packedkey(key)]
     * @output Const reference to the stored value.
     * @note   Requires a stored key inside the configured intervals.
     */
    template <typename Key_T>
    const Value_T& at(const Key_T& key_I1D_i) const {
        // packedkey(key) → lookup.
        const std::uint64_t packedkey_I = to_packedkey(key_I1D_i);
        const auto entry_ = value_T1D_key.find(packedkey_I);

        // packedkey(key) ∈ keys.
        assert(entry_ != value_T1D_key.end());

        // M[packedkey(key)] → value.
        return entry_->second;
    }

    /**
     * @brief  Read a stored value or the configured default value.
     * @math   value = M[packedkey(key)] if packedkey(key) ∈ keys, otherwise defaultvalue
     * @output Stored value or default value.
     * @note   Requires a key inside the configured intervals.
     */
    template <typename Key_T>
    Value_T read(const Key_T& key_I1D_i) const {
        // packedkey(key) → lookup.
        const std::uint64_t packedkey_I = to_packedkey(key_I1D_i);
        const auto entry_ = value_T1D_key.find(packedkey_I);

        // M[packedkey(key)] or defaultvalue → value.
        return (entry_ == value_T1D_key.end()) ? defaultvalue_T : entry_->second;
    }

    /**
     * @brief  Store a value by bounded integer-key tuple.
     * @math   M[packedkey(key)] ← value
     * @output Updated hash table.
     * @note   Requires a key inside the configured intervals; existing values are overwritten.
     */
    template <typename Key_T>
    void write(const Key_T& key_I1D_i, const Value_T& value_T) {
        // M[packedkey(key)] ← value.
        value_T1D_key[to_packedkey(key_I1D_i)] = value_T;
    }

    /**
     * @brief  Reserve storage for expected entries.
     * @math   N_reserved ≥ N
     * @output Hash table with reserved entry capacity.
     * @note   Existing values are unchanged.
     */
    void reserve(std::size_t count_I) {
        // N_reserved ← N.
        value_T1D_key.reserve(count_I);
    }

    /**
     * @brief  Remove all stored values.
     * @math   M → ∅
     * @output Empty hash table with unchanged bounds and default value.
     * @note   Reserved buckets may be retained.
     */
    void clear() {
        // M → ∅.
        value_T1D_key.clear();
    }

private:
    /**
     * @brief  Pack a bounded integer-key tuple using mixed-radix encoding.
     * @math   packedkey(key) = Σ_i(key_i - keymin_i)stride_i
     * @output Collision-free uint64_t packed key.
     * @note   Requires a matching key dimension and components inside the configured intervals.
     */
    template <typename Key_T>
    std::uint64_t to_packedkey(const Key_T& key_I1D_i) const {
        // dim(key) = dim(bounds).
        assert(key_I1D_i.size() == keymin_I1D_i.size());

        // packedkey(key) = Σ_i(key_i - keymin_i)stride_i.
        std::uint64_t packedkey_I = 0;
        for (std::size_t i_I = 0; i_I < key_I1D_i.size(); ++i_I) {
            // key_i ∈ [keymin_i, keymax_i].
            assert(key_I1D_i[i_I] >= keymin_I1D_i[i_I] && key_I1D_i[i_I] <= keymax_I1D_i[i_I]);

            // key_i - keymin_i → key offset.
            const std::uint64_t keyoffset_I = static_cast<std::uint64_t>(static_cast<long long>(key_I1D_i[i_I]) - static_cast<long long>(keymin_I1D_i[i_I]));
            packedkey_I += keyoffset_I * keystride_I1D_i[i_I];
        }

        // packedkey(key) → packed key.
        return packedkey_I;
    }
};
