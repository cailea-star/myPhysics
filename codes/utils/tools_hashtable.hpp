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
#include <istream>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <type_traits>
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

public:
    /**
     * @brief  Construct a mixed-radix hash table over closed key intervals.
     * @math   stride_i = ∏_{j<i}(keymax_j - keymin_j + 1)
     * @output Empty hash table with configured bounds.
     * @note   Requires valid nonempty bounds and a packed-key range within uint64_t.
     */
    MyHashTable(const std::vector<int>& keymin_I1D_i, const std::vector<int>& keymax_I1D_i)
    : keymin_I1D_i(keymin_I1D_i), keymax_I1D_i(keymax_I1D_i), keystride_I1D_i(keymin_I1D_i.size(), 1) {
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
    const Value_T& read(const Key_T& key_I1D_i) const {
        // packedkey(key) → lookup.
        const std::uint64_t packedkey_I = to_packedkey(key_I1D_i);
        const auto entry_ = value_T1D_key.find(packedkey_I);

        // packedkey(key) ∈ keys.
        assert(entry_ != value_T1D_key.end());

        // M[packedkey(key)] → value.
        return entry_->second;
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
     * @output Empty hash table with unchanged bounds.
     * @note   Reserved buckets may be retained.
     */
    void clear() {
        // M → ∅.
        value_T1D_key.clear();
    }

    /**
     * @brief  Serialize trivially copyable hash-table entries to a binary stream.
     * @math   stream ← N_entry ⊕ {(packedkey,value)}
     * @output Binary entry sequence written to output_.
     * @note   Bounds are external metadata; the binary representation is platform-dependent.
     */
    void to_stream(std::ostream& output_) const {
        static_assert(std::is_trivially_copyable_v<Value_T>, "to_stream requires a trivially copyable Value_T");

        // M → N_entry ⊕ {(packedkey,value)}.
        const std::uint64_t Nentry_I = static_cast<std::uint64_t>(value_T1D_key.size());
        output_.write(reinterpret_cast<const char*>(&Nentry_I), sizeof(Nentry_I));
        for (const auto& [packedkey_I, value_T] : value_T1D_key) {
            output_.write(reinterpret_cast<const char*>(&packedkey_I), sizeof(packedkey_I));
            output_.write(reinterpret_cast<const char*>(&value_T), sizeof(value_T));
        }

        // stream = valid.
        if (!output_) {throw std::runtime_error("[ERROR]: [MyHashTable::to_stream] binary write failed");}
    }

    /**
     * @brief  Deserialize trivially copyable hash-table entries from a binary stream.
     * @math   N_entry ⊕ {(packedkey,value)} → M
     * @output Hash table replaced by the entries read from input_.
     * @note   Bounds must match the external metadata used to create the stream.
     */
    void from_stream(std::istream& input_) {
        static_assert(std::is_trivially_copyable_v<Value_T>, "from_stream requires a trivially copyable Value_T");

        // N_entry ← stream.
        std::uint64_t Nentry_I = 0;
        input_.read(reinterpret_cast<char*>(&Nentry_I), sizeof(Nentry_I));
        if (!input_) {throw std::runtime_error("[ERROR]: [MyHashTable::from_stream] binary read failed");}

        // M ← {(packedkey,value)}.
        const std::uint64_t maxpackedkey_I = to_packedkey(keymax_I1D_i);
        value_T1D_key.clear();
        value_T1D_key.reserve(static_cast<std::size_t>(Nentry_I));
        for (std::uint64_t i_I = 0; i_I < Nentry_I; ++i_I) {
            std::uint64_t packedkey_I = 0;
            Value_T value_T{};
            input_.read(reinterpret_cast<char*>(&packedkey_I), sizeof(packedkey_I));
            input_.read(reinterpret_cast<char*>(&value_T), sizeof(value_T));
            if (!input_) {throw std::runtime_error("[ERROR]: [MyHashTable::from_stream] binary read failed");}
            if (packedkey_I > maxpackedkey_I) {throw std::runtime_error("[ERROR]: [MyHashTable::from_stream] packed key exceeds configured bounds");}
            value_T1D_key[packedkey_I] = value_T;
        }
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
