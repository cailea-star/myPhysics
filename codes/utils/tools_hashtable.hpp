/**
 * @file    tools_hashtable.hpp
 * @author  cailea
 * @date    2026-05-22
 * @brief   Sparse hash table for bounded integer-key tuples.
 */

#pragma once

#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <istream>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>

/**
 * @brief Store sparse values indexed by bounded integer-key tuples.
 * @note  Mixed-radix packing is collision-free inside the configured closed intervals.
 */
template <typename Value_T, std::size_t Rank_I>
class PackedHashTable {
private:
    std::array<int, Rank_I> keymin_I1D_i;
    std::array<int, Rank_I> keymax_I1D_i;
    std::array<std::uint64_t, Rank_I> keystride_I1D_i;
    std::unordered_map<std::uint64_t, Value_T> value_T1D_key;

public:
    /**
     * @brief  Construct a mixed-radix hash table over closed key intervals.
     * @math   stride_i = ∏_{j<i}(keymax_j - keymin_j + 1)
     * @output Empty hash table with configured bounds.
     * @note   Requires nonempty bounds; packed range must fit uint64_t.
     */
    PackedHashTable(const std::array<int, Rank_I>& keymin_I1D_i_, const std::array<int, Rank_I>& keymax_I1D_i_)
    : keymin_I1D_i(keymin_I1D_i_), keymax_I1D_i(keymax_I1D_i_), keystride_I1D_i{} {
        static_assert(Rank_I > 0, "PackedHashTable requires positive rank.");

        // {width_i} → ({stride_i},capacity).
        std::uint64_t capacity_I = 1;
        for (std::size_t i_I = 0; i_I < Rank_I; ++i_I) {
            assert(keymin_I1D_i_[i_I] <= keymax_I1D_i_[i_I]);

            const std::uint64_t width_I = static_cast<std::uint64_t>(static_cast<long long>(keymax_I1D_i_[i_I]) - static_cast<long long>(keymin_I1D_i_[i_I]) + 1LL);

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
    bool contains(const std::array<int, Rank_I>& key_I1D_i) const {
        return value_T1D_key.find(to_packedkey(key_I1D_i)) != value_T1D_key.end();
    }

    /**
     * @brief  Access a stored value by bounded integer-key tuple.
     * @math   value = M[packedkey(key)]
     * @output Const reference to the stored value.
     * @note   Requires a stored key inside the configured intervals.
     */
    const Value_T& read(const std::array<int, Rank_I>& key_I1D_i) const {
        // key → packedkey → entry → value.
        const std::uint64_t packedkey_I = to_packedkey(key_I1D_i);
        const auto entry_ = value_T1D_key.find(packedkey_I);

        assert(entry_ != value_T1D_key.end());

        return entry_->second;
    }

    /**
     * @brief  Store a value by bounded integer-key tuple.
     * @math   M[packedkey(key)] ← value
     * @output Updated hash table.
     * @note   Requires an in-range key; writes overwrite values.
     */
    void write(const std::array<int, Rank_I>& key_I1D_i, const Value_T& value_T) {
        value_T1D_key[to_packedkey(key_I1D_i)] = value_T;
    }

    /**
     * @brief  Reserve storage for expected entries.
     * @math   N_reserved ≥ N
     * @output Hash table with reserved entry capacity.
     * @note   Existing values are unchanged.
     */
    void reserve(std::size_t count_I) {
        value_T1D_key.reserve(count_I);
    }

    /**
     * @brief  Remove all stored values.
     * @math   M → ∅
     * @output Empty hash table with unchanged bounds.
     * @note   Reserved buckets may be retained.
     */
    void clear() {
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

        if (!output_) {throw std::runtime_error("[ERROR]: [PackedHashTable::to_stream] binary write failed");}
    }

    /**
     * @brief  Deserialize trivially copyable hash-table entries from a binary stream.
     * @math   N_entry ⊕ {(packedkey,value)} → M
     * @output Hash table replaced by the entries read from input_.
     * @note   Caller must reconstruct identical key bounds.
     */
    void from_stream(std::istream& input_) {
        static_assert(std::is_trivially_copyable_v<Value_T>, "from_stream requires a trivially copyable Value_T");

        std::uint64_t Nentry_I = 0;
        input_.read(reinterpret_cast<char*>(&Nentry_I), sizeof(Nentry_I));
        if (!input_) {throw std::runtime_error("[ERROR]: [PackedHashTable::from_stream] binary read failed");}

        // M ← {(packedkey,value)}.
        const std::uint64_t maxpackedkey_I = to_packedkey(keymax_I1D_i);
        value_T1D_key.clear();
        value_T1D_key.reserve(static_cast<std::size_t>(Nentry_I));
        for (std::uint64_t i_I = 0; i_I < Nentry_I; ++i_I) {
            std::uint64_t packedkey_I = 0;
            Value_T value_T{};
            input_.read(reinterpret_cast<char*>(&packedkey_I), sizeof(packedkey_I));
            input_.read(reinterpret_cast<char*>(&value_T), sizeof(value_T));
            if (!input_) {throw std::runtime_error("[ERROR]: [PackedHashTable::from_stream] binary read failed");}
            if (packedkey_I > maxpackedkey_I) {throw std::runtime_error("[ERROR]: [PackedHashTable::from_stream] packed key exceeds configured bounds");}
            value_T1D_key[packedkey_I] = value_T;
        }
    }

private:
    /**
     * @brief  Pack a bounded integer-key tuple using mixed-radix encoding.
     * @math   packedkey(key) = Σ_i(key_i - keymin_i)stride_i
     * @output Collision-free uint64_t packed key.
     * @note   Requires matching dimensions and in-range components.
     */
    std::uint64_t to_packedkey(const std::array<int, Rank_I>& key_I1D_i) const {
        // {key_i} → Σ_i(key_i-key_i^min)stride_i.
        std::uint64_t packedkey_I = 0;
        for (std::size_t i_I = 0; i_I < Rank_I; ++i_I) {
            assert(key_I1D_i[i_I] >= keymin_I1D_i[i_I] && key_I1D_i[i_I] <= keymax_I1D_i[i_I]);

            const std::uint64_t keyoffset_I = static_cast<std::uint64_t>(static_cast<long long>(key_I1D_i[i_I]) - static_cast<long long>(keymin_I1D_i[i_I]));
            packedkey_I += keyoffset_I * keystride_I1D_i[i_I];
        }

        return packedkey_I;
    }
};
