/**
 * @file    tools_string.hpp
 * @author  cailea
 * @date    2026-01-29
 * @brief   String replacement and splitting utilities.
 */

#pragma once

#include <cstddef>
#include <sstream>
#include <string>
#include <vector>

/**
 * @brief  Replace all non-overlapping matches from left to right.
 * @math   text_Str ← ReplaceAll(text_Str, from_Str, to_Str)
 * @output Input string modified in place.
 * @note   An empty from_Str leaves text_Str unchanged.
 */
inline void replace_all(std::string& text_Str, const std::string& from_Str, const std::string& to_Str) {
    if (from_Str.empty()) {return;}

    // text_Str[position_I:] : from_Str → to_Str.
    std::size_t position_I = 0;
    while ((position_I = text_Str.find(from_Str, position_I)) != std::string::npos) {
        text_Str.replace(position_I, from_Str.length(), to_Str);
        position_I += to_Str.length();
    }
}

/**
 * @brief  Split a string with std::getline delimiter semantics.
 * @math   tokens = Split(text_Str, delimiter)
 * @output Delimiter-separated tokens.
 * @note   Preserves leading and interior empty tokens.
 * @note   Omits trailing empties; empty input yields none.
 */
inline std::vector<std::string> split(const std::string& text_Str, char delimiter) {
    // text_Str → token_0, token_1, ….
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream token_stream(text_Str);
    while (std::getline(token_stream, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}
