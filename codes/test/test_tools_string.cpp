
#include <string>
#include <vector>
#include <iomanip>
#include <iostream>
#include "tools_string.hpp"

std::string head_line = " A, Z, N [Name] ";
std::string data_line = "  1, 1, 0 [H]    ";

int main() {
    std::cout << "Original head line: " << head_line << "\n";
    std::cout << "Original data line: " << data_line << "\n";

    // head_line → normalized comma-separated text.
    replace_all(head_line, " ", ""); // Remove spaces.
    replace_all(head_line, "\t", ""); // Remove tabs.
    replace_all(head_line, "[", ","); // [ → ,.
    replace_all(head_line, "]", ","); // ] → ,.
    // data_line → normalized comma-separated text.
    replace_all(data_line, " ", ""); // Remove spaces.
    replace_all(data_line, "\t", ""); // Remove tabs.
    replace_all(data_line, "[", ","); // [ → ,.
    replace_all(data_line, "]", ","); // ] → ,.

    std::cout << "Modified head line: " << head_line << "\n";
    std::cout << "Modified data line: " << data_line << "\n";

    // text → tokens.
    std::vector<std::string> head_tokens = split(head_line, ',');
    std::vector<std::string> data_tokens = split(data_line, ',');

    std::cout << "\nHead tokens:\n";
    for (const auto& token : head_tokens) {
        std::cout << "[" << std::setw(5) << token << "]";
    }

    std::cout << "\nData tokens:\n";
    for (const auto& token : data_tokens) {
        std::cout << "[" << std::setw(5) << token << "]";
    }
    std::cout << "\n";

    return 0;
}
