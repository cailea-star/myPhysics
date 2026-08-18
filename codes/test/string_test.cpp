
#include <string>
#include <vector>
#include <iomanip>
#include <iostream>
#include "string_tools.hpp"

std::string head_line = " A, Z, N [Name] ";
std::string data_line = "  1, 1, 0 [H]    ";

int main() {
    std::cout << "Original head line: " << head_line << "\n";
    std::cout << "Original data line: " << data_line << "\n";

    // 替换逗号为分号
    replace_all(head_line, " ", ""); // 去除空格
    replace_all(head_line, "\t", ""); // 去除制表符
    replace_all(head_line, "[", ","); // 替换左括号为逗号
    replace_all(head_line, "]", ","); // 替换右括号为逗号
    // 同样处理数据行
    replace_all(data_line, " ", ""); // 去除空格
    replace_all(data_line, "\t", ""); // 去除制表符
    replace_all(data_line, "[", ","); // 替换左括号为逗号
    replace_all(data_line, "]", ","); // 替换右括号为逗号

    std::cout << "Modified head line: " << head_line << "\n";
    std::cout << "Modified data line: " << data_line << "\n";

    // 分割字符串
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
