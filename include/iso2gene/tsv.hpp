#pragma once

#include <cstddef>
#include <fstream>
#include <string>
#include <vector>

namespace iso2gene {

std::vector<std::string> split_tsv_line(const std::string& line);
double parse_double_strict(
    const std::string& value,
    const std::string& path,
    std::size_t line_number,
    const std::string& column_name
);

class TsvReader {
public:
    explicit TsvReader(const std::string& path);

    bool read_row(std::vector<std::string>& fields);
    std::size_t line_number() const noexcept;
    const std::string& path() const noexcept;

private:
    std::string path_;
    std::ifstream input_;
    std::size_t line_number_ = 0;
};

int find_column(
    const std::vector<std::string>& header,
    const std::string& column_name,
    const std::string& path
);

} // namespace iso2gene
