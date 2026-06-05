#pragma once

#include <cstddef>
#include <memory>
#include <string>

namespace iso2gene {

class TextReader {
public:
    explicit TextReader(const std::string& path);
    ~TextReader();

    TextReader(const TextReader&) = delete;
    TextReader& operator=(const TextReader&) = delete;
    TextReader(TextReader&&) noexcept;
    TextReader& operator=(TextReader&&) noexcept;

    bool read_line(std::string& line);
    std::size_t line_number() const noexcept;
    const std::string& path() const noexcept;

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace iso2gene
