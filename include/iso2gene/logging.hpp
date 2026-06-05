#pragma once

#include <string>
#include <vector>

namespace iso2gene {

class Logger {
public:
    void info(const std::string& message);
    void warn(const std::string& message);

    bool has_warnings() const;
    const std::vector<std::string>& warnings() const;

    void write_warnings(const std::string& path) const;

private:
    std::vector<std::string> warnings_;
};

} // namespace iso2gene
