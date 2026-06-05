#pragma once

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace iso2gene {

template <typename T>
class Matrix {
public:
    Matrix() = default;

    Matrix(std::size_t rows, std::size_t cols, T init = T{})
        : rows_(rows), cols_(cols), data_(rows * cols, init) {}

    T& operator()(std::size_t row, std::size_t col) {
        return data_.at(index(row, col));
    }

    const T& operator()(std::size_t row, std::size_t col) const {
        return data_.at(index(row, col));
    }

    std::size_t rows() const noexcept {
        return rows_;
    }

    std::size_t cols() const noexcept {
        return cols_;
    }

private:
    std::size_t index(std::size_t row, std::size_t col) const {
        if (row >= rows_ || col >= cols_) {
            throw std::out_of_range("matrix index out of range");
        }
        return row * cols_ + col;
    }

    std::size_t rows_ = 0;
    std::size_t cols_ = 0;
    std::vector<T> data_;
};

} // namespace iso2gene
