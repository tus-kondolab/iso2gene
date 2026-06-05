#include "iso2gene/text_reader.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <memory>
#include <queue>
#include <sstream>
#include <string>
#include <vector>

#include "iso2gene/error.hpp"

#include "miniz.h"
#include "miniz_tinfl.h"

namespace iso2gene {

namespace {

constexpr std::size_t input_buffer_size = 64 * 1024;

bool ends_with_gz(const std::string& path) {
    if (path.size() < 3) {
        return false;
    }
    const char dot = path[path.size() - 3];
    const char g = path[path.size() - 2];
    const char z = path[path.size() - 1];
    return dot == '.'
        && (g == 'g' || g == 'G')
        && (z == 'z' || z == 'Z');
}

std::uint16_t little_u16(const unsigned char lo, const unsigned char hi) {
    return static_cast<std::uint16_t>(
        static_cast<std::uint16_t>(lo)
        | (static_cast<std::uint16_t>(hi) << 8U)
    );
}

std::uint32_t little_u32(const std::array<unsigned char, 4>& bytes) {
    return static_cast<std::uint32_t>(bytes[0])
        | (static_cast<std::uint32_t>(bytes[1]) << 8U)
        | (static_cast<std::uint32_t>(bytes[2]) << 16U)
        | (static_cast<std::uint32_t>(bytes[3]) << 24U);
}

class LineDecoder {
public:
    void feed(const unsigned char* data, const std::size_t size) {
        for (std::size_t i = 0; i < size; ++i) {
            process_char(static_cast<char>(data[i]));
        }
    }

    bool pop_line(std::string& line) {
        if (ready_.empty()) {
            return false;
        }
        line = std::move(ready_.front());
        ready_.pop();
        return true;
    }

    bool finish(std::string& line) {
        if (!finalized_) {
            if (pending_cr_) {
                push_current_line();
                pending_cr_ = false;
            } else if (!current_.empty()) {
                push_current_line();
            }
            finalized_ = true;
        }
        return pop_line(line);
    }

private:
    void process_char(const char c) {
        if (pending_cr_) {
            push_current_line();
            pending_cr_ = false;
            if (c == '\n') {
                return;
            }
        }

        if (c == '\n') {
            push_current_line();
        } else if (c == '\r') {
            pending_cr_ = true;
        } else {
            current_.push_back(c);
        }
    }

    void push_current_line() {
        if (first_line_) {
            remove_utf8_bom(current_);
            first_line_ = false;
        }
        ready_.push(std::move(current_));
        current_.clear();
    }

    static void remove_utf8_bom(std::string& line) {
        if (line.size() >= 3
            && static_cast<unsigned char>(line[0]) == 0xEFU
            && static_cast<unsigned char>(line[1]) == 0xBBU
            && static_cast<unsigned char>(line[2]) == 0xBFU) {
            line.erase(0, 3);
        }
    }

    std::queue<std::string> ready_;
    std::string current_;
    bool pending_cr_ = false;
    bool first_line_ = true;
    bool finalized_ = false;
};

class ByteSource {
public:
    virtual ~ByteSource() = default;
    virtual std::size_t read(unsigned char* buffer, std::size_t size) = 0;
};

class PlainFileSource final : public ByteSource {
public:
    explicit PlainFileSource(const std::string& path)
        : path_(path), input_(path, std::ios::binary) {
        if (!input_) {
            throw Iso2GeneError(ExitCode::io_error, "failed to open input file: " + path);
        }
    }

    std::size_t read(unsigned char* buffer, const std::size_t size) override {
        input_.read(reinterpret_cast<char*>(buffer), static_cast<std::streamsize>(size));
        const std::streamsize read_count = input_.gcount();
        if (read_count > 0) {
            return static_cast<std::size_t>(read_count);
        }
        if (input_.bad()) {
            throw Iso2GeneError(ExitCode::io_error, "failed to read input file: " + path_);
        }
        return 0;
    }

private:
    std::string path_;
    std::ifstream input_;
};

class GzipFileSource final : public ByteSource {
public:
    explicit GzipFileSource(const std::string& path)
        : path_(path), input_(path, std::ios::binary) {
        if (!input_) {
            throw Iso2GeneError(ExitCode::io_error, "failed to open gzip file: " + path);
        }
    }

    std::size_t read(unsigned char* buffer, const std::size_t size) override {
        while (true) {
            if (done_) {
                return 0;
            }
            if (!member_active_) {
                if (!start_next_member()) {
                    done_ = true;
                    return 0;
                }
            }

            const std::size_t produced = inflate_member(buffer, size);
            if (produced > 0) {
                return produced;
            }
        }
    }

private:
    bool refill_input() {
        if (input_pos_ < input_size_) {
            return true;
        }
        if (eof_) {
            return false;
        }

        input_.read(
            reinterpret_cast<char*>(input_buffer_.data()),
            static_cast<std::streamsize>(input_buffer_.size())
        );
        const std::streamsize read_count = input_.gcount();
        if (read_count > 0) {
            input_pos_ = 0;
            input_size_ = static_cast<std::size_t>(read_count);
            if (read_count < static_cast<std::streamsize>(input_buffer_.size())) {
                eof_ = true;
            }
            return true;
        }
        if (input_.bad()) {
            throw Iso2GeneError(ExitCode::io_error, "failed to read gzip file: " + path_);
        }
        eof_ = true;
        return false;
    }

    bool read_byte(unsigned char& byte) {
        if (!refill_input()) {
            return false;
        }
        byte = input_buffer_[input_pos_++];
        return true;
    }

    unsigned char require_byte(const std::string& context) {
        unsigned char byte = 0;
        if (!read_byte(byte)) {
            throw Iso2GeneError(
                ExitCode::io_error,
                "failed to read gzip file " + path_ + ": truncated " + context
            );
        }
        return byte;
    }

    void skip_bytes(const std::size_t count, const std::string& context) {
        for (std::size_t i = 0; i < count; ++i) {
            (void)require_byte(context);
        }
    }

    void skip_zero_terminated_string(const std::string& context) {
        while (true) {
            if (require_byte(context) == 0U) {
                return;
            }
        }
    }

    bool start_next_member() {
        unsigned char id1 = 0;
        if (!read_byte(id1)) {
            return false;
        }
        const unsigned char id2 = require_byte("gzip header");
        if (id1 != 0x1FU || id2 != 0x8BU) {
            throw Iso2GeneError(ExitCode::io_error, "invalid gzip header in " + path_);
        }

        const unsigned char compression_method = require_byte("gzip header");
        if (compression_method != 8U) {
            throw Iso2GeneError(
                ExitCode::io_error,
                "unsupported gzip compression method in " + path_
            );
        }

        const unsigned char flags = require_byte("gzip header");
        if ((flags & 0xE0U) != 0U) {
            throw Iso2GeneError(ExitCode::io_error, "invalid gzip header in " + path_);
        }

        skip_bytes(6, "gzip header");

        if ((flags & 0x04U) != 0U) {
            const unsigned char lo = require_byte("gzip extra field");
            const unsigned char hi = require_byte("gzip extra field");
            skip_bytes(little_u16(lo, hi), "gzip extra field");
        }
        if ((flags & 0x08U) != 0U) {
            skip_zero_terminated_string("gzip original file name");
        }
        if ((flags & 0x10U) != 0U) {
            skip_zero_terminated_string("gzip comment");
        }
        if ((flags & 0x02U) != 0U) {
            skip_bytes(2, "gzip header CRC");
        }

        tinfl_init(&inflater_);
        dict_offset_ = 0;
        crc32_ = static_cast<std::uint32_t>(MZ_CRC32_INIT);
        size_mod32_ = 0;
        member_active_ = true;
        return true;
    }

    std::size_t inflate_member(unsigned char* buffer, const std::size_t size) {
        if (size == 0) {
            return 0;
        }

        while (true) {
            if (input_pos_ == input_size_ && !refill_input()) {
                std::size_t input_bytes = 0;
                std::size_t output_bytes =
                    std::min(size, dictionary_.size() - dict_offset_);
                const tinfl_status status = tinfl_decompress(
                    &inflater_,
                    input_buffer_.data(),
                    &input_bytes,
                    dictionary_.data(),
                    dictionary_.data() + dict_offset_,
                    &output_bytes,
                    0
                );
                if (status < TINFL_STATUS_DONE) {
                    throw Iso2GeneError(
                        ExitCode::io_error,
                        "failed to read gzip file " + path_ + ": truncated gzip stream"
                    );
                }
                if (output_bytes > 0) {
                    copy_output(buffer, output_bytes);
                    return output_bytes;
                }
                if (status == TINFL_STATUS_DONE) {
                    finish_member();
                    return 0;
                }
                throw Iso2GeneError(
                    ExitCode::io_error,
                    "failed to read gzip file " + path_ + ": truncated gzip stream"
                );
            }

            std::size_t input_bytes = input_size_ - input_pos_;
            std::size_t output_bytes = std::min(size, dictionary_.size() - dict_offset_);
            const mz_uint32 flags = eof_
                ? 0U
                : static_cast<mz_uint32>(TINFL_FLAG_HAS_MORE_INPUT);
            const tinfl_status status = tinfl_decompress(
                &inflater_,
                input_buffer_.data() + input_pos_,
                &input_bytes,
                dictionary_.data(),
                dictionary_.data() + dict_offset_,
                &output_bytes,
                flags
            );
            input_pos_ += input_bytes;

            if (status < TINFL_STATUS_DONE) {
                throw Iso2GeneError(
                    ExitCode::io_error,
                    "failed to read gzip file " + path_ + ": invalid deflate stream"
                );
            }
            if (output_bytes > 0) {
                copy_output(buffer, output_bytes);
                if (status == TINFL_STATUS_DONE) {
                    finish_member();
                }
                return output_bytes;
            }
            if (status == TINFL_STATUS_DONE) {
                finish_member();
                return 0;
            }
            if (status == TINFL_STATUS_NEEDS_MORE_INPUT && eof_) {
                throw Iso2GeneError(
                    ExitCode::io_error,
                    "failed to read gzip file " + path_ + ": truncated gzip stream"
                );
            }
        }
    }

    void copy_output(unsigned char* buffer, const std::size_t output_bytes) {
        const unsigned char* produced = dictionary_.data() + dict_offset_;
        std::copy(produced, produced + output_bytes, buffer);
        crc32_ = static_cast<std::uint32_t>(mz_crc32(crc32_, produced, output_bytes));
        size_mod32_ += static_cast<std::uint32_t>(output_bytes);
        dict_offset_ = (dict_offset_ + output_bytes) & (TINFL_LZ_DICT_SIZE - 1U);
    }

    void finish_member() {
        const std::array<unsigned char, 4> expected_crc{
            require_byte("gzip trailer"),
            require_byte("gzip trailer"),
            require_byte("gzip trailer"),
            require_byte("gzip trailer")
        };
        const std::array<unsigned char, 4> expected_size{
            require_byte("gzip trailer"),
            require_byte("gzip trailer"),
            require_byte("gzip trailer"),
            require_byte("gzip trailer")
        };

        if (little_u32(expected_crc) != crc32_) {
            throw Iso2GeneError(
                ExitCode::io_error,
                "failed to read gzip file " + path_ + ": gzip CRC mismatch"
            );
        }
        if (little_u32(expected_size) != size_mod32_) {
            throw Iso2GeneError(
                ExitCode::io_error,
                "failed to read gzip file " + path_ + ": gzip size mismatch"
            );
        }

        member_active_ = false;
    }

    std::string path_;
    std::ifstream input_;
    std::array<unsigned char, input_buffer_size> input_buffer_{};
    std::size_t input_pos_ = 0;
    std::size_t input_size_ = 0;
    bool eof_ = false;
    bool done_ = false;
    bool member_active_ = false;

    tinfl_decompressor inflater_{};
    std::array<unsigned char, TINFL_LZ_DICT_SIZE> dictionary_{};
    std::size_t dict_offset_ = 0;
    std::uint32_t crc32_ = static_cast<std::uint32_t>(MZ_CRC32_INIT);
    std::uint32_t size_mod32_ = 0;
};

} // namespace

class TextReader::Impl {
public:
    explicit Impl(const std::string& path)
        : path_(path),
          source_(ends_with_gz(path)
              ? std::unique_ptr<ByteSource>(new GzipFileSource(path))
              : std::unique_ptr<ByteSource>(new PlainFileSource(path))) {}

    bool read_line(std::string& line) {
        while (true) {
            if (decoder_.pop_line(line)) {
                ++line_number_;
                return true;
            }
            if (source_done_) {
                if (decoder_.finish(line)) {
                    ++line_number_;
                    return true;
                }
                return false;
            }

            std::array<unsigned char, input_buffer_size> buffer;
            const std::size_t read_count = source_->read(buffer.data(), buffer.size());
            if (read_count == 0) {
                source_done_ = true;
            } else {
                decoder_.feed(buffer.data(), read_count);
            }
        }
    }

    std::size_t line_number() const noexcept {
        return line_number_;
    }

    const std::string& path() const noexcept {
        return path_;
    }

private:
    std::string path_;
    std::unique_ptr<ByteSource> source_;
    LineDecoder decoder_;
    std::size_t line_number_ = 0;
    bool source_done_ = false;
};

TextReader::TextReader(const std::string& path)
    : impl_(new Impl(path)) {}

TextReader::~TextReader() = default;

TextReader::TextReader(TextReader&&) noexcept = default;

TextReader& TextReader::operator=(TextReader&&) noexcept = default;

bool TextReader::read_line(std::string& line) {
    return impl_->read_line(line);
}

std::size_t TextReader::line_number() const noexcept {
    return impl_->line_number();
}

const std::string& TextReader::path() const noexcept {
    return impl_->path();
}

} // namespace iso2gene
