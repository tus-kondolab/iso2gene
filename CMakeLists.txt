cmake_minimum_required(VERSION 3.20)

project(iso2gene VERSION 0.1.1 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(iso2gene_miniz
    third_party/miniz/miniz.c
    third_party/miniz/miniz_tdef.c
    third_party/miniz/miniz_tinfl.c
)

target_include_directories(iso2gene_miniz PUBLIC third_party/miniz)
target_compile_definitions(iso2gene_miniz PUBLIC
    MINIZ_NO_ARCHIVE_APIS
    MINIZ_NO_ARCHIVE_WRITING_APIS
    MINIZ_NO_ZLIB_COMPATIBLE_NAMES
)

add_library(iso2gene_core
    src/cli.cpp
    src/error.cpp
    src/gtf.cpp
    src/id.cpp
    src/kallisto.cpp
    src/logging.cpp
    src/quant.cpp
    src/quant_reader.cpp
    src/rsem.cpp
    src/sample_sheet.cpp
    src/salmon.cpp
    src/summarize.cpp
    src/text_reader.cpp
    src/tsv.cpp
    src/tx2gene.cpp
    src/write_matrix.cpp
)

target_include_directories(iso2gene_core PUBLIC include)
target_link_libraries(iso2gene_core PRIVATE iso2gene_miniz)

add_executable(iso2gene src/main.cpp)
target_link_libraries(iso2gene PRIVATE iso2gene_core)

if(MSVC)
    target_compile_options(iso2gene_core PRIVATE /W4 /permissive-)
    target_compile_options(iso2gene PRIVATE /W4 /permissive-)
else()
    target_compile_options(iso2gene_core PRIVATE -Wall -Wextra -Wpedantic)
    target_compile_options(iso2gene PRIVATE -Wall -Wextra -Wpedantic)
endif()

enable_testing()

add_executable(iso2gene_tests tests/test_main.cpp)
target_link_libraries(iso2gene_tests PRIVATE iso2gene_core iso2gene_miniz)

if(MSVC)
    target_compile_options(iso2gene_tests PRIVATE /W4 /permissive-)
else()
    target_compile_options(iso2gene_tests PRIVATE -Wall -Wextra -Wpedantic)
endif()

add_test(NAME iso2gene_tests COMMAND iso2gene_tests)
set_tests_properties(iso2gene_tests PROPERTIES
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
)
