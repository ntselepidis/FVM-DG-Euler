#! /usr/bin/env python3

import glob
import os
import errno
import os.path

SUFFIXES = [".c", ".C", ".cpp", ".c++"]

def find_files(folder):
    files = sum((glob.glob("{}*{}".format(folder, s)) for s in SUFFIXES), [])
    return [os.path.basename(f) for f in sorted(files)]

def find_source_files(folder):
    return find_files(folder)

def find_subdirectories(folder):
    dirs = sorted(glob.glob(folder + "*/"))
    return [d for d in dirs if "CMake" not in d]

def format_sources(target, sources):
    ret = ""

    line_pattern = "  PRIVATE ${{CMAKE_CURRENT_LIST_DIR}}/{:s}\n"
    if sources:
        ret += "".join(["target_sources(" + target + "\n",
                        "".join(line_pattern.format(s) for s in sources),
                        ")\n\n"])

    return ret

def add_subdirectory(folder):
    line_pattern = "add_subdirectory({:s})\n"
    return  line_pattern.format(os.path.basename(folder[:-1]))

def remove_file(filename):
    try:
        os.remove(filename)
    except OSError as e:
        if e.errno != errno.ENOENT:
            raise

def append_to_file(filename, text):
    with open(filename, 'a') as f:
        f.write(text)

def recurse(base_directory, target):
    cmake_file = base_directory + "CMakeLists.txt"
    remove_file(cmake_file)

    source_files = find_source_files(base_directory)
    append_to_file(cmake_file, format_sources(target, source_files))

    for d in find_subdirectories(base_directory):
        recurse(d, target)
        append_to_file(cmake_file, add_subdirectory(base_directory + d))

def add_executable(cmake_file, target, source_file):
    line = """
target_sources({}
  PRIVATE ${{CMAKE_CURRENT_SOURCE_DIR}}/{}
)
""".format(target, source_file)
    append_to_file(cmake_file, line)


if __name__ == "__main__":

    cmake_file = "src/CMakeLists.txt"
    remove_file(cmake_file)

    base_directory = "src/"
    for d in find_subdirectories(base_directory):
        recurse(d, "fvm_scalar_obj")
        append_to_file(cmake_file, add_subdirectory(base_directory + d))

    add_executable(cmake_file, "fvm_scalar", "fvm_scalar.cpp")

    recurse("tests/", "unit_tests")
