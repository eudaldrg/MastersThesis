# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.15

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /home/eudald/Downloads/clion-2018.3.1/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /home/eudald/Downloads/clion-2018.3.1/bin/cmake/linux/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/eudald/Projects/MastersThesis/Programming2

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/eudald/Projects/MastersThesis/Programming2/cmake-build-debug

# Include any dependencies generated for this target.
include FFTW3/CMakeFiles/FFTW3.dir/depend.make

# Include the progress variables for this target.
include FFTW3/CMakeFiles/FFTW3.dir/progress.make

# Include the compile flags for this target's objects.
include FFTW3/CMakeFiles/FFTW3.dir/flags.make

# Object files for target FFTW3
FFTW3_OBJECTS =

# External object files for target FFTW3
FFTW3_EXTERNAL_OBJECTS =

../build_x64/Debug/lib/libFFTW3.a: FFTW3/CMakeFiles/FFTW3.dir/build.make
../build_x64/Debug/lib/libFFTW3.a: FFTW3/CMakeFiles/FFTW3.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/eudald/Projects/MastersThesis/Programming2/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Linking CXX static library ../../build_x64/Debug/lib/libFFTW3.a"
	cd /home/eudald/Projects/MastersThesis/Programming2/cmake-build-debug/FFTW3 && $(CMAKE_COMMAND) -P CMakeFiles/FFTW3.dir/cmake_clean_target.cmake
	cd /home/eudald/Projects/MastersThesis/Programming2/cmake-build-debug/FFTW3 && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/FFTW3.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
FFTW3/CMakeFiles/FFTW3.dir/build: ../build_x64/Debug/lib/libFFTW3.a

.PHONY : FFTW3/CMakeFiles/FFTW3.dir/build

FFTW3/CMakeFiles/FFTW3.dir/clean:
	cd /home/eudald/Projects/MastersThesis/Programming2/cmake-build-debug/FFTW3 && $(CMAKE_COMMAND) -P CMakeFiles/FFTW3.dir/cmake_clean.cmake
.PHONY : FFTW3/CMakeFiles/FFTW3.dir/clean

FFTW3/CMakeFiles/FFTW3.dir/depend:
	cd /home/eudald/Projects/MastersThesis/Programming2/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/eudald/Projects/MastersThesis/Programming2 /home/eudald/Projects/MastersThesis/Programming2/FFTW3 /home/eudald/Projects/MastersThesis/Programming2/cmake-build-debug /home/eudald/Projects/MastersThesis/Programming2/cmake-build-debug/FFTW3 /home/eudald/Projects/MastersThesis/Programming2/cmake-build-debug/FFTW3/CMakeFiles/FFTW3.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : FFTW3/CMakeFiles/FFTW3.dir/depend

