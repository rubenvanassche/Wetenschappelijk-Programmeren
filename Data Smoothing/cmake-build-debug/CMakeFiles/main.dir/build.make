# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.8

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
CMAKE_COMMAND = /Applications/CLion.app/Contents/bin/cmake/bin/cmake

# The command to remove a file.
RM = /Applications/CLion.app/Contents/bin/cmake/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/Users/ruben/Documents/School/Wetenschappelijk-Programmeren/Data Smoothing"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/Users/ruben/Documents/School/Wetenschappelijk-Programmeren/Data Smoothing/cmake-build-debug"

# Include any dependencies generated for this target.
include CMakeFiles/main.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/main.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/main.dir/flags.make

CMakeFiles/main.dir/main.cpp.o: CMakeFiles/main.dir/flags.make
CMakeFiles/main.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/Users/ruben/Documents/School/Wetenschappelijk-Programmeren/Data Smoothing/cmake-build-debug/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/main.dir/main.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/main.dir/main.cpp.o -c "/Users/ruben/Documents/School/Wetenschappelijk-Programmeren/Data Smoothing/main.cpp"

CMakeFiles/main.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main.dir/main.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/Users/ruben/Documents/School/Wetenschappelijk-Programmeren/Data Smoothing/main.cpp" > CMakeFiles/main.dir/main.cpp.i

CMakeFiles/main.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main.dir/main.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/Users/ruben/Documents/School/Wetenschappelijk-Programmeren/Data Smoothing/main.cpp" -o CMakeFiles/main.dir/main.cpp.s

CMakeFiles/main.dir/main.cpp.o.requires:

.PHONY : CMakeFiles/main.dir/main.cpp.o.requires

CMakeFiles/main.dir/main.cpp.o.provides: CMakeFiles/main.dir/main.cpp.o.requires
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/main.cpp.o.provides.build
.PHONY : CMakeFiles/main.dir/main.cpp.o.provides

CMakeFiles/main.dir/main.cpp.o.provides.build: CMakeFiles/main.dir/main.cpp.o


CMakeFiles/main.dir/PointsWriter.cpp.o: CMakeFiles/main.dir/flags.make
CMakeFiles/main.dir/PointsWriter.cpp.o: ../PointsWriter.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/Users/ruben/Documents/School/Wetenschappelijk-Programmeren/Data Smoothing/cmake-build-debug/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/main.dir/PointsWriter.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/main.dir/PointsWriter.cpp.o -c "/Users/ruben/Documents/School/Wetenschappelijk-Programmeren/Data Smoothing/PointsWriter.cpp"

CMakeFiles/main.dir/PointsWriter.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main.dir/PointsWriter.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/Users/ruben/Documents/School/Wetenschappelijk-Programmeren/Data Smoothing/PointsWriter.cpp" > CMakeFiles/main.dir/PointsWriter.cpp.i

CMakeFiles/main.dir/PointsWriter.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main.dir/PointsWriter.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/Users/ruben/Documents/School/Wetenschappelijk-Programmeren/Data Smoothing/PointsWriter.cpp" -o CMakeFiles/main.dir/PointsWriter.cpp.s

CMakeFiles/main.dir/PointsWriter.cpp.o.requires:

.PHONY : CMakeFiles/main.dir/PointsWriter.cpp.o.requires

CMakeFiles/main.dir/PointsWriter.cpp.o.provides: CMakeFiles/main.dir/PointsWriter.cpp.o.requires
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/PointsWriter.cpp.o.provides.build
.PHONY : CMakeFiles/main.dir/PointsWriter.cpp.o.provides

CMakeFiles/main.dir/PointsWriter.cpp.o.provides.build: CMakeFiles/main.dir/PointsWriter.cpp.o


# Object files for target main
main_OBJECTS = \
"CMakeFiles/main.dir/main.cpp.o" \
"CMakeFiles/main.dir/PointsWriter.cpp.o"

# External object files for target main
main_EXTERNAL_OBJECTS =

main: CMakeFiles/main.dir/main.cpp.o
main: CMakeFiles/main.dir/PointsWriter.cpp.o
main: CMakeFiles/main.dir/build.make
main: /usr/local/Cellar/gsl/2.4/lib/libgsl.dylib
main: /usr/local/Cellar/gsl/2.4/lib/libgslcblas.dylib
main: CMakeFiles/main.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/Users/ruben/Documents/School/Wetenschappelijk-Programmeren/Data Smoothing/cmake-build-debug/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable main"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/main.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/main.dir/build: main

.PHONY : CMakeFiles/main.dir/build

CMakeFiles/main.dir/requires: CMakeFiles/main.dir/main.cpp.o.requires
CMakeFiles/main.dir/requires: CMakeFiles/main.dir/PointsWriter.cpp.o.requires

.PHONY : CMakeFiles/main.dir/requires

CMakeFiles/main.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/main.dir/cmake_clean.cmake
.PHONY : CMakeFiles/main.dir/clean

CMakeFiles/main.dir/depend:
	cd "/Users/ruben/Documents/School/Wetenschappelijk-Programmeren/Data Smoothing/cmake-build-debug" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/Users/ruben/Documents/School/Wetenschappelijk-Programmeren/Data Smoothing" "/Users/ruben/Documents/School/Wetenschappelijk-Programmeren/Data Smoothing" "/Users/ruben/Documents/School/Wetenschappelijk-Programmeren/Data Smoothing/cmake-build-debug" "/Users/ruben/Documents/School/Wetenschappelijk-Programmeren/Data Smoothing/cmake-build-debug" "/Users/ruben/Documents/School/Wetenschappelijk-Programmeren/Data Smoothing/cmake-build-debug/CMakeFiles/main.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/main.dir/depend

