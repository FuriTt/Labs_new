# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.16

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

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "C:\Program Files\JetBrains\CLion 2020.1.1\bin\cmake\win\bin\cmake.exe"

# The command to remove a file.
RM = "C:\Program Files\JetBrains\CLion 2020.1.1\bin\cmake\win\bin\cmake.exe" -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = C:\Users\Brux\Documents\GitHub\Labs_new\cpp_code

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = C:\Users\Brux\Documents\GitHub\Labs_new\cpp_code\cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/teylor.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/teylor.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/teylor.dir/flags.make

CMakeFiles/teylor.dir/teylor.cpp.obj: CMakeFiles/teylor.dir/flags.make
CMakeFiles/teylor.dir/teylor.cpp.obj: ../teylor.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\Brux\Documents\GitHub\Labs_new\cpp_code\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/teylor.dir/teylor.cpp.obj"
	C:\PROGRA~2\MINGW-~1\I686-8~1.0-P\mingw32\bin\G__~1.EXE  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\teylor.dir\teylor.cpp.obj -c C:\Users\Brux\Documents\GitHub\Labs_new\cpp_code\teylor.cpp

CMakeFiles/teylor.dir/teylor.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/teylor.dir/teylor.cpp.i"
	C:\PROGRA~2\MINGW-~1\I686-8~1.0-P\mingw32\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\Users\Brux\Documents\GitHub\Labs_new\cpp_code\teylor.cpp > CMakeFiles\teylor.dir\teylor.cpp.i

CMakeFiles/teylor.dir/teylor.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/teylor.dir/teylor.cpp.s"
	C:\PROGRA~2\MINGW-~1\I686-8~1.0-P\mingw32\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S C:\Users\Brux\Documents\GitHub\Labs_new\cpp_code\teylor.cpp -o CMakeFiles\teylor.dir\teylor.cpp.s

# Object files for target teylor
teylor_OBJECTS = \
"CMakeFiles/teylor.dir/teylor.cpp.obj"

# External object files for target teylor
teylor_EXTERNAL_OBJECTS =

teylor.exe: CMakeFiles/teylor.dir/teylor.cpp.obj
teylor.exe: CMakeFiles/teylor.dir/build.make
teylor.exe: CMakeFiles/teylor.dir/linklibs.rsp
teylor.exe: CMakeFiles/teylor.dir/objects1.rsp
teylor.exe: CMakeFiles/teylor.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=C:\Users\Brux\Documents\GitHub\Labs_new\cpp_code\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable teylor.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\teylor.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/teylor.dir/build: teylor.exe

.PHONY : CMakeFiles/teylor.dir/build

CMakeFiles/teylor.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\teylor.dir\cmake_clean.cmake
.PHONY : CMakeFiles/teylor.dir/clean

CMakeFiles/teylor.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" C:\Users\Brux\Documents\GitHub\Labs_new\cpp_code C:\Users\Brux\Documents\GitHub\Labs_new\cpp_code C:\Users\Brux\Documents\GitHub\Labs_new\cpp_code\cmake-build-debug C:\Users\Brux\Documents\GitHub\Labs_new\cpp_code\cmake-build-debug C:\Users\Brux\Documents\GitHub\Labs_new\cpp_code\cmake-build-debug\CMakeFiles\teylor.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/teylor.dir/depend

