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
CMAKE_SOURCE_DIR = C:\Users\Brux\Documents\GitHub\Labs_new\untitled2

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = C:\Users\Brux\Documents\GitHub\Labs_new\untitled2\cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/1.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/1.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/1.dir/flags.make

CMakeFiles/1.dir/main.cpp.obj: CMakeFiles/1.dir/flags.make
CMakeFiles/1.dir/main.cpp.obj: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\Brux\Documents\GitHub\Labs_new\untitled2\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/1.dir/main.cpp.obj"
	C:\PROGRA~2\MINGW-~1\I686-8~1.0-P\mingw32\bin\G__~1.EXE  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\1.dir\main.cpp.obj -c C:\Users\Brux\Documents\GitHub\Labs_new\untitled2\main.cpp

CMakeFiles/1.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/1.dir/main.cpp.i"
	C:\PROGRA~2\MINGW-~1\I686-8~1.0-P\mingw32\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\Users\Brux\Documents\GitHub\Labs_new\untitled2\main.cpp > CMakeFiles\1.dir\main.cpp.i

CMakeFiles/1.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/1.dir/main.cpp.s"
	C:\PROGRA~2\MINGW-~1\I686-8~1.0-P\mingw32\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S C:\Users\Brux\Documents\GitHub\Labs_new\untitled2\main.cpp -o CMakeFiles\1.dir\main.cpp.s

# Object files for target 1
1_OBJECTS = \
"CMakeFiles/1.dir/main.cpp.obj"

# External object files for target 1
1_EXTERNAL_OBJECTS =

1.exe: CMakeFiles/1.dir/main.cpp.obj
1.exe: CMakeFiles/1.dir/build.make
1.exe: CMakeFiles/1.dir/linklibs.rsp
1.exe: CMakeFiles/1.dir/objects1.rsp
1.exe: CMakeFiles/1.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=C:\Users\Brux\Documents\GitHub\Labs_new\untitled2\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable 1.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\1.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/1.dir/build: 1.exe

.PHONY : CMakeFiles/1.dir/build

CMakeFiles/1.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\1.dir\cmake_clean.cmake
.PHONY : CMakeFiles/1.dir/clean

CMakeFiles/1.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" C:\Users\Brux\Documents\GitHub\Labs_new\untitled2 C:\Users\Brux\Documents\GitHub\Labs_new\untitled2 C:\Users\Brux\Documents\GitHub\Labs_new\untitled2\cmake-build-debug C:\Users\Brux\Documents\GitHub\Labs_new\untitled2\cmake-build-debug C:\Users\Brux\Documents\GitHub\Labs_new\untitled2\cmake-build-debug\CMakeFiles\1.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/1.dir/depend
