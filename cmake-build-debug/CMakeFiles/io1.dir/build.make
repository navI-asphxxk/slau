# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.20

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "C:\Program Files\JetBrains\CLion 2021.2.2\bin\cmake\win\bin\cmake.exe"

# The command to remove a file.
RM = "C:\Program Files\JetBrains\CLion 2021.2.2\bin\cmake\win\bin\cmake.exe" -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = C:\Users\maksi\CLionProjects\io1

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = C:\Users\maksi\CLionProjects\io1\cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/io1.dir/depend.make
# Include the progress variables for this target.
include CMakeFiles/io1.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/io1.dir/flags.make

CMakeFiles/io1.dir/main.cpp.obj: CMakeFiles/io1.dir/flags.make
CMakeFiles/io1.dir/main.cpp.obj: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\maksi\CLionProjects\io1\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/io1.dir/main.cpp.obj"
	C:\PROGRA~1\MINGW-~1\X86_64~1.0-P\mingw64\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\io1.dir\main.cpp.obj -c C:\Users\maksi\CLionProjects\io1\main.cpp

CMakeFiles/io1.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/io1.dir/main.cpp.i"
	C:\PROGRA~1\MINGW-~1\X86_64~1.0-P\mingw64\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\Users\maksi\CLionProjects\io1\main.cpp > CMakeFiles\io1.dir\main.cpp.i

CMakeFiles/io1.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/io1.dir/main.cpp.s"
	C:\PROGRA~1\MINGW-~1\X86_64~1.0-P\mingw64\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S C:\Users\maksi\CLionProjects\io1\main.cpp -o CMakeFiles\io1.dir\main.cpp.s

# Object files for target io1
io1_OBJECTS = \
"CMakeFiles/io1.dir/main.cpp.obj"

# External object files for target io1
io1_EXTERNAL_OBJECTS =

io1.exe: CMakeFiles/io1.dir/main.cpp.obj
io1.exe: CMakeFiles/io1.dir/build.make
io1.exe: CMakeFiles/io1.dir/linklibs.rsp
io1.exe: CMakeFiles/io1.dir/objects1.rsp
io1.exe: CMakeFiles/io1.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=C:\Users\maksi\CLionProjects\io1\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable io1.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\io1.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/io1.dir/build: io1.exe
.PHONY : CMakeFiles/io1.dir/build

CMakeFiles/io1.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\io1.dir\cmake_clean.cmake
.PHONY : CMakeFiles/io1.dir/clean

CMakeFiles/io1.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" C:\Users\maksi\CLionProjects\io1 C:\Users\maksi\CLionProjects\io1 C:\Users\maksi\CLionProjects\io1\cmake-build-debug C:\Users\maksi\CLionProjects\io1\cmake-build-debug C:\Users\maksi\CLionProjects\io1\cmake-build-debug\CMakeFiles\io1.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/io1.dir/depend

