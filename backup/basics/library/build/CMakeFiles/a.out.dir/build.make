# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.7

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/ratnesh/Desktop/Computational-Gas-Dynamics/c++/2D

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/ratnesh/Desktop/Computational-Gas-Dynamics/c++/2D/build

# Include any dependencies generated for this target.
include CMakeFiles/a.out.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/a.out.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/a.out.dir/flags.make

CMakeFiles/a.out.dir/source/advection.cc.o: CMakeFiles/a.out.dir/flags.make
CMakeFiles/a.out.dir/source/advection.cc.o: ../source/advection.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ratnesh/Desktop/Computational-Gas-Dynamics/c++/2D/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/a.out.dir/source/advection.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/a.out.dir/source/advection.cc.o -c /home/ratnesh/Desktop/Computational-Gas-Dynamics/c++/2D/source/advection.cc

CMakeFiles/a.out.dir/source/advection.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/a.out.dir/source/advection.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ratnesh/Desktop/Computational-Gas-Dynamics/c++/2D/source/advection.cc > CMakeFiles/a.out.dir/source/advection.cc.i

CMakeFiles/a.out.dir/source/advection.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/a.out.dir/source/advection.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ratnesh/Desktop/Computational-Gas-Dynamics/c++/2D/source/advection.cc -o CMakeFiles/a.out.dir/source/advection.cc.s

CMakeFiles/a.out.dir/source/advection.cc.o.requires:

.PHONY : CMakeFiles/a.out.dir/source/advection.cc.o.requires

CMakeFiles/a.out.dir/source/advection.cc.o.provides: CMakeFiles/a.out.dir/source/advection.cc.o.requires
	$(MAKE) -f CMakeFiles/a.out.dir/build.make CMakeFiles/a.out.dir/source/advection.cc.o.provides.build
.PHONY : CMakeFiles/a.out.dir/source/advection.cc.o.provides

CMakeFiles/a.out.dir/source/advection.cc.o.provides.build: CMakeFiles/a.out.dir/source/advection.cc.o


CMakeFiles/a.out.dir/source/field.cc.o: CMakeFiles/a.out.dir/flags.make
CMakeFiles/a.out.dir/source/field.cc.o: ../source/field.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ratnesh/Desktop/Computational-Gas-Dynamics/c++/2D/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/a.out.dir/source/field.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/a.out.dir/source/field.cc.o -c /home/ratnesh/Desktop/Computational-Gas-Dynamics/c++/2D/source/field.cc

CMakeFiles/a.out.dir/source/field.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/a.out.dir/source/field.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ratnesh/Desktop/Computational-Gas-Dynamics/c++/2D/source/field.cc > CMakeFiles/a.out.dir/source/field.cc.i

CMakeFiles/a.out.dir/source/field.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/a.out.dir/source/field.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ratnesh/Desktop/Computational-Gas-Dynamics/c++/2D/source/field.cc -o CMakeFiles/a.out.dir/source/field.cc.s

CMakeFiles/a.out.dir/source/field.cc.o.requires:

.PHONY : CMakeFiles/a.out.dir/source/field.cc.o.requires

CMakeFiles/a.out.dir/source/field.cc.o.provides: CMakeFiles/a.out.dir/source/field.cc.o.requires
	$(MAKE) -f CMakeFiles/a.out.dir/build.make CMakeFiles/a.out.dir/source/field.cc.o.provides.build
.PHONY : CMakeFiles/a.out.dir/source/field.cc.o.provides

CMakeFiles/a.out.dir/source/field.cc.o.provides.build: CMakeFiles/a.out.dir/source/field.cc.o


# Object files for target a.out
a_out_OBJECTS = \
"CMakeFiles/a.out.dir/source/advection.cc.o" \
"CMakeFiles/a.out.dir/source/field.cc.o"

# External object files for target a.out
a_out_EXTERNAL_OBJECTS =

a.out: CMakeFiles/a.out.dir/source/advection.cc.o
a.out: CMakeFiles/a.out.dir/source/field.cc.o
a.out: CMakeFiles/a.out.dir/build.make
a.out: CMakeFiles/a.out.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/ratnesh/Desktop/Computational-Gas-Dynamics/c++/2D/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable a.out"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/a.out.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/a.out.dir/build: a.out

.PHONY : CMakeFiles/a.out.dir/build

CMakeFiles/a.out.dir/requires: CMakeFiles/a.out.dir/source/advection.cc.o.requires
CMakeFiles/a.out.dir/requires: CMakeFiles/a.out.dir/source/field.cc.o.requires

.PHONY : CMakeFiles/a.out.dir/requires

CMakeFiles/a.out.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/a.out.dir/cmake_clean.cmake
.PHONY : CMakeFiles/a.out.dir/clean

CMakeFiles/a.out.dir/depend:
	cd /home/ratnesh/Desktop/Computational-Gas-Dynamics/c++/2D/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ratnesh/Desktop/Computational-Gas-Dynamics/c++/2D /home/ratnesh/Desktop/Computational-Gas-Dynamics/c++/2D /home/ratnesh/Desktop/Computational-Gas-Dynamics/c++/2D/build /home/ratnesh/Desktop/Computational-Gas-Dynamics/c++/2D/build /home/ratnesh/Desktop/Computational-Gas-Dynamics/c++/2D/build/CMakeFiles/a.out.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/a.out.dir/depend

