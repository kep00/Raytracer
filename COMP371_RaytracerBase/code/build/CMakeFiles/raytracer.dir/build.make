# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

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
CMAKE_SOURCE_DIR = /COMP371/COMP371_RaytracerBase/code

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /COMP371/COMP371_RaytracerBase/code/build

# Include any dependencies generated for this target.
include CMakeFiles/raytracer.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/raytracer.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/raytracer.dir/flags.make

CMakeFiles/raytracer.dir/main.cpp.o: CMakeFiles/raytracer.dir/flags.make
CMakeFiles/raytracer.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/COMP371/COMP371_RaytracerBase/code/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/raytracer.dir/main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/raytracer.dir/main.cpp.o -c /COMP371/COMP371_RaytracerBase/code/main.cpp

CMakeFiles/raytracer.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/raytracer.dir/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /COMP371/COMP371_RaytracerBase/code/main.cpp > CMakeFiles/raytracer.dir/main.cpp.i

CMakeFiles/raytracer.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/raytracer.dir/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /COMP371/COMP371_RaytracerBase/code/main.cpp -o CMakeFiles/raytracer.dir/main.cpp.s

CMakeFiles/raytracer.dir/external/simpleppm.cpp.o: CMakeFiles/raytracer.dir/flags.make
CMakeFiles/raytracer.dir/external/simpleppm.cpp.o: ../external/simpleppm.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/COMP371/COMP371_RaytracerBase/code/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/raytracer.dir/external/simpleppm.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/raytracer.dir/external/simpleppm.cpp.o -c /COMP371/COMP371_RaytracerBase/code/external/simpleppm.cpp

CMakeFiles/raytracer.dir/external/simpleppm.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/raytracer.dir/external/simpleppm.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /COMP371/COMP371_RaytracerBase/code/external/simpleppm.cpp > CMakeFiles/raytracer.dir/external/simpleppm.cpp.i

CMakeFiles/raytracer.dir/external/simpleppm.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/raytracer.dir/external/simpleppm.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /COMP371/COMP371_RaytracerBase/code/external/simpleppm.cpp -o CMakeFiles/raytracer.dir/external/simpleppm.cpp.s

CMakeFiles/raytracer.dir/external/test_eigen.cpp.o: CMakeFiles/raytracer.dir/flags.make
CMakeFiles/raytracer.dir/external/test_eigen.cpp.o: ../external/test_eigen.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/COMP371/COMP371_RaytracerBase/code/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/raytracer.dir/external/test_eigen.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/raytracer.dir/external/test_eigen.cpp.o -c /COMP371/COMP371_RaytracerBase/code/external/test_eigen.cpp

CMakeFiles/raytracer.dir/external/test_eigen.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/raytracer.dir/external/test_eigen.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /COMP371/COMP371_RaytracerBase/code/external/test_eigen.cpp > CMakeFiles/raytracer.dir/external/test_eigen.cpp.i

CMakeFiles/raytracer.dir/external/test_eigen.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/raytracer.dir/external/test_eigen.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /COMP371/COMP371_RaytracerBase/code/external/test_eigen.cpp -o CMakeFiles/raytracer.dir/external/test_eigen.cpp.s

CMakeFiles/raytracer.dir/external/test_json.cpp.o: CMakeFiles/raytracer.dir/flags.make
CMakeFiles/raytracer.dir/external/test_json.cpp.o: ../external/test_json.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/COMP371/COMP371_RaytracerBase/code/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/raytracer.dir/external/test_json.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/raytracer.dir/external/test_json.cpp.o -c /COMP371/COMP371_RaytracerBase/code/external/test_json.cpp

CMakeFiles/raytracer.dir/external/test_json.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/raytracer.dir/external/test_json.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /COMP371/COMP371_RaytracerBase/code/external/test_json.cpp > CMakeFiles/raytracer.dir/external/test_json.cpp.i

CMakeFiles/raytracer.dir/external/test_json.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/raytracer.dir/external/test_json.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /COMP371/COMP371_RaytracerBase/code/external/test_json.cpp -o CMakeFiles/raytracer.dir/external/test_json.cpp.s

CMakeFiles/raytracer.dir/external/test_ppm.cpp.o: CMakeFiles/raytracer.dir/flags.make
CMakeFiles/raytracer.dir/external/test_ppm.cpp.o: ../external/test_ppm.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/COMP371/COMP371_RaytracerBase/code/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/raytracer.dir/external/test_ppm.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/raytracer.dir/external/test_ppm.cpp.o -c /COMP371/COMP371_RaytracerBase/code/external/test_ppm.cpp

CMakeFiles/raytracer.dir/external/test_ppm.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/raytracer.dir/external/test_ppm.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /COMP371/COMP371_RaytracerBase/code/external/test_ppm.cpp > CMakeFiles/raytracer.dir/external/test_ppm.cpp.i

CMakeFiles/raytracer.dir/external/test_ppm.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/raytracer.dir/external/test_ppm.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /COMP371/COMP371_RaytracerBase/code/external/test_ppm.cpp -o CMakeFiles/raytracer.dir/external/test_ppm.cpp.s

CMakeFiles/raytracer.dir/src/RayTracer.cpp.o: CMakeFiles/raytracer.dir/flags.make
CMakeFiles/raytracer.dir/src/RayTracer.cpp.o: ../src/RayTracer.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/COMP371/COMP371_RaytracerBase/code/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/raytracer.dir/src/RayTracer.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/raytracer.dir/src/RayTracer.cpp.o -c /COMP371/COMP371_RaytracerBase/code/src/RayTracer.cpp

CMakeFiles/raytracer.dir/src/RayTracer.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/raytracer.dir/src/RayTracer.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /COMP371/COMP371_RaytracerBase/code/src/RayTracer.cpp > CMakeFiles/raytracer.dir/src/RayTracer.cpp.i

CMakeFiles/raytracer.dir/src/RayTracer.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/raytracer.dir/src/RayTracer.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /COMP371/COMP371_RaytracerBase/code/src/RayTracer.cpp -o CMakeFiles/raytracer.dir/src/RayTracer.cpp.s

CMakeFiles/raytracer.dir/src/arealight.cpp.o: CMakeFiles/raytracer.dir/flags.make
CMakeFiles/raytracer.dir/src/arealight.cpp.o: ../src/arealight.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/COMP371/COMP371_RaytracerBase/code/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/raytracer.dir/src/arealight.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/raytracer.dir/src/arealight.cpp.o -c /COMP371/COMP371_RaytracerBase/code/src/arealight.cpp

CMakeFiles/raytracer.dir/src/arealight.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/raytracer.dir/src/arealight.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /COMP371/COMP371_RaytracerBase/code/src/arealight.cpp > CMakeFiles/raytracer.dir/src/arealight.cpp.i

CMakeFiles/raytracer.dir/src/arealight.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/raytracer.dir/src/arealight.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /COMP371/COMP371_RaytracerBase/code/src/arealight.cpp -o CMakeFiles/raytracer.dir/src/arealight.cpp.s

CMakeFiles/raytracer.dir/src/obj.cpp.o: CMakeFiles/raytracer.dir/flags.make
CMakeFiles/raytracer.dir/src/obj.cpp.o: ../src/obj.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/COMP371/COMP371_RaytracerBase/code/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CMakeFiles/raytracer.dir/src/obj.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/raytracer.dir/src/obj.cpp.o -c /COMP371/COMP371_RaytracerBase/code/src/obj.cpp

CMakeFiles/raytracer.dir/src/obj.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/raytracer.dir/src/obj.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /COMP371/COMP371_RaytracerBase/code/src/obj.cpp > CMakeFiles/raytracer.dir/src/obj.cpp.i

CMakeFiles/raytracer.dir/src/obj.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/raytracer.dir/src/obj.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /COMP371/COMP371_RaytracerBase/code/src/obj.cpp -o CMakeFiles/raytracer.dir/src/obj.cpp.s

CMakeFiles/raytracer.dir/src/pointlight.cpp.o: CMakeFiles/raytracer.dir/flags.make
CMakeFiles/raytracer.dir/src/pointlight.cpp.o: ../src/pointlight.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/COMP371/COMP371_RaytracerBase/code/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object CMakeFiles/raytracer.dir/src/pointlight.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/raytracer.dir/src/pointlight.cpp.o -c /COMP371/COMP371_RaytracerBase/code/src/pointlight.cpp

CMakeFiles/raytracer.dir/src/pointlight.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/raytracer.dir/src/pointlight.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /COMP371/COMP371_RaytracerBase/code/src/pointlight.cpp > CMakeFiles/raytracer.dir/src/pointlight.cpp.i

CMakeFiles/raytracer.dir/src/pointlight.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/raytracer.dir/src/pointlight.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /COMP371/COMP371_RaytracerBase/code/src/pointlight.cpp -o CMakeFiles/raytracer.dir/src/pointlight.cpp.s

CMakeFiles/raytracer.dir/src/rectangle.cpp.o: CMakeFiles/raytracer.dir/flags.make
CMakeFiles/raytracer.dir/src/rectangle.cpp.o: ../src/rectangle.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/COMP371/COMP371_RaytracerBase/code/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object CMakeFiles/raytracer.dir/src/rectangle.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/raytracer.dir/src/rectangle.cpp.o -c /COMP371/COMP371_RaytracerBase/code/src/rectangle.cpp

CMakeFiles/raytracer.dir/src/rectangle.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/raytracer.dir/src/rectangle.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /COMP371/COMP371_RaytracerBase/code/src/rectangle.cpp > CMakeFiles/raytracer.dir/src/rectangle.cpp.i

CMakeFiles/raytracer.dir/src/rectangle.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/raytracer.dir/src/rectangle.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /COMP371/COMP371_RaytracerBase/code/src/rectangle.cpp -o CMakeFiles/raytracer.dir/src/rectangle.cpp.s

CMakeFiles/raytracer.dir/src/sphere.cpp.o: CMakeFiles/raytracer.dir/flags.make
CMakeFiles/raytracer.dir/src/sphere.cpp.o: ../src/sphere.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/COMP371/COMP371_RaytracerBase/code/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building CXX object CMakeFiles/raytracer.dir/src/sphere.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/raytracer.dir/src/sphere.cpp.o -c /COMP371/COMP371_RaytracerBase/code/src/sphere.cpp

CMakeFiles/raytracer.dir/src/sphere.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/raytracer.dir/src/sphere.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /COMP371/COMP371_RaytracerBase/code/src/sphere.cpp > CMakeFiles/raytracer.dir/src/sphere.cpp.i

CMakeFiles/raytracer.dir/src/sphere.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/raytracer.dir/src/sphere.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /COMP371/COMP371_RaytracerBase/code/src/sphere.cpp -o CMakeFiles/raytracer.dir/src/sphere.cpp.s

# Object files for target raytracer
raytracer_OBJECTS = \
"CMakeFiles/raytracer.dir/main.cpp.o" \
"CMakeFiles/raytracer.dir/external/simpleppm.cpp.o" \
"CMakeFiles/raytracer.dir/external/test_eigen.cpp.o" \
"CMakeFiles/raytracer.dir/external/test_json.cpp.o" \
"CMakeFiles/raytracer.dir/external/test_ppm.cpp.o" \
"CMakeFiles/raytracer.dir/src/RayTracer.cpp.o" \
"CMakeFiles/raytracer.dir/src/arealight.cpp.o" \
"CMakeFiles/raytracer.dir/src/obj.cpp.o" \
"CMakeFiles/raytracer.dir/src/pointlight.cpp.o" \
"CMakeFiles/raytracer.dir/src/rectangle.cpp.o" \
"CMakeFiles/raytracer.dir/src/sphere.cpp.o"

# External object files for target raytracer
raytracer_EXTERNAL_OBJECTS =

raytracer: CMakeFiles/raytracer.dir/main.cpp.o
raytracer: CMakeFiles/raytracer.dir/external/simpleppm.cpp.o
raytracer: CMakeFiles/raytracer.dir/external/test_eigen.cpp.o
raytracer: CMakeFiles/raytracer.dir/external/test_json.cpp.o
raytracer: CMakeFiles/raytracer.dir/external/test_ppm.cpp.o
raytracer: CMakeFiles/raytracer.dir/src/RayTracer.cpp.o
raytracer: CMakeFiles/raytracer.dir/src/arealight.cpp.o
raytracer: CMakeFiles/raytracer.dir/src/obj.cpp.o
raytracer: CMakeFiles/raytracer.dir/src/pointlight.cpp.o
raytracer: CMakeFiles/raytracer.dir/src/rectangle.cpp.o
raytracer: CMakeFiles/raytracer.dir/src/sphere.cpp.o
raytracer: CMakeFiles/raytracer.dir/build.make
raytracer: CMakeFiles/raytracer.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/COMP371/COMP371_RaytracerBase/code/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Linking CXX executable raytracer"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/raytracer.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/raytracer.dir/build: raytracer

.PHONY : CMakeFiles/raytracer.dir/build

CMakeFiles/raytracer.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/raytracer.dir/cmake_clean.cmake
.PHONY : CMakeFiles/raytracer.dir/clean

CMakeFiles/raytracer.dir/depend:
	cd /COMP371/COMP371_RaytracerBase/code/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /COMP371/COMP371_RaytracerBase/code /COMP371/COMP371_RaytracerBase/code /COMP371/COMP371_RaytracerBase/code/build /COMP371/COMP371_RaytracerBase/code/build /COMP371/COMP371_RaytracerBase/code/build/CMakeFiles/raytracer.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/raytracer.dir/depend

