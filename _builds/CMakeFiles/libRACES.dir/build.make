# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.27

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

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.27.7/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.27.7/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/gcaravagna/Documents/GitHub/rRACES/RACES

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/gcaravagna/Documents/GitHub/rRACES/_builds

# Include any dependencies generated for this target.
include CMakeFiles/libRACES.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/libRACES.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/libRACES.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/libRACES.dir/flags.make

# Object files for target libRACES
libRACES_OBJECTS =

# External object files for target libRACES
libRACES_EXTERNAL_OBJECTS = \
"/Users/gcaravagna/Documents/GitHub/rRACES/_builds/CMakeFiles/objlib.dir/src/basic_IO.cpp.o" \
"/Users/gcaravagna/Documents/GitHub/rRACES/_builds/CMakeFiles/objlib.dir/src/archive.cpp.o" \
"/Users/gcaravagna/Documents/GitHub/rRACES/_builds/CMakeFiles/objlib.dir/src/driver_genotype.cpp.o" \
"/Users/gcaravagna/Documents/GitHub/rRACES/_builds/CMakeFiles/objlib.dir/src/species.cpp.o" \
"/Users/gcaravagna/Documents/GitHub/rRACES/_builds/CMakeFiles/objlib.dir/src/position.cpp.o" \
"/Users/gcaravagna/Documents/GitHub/rRACES/_builds/CMakeFiles/objlib.dir/src/simulation.cpp.o" \
"/Users/gcaravagna/Documents/GitHub/rRACES/_builds/CMakeFiles/objlib.dir/src/tissue.cpp.o" \
"/Users/gcaravagna/Documents/GitHub/rRACES/_builds/CMakeFiles/objlib.dir/src/cell.cpp.o" \
"/Users/gcaravagna/Documents/GitHub/rRACES/_builds/CMakeFiles/objlib.dir/src/cell_event.cpp.o" \
"/Users/gcaravagna/Documents/GitHub/rRACES/_builds/CMakeFiles/objlib.dir/src/simulation_event.cpp.o" \
"/Users/gcaravagna/Documents/GitHub/rRACES/_builds/CMakeFiles/objlib.dir/src/driver_mutation.cpp.o" \
"/Users/gcaravagna/Documents/GitHub/rRACES/_builds/CMakeFiles/objlib.dir/src/rate_update.cpp.o" \
"/Users/gcaravagna/Documents/GitHub/rRACES/_builds/CMakeFiles/objlib.dir/src/tissue_sample.cpp.o" \
"/Users/gcaravagna/Documents/GitHub/rRACES/_builds/CMakeFiles/objlib.dir/src/sample_specification.cpp.o" \
"/Users/gcaravagna/Documents/GitHub/rRACES/_builds/CMakeFiles/objlib.dir/src/sampling.cpp.o" \
"/Users/gcaravagna/Documents/GitHub/rRACES/_builds/CMakeFiles/objlib.dir/src/event_wrapper.cpp.o" \
"/Users/gcaravagna/Documents/GitHub/rRACES/_builds/CMakeFiles/objlib.dir/src/timed_event.cpp.o" \
"/Users/gcaravagna/Documents/GitHub/rRACES/_builds/CMakeFiles/objlib.dir/src/logger.cpp.o" \
"/Users/gcaravagna/Documents/GitHub/rRACES/_builds/CMakeFiles/objlib.dir/src/binary_logger.cpp.o" \
"/Users/gcaravagna/Documents/GitHub/rRACES/_builds/CMakeFiles/objlib.dir/src/statistics.cpp.o" \
"/Users/gcaravagna/Documents/GitHub/rRACES/_builds/CMakeFiles/objlib.dir/src/palette.cpp.o" \
"/Users/gcaravagna/Documents/GitHub/rRACES/_builds/CMakeFiles/objlib.dir/src/plot_2D.cpp.o" \
"/Users/gcaravagna/Documents/GitHub/rRACES/_builds/CMakeFiles/objlib.dir/src/progress_bar.cpp.o" \
"/Users/gcaravagna/Documents/GitHub/rRACES/_builds/CMakeFiles/objlib.dir/src/context.cpp.o" \
"/Users/gcaravagna/Documents/GitHub/rRACES/_builds/CMakeFiles/objlib.dir/src/snv_signature.cpp.o" \
"/Users/gcaravagna/Documents/GitHub/rRACES/_builds/CMakeFiles/objlib.dir/src/snv.cpp.o" \
"/Users/gcaravagna/Documents/GitHub/rRACES/_builds/CMakeFiles/objlib.dir/src/allele.cpp.o" \
"/Users/gcaravagna/Documents/GitHub/rRACES/_builds/CMakeFiles/objlib.dir/src/genomic_region.cpp.o" \
"/Users/gcaravagna/Documents/GitHub/rRACES/_builds/CMakeFiles/objlib.dir/src/genomic_position.cpp.o" \
"/Users/gcaravagna/Documents/GitHub/rRACES/_builds/CMakeFiles/objlib.dir/src/fasta_reader.cpp.o" \
"/Users/gcaravagna/Documents/GitHub/rRACES/_builds/CMakeFiles/objlib.dir/src/fasta_utils.cpp.o" \
"/Users/gcaravagna/Documents/GitHub/rRACES/_builds/CMakeFiles/objlib.dir/src/genome_mutations.cpp.o" \
"/Users/gcaravagna/Documents/GitHub/rRACES/_builds/CMakeFiles/objlib.dir/src/cna.cpp.o" \
"/Users/gcaravagna/Documents/GitHub/rRACES/_builds/CMakeFiles/objlib.dir/src/mutation_engine.cpp.o" \
"/Users/gcaravagna/Documents/GitHub/rRACES/_builds/CMakeFiles/objlib.dir/src/mutational_properties.cpp.o" \
"/Users/gcaravagna/Documents/GitHub/rRACES/_builds/CMakeFiles/objlib.dir/src/read_simulator.cpp.o" \
"/Users/gcaravagna/Documents/GitHub/rRACES/_builds/CMakeFiles/objlib.dir/src/position_set.cpp.o" \
"/Users/gcaravagna/Documents/GitHub/rRACES/_builds/CMakeFiles/objlib.dir/src/phyloXML.cpp.o" \
"/Users/gcaravagna/Documents/GitHub/rRACES/_builds/CMakeFiles/objlib.dir/src/phylogenetic_forest.cpp.o"

RACES.0.1.0.dylib: CMakeFiles/objlib.dir/src/basic_IO.cpp.o
RACES.0.1.0.dylib: CMakeFiles/objlib.dir/src/archive.cpp.o
RACES.0.1.0.dylib: CMakeFiles/objlib.dir/src/driver_genotype.cpp.o
RACES.0.1.0.dylib: CMakeFiles/objlib.dir/src/species.cpp.o
RACES.0.1.0.dylib: CMakeFiles/objlib.dir/src/position.cpp.o
RACES.0.1.0.dylib: CMakeFiles/objlib.dir/src/simulation.cpp.o
RACES.0.1.0.dylib: CMakeFiles/objlib.dir/src/tissue.cpp.o
RACES.0.1.0.dylib: CMakeFiles/objlib.dir/src/cell.cpp.o
RACES.0.1.0.dylib: CMakeFiles/objlib.dir/src/cell_event.cpp.o
RACES.0.1.0.dylib: CMakeFiles/objlib.dir/src/simulation_event.cpp.o
RACES.0.1.0.dylib: CMakeFiles/objlib.dir/src/driver_mutation.cpp.o
RACES.0.1.0.dylib: CMakeFiles/objlib.dir/src/rate_update.cpp.o
RACES.0.1.0.dylib: CMakeFiles/objlib.dir/src/tissue_sample.cpp.o
RACES.0.1.0.dylib: CMakeFiles/objlib.dir/src/sample_specification.cpp.o
RACES.0.1.0.dylib: CMakeFiles/objlib.dir/src/sampling.cpp.o
RACES.0.1.0.dylib: CMakeFiles/objlib.dir/src/event_wrapper.cpp.o
RACES.0.1.0.dylib: CMakeFiles/objlib.dir/src/timed_event.cpp.o
RACES.0.1.0.dylib: CMakeFiles/objlib.dir/src/logger.cpp.o
RACES.0.1.0.dylib: CMakeFiles/objlib.dir/src/binary_logger.cpp.o
RACES.0.1.0.dylib: CMakeFiles/objlib.dir/src/statistics.cpp.o
RACES.0.1.0.dylib: CMakeFiles/objlib.dir/src/palette.cpp.o
RACES.0.1.0.dylib: CMakeFiles/objlib.dir/src/plot_2D.cpp.o
RACES.0.1.0.dylib: CMakeFiles/objlib.dir/src/progress_bar.cpp.o
RACES.0.1.0.dylib: CMakeFiles/objlib.dir/src/context.cpp.o
RACES.0.1.0.dylib: CMakeFiles/objlib.dir/src/snv_signature.cpp.o
RACES.0.1.0.dylib: CMakeFiles/objlib.dir/src/snv.cpp.o
RACES.0.1.0.dylib: CMakeFiles/objlib.dir/src/allele.cpp.o
RACES.0.1.0.dylib: CMakeFiles/objlib.dir/src/genomic_region.cpp.o
RACES.0.1.0.dylib: CMakeFiles/objlib.dir/src/genomic_position.cpp.o
RACES.0.1.0.dylib: CMakeFiles/objlib.dir/src/fasta_reader.cpp.o
RACES.0.1.0.dylib: CMakeFiles/objlib.dir/src/fasta_utils.cpp.o
RACES.0.1.0.dylib: CMakeFiles/objlib.dir/src/genome_mutations.cpp.o
RACES.0.1.0.dylib: CMakeFiles/objlib.dir/src/cna.cpp.o
RACES.0.1.0.dylib: CMakeFiles/objlib.dir/src/mutation_engine.cpp.o
RACES.0.1.0.dylib: CMakeFiles/objlib.dir/src/mutational_properties.cpp.o
RACES.0.1.0.dylib: CMakeFiles/objlib.dir/src/read_simulator.cpp.o
RACES.0.1.0.dylib: CMakeFiles/objlib.dir/src/position_set.cpp.o
RACES.0.1.0.dylib: CMakeFiles/objlib.dir/src/phyloXML.cpp.o
RACES.0.1.0.dylib: CMakeFiles/objlib.dir/src/phylogenetic_forest.cpp.o
RACES.0.1.0.dylib: CMakeFiles/libRACES.dir/build.make
RACES.0.1.0.dylib: /usr/local/lib/libboost_program_options-mt.dylib
RACES.0.1.0.dylib: CMakeFiles/libRACES.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/Users/gcaravagna/Documents/GitHub/rRACES/_builds/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Linking CXX shared library RACES.dylib"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/libRACES.dir/link.txt --verbose=$(VERBOSE)
	$(CMAKE_COMMAND) -E cmake_symlink_library RACES.0.1.0.dylib RACES.0.1.0.dylib RACES.dylib

RACES.dylib: RACES.0.1.0.dylib
	@$(CMAKE_COMMAND) -E touch_nocreate RACES.dylib

# Rule to build all files generated by this target.
CMakeFiles/libRACES.dir/build: RACES.dylib
.PHONY : CMakeFiles/libRACES.dir/build

CMakeFiles/libRACES.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/libRACES.dir/cmake_clean.cmake
.PHONY : CMakeFiles/libRACES.dir/clean

CMakeFiles/libRACES.dir/depend:
	cd /Users/gcaravagna/Documents/GitHub/rRACES/_builds && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/gcaravagna/Documents/GitHub/rRACES/RACES /Users/gcaravagna/Documents/GitHub/rRACES/RACES /Users/gcaravagna/Documents/GitHub/rRACES/_builds /Users/gcaravagna/Documents/GitHub/rRACES/_builds /Users/gcaravagna/Documents/GitHub/rRACES/_builds/CMakeFiles/libRACES.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/libRACES.dir/depend

