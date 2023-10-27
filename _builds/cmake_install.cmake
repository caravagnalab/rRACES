# Install script for directory: /Users/gcaravagna/Documents/GitHub/rRACES/RACES

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/Users/gcaravagna/Documents/GitHub/rRACES/_install")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/Library/Developer/CommandLineTools/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/Users/gcaravagna/Documents/GitHub/rRACES/_builds/RACES.0.1.0.dylib")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/RACES.0.1.0.dylib" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/RACES.0.1.0.dylib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Library/Developer/CommandLineTools/usr/bin/strip" -x "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/RACES.0.1.0.dylib")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/Users/gcaravagna/Documents/GitHub/rRACES/_builds/RACES.dylib")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/variables.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/basic_IO.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/archive.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/driver_genotype.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/driver_genotype_id.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/species.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/position.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/simulation.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/tissue.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/cell.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/cell_event.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/simulation_event.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/driver_mutation.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/rate_update.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/tissue_sample.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/sample_specification.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/sampling.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/event_wrapper.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/timed_event.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/logger.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/binary_logger.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/statistics.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/palette.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/plot_2D.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/progress_bar.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/phyloXML.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/phylogenetic_forest.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/palette.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/statistics.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/tissue_plotter.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/time.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/simulation.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/context.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/context_index.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/snv_signature.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/snv.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/allele.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/genomic_region.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/genomic_position.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/fasta_reader.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/fasta_utils.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/cna.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/mutation_engine.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/mutational_properties.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/read_simulator.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/position_set.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/phyloXML.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/phylogenetic_forest.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/indicators.hpp"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/Users/gcaravagna/Documents/GitHub/rRACES/_builds/libRACES.a")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libRACES.a" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libRACES.a")
    execute_process(COMMAND "/Library/Developer/CommandLineTools/usr/bin/ranlib" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libRACES.a")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/variables.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/basic_IO.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/archive.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/driver_genotype.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/driver_genotype_id.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/species.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/position.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/simulation.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/tissue.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/cell.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/cell_event.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/simulation_event.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/driver_mutation.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/rate_update.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/tissue_sample.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/sample_specification.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/sampling.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/event_wrapper.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/timed_event.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/logger.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/binary_logger.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/statistics.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/palette.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/plot_2D.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/progress_bar.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/phyloXML.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/phylogenetic_forest.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/palette.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/statistics.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/tissue_plotter.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/time.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/simulation.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/context.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/context_index.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/snv_signature.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/snv.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/allele.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/genomic_region.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/genomic_position.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/fasta_reader.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/fasta_utils.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/cna.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/mutation_engine.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/mutational_properties.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/read_simulator.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/position_set.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/phyloXML.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/phylogenetic_forest.hpp"
    "/Users/gcaravagna/Documents/GitHub/rRACES/RACES/include/indicators.hpp"
    )
endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/Users/gcaravagna/Documents/GitHub/rRACES/_builds/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
