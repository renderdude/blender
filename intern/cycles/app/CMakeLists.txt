# SPDX-FileCopyrightText: 2011-2022 Blender Foundation
#
# SPDX-License-Identifier: Apache-2.0

#####################################################################
# Cycles standalone executable
#####################################################################
function(print_target_properties tgt)
    if(NOT TARGET ${tgt})
      message("There is no target named '${tgt}'")
      return()
    endif()

    foreach (prop ${CMAKE_PROPERTY_LIST})
        string(REPLACE "<CONFIG>" "${CMAKE_BUILD_TYPE}" prop ${prop})
        get_target_property(propval ${tgt} ${prop})
        if (propval)
            message ("${tgt} ${prop} = ${propval}")
        endif()
    endforeach(prop)
endfunction(print_target_properties)

set(INC
  ..
)
set(INC_SYS
)

set(LIB
  cycles_device
  cycles_kernel
  cycles_scene
  cycles_session
  cycles_bvh
  cycles_subd
  cycles_graph
  cycles_util
)

if(WITH_CYCLES_DISTRIBUTED)
  list(APPEND LIB cycles_distributed)
endif()

include (CheckCXXSourceCompiles)

check_cxx_source_compiles ("
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
int main() {
   int fd = open(\"foo\", O_RDONLY);
   struct stat s;
   fstat(fd, &s);
   size_t len = s.st_size;
   void *ptr = mmap(0, len, PROT_READ, MAP_FILE | MAP_SHARED, fd, 0);
   munmap(ptr, len);   
}
" HAVE_MMAP)

if (HAVE_MMAP)
  add_definitions(-DHAVE_MMAP)
ENDIF ()

if(WITH_ALEMBIC)
  add_definitions(-DWITH_ALEMBIC)
  list(APPEND INC_SYS
    ${ALEMBIC_INCLUDE_DIRS}
  )
  list(APPEND LIB
    ${ALEMBIC_LIBRARIES}
  )
endif()

if(WITH_CYCLES_OSL)
  list(APPEND LIB cycles_kernel_osl)
endif()

if(CYCLES_STANDALONE_REPOSITORY)
  list(APPEND LIB extern_sky)
else()
  list(APPEND LIB bf_intern_sky)
endif()

if(WITH_CYCLES_STANDALONE)
  list(APPEND INC_SYS ${MATERIALX_INCLUDE_DIRS})
  list(APPEND LIB ${MATERIALX_LIBRARIES} double-conversion::double-conversion)
endif()

if(WITH_CYCLES_STANDALONE AND WITH_CYCLES_STANDALONE_GUI)
  list(APPEND INC_SYS
    ${Epoxy_INCLUDE_DIRS}
    ${SDL2_INCLUDE_DIRS}
  )
  list(APPEND LIB ${Epoxy_LIBRARIES} ${SDL2_LIBRARIES})
endif()

if(WITH_USD)
  # Silence warning from USD headers using deprecated TBB header.
  add_definitions(
    -D__TBB_show_deprecation_message_atomic_H
    -D__TBB_show_deprecation_message_task_H
  )

  list(APPEND INC_SYS
    ${USD_INCLUDE_DIRS}
  )
  list(APPEND LIB
    cycles_hydra
    ${USD_LIBRARIES}
    ${PYTHON_LIBRARIES}
  )
endif()

cycles_external_libraries_append(LIB)

# Common configuration.

include_directories(${INC})
include_directories(SYSTEM ${INC_SYS})

# Application build targets

if(WITH_CYCLES_STANDALONE)
  set(SRC
    cycles_standalone.cpp
    cycles_xml.cpp
    cycles_xml.h
    oiio_output_driver.cpp
    oiio_output_driver.h
    rib_parser/error.h
    rib_parser/intern_cache.h
    rib_parser/parsed_parameter.cpp
    rib_parser/parsed_parameter.h
    rib_parser/param_dict.cpp
    rib_parser/param_dict.h
    rib_parser/parser.cpp
    rib_parser/parser.h
    rib_parser/parallel.cpp
    rib_parser/parallel.h
    rib_parser/ri_api.cpp
    rib_parser/ri_api.h
    rib_parser/exporters/attribute.cpp
    rib_parser/exporters/attribute.h
    rib_parser/exporters/curves.cpp
    rib_parser/exporters/curves.h
    rib_parser/exporters/geometry.cpp
    rib_parser/exporters/geometry.h
    rib_parser/exporters/lights.cpp
    rib_parser/exporters/lights.h
    rib_parser/exporters/node_util.cpp
    rib_parser/exporters/node_util.h
    rib_parser/exporters/materials/convert_lama_network.h
    rib_parser/exporters/materials/convert_lama_network.cpp
    rib_parser/exporters/materials/materials.cpp
    rib_parser/exporters/materials/materials.h
    rib_parser/exporters/materials/rib_to_cycles.cpp
    rib_parser/exporters/materials/rib_to_cycles.h
    rib_parser/exporters/materials/shader_defaults.h
    tev/display_driver.cpp
    tev/display_driver.h
  )

  if(WITH_CYCLES_STANDALONE_GUI)
    list(APPEND SRC
      opengl/display_driver.cpp
      opengl/display_driver.h
      opengl/shader.cpp
      opengl/shader.h
      opengl/window.cpp
      opengl/window.h
    )
  endif()

  add_executable(cycles ${SRC} ${INC} ${INC_SYS})
  unset(SRC)

  target_link_libraries(cycles PRIVATE ${LIB})
  message(MATERIALX_FOUND=${MATERIALX_FOUND})
  if (MATERIALX_FOUND)
    message(MATERIALX_BASE_DIR=${MATERIALX_BASE_DIR})
    add_definitions(-DMATERIALX_BASE_DIR="/Users/mbolsta/test/blender/lib/darwin/materialx_1.38.7")
    add_definitions(-DMATERIALX_INCLUDE_DIR="/Users/mbolsta/test/blender/lib/darwin/materialx_1.38.7/include")
  endif()

  if(APPLE)
    if(WITH_CYCLES_STANDALONE_GUI)
      set_property(TARGET cycles APPEND_STRING PROPERTY LINK_FLAGS
        " -framework Cocoa -framework CoreAudio -framework AudioUnit -framework AudioToolbox -framework ForceFeedback -framework CoreVideo")
    endif()
  set_target_properties(cycles PROPERTIES
    XCODE_GENERATE_SCHEME TRUE
    XCODE_SCHEME_WORKING_DIRECTORY "/Users/mbolsta"
    XCODE_SCHEME_ARGUMENTS "--display-type tev;--display-server 127.0.0.1:14158;--shadingsys osl;--samples 16;--tile-size 64;--threads 0;/Users/mbolsta/osl_test.ris1.0001.rib"
  )
  endif()

  if(UNIX AND NOT APPLE)
    set_target_properties(cycles PROPERTIES INSTALL_RPATH $ORIGIN/lib)
  endif()

  if(CYCLES_STANDALONE_REPOSITORY)
    cycles_install_libraries(cycles)
  endif()

  if(WITH_CYCLES_OSL)
    add_subdirectory(rib_parser/shaders)
  endif()

  install(PROGRAMS
    $<TARGET_FILE:cycles>
    DESTINATION ${CMAKE_INSTALL_PREFIX})
endif()

if(WITH_CYCLES_PRECOMPUTE)
  set(SRC
    cycles_precompute.cpp
  )

  add_executable(cycles_precompute ${SRC} ${INC} ${INC_SYS})
  unset(SRC)

  target_link_libraries(cycles_precompute PRIVATE ${LIB})

  if(UNIX AND NOT APPLE)
    set_target_properties(cycles_precompute PROPERTIES INSTALL_RPATH $ORIGIN/lib)
  endif()

  install(PROGRAMS
    $<TARGET_FILE:cycles_precompute>
    DESTINATION ${CMAKE_INSTALL_PREFIX})
endif()
