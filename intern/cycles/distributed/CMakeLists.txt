# SPDX-FileCopyrightText: 2011-2022 Blender Foundation
#
# SPDX-License-Identifier: Apache-2.0


set(INC
  ..
)

set(SRC
  distributed.cpp
)

set(SRC_HEADERS
  distributed.h
)

set(LIB
)

list(APPEND INC_SYS
  mpl
)

list(APPEND INC_SYS
  ${ZeroMQ_INCLUDE_DIR}
)

list(APPEND LIB
  ${ZeroMQ_LIBRARIES}
)

list(APPEND INC_SYS
  ${MPI_C_INCLUDE_DIRS}
)

list(APPEND LIB
  ${MPI_C_LIBRARIES}
)

include_directories(${INC})
include_directories(SYSTEM ${INC_SYS})

cycles_add_library(cycles_distributed "${LIB}" ${SRC} ${SRC_HEADERS})
