# SPDX-FileCopyrightText: 2023 Blender Authors
#
# SPDX-License-Identifier: GPL-2.0-or-later

set(INC
  ../include
  ../../blenkernel
  ../../blenloader
  ../../blentranslation
  ../../bmesh
  ../../gpu
  ../../imbuf
  ../../makesrna
  ../../windowmanager

  # RNA_prototypes.hh
  ${CMAKE_BINARY_DIR}/source/blender/makesrna
)

set(INC_SYS
)

set(SRC
  info_draw.cc
  info_ops.cc
  info_report.cc
  info_stats.cc
  space_info.cc
  textview.cc

  info_intern.hh
  textview.hh
)

set(LIB
  PRIVATE bf::blenfont
  PRIVATE bf::blenlib
  PRIVATE bf::depsgraph
  PRIVATE bf::dna
  PRIVATE bf::intern::guardedalloc
  PRIVATE bf::extern::fmtlib
)


blender_add_lib(bf_editor_space_info "${SRC}" "${INC}" "${INC_SYS}" "${LIB}")

# RNA_prototypes.hh dna_type_offsets.h
add_dependencies(bf_editor_space_info bf_rna)
