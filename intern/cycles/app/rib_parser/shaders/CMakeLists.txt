# SPDX-License-Identifier: Apache-2.0
# Copyright 2011-2022 Blender Foundation

# OSL node shaders

set(SRC_OSO
  PxrAdjustNormal.oso
  PxrArithmetic.oso
  PxrAttribute.oso
  PxrBlackBody.oso
  PxrBlend.oso
  PxrBlenderPrincipledInputs.oso
  PxrBump.oso
  PxrBumpManifold2D.oso
  PxrBumpMixer.oso
  PxrBumpRoughness.oso
  PxrChecker.oso
  PxrClamp.oso
  PxrColorCorrect.oso
  PxrColorGrade.oso
  PxrColorSpace.oso
  PxrCross.oso
  PxrCurvature.oso
  PxrDirt.oso
  PxrDispScalarLayer.oso
  PxrDispTransform.oso
  PxrDispVectorLayer.oso
  PxrDisplace.oso
  PxrDot.oso
  PxrEnvGround.oso
  PxrExposure.oso
  PxrFacingRatio.oso
  PxrFlakes.oso
  PxrFractal.oso
  PxrGamma.oso
  PxrGrid.oso
  PxrHSL.oso
  PxrHairColor.oso
  PxrHexTileManifold.oso
  PxrInvert.oso
  PxrLayer.oso
  PxrLayerMixer.oso
  PxrLayeredBlend.oso
  PxrLayeredTexture.oso
  PxrManifold2D.oso
  PxrManifold3D.oso
  PxrMatteID.oso
  PxrMetallicWorkflow.oso
  PxrMix.oso
  PxrMultiTexture.oso
  PxrPhasorNoise.oso
  PxrPrimvar.oso
  PxrProjectionLayer.oso
  PxrProjectionStack.oso
  PxrProjector.oso
  PxrPtexture.oso
  PxrRadialDensity.oso
  PxrRamp.oso
  PxrRandomTextureManifold.oso
  PxrRemap.oso
  PxrRoundCube.oso
  PxrSetRange.oso
  PxrShadedSide.oso
  PxrStylizedControl.oso
  PxrSwitch.oso
  PxrTangentField.oso
  PxrTee.oso
  PxrTexture.oso
  PxrThinFilm.oso
  PxrThreshold.oso
  PxrTileManifold.oso
  PxrToFloat.oso
  PxrToFloat3.oso
  PxrVariable.oso
  PxrVary.oso
  PxrVoronoise.oso
  PxrWireframe.oso
  PxrWorley.oso
)

set(TGT_OSO

)

# TODO, add a module to compile OSL
foreach(_file ${SRC_OSO})
  set(_OSO_SOURCE_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${_file})
  set_source_files_properties(${_file} PROPERTIES HEADER_FILE_ONLY TRUE)
  string(REPLACE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_BINARY_DIR}/intern/cycles/kernel/osl/shaders _OSO_TARGET_FILE ${_OSO_SOURCE_FILE})
  add_custom_command(
          OUTPUT ${_OSO_TARGET_FILE}
          COMMAND ${CMAKE_COMMAND} -E copy
                  ${_OSO_SOURCE_FILE}
                  ${_OSO_TARGET_FILE})
  list(APPEND TGT_OSO
    ${_OSO_TARGET_FILE}
  )

  unset(_OSO_TARGET_FILE)
  unset(_OSO_SOURCE_FILE)
endforeach()

add_custom_target(cycles_pxr_shaders ALL DEPENDS ${TGT_OSO} SOURCES ${SRC_OSL})
add_dependencies(cycles cycles_pxr_shaders)
cycles_set_solution_folder(cycles_pxr_shaders)

