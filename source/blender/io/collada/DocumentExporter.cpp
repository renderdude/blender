/* SPDX-FileCopyrightText: 2009-2023 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup collada
 */

#include <cstdio>
#include <cstdlib>

#include "COLLADASWAsset.h"
#include "COLLADASWCamera.h"
#include "COLLADASWException.h"
#include "COLLADASWScene.h"

#include "DNA_object_types.h"
#include "DNA_scene_types.h"

#include "BLI_fileops.h"
#include "BLI_path_utils.hh"
#include "BLI_string.h"

#include "BKE_animsys.h"
#include "BKE_appdir.hh"
#include "BKE_blender_version.h"
#include "BKE_customdata.hh"
#include "BKE_fcurve.hh"

#include "ED_keyframing.hh"
#ifdef WITH_BUILDINFO
extern "C" char build_commit_date[];
extern "C" char build_commit_time[];
extern "C" char build_hash[];
#endif

#include "RNA_access.hh"

#include "DocumentExporter.h"
#include "collada_internal.h"
#include "collada_utils.h"

/* can probably go after refactor is complete */

#include "AnimationExporter.h"
#include "ArmatureExporter.h"
#include "CameraExporter.h"
#include "ControllerExporter.h"
#include "EffectExporter.h"
#include "GeometryExporter.h"
#include "ImageExporter.h"
#include "LightExporter.h"
#include "MaterialExporter.h"
#include "SceneExporter.h"

#include <cerrno>

const char *bc_CustomData_get_layer_name(const CustomData *data, const eCustomDataType type, int n)
{
  int layer_index = CustomData_get_layer_index(data, type);
  if (layer_index < 0) {
    return nullptr;
  }

  return data->layers[layer_index + n].name;
}

const char *bc_CustomData_get_active_layer_name(const CustomData *data, const eCustomDataType type)
{
  /* get the layer index of the active layer of type */
  int layer_index = CustomData_get_active_layer_index(data, type);
  if (layer_index < 0) {
    return nullptr;
  }

  return data->layers[layer_index].name;
}

DocumentExporter::DocumentExporter(BlenderContext &blender_context,
                                   ExportSettings *export_settings)
    : blender_context(blender_context),
      export_settings(BCExportSettings(export_settings, blender_context))
{
}

static COLLADABU::NativeString make_temp_filepath(const char *name, const char *extension)
{
  char tempfile[FILE_MAX];

  if (name == nullptr) {
    name = "untitled";
  }

  BLI_path_join(tempfile, sizeof(tempfile), BKE_tempdir_session(), name);

  if (extension) {
    BLI_path_extension_ensure(tempfile, FILE_MAX, extension);
  }

  COLLADABU::NativeString native_filename = COLLADABU::NativeString(
      tempfile, COLLADABU::NativeString::ENCODING_UTF8);
  return native_filename;
}

/* TODO: it would be better to instantiate animations rather than create a new one per object
 * COLLADA allows this through multiple <channel>s in <animation>.
 * For this to work, we need to know objects that use a certain action. */

int DocumentExporter::exportCurrentScene()
{
  Scene *sce = blender_context.get_scene();
  bContext *C = blender_context.get_context();

  PointerRNA unit_settings;
  PropertyRNA *system; /* unused, *scale; */

  clear_global_id_map();

  COLLADABU::NativeString native_filename = make_temp_filepath(nullptr, ".dae");
  COLLADASW::StreamWriter *writer;
  try {
    writer = new COLLADASW::StreamWriter(native_filename);
  }
  catch (COLLADASW::StreamWriterException &e) {
    e.printMessage();
    fprintf(stderr, "Collada: No Objects will be exported.\n");
    return 1;
  }

  /* open <collada> */
  writer->startDocument();

  /* <asset> */
  COLLADASW::Asset asset(writer);

  PointerRNA sceneptr = RNA_id_pointer_create(&sce->id);
  unit_settings = RNA_pointer_get(&sceneptr, "unit_settings");
  system = RNA_struct_find_property(&unit_settings, "system");
  // scale = RNA_struct_find_property(&unit_settings, "scale_length");

  std::string unitname = "meter";
  float linearmeasure = RNA_float_get(&unit_settings, "scale_length");

  switch (RNA_property_enum_get(&unit_settings, system)) {
    case USER_UNIT_NONE:
    case USER_UNIT_METRIC:
      if (linearmeasure == 0.001f) {
        unitname = "millimeter";
      }
      else if (linearmeasure == 0.01f) {
        unitname = "centimeter";
      }
      else if (linearmeasure == 0.1f) {
        unitname = "decimeter";
      }
      else if (linearmeasure == 1.0f) {
        unitname = "meter";
      }
      else if (linearmeasure == 1000.0f) {
        unitname = "kilometer";
      }
      break;
    case USER_UNIT_IMPERIAL:
      if (linearmeasure == 0.0254f) {
        unitname = "inch";
      }
      else if (linearmeasure == 0.3048f) {
        unitname = "foot";
      }
      else if (linearmeasure == 0.9144f) {
        unitname = "yard";
      }
      break;
    default:
      break;
  }

  asset.setUnit(unitname, linearmeasure);
  asset.setUpAxisType(COLLADASW::Asset::Z_UP);
  asset.getContributor().mAuthor = "Blender User";
  char version_buf[128];
#ifdef WITH_BUILDINFO
  SNPRINTF(version_buf,
           "Blender %s commit date:%s, commit time:%s, hash:%s",
           BKE_blender_version_string(),
           build_commit_date,
           build_commit_time,
           build_hash);
#else
  SNPRINTF(version_buf, "Blender %s", BKE_blender_version_string());
#endif
  asset.getContributor().mAuthoringTool = version_buf;
  asset.add();

  LinkNode *export_set = this->export_settings.get_export_set();
  /* <library_cameras> */
  if (bc_has_object_type(export_set, OB_CAMERA)) {
    CamerasExporter ce(writer, this->export_settings);
    ce.exportCameras(sce);
  }

  /* <library_lights> */
  if (bc_has_object_type(export_set, OB_LAMP)) {
    LightsExporter le(writer, this->export_settings);
    le.exportLights(sce);
  }

  /* <library_effects> */
  EffectsExporter ee(writer, this->export_settings, key_image_map);
  ee.exportEffects(C, sce);

  /* <library_images> */
  ImagesExporter ie(writer, this->export_settings, key_image_map);
  ie.exportImages(sce);

  /* <library_materials> */
  MaterialsExporter me(writer, this->export_settings);
  me.exportMaterials(sce);

  /* <library_geometries> */
  if (bc_has_object_type(export_set, OB_MESH)) {
    GeometryExporter ge(blender_context, writer, this->export_settings);
    ge.exportGeom();
  }

  /* <library_controllers> */
  ArmatureExporter arm_exporter(blender_context, writer, this->export_settings);
  ControllerExporter controller_exporter(blender_context, writer, this->export_settings);
  if (bc_has_object_type(export_set, OB_ARMATURE) || this->export_settings.get_include_shapekeys())
  {
    controller_exporter.export_controllers();
  }

  /* <library_visual_scenes> */

  SceneExporter se(blender_context, writer, &arm_exporter, this->export_settings);

  if (this->export_settings.get_include_animations()) {
    /* <library_animations> */
    AnimationExporter ae(writer, this->export_settings);
    ae.exportAnimations();
  }

  se.exportScene();

  /* <scene> */
  std::string scene_name(translate_id(id_name(sce)));
  COLLADASW::Scene scene(writer, COLLADASW::URI(COLLADABU::Utils::EMPTY_STRING, scene_name));
  scene.add();

  /* close <Collada> */
  writer->endDocument();
  delete writer;

  /* Finally move the created document into place */
  fprintf(stdout, "Collada export to: %s\n", this->export_settings.get_filepath());
  int status = BLI_rename_overwrite(native_filename.c_str(), this->export_settings.get_filepath());
  if (status != 0) {
    status = BLI_copy(native_filename.c_str(), this->export_settings.get_filepath());
    BLI_delete(native_filename.c_str(), false, false);
  }
  return status;
}

/*
 * NOTES:
 *
 * AnimationExporter::sample_animation enables all curves on armature, this is undesirable for a
 * user
 */
