/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2001-2002 NaN Holding BV. All rights reserved. */

/** \file
 * \ingroup bke
 */
#include <string.h> /* for memcpy */

#include "MEM_guardedalloc.h"

#include "DNA_mesh_types.h"
#include "DNA_meshdata_types.h"

#include "BLI_bitmap.h"
#include "BLI_edgehash.h"
#include "BLI_ghash.h"
#include "BLI_math_vector.h"
#include "BLI_utildefines.h"
#include "BLI_utildefines_stack.h"

#include "BKE_customdata.h"
#include "BKE_lib_id.h"
#include "BKE_mesh.h"
#include "BKE_mesh_mapping.h"

#if 0

/**
 * Poly compare with vtargetmap
 * Function used by #BKE_mesh_merge_verts.
 * The function compares poly_source after applying vtargetmap, with poly_target.
 * The two polys are identical if they share the same vertices in the same order,
 * or in reverse order, but starting position loopstart may be different.
 * The function is called with direct_reverse=1 for same order (i.e. same normal),
 * and may be called again with direct_reverse=-1 for reverse order.
 * \return 1 if polys are identical,  0 if polys are different.
 */
static int cddm_poly_compare(const int *corner_verts,
                             const MPoly *mpoly_source,
                             const MPoly *mpoly_target,
                             const int *vtargetmap,
                             const int direct_reverse)
{
  int vert_source, first_vert_source, vert_target;
  int i_loop_source;
  int i_loop_target, i_loop_target_start, i_loop_target_offset, i_loop_target_adjusted;
  bool compare_completed = false;
  bool same_loops = false;

  const int *corner_vert_source, *corner_vert_target;

  BLI_assert(ELEM(direct_reverse, 1, -1));

  i_loop_source = 0;
  corner_vert_source = corner_verts + mpoly_source->loopstart;
  vert_source = *corner_vert_source;

  if (vtargetmap[vert_source] != -1) {
    vert_source = vtargetmap[vert_source];
  }
  else {
    /* All source loop vertices should be mapped */
    BLI_assert(false);
  }

  /* Find same vertex within mpoly_target's loops */
  corner_vert_target = corner_verts + mpoly_target->loopstart;
  for (i_loop_target = 0; i_loop_target < mpoly_target->totloop;
       i_loop_target++, corner_vert_target++) {
    if (*corner_vert_target == vert_source) {
      break;
    }
  }

  /* If same vertex not found, then polys cannot be equal */
  if (i_loop_target >= mpoly_target->totloop) {
    return false;
  }

  /* Now corner_vert_source and m_loop_target have one identical vertex */
  /* corner_vert_source is at position 0, while m_loop_target has advanced to find identical vertex
   */
  /* Go around the loop and check that all vertices match in same order */
  /* Skipping source loops when consecutive source vertices are mapped to same target vertex */

  i_loop_target_start = i_loop_target;
  i_loop_target_offset = 0;
  first_vert_source = vert_source;

  compare_completed = false;
  same_loops = false;

  while (!compare_completed) {

    vert_target = *corner_vert_target;

    /* First advance i_loop_source, until it points to different vertex, after mapping applied */
    do {
      i_loop_source++;

      if (i_loop_source == mpoly_source->totloop) {
        /* End of loops for source, must match end of loop for target. */
        if (i_loop_target_offset == mpoly_target->totloop - 1) {
          compare_completed = true;
          same_loops = true;
          break; /* Polys are identical */
        }

        compare_completed = true;
        same_loops = false;
        break; /* Polys are different */
      }

      corner_vert_source++;
      vert_source = *corner_vert_source;

      if (vtargetmap[vert_source] != -1) {
        vert_source = vtargetmap[vert_source];
      }
      else {
        /* All source loop vertices should be mapped */
        BLI_assert(false);
      }

    } while (vert_source == vert_target);

    if (compare_completed) {
      break;
    }

    /* Now advance i_loop_target as well */
    i_loop_target_offset++;

    if (i_loop_target_offset == mpoly_target->totloop) {
      /* End of loops for target only, that means no match */
      /* except if all remaining source vertices are mapped to first target */
      for (; i_loop_source < mpoly_source->totloop; i_loop_source++, corner_vert_source++) {
        vert_source = vtargetmap[*corner_vert_source];
        if (vert_source != first_vert_source) {
          compare_completed = true;
          same_loops = false;
          break;
        }
      }
      if (!compare_completed) {
        same_loops = true;
      }
      break;
    }

    /* Adjust i_loop_target for cycling around and for direct/reverse order
     * defined by delta = +1 or -1 */
    i_loop_target_adjusted = (i_loop_target_start + direct_reverse * i_loop_target_offset) %
                             mpoly_target->totloop;
    if (i_loop_target_adjusted < 0) {
      i_loop_target_adjusted += mpoly_target->totloop;
    }
    corner_vert_target = corner_verts + mpoly_target->loopstart + i_loop_target_adjusted;
    vert_target = *corner_vert_target;

    if (vert_target != vert_source) {
      same_loops = false; /* Polys are different */
      break;
    }
  }
  return same_loops;
}

/* Utility stuff for using GHash with polys, used by vertex merging. */

typedef struct PolyKey {
  int poly_index; /* index of the MPoly within the derived mesh */
  int totloops;   /* number of loops in the poly */
  uint hash_sum;  /* Sum of all vertices indices */
  uint hash_xor;  /* Xor of all vertices indices */
} PolyKey;

static uint poly_gset_hash_fn(const void *key)
{
  const PolyKey *pk = key;
  return pk->hash_sum;
}

static bool poly_gset_compare_fn(const void *k1, const void *k2)
{
  const PolyKey *pk1 = k1;
  const PolyKey *pk2 = k2;
  if ((pk1->hash_sum == pk2->hash_sum) && (pk1->hash_xor == pk2->hash_xor) &&
      (pk1->totloops == pk2->totloops)) {
    /* Equality - note that this does not mean equality of polys */
    return false;
  }

  return true;
}

#endif

Mesh *BKE_mesh_merge_verts(Mesh *mesh,
                           const int *vtargetmap,
                           const int tot_vtargetmap,
                           const int merge_mode)
{
  return NULL;
#if 0
  /* This was commented out back in 2013, see commit f45d8827bafe6b9eaf9de42f4054e9d84a21955d. */
  // #define USE_LOOPS

  Mesh *result = NULL;

  const int totvert = mesh->totvert;
  const int totedge = mesh->totedge;
  const int totloop = mesh->totloop;
  const int totpoly = mesh->totpoly;
  const MEdge *src_edges = BKE_mesh_edges(mesh);
  const MPoly *src_polys = BKE_mesh_poly_offsets(mesh);
  const int *src_corner_verts = BKE_mesh_corner_verts(mesh);
  const int *src_corner_edges = BKE_mesh_corner_edges(mesh);

  const int totvert_final = totvert - tot_vtargetmap;

  int *oldv = MEM_malloc_arrayN(totvert_final, sizeof(*oldv), __func__);
  int *newv = MEM_malloc_arrayN(totvert, sizeof(*newv), __func__);
  STACK_DECLARE(oldv);

  /* NOTE: create (totedge + totloop) elements because partially invalid polys due to merge may
   * require generating new edges, and while in 99% cases we'll still end with less final edges
   * than totedge, cases can be forged that would end requiring more. */
  const MEdge *med;
  MEdge *medge = MEM_malloc_arrayN((totedge + totloop), sizeof(*medge), __func__);
  int *olde = MEM_malloc_arrayN((totedge + totloop), sizeof(*olde), __func__);
  int *newe = MEM_malloc_arrayN((totedge + totloop), sizeof(*newe), __func__);
  STACK_DECLARE(medge);
  STACK_DECLARE(olde);

  int *corner_verts = MEM_malloc_arrayN(totloop, sizeof(int), __func__);
  int *corner_edges = MEM_malloc_arrayN(totloop, sizeof(int), __func__);
  int *oldl = MEM_malloc_arrayN(totloop, sizeof(*oldl), __func__);
#  ifdef USE_LOOPS
  int *newl = MEM_malloc_arrayN(totloop, sizeof(*newl), __func__);
#  endif
  STACK_DECLARE(corner_verts);
  STACK_DECLARE(corner_edges);
  STACK_DECLARE(oldl);

  const MPoly *mp;
  MPoly *mpoly = MEM_malloc_arrayN(totpoly, sizeof(*medge), __func__);
  int *oldp = MEM_malloc_arrayN(totpoly, sizeof(*oldp), __func__);
  STACK_DECLARE(mpoly);
  STACK_DECLARE(oldp);

  EdgeHash *ehash = BLI_edgehash_new_ex(__func__, totedge);

  int i, j, c;

  PolyKey *poly_keys;
  GSet *poly_gset = NULL;
  MeshElemMap *poly_map = NULL;
  int *poly_map_mem = NULL;

  STACK_INIT(oldv, totvert_final);
  STACK_INIT(olde, totedge);
  STACK_INIT(oldl, totloop);
  STACK_INIT(oldp, totpoly);

  STACK_INIT(medge, totedge);
  STACK_INIT(corner_verts, totloop);
  STACK_INIT(corner_edges, totloop);
  STACK_INIT(mpoly, totpoly);

  /* fill newv with destination vertex indices */
  c = 0;
  for (i = 0; i < totvert; i++) {
    if (vtargetmap[i] == -1) {
      STACK_PUSH(oldv, i);
      newv[i] = c++;
    }
    else {
      /* dummy value */
      newv[i] = 0;
    }
  }

  /* now link target vertices to destination indices */
  for (i = 0; i < totvert; i++) {
    if (vtargetmap[i] != -1) {
      newv[i] = newv[vtargetmap[i]];
    }
  }

  /* Don't remap vertices in cddm->mloop, because we need to know the original
   * indices in order to skip faces with all vertices merged.
   * The "update loop indices..." section further down remaps vertices in mloop.
   */

  /* now go through and fix edges and faces */
  med = src_edges;
  c = 0;
  for (i = 0; i < totedge; i++, med++) {
    const uint v1 = (vtargetmap[med->v1] != -1) ? vtargetmap[med->v1] : med->v1;
    const uint v2 = (vtargetmap[med->v2] != -1) ? vtargetmap[med->v2] : med->v2;
    if (LIKELY(v1 != v2)) {
      void **val_p;

      if (BLI_edgehash_ensure_p(ehash, v1, v2, &val_p)) {
        newe[i] = POINTER_AS_INT(*val_p);
      }
      else {
        STACK_PUSH(olde, i);
        STACK_PUSH(medge, *med);
        newe[i] = c;
        *val_p = POINTER_FROM_INT(c);
        c++;
      }
    }
    else {
      newe[i] = -1;
    }
  }

  if (merge_mode == MESH_MERGE_VERTS_DUMP_IF_EQUAL) {
    /* In this mode, we need to determine,  whenever a poly' vertices are all mapped */
    /* if the targets already make up a poly, in which case the new poly is dropped */
    /* This poly equality check is rather complex.
     * We use a BLI_ghash to speed it up with a first level check */
    PolyKey *mpgh;
    poly_keys = MEM_malloc_arrayN(totpoly, sizeof(PolyKey), __func__);
    poly_gset = BLI_gset_new_ex(poly_gset_hash_fn, poly_gset_compare_fn, __func__, totpoly);
    /* Duplicates allowed because our compare function is not pure equality */
    BLI_gset_flag_set(poly_gset, GHASH_FLAG_ALLOW_DUPES);

    mp = src_polys;
    mpgh = poly_keys;
    for (i = 0; i < totpoly; i++, mp++, mpgh++) {
      mpgh->poly_index = i;
      mpgh->totloops = mp->totloop;
      mpgh->hash_sum = mpgh->hash_xor = 0;
      for (j = 0; j < mp->totloop; j++) {
        const int vert_i = src_corner_verts[mp->loopstart + j];
        mpgh->hash_sum += vert_i;
        mpgh->hash_xor ^= vert_i;
      }
      BLI_gset_insert(poly_gset, mpgh);
    }

    /* Can we optimize by reusing an old `pmap`? How do we know an old `pmap` is stale? */
    /* When called by `MOD_array.c` the `cddm` has just been created, so it has no valid `pmap`. */
    BKE_mesh_vert_poly_map_create(
        &poly_map, &poly_map_mem, src_polys, src_corner_verts, totvert, totpoly, totloop);
  } /* done preparing for fast poly compare */

  BLI_bitmap *vert_tag = BLI_BITMAP_NEW(mesh->totvert, __func__);

  mp = src_polys;
  for (i = 0; i < totpoly; i++, mp++) {
    MPoly *mp_new;

    /* check faces with all vertices merged */
    bool all_verts_merged = true;

    for (j = 0; j < mp->totloop; j++) {
      const int vert_i = src_corner_verts[mp->loopstart + j];
      if (vtargetmap[vert_i] == -1) {
        all_verts_merged = false;
        /* This will be used to check for poly using several time the same vert. */
        BLI_BITMAP_DISABLE(vert_tag, vert_i);
      }
      else {
        /* This will be used to check for poly using several time the same vert. */
        BLI_BITMAP_DISABLE(vert_tag, vtargetmap[vert_i]);
      }
    }

    if (UNLIKELY(all_verts_merged)) {
      if (merge_mode == MESH_MERGE_VERTS_DUMP_IF_MAPPED) {
        /* In this mode, all vertices merged is enough to dump face */
        continue;
      }
      if (merge_mode == MESH_MERGE_VERTS_DUMP_IF_EQUAL) {
        /* Additional condition for face dump:  target vertices must make up an identical face.
         * The test has 2 steps:
         * 1) first step is fast `ghash` lookup, but not fail-proof.
         * 2) second step is thorough but more costly poly compare. */
        int i_poly, v_target;
        bool found = false;
        PolyKey pkey;

        /* Use poly_gset for fast (although not 100% certain) identification of same poly */
        /* First, make up a poly_summary structure */
        pkey.hash_sum = pkey.hash_xor = 0;
        pkey.totloops = 0;
        for (j = 0; j < mp->totloop; j++) {
          const int vert_i = src_corner_verts[mp->loopstart + j];
          v_target = vtargetmap[vert_i]; /* Cannot be -1, they are all mapped */
          pkey.hash_sum += v_target;
          pkey.hash_xor ^= v_target;
          pkey.totloops++;
        }
        if (BLI_gset_haskey(poly_gset, &pkey)) {

          /* There might be a poly that matches this one.
           * We could just leave it there and say there is, and do a "continue".
           * ... but we are checking whether there is an exact poly match.
           * It's not so costly in terms of CPU since it's very rare, just a lot of complex code.
           */

          /* Consider current loop again */
          /* Consider the target of the loop's first vert */
          v_target = vtargetmap[src_corner_verts[mp->loopstart]];
          /* Now see if v_target belongs to a poly that shares all vertices with source poly,
           * in same order, or reverse order */

          for (i_poly = 0; i_poly < poly_map[v_target].count; i_poly++) {
            const MPoly *target_poly = src_polys + *(poly_map[v_target].indices + i_poly);

            if (cddm_poly_compare(corner_verts, mp, target_poly, vtargetmap, +1) ||
                cddm_poly_compare(corner_verts, mp, target_poly, vtargetmap, -1)) {
              found = true;
              break;
            }
          }
          if (found) {
            /* Current poly's vertices are mapped to a poly that is strictly identical */
            /* Current poly is dumped */
            continue;
          }
        }
      }
    }

    /* Here either the poly's vertices were not all merged
     * or they were all merged, but targets do not make up an identical poly,
     * the poly is retained.
     */
    c = 0;
    int *last_valid_corner_vert = NULL;
    int *last_valid_corner_edge = NULL;
    int *first_valid_corner_vert = NULL;
    int *first_valid_corner_edge = NULL;
    bool need_edge_from_last_valid_ml = false;
    bool need_edge_to_first_valid_ml = false;
    int created_edges = 0;
    for (j = 0; j < mp->totloop; j++) {
      const int orig_vert_i = src_corner_verts[mp->loopstart + j];
      const int orig_edge_i = src_corner_edges[mp->loopstart + j];
      const uint mlv = (vtargetmap[orig_vert_i] != -1) ? vtargetmap[orig_vert_i] : orig_vert_i;
#  ifndef NDEBUG
      {
        const int next_corner_vert = src_corner_verts[mp->loopstart + ((j + 1) % mp->totloop)];
        uint next_mlv = (vtargetmap[next_corner_vert] != -1) ? vtargetmap[next_corner_vert] :
                                                               next_corner_vert;
        med = src_edges + orig_edge_i;
        uint v1 = (vtargetmap[med->v1] != -1) ? vtargetmap[med->v1] : med->v1;
        uint v2 = (vtargetmap[med->v2] != -1) ? vtargetmap[med->v2] : med->v2;
        BLI_assert((mlv == v1 && next_mlv == v2) || (mlv == v2 && next_mlv == v1));
      }
#  endif
      /* A loop is only valid if its matching edge is,
       * and it's not reusing a vertex already used by this poly. */
      if (LIKELY((newe[orig_edge_i] != -1) && !BLI_BITMAP_TEST(vert_tag, mlv))) {
        BLI_BITMAP_ENABLE(vert_tag, mlv);

        if (UNLIKELY(last_valid_corner_vert != NULL && need_edge_from_last_valid_ml)) {
          /* We need to create a new edge between last valid loop and this one! */
          void **val_p;

          uint v1 = (vtargetmap[*last_valid_corner_vert] != -1) ?
                        vtargetmap[*last_valid_corner_vert] :
                        *last_valid_corner_vert;
          uint v2 = mlv;
          BLI_assert(v1 != v2);
          if (BLI_edgehash_ensure_p(ehash, v1, v2, &val_p)) {
            *last_valid_corner_edge = POINTER_AS_INT(*val_p);
          }
          else {
            const int new_eidx = STACK_SIZE(medge);
            STACK_PUSH(olde, olde[*last_valid_corner_edge]);
            STACK_PUSH(medge, src_edges[*last_valid_corner_edge]);
            medge[new_eidx].v1 = *last_valid_corner_vert;
            medge[new_eidx].v2 = orig_vert_i;
            /* DO NOT change newe mapping,
             * could break actual values due to some deleted original edges. */
            *val_p = POINTER_FROM_INT(new_eidx);
            created_edges++;

            *last_valid_corner_edge = new_eidx;
          }
          need_edge_from_last_valid_ml = false;
        }

#  ifdef USE_LOOPS
        newl[j + mp->loopstart] = STACK_SIZE(corner_verts);
#  endif
        STACK_PUSH(oldl, j + mp->loopstart);
        last_valid_corner_vert = STACK_PUSH_RET_PTR(corner_verts);
        last_valid_corner_edge = STACK_PUSH_RET_PTR(corner_edges);
        *last_valid_corner_vert = orig_vert_i;
        *last_valid_corner_edge = orig_edge_i;
        if (first_valid_corner_vert == NULL) {
          first_valid_corner_vert = last_valid_corner_vert;
        }
        if (first_valid_corner_edge == NULL) {
          first_valid_corner_edge = last_valid_corner_edge;
        }
        c++;

        /* We absolutely HAVE to handle edge index remapping here, otherwise potential newly
         * created edges in that part of code make remapping later totally unreliable. */
        BLI_assert(newe[orig_edge_i] != -1);
        *last_valid_corner_edge = newe[orig_edge_i];
      }
      else {
        if (last_valid_corner_vert != NULL) {
          need_edge_from_last_valid_ml = true;
        }
        else {
          need_edge_to_first_valid_ml = true;
        }
      }
    }
    if (UNLIKELY(last_valid_corner_vert != NULL &&
                 !ELEM(first_valid_corner_vert, NULL, last_valid_corner_vert) &&
                 (need_edge_to_first_valid_ml || need_edge_from_last_valid_ml))) {
      /* We need to create a new edge between last valid loop and first valid one! */
      void **val_p;

      uint v1 = (vtargetmap[*last_valid_corner_vert] != -1) ? vtargetmap[*last_valid_corner_vert] :
                                                              *last_valid_corner_vert;
      uint v2 = (vtargetmap[*first_valid_corner_vert] != -1) ?
                    vtargetmap[*first_valid_corner_vert] :
                    *first_valid_corner_vert;
      BLI_assert(v1 != v2);
      if (BLI_edgehash_ensure_p(ehash, v1, v2, &val_p)) {
        *last_valid_corner_edge = POINTER_AS_INT(*val_p);
      }
      else {
        const int new_eidx = STACK_SIZE(medge);
        STACK_PUSH(olde, olde[*last_valid_corner_edge]);
        STACK_PUSH(medge, src_edges[*last_valid_corner_edge]);
        medge[new_eidx].v1 = *last_valid_corner_vert;
        medge[new_eidx].v2 = *first_valid_corner_vert;
        /* DO NOT change newe mapping,
         * could break actual values due to some deleted original edges. */
        *val_p = POINTER_FROM_INT(new_eidx);
        created_edges++;

        *last_valid_corner_edge = new_eidx;
      }
      need_edge_to_first_valid_ml = need_edge_from_last_valid_ml = false;
    }

    if (UNLIKELY(c == 0)) {
      BLI_assert(created_edges == 0);
      continue;
    }
    if (UNLIKELY(c < 3)) {
      STACK_DISCARD(oldl, c);
      STACK_DISCARD(corner_verts, c);
      STACK_DISCARD(corner_edges, c);
      if (created_edges > 0) {
        for (j = STACK_SIZE(medge) - created_edges; j < STACK_SIZE(medge); j++) {
          BLI_edgehash_remove(ehash, medge[j].v1, medge[j].v2, NULL);
        }
        STACK_DISCARD(olde, created_edges);
        STACK_DISCARD(medge, created_edges);
      }
      continue;
    }

    mp_new = STACK_PUSH_RET_PTR(mpoly);
    *mp_new = *mp;
    mp_new->totloop = c;
    BLI_assert(mp_new->totloop >= 3);
    mp_new->loopstart = STACK_SIZE(corner_verts) - c;

    STACK_PUSH(oldp, i);
  } /* End of the loop that tests polys. */

  if (poly_gset) {
    // printf("hash quality %.6f\n", BLI_gset_calc_quality(poly_gset));

    BLI_gset_free(poly_gset, NULL);
    MEM_freeN(poly_keys);
  }

  /* Create new cddm. */
  result = BKE_mesh_new_nomain_from_template(
      mesh, totvert_final, STACK_SIZE(medge), 0, STACK_SIZE(corner_verts), STACK_SIZE(mpoly));

  /* Update edge indices and copy customdata. */
  MEdge *new_med = medge;
  for (i = 0; i < result->totedge; i++, new_med++) {
    BLI_assert(newv[new_med->v1] != -1);
    new_med->v1 = newv[new_med->v1];
    BLI_assert(newv[new_med->v2] != -1);
    new_med->v2 = newv[new_med->v2];

    /* Can happen in case vtargetmap contains some double chains, we do not support that. */
    BLI_assert(new_med->v1 != new_med->v2);

    CustomData_copy_data(&mesh->edata, &result->edata, olde[i], i, 1);
  }

  /* Update loop indices and copy customdata. */
  for (i = 0; i < result->totloop; i++) {
    /* Edge remapping has already be done in main loop handling part above. */
    BLI_assert(newv[corner_verts[i]] != -1);
    corner_verts[i] = newv[corner_verts[i]];

    CustomData_copy_data(&mesh->ldata, &result->ldata, oldl[i], i, 1);
  }

  /* Copy vertex customdata. */
  for (i = 0; i < result->totvert; i++) {
    CustomData_copy_data(&mesh->vdata, &result->vdata, oldv[i], i, 1);
  }

  /* Copy poly customdata. */
  mp = mpoly;
  for (i = 0; i < result->totpoly; i++, mp++) {
    CustomData_copy_data(&mesh->pdata, &result->pdata, oldp[i], i, 1);
  }

  /* Copy over data. #CustomData_add_layer can do this, need to look it up. */
  if (STACK_SIZE(medge)) {
    memcpy(BKE_mesh_edges_for_write(result), medge, sizeof(MEdge) * STACK_SIZE(medge));
  }
  if (STACK_SIZE(corner_verts)) {
    memcpy(BKE_mesh_corner_verts_for_write(result),
           corner_verts,
           sizeof(int) * STACK_SIZE(corner_verts));
  }
  if (STACK_SIZE(corner_edges)) {
    memcpy(BKE_mesh_corner_edges_for_write(result),
           corner_edges,
           sizeof(int) * STACK_SIZE(corner_edges));
  }
  if (STACK_SIZE(mpoly)) {
    memcpy(BKE_mesh_poly_offsets_for_write(result), mpoly, sizeof(MPoly) * STACK_SIZE(mpoly));
  }

  MEM_freeN(medge);
  MEM_freeN(corner_verts);
  MEM_freeN(corner_edges);
  MEM_freeN(mpoly);

  MEM_freeN(newv);
  MEM_freeN(newe);
#  ifdef USE_LOOPS
  MEM_freeN(newl);
#  endif

  MEM_freeN(oldv);
  MEM_freeN(olde);
  MEM_freeN(oldl);
  MEM_freeN(oldp);

  MEM_freeN(vert_tag);

  BLI_edgehash_free(ehash, NULL);

  if (poly_map != NULL) {
    MEM_freeN(poly_map);
  }
  if (poly_map_mem != NULL) {
    MEM_freeN(poly_map_mem);
  }

  BKE_id_free(NULL, mesh);

  return result;
#endif
}
