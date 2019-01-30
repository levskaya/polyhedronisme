// Polyhédronisme
// ===================================================================================================
//
// A toy for constructing and manipulating polyhedra.
//
// Copyright 2019, Anselm Levskaya
// Released under the MIT License
//

/* eslint-disable camelcase */
/* eslint-disable standard/array-bracket-even-spacing */
/* eslint-disable no-multi-spaces */

import { _ } from 'underscore';
import * as martinez from 'martinez-polygon-clipping';

import { abs, add, sub, mag2, __range__,
  calcCentroid, calcCentroid2d, dot2d, sqrt, sin, planararea2d,
  cos, PI, pow, dot, tween, tween2d, perspT } from './geo';
import { VISIBILITY_ZNORMAL_CUTOFF, persp_z_max, persp_z_min, persp_ratio, perspective_scale } from './main';
import { Polyhedron }  from './polyhedron';
import { TopoSorter } from './toposort';

export const testShape = function () {
  const poly = new Polyhedron();
  poly.name = 'X';
  poly.faces = [];
  poly.vertices = [];
  const h0 = 1.0;
  const r0 = 1.0;
  const theta = PI / 32;
  const n = 5;
  for (let i = 0; i < n; i++) {
    poly.faces.push([3 * i, 3 * i + 1, 3 * i + 2]);
    let h = h0 / pow(i + 1, 0.8);
    let r = r0 / sqrt(i + 1);
    poly.vertices.push(add([0, 0, 0], [i*0.1, 0, 0]));
    poly.vertices.push(add([h, r * cos(i * theta), r * sin(i * theta)], [i*0.1, 0, 0]));
    poly.vertices.push(add([0, r * cos(i * theta), r * sin(i * theta)], [i*0.1, 0, 0]));
  }
  return poly;
};


const simple_intersect = function (A_Vs2d, B_Vs2d) {
  // find 2d intersection of two polygons each specified as list of 2d vertices.
  try {
    // martinez lib expects geojson multi-path polygon which uses repeated final vertices.
    const A_Vs2d_ = A_Vs2d.concat([A_Vs2d[0]]);
    const B_Vs2d_ = B_Vs2d.concat([B_Vs2d[0]]);
    const result = martinez.intersection([A_Vs2d_], [B_Vs2d_]);
    // assume simple intersection structure
    // take 1st polygon, 1st path, drop last repeated vertex
    if (result.length > 0 && result[0].length > 0) {
      return result[0][0].slice(0, result[0][0].length);
    } else {
      return [];
    }
  } catch (e) { // catch rare errors
    console.log(e);
    return [];
  }
}

const simple_difference = function (A_Vs2d, B_Vs2d) {
  // find 2d intersection of two polygons each specified as list of 2d vertices.
  try {
    // martinez lib expects geojson multi-path polygon which uses repeated final vertices.
    const A_Vs2d_ = A_Vs2d.concat([A_Vs2d[0]]);
    const B_Vs2d_ = B_Vs2d.concat([B_Vs2d[0]]);
    const result = martinez.diff([A_Vs2d_], [B_Vs2d_]);
    // assume simple intersection structure
    // take 1st polygon, 1st path, drop last repeated vertex
    if (result.length > 0 && result[0].length > 0) {
      //return result[0][0].slice(0, result[0][0].length);
      console.log(result);
    } else {
      return [];
    }
  } catch (e) { // catch rare errors
    console.log(e);
    return [];
  }
}
window.simple_difference = simple_difference;

const get_inset_verts = function (verts, alpha) {
  let c = calcCentroid(verts);
  return verts.map(v => tween(v, c, alpha))
}
const get_inset_verts2d = function (verts, alpha) {
  let c = calcCentroid2d(verts);
  return verts.map(v => tween2d(v, c, alpha))
}

const get_extents = function (verts) {
  let xs = verts.map(x => x[0]);
  let ys = verts.map(x => x[1]);
  let zs = verts.map(x => x[2]);
  return [
    [Math.min.apply(null, xs), Math.max.apply(null, xs)],
    [Math.min.apply(null, ys), Math.max.apply(null, ys)],
    [Math.min.apply(null, zs), Math.max.apply(null, zs)]
  ];
}

const get_2dextents = function (verts) {
  let xs = verts.map(x => x[0]);
  let ys = verts.map(x => x[1]);
  return [
    [Math.min.apply(null, xs), Math.max.apply(null, xs)],
    [Math.min.apply(null, ys), Math.max.apply(null, ys)]
  ];
};

const get_2d_vertices = function (verts) {
  return verts.map(v => perspT(v, persp_z_max, persp_z_min, persp_ratio, perspective_scale));
};

// Project screen ray onto 3d plane with centroid at plane_c and normal vector plane_n.
const project_screen_ray_onto_plane = function (sigma, z0, screen_coord, plane_c, plane_n) {
  const numerator = sigma * dot(plane_c, plane_n) - z0 * dot2d(screen_coord, [plane_n[0], plane_n[1]]);
  const denom = dot(plane_n, [screen_coord[0], screen_coord[1], sigma])
  const z = abs(denom) < 1e-6 ? 1e6 : numerator / denom;
  return [(z + z0) / sigma * screen_coord[0], (z + z0) / sigma * screen_coord[1], z]
};

// Determine whether polygon P obscures polygon Q.
export const P_obscures_Q = function (P, Q) {
  const AREA_CUTOFF = 0.01;
  // perspective transform constant parameters
  const origin = [0, 0, ((persp_z_max * persp_ratio) - persp_z_min) / (1 - persp_ratio)];
  const sigma = (perspective_scale * persp_ratio) / (1 - persp_ratio);
  const z0 = ((persp_z_max * persp_ratio) - persp_z_min) / (1 - persp_ratio);
  // pull out relevant polygon data - centroids, normals, 2d coords, 2d extents
  const [ Pc, Pn, Pvs2d, Pes2d, Pi ] = P;
  const [ Qc, Qn, Qvs2d, Qes2d, Qi ] = Q;

  // fast check if projected 2d extents overlap
  const [[Pxmin, Pxmax], [Pymin, Pymax]] = Pes2d;
  const [[Qxmin, Qxmax], [Qymin, Qymax]] = Qes2d;
  if (Qxmax < Pxmin || Pxmax < Qxmin || Qymax < Pymin || Pymax < Qymin) {
    return 0;
  }
  // perform more expensive full polygon intersection test
  const isect = simple_intersect(Pvs2d, Qvs2d);
  if (isect.length === 0) {  // no intersection
    return 0;
  }
  // window._traces.push(isect); // debug

  const c2d = calcCentroid2d(isect);  // get 2d centroid of projected 2d intersection area
  const area2d = planararea2d(isect); // measure intersection area in pixel coordinates
  if (area2d < AREA_CUTOFF) {  // if intersection is tiny, ignore
    return 0;
  }
  // measure depth of reverse-projected intersection centroid on both 3d polygonal planes
  const Pdepth = mag2(sub(project_screen_ray_onto_plane(sigma, z0, c2d, Pc, Pn), origin));
  const Qdepth = mag2(sub(project_screen_ray_onto_plane(sigma, z0, c2d, Qc, Qn), origin));
  // return depth ordering as sign
  return (Pdepth > Qdepth) ? -1 : 1;
};


export const sortFacesCarefully = function (poly) {
  const centroids  = poly.centers();
  const normals    = poly.normals();
  const faceverts  = poly.faces.map(f => f.map(v => poly.vertices[v]));
  const faceverts2d  = faceverts.map(get_2d_vertices);
  const extents2d = faceverts2d.map(get_2dextents);
  const idxs = __range__(0, poly.faces.length, false);
  let alldata = _.zip(centroids, normals, faceverts2d, extents2d, idxs);
  const znormal = x => x[1][2];
  const original_face_idx = x => x[4];
  let visible_faces = [];
  let zsortIndex = [];
  // debug
  // window._pts = [];
  // window._traces = [];

  // pre-sort by z-centroid
  alldata.sort((a, b) => a[0][2] - b[0][2]);

  // remove hidden faces from consideration
  for (let i = 0; i < alldata.length; i++) {
    if (znormal(alldata[i]) < 0) { // hidden face
      zsortIndex.push(original_face_idx(alldata[i]));
    } else {
      visible_faces.push(i);
    }
  }

  // prepare obscuration / visibility graph
  let toposorter = new TopoSorter(visible_faces.length);
  for (let i = 0; i < visible_faces.length; i++) {
    for (let j = i + 1; j < visible_faces.length; j++) {
      let fi = visible_faces[i];
      let fj = visible_faces[j];
      let res = P_obscures_Q(alldata[fi], alldata[fj]);
      if (res === 1) {
        toposorter.addEdge(j, i);
      } else if (res === -1) {
        toposorter.addEdge(i, j);
      }
    }
  }
  // get topological sort indices
  let sort_indices = toposorter.sort();
  // the following looks confusing because we have to chase back through
  // two index sets to get the _original_ face indices for the polyhedron
  let original_face_indices = alldata.map(original_face_idx);
  for (let sorted_idx of sort_indices) {
    zsortIndex.push(original_face_indices[visible_faces[sorted_idx]]);
  }

  // sort all face-associated properties according to the topological sort
  poly.faces = zsortIndex.map(idx => poly.faces[idx]);
  poly.face_classes = zsortIndex.map(idx => poly.face_classes[idx]);
  poly.visibility = zsortIndex.map(idx => normals[idx][2] > VISIBILITY_ZNORMAL_CUTOFF); 
};
