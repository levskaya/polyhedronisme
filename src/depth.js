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

import { abs, add, sub, mult, mag2, linePointDist2, normal, __range__, 
  calcCentroid, calcCentroid2d, round, pointInPolygon, dot2d,
  random, sigfigs, planararea, faceSignature, sqrt, sin, planararea2d,
  cos, PI, pow, unit, mag, dot, sign, tween, tween2d, perspT, magic_depth, cross2d } from './geo';

import { PALETTE, COLOR_METHOD, COLOR_SENSITIVITY, persp_z_max, persp_z_min, persp_ratio, perspective_scale } from './main';

import { Polyhedron }  from './polyhedron';

const polygonClipping = require('polygon-clipping').default;

export const testShape = function () {
  const poly = new Polyhedron();
  poly.name = 'X';
  poly.faces = [];
  poly.vertices = [];
  const h0 = 1.0;
  const r0 = 1.0;
  const theta = PI / 32;
  const n = 3;
  for (let i = 0; i < n; i++) {
    poly.faces.push([3 * i, 3 * i + 1, 3 * i + 2]);
    let h = h0 / pow(i + 1, 0.8);
    let r = r0 / sqrt(i + 1);
    poly.vertices.push(add([0, 0, 0], [i*0.2, 0, 0]));
    poly.vertices.push(add([h, r * cos(i * theta), r * sin(i * theta)], [i*0.2, 0, 0]));
    poly.vertices.push(add([0, r * cos(i * theta), r * sin(i * theta)], [i*0.2, 0, 0]));
  }
  return poly;
};


const simple_intersect = function (A_Vs2d, B_Vs2d) {
  // polygon-clipping expects a multi-polygon for A, B
  try {
    const result = polygonClipping.intersection([A_Vs2d], [B_Vs2d]);
    // multipolygon - 1st polygon, 1st sub-polygon, drop last repeated vertex
    if (result.length > 0 && result[0].length > 0) {
      return result[0][0].slice(0, result[0][0].length);
    } else {
      return [];
    }
  } catch (e) { // library does error-out sometimes... :(
    console.log(e);
    return [];
  }
}


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
    [Math.min.apply(null, ys), Math.max.apply(null, ys)],
  ];
}

const get_2d_vertices = function (tvec) {
  const mapfn = function (vs) {
    return vs.map(v => perspT(add(tvec, v), persp_z_max, persp_z_min, persp_ratio, perspective_scale));
  }
  return mapfn;
}


const project_screen_ray_onto_plane = function (origin, sigma, z0, v2d, bc, bn) {
  // const numerator = sigma * dot(bc, bn) - (z0 + 3) * cross2d(v2d, [bn[0], bn[1]]) - sigma*bn[2]*3;
  const numerator = sigma * dot(bc, bn) - z0 * dot2d(v2d, [bn[0], bn[1]]);
  const denom = dot(bn, [v2d[0], v2d[1], sigma])
  const z = abs(denom) < 1e-6 ? 1e6 : numerator / denom;
  return [(z + z0) / sigma * v2d[0], (z + z0) / sigma * v2d[1], z]
}

export const P_obscures_Q = function (P, Q) {
  const origin = [0, 0, ((persp_z_max * persp_ratio) - persp_z_min) / (1 - persp_ratio)];
  // const magicZ = magic_depth(persp_z_max, persp_z_min, persp_ratio, perspective_scale);
  const sigma = (perspective_scale * persp_ratio) / (1 - persp_ratio);
  const z0 = ((persp_z_max * persp_ratio) - persp_z_min) / (1 - persp_ratio);
  const [ Pc, Pn, Pvs, Pes, Pvs2d, Pes2d, Pm, Pi ] = P;
  const [ Qc, Qn, Qvs, Qes, Qvs2d, Qes2d, Qm, Qi ] = Q;
  const [[Pxmin, Pxmax], [Pymin, Pymax]] = Pes2d;
  const [[Qxmin, Qxmax], [Qymin, Qymax]] = Qes2d;

  if (Qxmax < Pxmin || Pxmax < Qxmin || Qymax < Pymin || Pymax < Qymin) {
    // projected 2d extents don't overlap
    return 0;
  }
  const isect = simple_intersect(Pvs2d, Qvs2d);
  if (isect.length === 0) { // No intersection.
    return 0
  }
  window._traces.push(isect);
  const c2d = calcCentroid2d(isect);
  const A2d = planararea2d(isect);
  const Pdepth = mag2(sub(project_screen_ray_onto_plane(origin, sigma, z0, c2d, Pc, Pn), origin));
  const Qdepth = mag2(sub(project_screen_ray_onto_plane(origin, sigma, z0, c2d, Qc, Qn), origin));
  //console.log('Pi, depth', Pi, Pdepth, 'Qi, depth', Qi, Qdepth, 'A2d', A2d);
  if (A2d < 0.01) {  // if intersection is tiny, ignore...
    return 0;
  }
  // remember that deeper points are more negative, so Pz > Qz if P obscures Q.
  return (Pdepth > Qdepth) ? -1 : 1;
};
// ------------------------------------------------------------------------------------------------------------


export const hardsortfaces = function (poly, tvec) {
  const centroids  = poly.centers();
  const normals    = poly.normals();
  const faceverts  = poly.faces.map(f => f.map(v => poly.vertices[v]));
  const extents = faceverts.map(get_extents);
  const faceverts2d  = faceverts.map(get_2d_vertices(tvec));
  const extents2d = faceverts2d.map(get_2dextents);
  const marks = new Array(poly.faces.length).fill(0);
  const idxs = __range__(0, poly.faces.length, false);
  let alldata = _.zip(centroids, normals, faceverts, extents, faceverts2d, extents2d, marks, idxs);
  let sortdata = [];
  let outdata = [];

  const sigma = (perspective_scale * persp_ratio) / (1 - persp_ratio);
  const z0 = ((persp_z_max * persp_ratio) - persp_z_min) / (1 - persp_ratio);
  const origin = [0, 0, ((persp_z_max * persp_ratio) - persp_z_min) / (1 - persp_ratio)];
  // const magicZ = magic_depth(persp_z_max, persp_z_min, persp_ratio, perspective_scale);
  // let rvec = unit([faceverts2d[0][0][0], faceverts2d[0][0][1], magicZ]);
  // let pt = add(origin, mult(d, rvec));
  //const tvec = [0,0,3];
  let pt = project_screen_ray_onto_plane(origin, sigma, z0, faceverts2d[0][0], centroids[0], normals[0]);
  //console.log('check', faceverts[0][0]);
  //console.log('check', pt);
  let pt2 = project_screen_ray_onto_plane(origin, sigma, z0, faceverts2d[0][1], centroids[0], normals[0]);
  //console.log('check', faceverts[0][1]);
  //console.log('check', pt2);

  // debug  
  window._pts = [];
  window._traces = [];

  // pre-sort
  alldata.sort((a, b) => a[0][2] - b[0][2]);

  // remove hidden faces from consideration
  for (let i = 0; i < alldata.length; i++) {
    if (alldata[i][1][2] < 0) { // hidden face
      outdata.push(alldata[i]);
    } else {
      sortdata.push(alldata[i]);
    }
  }

  // prepare obscuration / visibility graph
  let visgraph = __range__(0, sortdata.length, false).map(i => []);
  let invgraph = __range__(0, sortdata.length, false).map(i => []);
  let rstr=""
  for (let i = 0; i < sortdata.length; i++) {
    for (let j = i + 1; j < sortdata.length; j++) {
      let res = P_obscures_Q(sortdata[i], sortdata[j]);
      if (res === 1) {
        rstr += `${i}~|${j} `;
        visgraph[j].push(i);
        invgraph[i].push(j);
      } else if (res === -1) {
        rstr += `${j}~|${i} `;
        visgraph[i].push(j);
        invgraph[j].push(i);
      }
    }
  }
  //console.log(rstr);

  //console.log('visgraph', visgraph[0], visgraph[1], visgraph[2]);
  //console.log('invgraph', invgraph[0], invgraph[1], invgraph[2]);
  const tmp1 = Object.assign({}, visgraph);
  const tmp2 = Object.assign({}, invgraph);

  // initialize topological sort
  let L = [];
  let S = [];
  for (let i = 0; i < sortdata.length; i++) {
    if (invgraph[i].length === 0) {
      S.push(i);
    }
  }
  const remove_edge = function (n, m) {
    //console.log('remove ', n, m);
    const idxm = visgraph[n].indexOf(m);
    const idxn = invgraph[m].indexOf(n);
    if (idxm > -1) {
      visgraph[n].splice(idxm, 1);
    }
    if (idxn > -1) {
      invgraph[m].splice(idxn, 1);
    }
  };
  // topological sort
  while (S.length > 0) {
    //console.log('S', S);
    let n = S.shift();
    L.push(n);
    //console.log('visgraph[n]', visgraph[n]);
    for (let m of visgraph[n].slice()) { // must work with a copy of data
      //console.log('m', m);
      remove_edge(n, m);
      if (invgraph[m].length === 0) {
        S.push(m);
      }
    }
  }
  //console.log('post-visgraph', visgraph[0], visgraph[1], visgraph[2]);
  //console.log('post-invgraph', invgraph[0], invgraph[1], invgraph[2]);
  // check for cycles
  for (let i = 0; i < sortdata.length; i++) {
    if (visgraph[i].length !== 0 || invgraph[i].length !== 0) {
      console.log('TOPOLOGICAL SORT FAILED, CYCLE DETECTED.');
      //console.log(tmp1);
      //console.log(tmp2);
      return;
    }
  }
  // sortdata = L.map(idx => sortdata[idx]);
  outdata = outdata.concat(L.map(idx => sortdata[idx]));
  const zsortIndex = outdata.map(x => x[7]);
  //console.log(zsortIndex);
  // let fmap = {}
  // for (let i = 0; i < zsortIndex.length; i++) {
  //   fmap[zsortIndex[i]] = i;
  // }
  // window._pts = window._pts.map(x => [x[0], fmap[x[1]], fmap[x[2]], x[3]])

  // sort all face-associated properties
  poly.faces = zsortIndex.map(idx => poly.faces[idx]);
  poly.face_classes = zsortIndex.map(idx => poly.face_classes[idx]);
};





export const __hardsortfaces = function (poly, tvec) {
  const centroids  = poly.centers();
  const normals    = poly.normals();
  const faceverts  = poly.faces.map(f => f.map(v => poly.vertices[v]));
  const extents = faceverts.map(get_extents);
  const faceverts2d  = faceverts.map(get_2d_vertices(tvec));
  const extents2d = faceverts2d.map(get_2dextents);
  const marks = new Array(poly.faces.length).fill(0);
  const idxs = __range__(0, poly.faces.length, false);
  let alldata = _.zip(centroids, normals, faceverts, extents, faceverts2d, extents2d, marks, idxs);
  let indata = [];
  let outdata = [];
  // window._pts = [];
  // pre-sort by minimum (furthest) z-extent of polygons
  // alldata.sort((a, b) => a[3][2][0] - b[3][2][0]);
  alldata.sort((a, b) => a[0][2] - b[0][2]);

  for (let i = 0; i < alldata.length; i++) {
    if (alldata[i][1][2] < 0) { // hidden face
      outdata.push(alldata[i]);
    } else {
      indata.push(alldata[i]);
    }
  }

  outdata = outdata.concat(indata);
  const zsortIndex = outdata.map(x => x[7]);

  // let fmap = {}
  // for (let i = 0; i < zsortIndex.length; i++) {
  //   fmap[zsortIndex[i]] = i;
  // }
  // window._pts = window._pts.map(x => [x[0], fmap[x[1]], fmap[x[2]], x[3]])

  // sort all face-associated properties
  poly.faces = zsortIndex.map(idx => poly.faces[idx]);
  poly.face_classes = zsortIndex.map(idx => poly.face_classes[idx]);
};


export const __P_obscures_Q = function (P, Q) {
  const ray_origin = [0, 0, ((persp_z_max * persp_ratio) - persp_z_min) / (1 - persp_ratio)];
  const [ Pc, Pn, Pvs, Pes, Pvs2d, Pes2d, Pm, Pi ] = P;
  const [ Qc, Qn, Qvs, Qes, Qvs2d, Qes2d, Qm, Qi ] = Q;
  const [[Pxmin, Pxmax], [Pymin, Pymax]] = Pes2d;
  const [[Qxmin, Qxmax], [Qymin, Qymax]] = Qes2d;
  if (false) {
    0;
    // } else if (Qxmax < Pxmin || Pxmax < Qxmin || Qymax < Pymin || Pymax < Qymin) {
    //   // projected 2d extents don't overlap
    //   return false;
  } else {
    let Pvs3d_inset = get_inset_verts(Pvs, 0.01);
    let Pvs2d_inset = get_inset_verts2d(Pvs2d, 0.01);
    let Qvs3d_inset = get_inset_verts(Qvs, 0.01);
    let Qvs2d_inset = get_inset_verts2d(Qvs2d, 0.01);
  
    const PpQds = a_onto_b(ray_origin, Pvs3d_inset, Qc, Qn);
    const QpPds = a_onto_b(ray_origin, Qvs3d_inset, Pc, Pn);
  
    let test = false;
    for (let i = 0; i < Pvs2d.length; i++) {
      // P pt above Q plane, do we intersect in projection?
      if (pointInPolygon(Pvs2d_inset[i], Qvs2d)) {
        window._pts.push([Pvs2d_inset[i], Pi, Qi, PpQds[i]]);
        // return true;
        if (PpQds[i] > 0) test = true;
      }
    }
    for (let i = 0; i < Qvs2d.length; i++) {
      if (pointInPolygon(Qvs2d_inset[i], Pvs2d)) {
        window._pts.push([Qvs2d_inset[i], Qi, Pi, QpPds[i]]);
        // return true;
        if (QpPds[i] < 0) test = true;
      }
    }
    if (test) { return true; }
  }
  return false;
}

const __a_onto_b = function (origin, avs, bc, bn) {
  const numerator = dot(sub(bc, origin), bn);
  let vds = []
  for (let v of avs) {
    const o_v = sub(v, origin);
    const len_o_v = mag(o_v);
    const denom = dot(unit(o_v), bn);
    const len_o_v_to_b_plane = abs(denom) < 1e-6 ? 1e6 : numerator / denom;
    vds.push(len_o_v_to_b_plane - len_o_v);
  }
  return vds;
}
  