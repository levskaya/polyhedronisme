// Polyhédronisme
// ===================================================================================================
//
// A toy for constructing and manipulating polyhedra and other meshes
//
// Copyright 2026, Anselm Levskaya
// Released under the MIT License
//

// Polyhedra Triangulation
// ===================================================================================================
//
// Routines for transforming N-face meshes into triangular meshes (via O'Rourke's 2D ear-clipping
// algorithm), necessary for exporting STL or VRML for 3D printing.
//

import { copyVecArray, project2dface } from './geo.js';
import { Polyhedron } from './polyhedron.js';

const xor = (x, y) => Boolean(x) !== Boolean(y);

// Signed twice-area of triangle (va, vb, vc)
const area2 = (va, vb, vc) =>
  ((vb[0] - va[0]) * (vc[1] - va[1])) - ((vc[0] - va[0]) * (vb[1] - va[1]));

const isLeft = (va, vb, vc) => area2(va, vb, vc) > 0;
const isLeftOn = (va, vb, vc) => area2(va, vb, vc) >= 0;
const isCollinear = (va, vb, vc) => area2(va, vb, vc) === 0;

// Proper intersection of segments (va, vb) and (vc, vd)
const intersectProp = function (va, vb, vc, vd) {
  if (isCollinear(va, vb, vc) || isCollinear(va, vb, vd) ||
      isCollinear(vc, vd, va) || isCollinear(vc, vd, vb)) {
    return false;
  }
  return xor(isLeft(va, vb, vc), isLeft(va, vb, vd)) &&
         xor(isLeft(vc, vd, va), isLeft(vc, vd, vb));
};

// Ear-clipping triangulation of a 2D polygon.
// Takes an array of 2D coordinates in face ordering and returns an array of
// local vertex index triplets [i, j, k] preserving the polygon's winding order.
const earClipTriangulate = function (verts) {
  const MAX_ITERS = 999;
  const tris = [];
  const ear = [];
  let facelen = verts.length;

  const inCone = function (a, b) {
    const a1 = (a + 1) % facelen;
    const a0 = (a - 1 + facelen) % facelen;
    if (isLeftOn(verts[a], verts[a1], verts[a0])) {
      return isLeft(verts[a], verts[b], verts[a0]) && isLeft(verts[b], verts[a], verts[a1]);
    }
    return !(isLeftOn(verts[a], verts[b], verts[a1]) && isLeftOn(verts[b], verts[a], verts[a0]));
  };

  const diagonalie = function (a, b) {
    for (let c = 0; c < facelen; c++) {
      const c1 = (c + 1) % facelen;
      if ((c !== a) && (c1 !== a) && (c !== b) && (c1 !== b) &&
          intersectProp(verts[a], verts[b], verts[c], verts[c1])) {
        return false;
      }
    }
    return true;
  };

  const isDiagonal = (a, b) => inCone(a, b) && inCone(b, a) && diagonalie(a, b);

  verts = verts.slice();
  for (let v1 = 0; v1 < facelen; v1++) {
    const v2 = (v1 + 1) % facelen;
    const v0 = (v1 - 1 + facelen) % facelen;
    ear[v1] = isDiagonal(v0, v2);
  }

  const origIdx = Array.from({ length: facelen }, (_, i) => i);
  let n = facelen;
  let outerLimit = MAX_ITERS;
  let head = 0;

  while (outerLimit > 0 && n > 3) {
    outerLimit -= 1;
    let v2 = head;
    let innerLimit = MAX_ITERS;
    while (true) {
      innerLimit -= 1;
      let clipped = false;
      if (ear[v2]) {
        let v3 = (v2 + 1) % facelen;
        let v4 = (v3 + 1) % facelen;
        let v1 = (v2 - 1 + facelen) % facelen;
        let v0 = (v1 - 1 + facelen) % facelen;
        tris.push([origIdx[v1], origIdx[v2], origIdx[v3]]);

        // Remove clipped ear vertex v2 from polygon ring
        verts.splice(v2, 1);
        origIdx.splice(v2, 1);
        ear.splice(v2, 1);
        if (v0 > v2) v0 -= 1;
        if (v1 > v2) v1 -= 1;
        if (v3 > v2) v3 -= 1;
        if (v4 > v2) v4 -= 1;
        facelen--;
        n--;

        // Update ear status of neighboring vertices in the reduced polygon
        ear[v1] = isDiagonal(v0, v3);
        ear[v3] = isDiagonal(v1, v4);
        head = v3;
        clipped = true;
      }
      v2 = (v2 + 1) % facelen;
      if (innerLimit <= 0 || clipped || v2 === head) break;
    }
  }

  // Emit final remaining triangle (or fallback fan if ear-clipping hit limiter on degenerate face)
  for (let k = 1; k < origIdx.length - 1; k++) {
    tris.push([origIdx[0], origIdx[k], origIdx[k + 1]]);
  }

  return tris;
};

// Driver routine: projects each 3D face to 2D, ear-clips into triangles,
// and builds a new Polyhedron preserving original face colors.
export const triangulate = function (poly, colors = false) {
  console.log(`Triangulating faces of ${poly.name}...`);

  const newpoly = new Polyhedron();
  newpoly.vertices = copyVecArray(poly.vertices);
  newpoly.face_classes = [];
  const copyColors = Boolean(colors && poly.face_classes);

  for (let i = 0; i < poly.faces.length; i++) {
    const f = poly.faces[i];
    if (f.length > 3) {
      const face2d = project2dface(f.map(v => poly.vertices[v]));
      const tris = earClipTriangulate(face2d);
      for (let j = 0; j < tris.length; j++) {
        const tri = tris[j];
        newpoly.faces.push([f[tri[0]], f[tri[1]], f[tri[2]]]);
        if (copyColors) newpoly.face_classes.push(poly.face_classes[i]);
      }
    } else {
      newpoly.faces.push([f[0], f[1], f[2]]);
      if (copyColors) newpoly.face_classes.push(poly.face_classes[i]);
    }
  }
  newpoly.name = poly.name; // don't change the name for export
  return newpoly;
};
