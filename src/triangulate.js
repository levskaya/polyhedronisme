// Polyhédronisme
// ===================================================================================================
//
// A toy for constructing and manipulating polyhedra and other meshes
//
// Copyright 2019, Anselm Levskaya
// Released under the MIT License
//

/* eslint-disable camelcase */
/* eslint-disable standard/array-bracket-even-spacing */
/* eslint-disable no-multi-spaces */

// Polyhedra Functions
// ===================================================================================================
//
// Set of routines for transforming N-face meshes into triangular meshes, necessary for exporting
// STL or VRML for 3D Printing.
//

import { copyVecArray, project2dface, __range__ } from './geo.js';
import { Polyhedron } from './polyhedron.js';

// Ear-clipping triangulation of a 2D polygon.
// Takes an array of 2D coordinates in face ordering and returns an array of
// local vertex index triplets [i, j, k] preserving the polygon's winding order.
const earClipTriangulate = function (verts) {
  let v0, v2;
  const limiter = 999;
  const tris = [];
  const ear = [];
  let facelen = verts.length;

  const XOR = (x, y) => (x || y) && !(x && y);
  const Area2 = (Va, Vb, Vc) => ((Vb[0] - Va[0]) * (Vc[1] - Va[1])) - ((Vc[0] - Va[0]) * (Vb[1] - Va[1]));
  const Left = (Va, Vb, Vc) => Area2(Va, Vb, Vc) > 0;
  const LeftOn = (Va, Vb, Vc) => Area2(Va, Vb, Vc) >= 0;
  const Collinear = (Va, Vb, Vc) => Area2(Va, Vb, Vc) === 0;

  const IntersectProp = function (Va, Vb, Vc, Vd) {
    if (Collinear(Va, Vb, Vc) || Collinear(Va, Vb, Vd) ||
        Collinear(Vc, Vd, Va) || Collinear(Vc, Vd, Vb)) {
      return false;
    }
    return XOR(Left(Va, Vb, Vc), Left(Va, Vb, Vd)) && XOR(Left(Vc, Vd, Va), Left(Vc, Vd, Vb));
  };

  const InCone = function (a, b) {
    const a1 = (a + 1 + facelen) % facelen;
    const a0 = ((a - 1) + facelen) % facelen;
    if (LeftOn(verts[a], verts[a1], verts[a0])) {
      return (Left(verts[a], verts[b], verts[a0]) && Left(verts[b], verts[a], verts[a1]));
    }
    return !(LeftOn(verts[a], verts[b], verts[a1]) && LeftOn(verts[b], verts[a], verts[a0]));
  };

  const Diagonalie = function (a, b) {
    let c = 0;
    while (true) {
      const c1 = (c + 1 + facelen) % facelen;
      if ((c !== a) && (c1 !== a) && (c !== b) && (c1 !== b) &&
          IntersectProp(verts[a], verts[b], verts[c], verts[c1])) {
        return false;
      }
      c = (c + 1 + facelen) % facelen;
      if (c === 0) { break; }
    }
    return true;
  };

  const Diagonal = (a, b) => InCone(a, b) && InCone(b, a) && Diagonalie(a, b);

  verts = verts.slice();
  let v1 = 0;
  while (true) {
    v2 = (v1 + 1 + facelen) % facelen; // v1.next
    v0 = ((v1 - 1) + facelen) % facelen; // v1.prev
    ear[v1] = Diagonal(v0, v2);
    v1 = (v1 + 1 + facelen) % facelen;
    if (v1 === 0) { break; }
  }

  const origIdx = __range__(0, facelen - 1, true);
  let n = facelen;
  let z = limiter;
  let head = 0;
  while ((z > 0) && (n > 3)) {
    z -= 1;
    v2 = head;
    let y = limiter;
    while (true) {
      y -= 1;
      let broke = false;
      if (ear[v2]) {
        let v3 = (v2 + 1 + facelen) % facelen; // v2.next
        let v4 = (v3 + 1 + facelen) % facelen; // v3.next
        v1 = ((v2 - 1) + facelen) % facelen; // v2.prev
        v0 = ((v1 - 1) + facelen) % facelen; // v1.prev
        tris.push([ origIdx[v1], origIdx[v2], origIdx[v3] ]);

        // Remove clipped ear vertex v2 from polygon ring
        verts.splice(v2, 1);
        origIdx.splice(v2, 1);
        ear.splice(v2, 1);
        if (v0 > v2) { v0 -= 1; }
        if (v1 > v2) { v1 -= 1; }
        if (v3 > v2) { v3 -= 1; }
        if (v4 > v2) { v4 -= 1; }
        facelen--;
        n--;

        // Update ear status of neighboring vertices in the reduced polygon
        ear[v1] = Diagonal(v0, v3);
        ear[v3] = Diagonal(v1, v4);
        head = v3;
        broke = true;
      }
      v2 = (v2 + 1 + facelen) % facelen; // v2.next
      if ((y <= 0) || broke || (v2 === head)) { break; }
    }
  }

  // Emit final remaining triangle (or fallback fan if ear-clipping hit limiter on degenerate face)
  for (let k = 1; k < origIdx.length - 1; k++) {
    tris.push([ origIdx[0], origIdx[k], origIdx[k + 1] ]);
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
      const twoDface = project2dface(f.map(v => poly.vertices[v]));
      const tris = earClipTriangulate(twoDface);
      for (let j = 0; j < tris.length; j++) {
        const tri = tris[j];
        newpoly.faces.push([ f[tri[0]], f[tri[1]], f[tri[2]] ]);
        if (copyColors) { newpoly.face_classes.push(poly.face_classes[i]); }
      }
    } else {
      newpoly.faces.push([ f[0], f[1], f[2] ]);
      if (copyColors) { newpoly.face_classes.push(poly.face_classes[i]); }
    }
  }
  newpoly.name = poly.name; // don't change the name for export
  return newpoly;
};
