// Polyhédronisme
// ===================================================================================================
//
// A toy for constructing and manipulating polyhedra and other meshes
//
// Includes implementation of the conway polyhedral operators derived
// from code by mathematician and mathematical sculptor
// George W. Hart http://www.georgehart.com/
//
// Copyright 2026, Anselm Levskaya
// Released under the MIT License

/* eslint-disable camelcase */
/* eslint-disable standard/array-bracket-even-spacing */
/* eslint-disable no-multi-spaces */

import { add, sub, mult, unit, midpoint, oneThird, tween } from './geo.js';
import { Polyhedron } from './polyhedron.js';

// ===================================================================================================
// Polyhedron Flagset Construct
//
// A Flag is an associative triple of a face index and two adjacent vertex vertidxs,
// listed in geometric clockwise order (staring into the normal)
//
// Face_i -> V_i -> V_j
//
// They are a useful abstraction for defining topological transformations of the polyhedral mesh, as
// one can refer to vertices and faces that don't yet exist or haven't been traversed yet in the
// transformation code.
//
// A flag is similar in concept to a directed halfedge in halfedge data structures.
//
export const MAX_FACE_SIDEDNESS = 1000; // GLOBAL

export class Polyflag {
  constructor () {
    this.flags = new Map();    // flags.get(faceName).get(vertName1) = vertName2
    this.vertidxs = new Map(); // vertidxs.get(vertName) = vertex index in output Polyhedron
    this.vertices = new Map(); // vertices.get(vertName) = [x, y, z] coordinates
  }

  // Add a new named vertex if not already present.
  newV (vertName, coordinates) {
    if (!this.vertices.has(vertName)) {
      this.vertices.set(vertName, coordinates);
    }
  }

  // Add a single directed half-edge (v1 -> v2) to faceName.
  newFlag (faceName, vertName1, vertName2) {
    let face = this.flags.get(faceName);
    if (face === undefined) {
      face = new Map();
      this.flags.set(faceName, face);
    }
    face.set(vertName1, vertName2);
  }

  // Add a closed polygon cycle [v0, v1, ..., vk-1] (v0 -> v1 -> ... -> vk-1 -> v0) to faceName.
  newFace (faceName, vertNames) {
    let face = this.flags.get(faceName);
    if (face === undefined) {
      face = new Map();
      this.flags.set(faceName, face);
    }
    const len = vertNames.length;
    let v1 = vertNames[len - 1];
    for (let i = 0; i < len; i++) {
      const v2 = vertNames[i];
      face.set(v1, v2);
      v1 = v2;
    }
  }

  // Add an open directed chain [v0, v1, ..., vk-1] (v0 -> v1 -> ... -> vk-1) to faceName.
  newChain (faceName, vertNames) {
    let face = this.flags.get(faceName);
    if (face === undefined) {
      face = new Map();
      this.flags.set(faceName, face);
    }
    for (let i = 0; i < vertNames.length - 1; i++) {
      face.set(vertNames[i], vertNames[i + 1]);
    }
  }

  // Validate that all vertices are used and every face forms balanced closed cycles in O(V + E).
  check () {
    let bad = false;
    const degrees = new Map();
    for (const vName of this.vertices.keys()) {
      degrees.set(vName, [0, 0]);
    }

    for (const [facename, face] of this.flags) {
      for (const [fromV, toV] of face) {
        const dFrom = degrees.get(fromV);
        const dTo = degrees.get(toV);
        if (!dFrom) {
          console.log(`face ${facename} references undeclared vertex ${fromV}!`);
          bad = true;
        } else {
          dFrom[0] += 1;
        }
        if (!dTo) {
          console.log(`face ${facename} references undeclared vertex ${toV}!`);
          bad = true;
        } else {
          dTo[1] += 1;
        }
        // Within each face map, every target vertex must also be a source vertex
        if (!face.has(toV)) {
          console.log(`face ${facename} vertex named ${toV} unbalanced!`);
          bad = true;
        }
      }
    }

    for (const [vName, deg] of degrees) {
      if (deg[0] === 0 && deg[1] === 0) {
        console.log(`vertex named ${vName} not used!`);
        bad = true;
      } else if (deg[0] !== deg[1]) {
        console.log(`vertex named ${vName} unbalanced! ${deg}`);
        bad = true;
      }
    }

    if (bad) {
      throw new Error('Malformed polyflag!');
    }
  }

  topoly () {
    const poly = new Polyhedron();
    poly.vertices = new Array(this.vertices.size);

    let ctr = 0; // first number the vertices
    for (const [vname, coords] of this.vertices) {
      poly.vertices[ctr] = coords;
      this.vertidxs.set(vname, ctr++);
    }

    poly.faces = new Array(this.flags.size);
    ctr = 0;
    for (const [fname, face] of this.flags) {
      const first = face.values().next();
      if (first.done) { continue; }
      const v0 = first.value; // grab first vertex as starting point
      const polyFace = [];
      // build face out of all the edge relations in the flag map
      let v = v0;
      let faceCTR = 0;
      do {
        const vIdx = this.vertidxs.get(v);
        if (vIdx === undefined) {
          throw new Error(`Undefined vertex "${v}" in face "${fname}"`);
        }
        polyFace.push(vIdx);
        v = face.get(v);
        faceCTR++;
        // necessary during development to prevent browser hangs on badly formed flagsets
        if (faceCTR > MAX_FACE_SIDEDNESS) {
          throw new Error(`Bad flag spec, have a neverending face: ${fname}`);
        }
      } while (v !== v0);
      poly.faces[ctr++] = polyFace;
    }
    if (ctr < poly.faces.length) {
      poly.faces.length = ctr;
    }
    poly.name = 'unknown polyhedron';
    return poly;
  }
}


// ===================================================================================================
// Polyhedron Operators
// ===================================================================================================

// for each vertex of new polyhedron:
//     call newV(Vname, xyz) with a symbolic name and coordinates
// for each flag of new polyhedron:
//     call newFlag(Fname, Vname1, Vname2) with a symbolic name for the new face
//     and the symbolic names for the two vertices forming an oriented edge,
//     or newFace(Fname, [Vname1, Vname2, ...]) for a closed polygon cycle.
// ORIENTATION -must- be dealt with properly to make a manifold (correct) mesh.
// Specifically, no edge v1->v2 can ever be crossed in the -same direction- by
// two different faces - the convention is to define each face with a clockwise set of
// vertices.
//
// call topoly() to assemble flags into polyhedron structure by following the orbits
// of the vertex mapping stored in the flagset for each new face
//
// set name as appropriate

// helper func to ensure unique names of midpoints
export const midName = (v1, v2) => (v1 < v2 ? `${v1}_${v2}` : `${v2}_${v1}`);

// Iterate over directed edges (v1, v2) of a face in clockwise order
const forEachEdge = (face, fn) => {
  const len = face.length;
  let v1 = face[len - 1];
  for (let i = 0; i < len; i++) {
    const v2 = face[i];
    fn(v1, v2);
    v1 = v2;
  }
};

// Iterate over consecutive vertex triplets (v1, v2, v3) of a face in clockwise order
const forEachCorner = (face, fn) => {
  const len = face.length;
  let v1 = face[len - 2];
  let v2 = face[len - 1];
  for (let i = 0; i < len; i++) {
    const v3 = face[i];
    fn(v1, v2, v3);
    v1 = v2;
    v2 = v3;
  }
};


// Kis
// ------------------------------------------------------------------------------------------
// Kis (abbreviated from triakis) transforms an N-sided face into an N-pyramid rooted at the
// same base vertices.
// only kis n-sided faces, but n==0 means kis all.
//
export const kisN = function (poly, n = 0, apexdist = 0.1) {
  if (!n) { n = 0; }
  console.log(`Taking kis of ${n === 0 ? '' : n}-sided faces of ${poly.name}...`);

  const flag = new Polyflag();
  for (let i = 0; i < poly.vertices.length; i++) {
    // each old vertex is a new vertex
    flag.newV(`v${i}`, poly.vertices[i]);
  }

  const normals = poly.normals();
  const centers = poly.centers();
  let foundAny = false;
  for (let i = 0; i < poly.faces.length; i++) {
    const f = poly.faces[i];
    if (n === 0 || f.length === n) {
      foundAny = true;
      const apex = `apex${i}`;
      // new vertex in center of n-sided face
      flag.newV(apex, add(centers[i], mult(apexdist, normals[i])));
      forEachEdge(f, (v1, v2) => {
        flag.newFace(`${i}_v${v1}`, [`v${v1}`, `v${v2}`, apex]);
      });
    } else {
      flag.newFace(`f${i}`, f.map(v => `v${v}`));
    }
  }

  if (!foundAny) {
    console.log(`No ${n}-fold components were found.`);
  }

  const newpoly = flag.topoly();
  newpoly.name = `k${n === 0 ? '' : n}${poly.name}`;
  return newpoly;
};


// Ambo
// ------------------------------------------------------------------------------------------
// The best way to think of the ambo operator is as a topological "tween" between a polyhedron
// and its dual polyhedron.  Thus the ambo of a dual polyhedron is the same as the ambo of the
// original. Also called "Rectify".
//
export const ambo = function (poly) {
  console.log(`Taking ambo of ${poly.name}...`);
  const flag = new Polyflag();

  // For each face f in the original poly
  for (let i = 0; i < poly.faces.length; i++) {
    forEachCorner(poly.faces[i], (v1, v2, v3) => {
      const m12 = midName(v1, v2);
      const m23 = midName(v2, v3);
      if (v1 < v2) { // vertices are the midpoints of all edges of original poly
        flag.newV(m12, midpoint(poly.vertices[v1], poly.vertices[v2]));
      }
      // two new flags:
      // One whose face corresponds to the original f:
      flag.newFlag(`orig${i}`, m12, m23);
      // Another flag whose face corresponds to (the truncated) v2:
      flag.newFlag(`dual${v2}`, m23, m12);
    });
  }

  const newpoly = flag.topoly();
  newpoly.name = `a${poly.name}`;
  return newpoly;
};


// Gyro
// ----------------------------------------------------------------------------------------------
// This is the dual operator to "snub", i.e dual*Gyro = Snub.  It is a bit easier to implement
// this way.
//
// Snub creates at each vertex a new face, expands and twists it, and adds two new triangles to
// replace each edge.
//
export const gyro = function (poly) {
  console.log(`Taking gyro of ${poly.name}...`);

  const flag = new Polyflag();

  // each old vertex is a new vertex
  for (let i = 0; i < poly.vertices.length; i++) {
    flag.newV(`v${i}`, unit(poly.vertices[i]));
  }

  // new vertices in center of each face
  const centers = poly.centers();
  for (let i = 0; i < poly.faces.length; i++) {
    flag.newV(`center${i}`, unit(centers[i]));
  }

  for (let i = 0; i < poly.faces.length; i++) {
    forEachCorner(poly.faces[i], (v1, v2, v3) => {
      flag.newV(`${v1}~${v2}`, oneThird(poly.vertices[v1], poly.vertices[v2])); // new v in face
      // 5-sided kite/pentagon for each directed edge of face i
      flag.newFace(`${i}f${v1}`, [
        `center${i}`,
        `${v1}~${v2}`,
        `${v2}~${v1}`,
        `v${v2}`,
        `${v2}~${v3}`
      ]);
    });
  }

  const newpoly = flag.topoly();
  newpoly.name = `g${poly.name}`;
  return newpoly;
};


// Propellor
// ------------------------------------------------------------------------------------------
// builds a new 'skew face' by making new points along edges, 1/3rd the distance from v1->v2,
// then connecting these into a new inset face.  This breaks rotational symmetry about the
// faces, whirling them into gyres
//
export const propellor = function (poly) {
  console.log(`Taking propellor of ${poly.name}...`);

  const flag = new Polyflag();

  // each old vertex is a new vertex
  for (let i = 0; i < poly.vertices.length; i++) {
    flag.newV(`v${i}`, unit(poly.vertices[i]));
  }

  for (let i = 0; i < poly.faces.length; i++) {
    forEachCorner(poly.faces[i], (v1, v2, v3) => {
      // new v in face, 1/3rd along edge
      flag.newV(`${v1}~${v2}`, oneThird(poly.vertices[v1], poly.vertices[v2]));
      // inset skew face of same sidedness as original face i
      flag.newFlag(`in${i}`, `${v1}~${v2}`, `${v2}~${v3}`);
      // quadrilateral around corner v2
      flag.newFace(`${i}f${v2}`, [
        `${v1}~${v2}`,
        `${v2}~${v1}`,
        `v${v2}`,
        `${v2}~${v3}`
      ]);
    });
  }

  const newpoly = flag.topoly();
  newpoly.name = `p${poly.name}`;
  return newpoly;
};


// Reflection
// ------------------------------------------------------------------------------------------
// geometric reflection through origin
export const reflect = function (poly) {
  console.log(`Taking reflection of ${poly.name}...`);
  const newVertices = poly.vertices.map(v => mult(-1, v));
  const newFaces = poly.faces.map(f => f.slice().reverse());
  return new Polyhedron(newVertices, newFaces, `r${poly.name}`);
};


// Dual
// ------------------------------------------------------------------------------------------------
// The dual of a polyhedron is another mesh wherein:
// - every face in the original becomes a vertex in the dual
// - every vertex in the original becomes a face in the dual
//
// So N_faces, N_vertices = N_dualfaces, N_dualvertices
//
// The new vertex coordinates are convenient to set to the original face centroids.
//
export const dual = function (poly) {
  console.log(`Taking dual of ${poly.name}...`);

  const flag = new Polyflag();

  // Pre-register dual faces in primal vertex order (0 .. V-1) so dpoly.faces[v]
  // corresponds directly to poly.vertices[v] in O(1) without a post-hoc intersection search.
  const edgeToFace = new Array(poly.vertices.length);
  for (let v = 0; v < poly.vertices.length; v++) {
    edgeToFace[v] = new Map();
    flag.flags.set(v, new Map());
  }

  for (let i = 0; i < poly.faces.length; i++) {
    forEachEdge(poly.faces[i], (v1, v2) => {
      // THIS ASSUMES that no 2 faces that share an edge share it in the same orientation!
      // which of course never happens for proper manifold meshes, so get your meshes right.
      edgeToFace[v1].set(v2, i);
    });
  }

  const centers = poly.centers();
  for (let i = 0; i < poly.faces.length; i++) {
    flag.newV(i, centers[i]);
  }

  for (let i = 0; i < poly.faces.length; i++) {
    forEachEdge(poly.faces[i], (v1, v2) => {
      flag.newFlag(v1, edgeToFace[v2].get(v1), i);
    });
  }

  const dpoly = flag.topoly(); // build topological dual from flags
  dpoly.name = poly.name[0] !== 'd' ? `d${poly.name}` : poly.name.slice(1);
  return dpoly;
};


// Chamfer
// ----------------------------------------------------------------------------------------
// A truncation along a polyhedron's edges.
// Chamfering or edge-truncation is similar to expansion, moving faces apart and outward,
// but also maintains the original vertices. Adds a new hexagonal face in place of each
// original edge.
// A polyhedron with e edges will have a chamfered form containing 2e new vertices,
// 3e new edges, and e new hexagonal faces. -- Wikipedia
// See also http://dmccooey.com/polyhedra/Chamfer.html
//
// Q: what is the dual operation of chamfering? I.e.
// if cX = dxdX, and xX = dcdX, what operation is x?
// A: subdivision
//
// We could "almost" do this in terms of already-implemented operations:
// cC = t4daC = t4jC, cO = t3daO, cD = t5daD, cI = t3daI
// But it doesn't work for cases like T.

export const chamfer = function (poly, dist = 0.5) {
  console.log(`Taking chamfer of ${poly.name}...`);
  if (!dist) { dist = 0.5; }

  const flag = new Polyflag();
  const normals = poly.normals();

  // Move each old vertex further from the origin once.
  for (let v = 0; v < poly.vertices.length; v++) {
    flag.newV(`v${v}`, mult(1.0 + dist, poly.vertices[v]));
  }

  // For each face f in the original poly
  for (let i = 0; i < poly.faces.length; i++) {
    const f = poly.faces[i];
    const normalOffset = mult(dist * 1.5, normals[i]);

    for (let v of f) {
      // Add a new vertex, moved parallel to normal.
      flag.newV(`${i}_${v}`, add(poly.vertices[v], normalOffset));
    }

    forEachEdge(f, (v1, v2) => {
      const v1new = `${i}_${v1}`;
      const v2new = `${i}_${v2}`;
      // One flag whose face corresponds to the original face:
      flag.newFlag(`orig${i}`, v1new, v2new);
      // And three flags for the edges of the new hexagon replacing edge (v1, v2):
      flag.newChain(`hex${midName(v1, v2)}`, [`v${v2}`, v2new, v1new, `v${v1}`]);
    });
  }

  const newpoly = flag.topoly();
  newpoly.name = `c${poly.name}`;
  return newpoly;
};


// Whirl
// ----------------------------------------------------------------------------------------------
// Gyro followed by truncation of vertices centered on original faces.
// This creates 2 new hexagons for every original edge.
// (https://en.wikipedia.org/wiki/Conway_polyhedron_notation#Operations_on_polyhedra)
//
export const whirl = function (poly) {
  console.log(`Taking whirl of ${poly.name}...`);

  const flag = new Polyflag();

  // each old vertex is a new vertex
  for (let i = 0; i < poly.vertices.length; i++) {
    flag.newV(`v${i}`, unit(poly.vertices[i]));
  }

  // new vertices around center of each face
  const centers = poly.centers();

  for (let i = 0; i < poly.faces.length; i++) {
    forEachCorner(poly.faces[i], (v1, v2, v3) => {
      // New vertex along edge
      const v1_2 = oneThird(poly.vertices[v1], poly.vertices[v2]);
      flag.newV(`${v1}~${v2}`, v1_2);
      // New vertices near center of face
      const cv1name = `center${i}~${v1}`;
      const cv2name = `center${i}~${v2}`;
      flag.newV(cv1name, unit(oneThird(centers[i], v1_2)));

      // New hexagon for each original directed edge
      flag.newFace(`${i}f${v1}`, [
        cv1name,
        `${v1}~${v2}`,
        `${v2}~${v1}`,
        `v${v2}`,
        `${v2}~${v3}`,
        cv2name
      ]);
      // New face in center of each old face
      flag.newFlag(`c${i}`, cv1name, cv2name);
    });
  }

  const newpoly = flag.topoly();
  newpoly.name = `w${poly.name}`;
  return newpoly;
};


// Quinto
// ----------------------------------------------------------------------------------------------
// This creates a pentagon for every point in the original face, as well as one new inset face.
export const quinto = function (poly) {
  console.log(`Taking quinto of ${poly.name}...`);
  const flag = new Polyflag();
  const centers = poly.centers();

  // For each face f in the original poly
  for (let i = 0; i < poly.faces.length; i++) {
    const centroid = centers[i];
    forEachCorner(poly.faces[i], (v1, v2, v3) => {
      const t12 = midName(v1, v2);
      const t23 = midName(v2, v3);
      const in12 = `inner_${i}_${t12}`;
      const in23 = `inner_${i}_${t23}`;

      // for each face-corner, we make two new points:
      const midpt = midpoint(poly.vertices[v1], poly.vertices[v2]);
      flag.newV(t12, midpt);
      flag.newV(in12, midpoint(midpt, centroid));
      // and add the old corner-vertex
      flag.newV(`v${v2}`, poly.vertices[v2]);

      // pentagon for each vertex in original face
      flag.newFace(`f${i}_${v2}`, [in12, t12, `v${v2}`, t23, in23]);

      // inner rotated face of same vertex-number as original
      flag.newFlag(`f_in_${i}`, in12, in23);
    });
  }

  const newpoly = flag.topoly();
  newpoly.name = `q${poly.name}`;
  return newpoly;
};


// "Loft" operator / inset / extrude
// ------------------------------------------------------------------------------------------
export const insetN = function (poly, n = 0, inset_dist = 0.5, popout_dist = -0.2) {
  if (!n) { n = 0; }
  if (inset_dist === undefined) { inset_dist = 0.5; }
  if (popout_dist === undefined) { popout_dist = -0.2; }

  console.log(`Taking inset of ${n === 0 ? '' : n}-sided faces of ${poly.name}...`);

  const flag = new Polyflag();
  for (let i = 0; i < poly.vertices.length; i++) {
    // each old vertex is a new vertex
    flag.newV(`v${i}`, poly.vertices[i]);
  }

  const normals = poly.normals();
  const centers = poly.centers();
  let foundAny = false; // alert if don't find any

  for (let i = 0; i < poly.faces.length; i++) {
    const f = poly.faces[i];
    if (n === 0 || f.length === n) {
      foundAny = true;
      const popVec = mult(popout_dist, normals[i]);
      for (let v of f) {
        flag.newV(`f${i}v${v}`, add(tween(poly.vertices[v], centers[i], inset_dist), popVec));
      }
      forEachEdge(f, (v1, v2) => {
        flag.newFace(`${i}v${v1}`, [`v${v1}`, `v${v2}`, `f${i}v${v2}`, `f${i}v${v1}`]);
        // new inset, extruded face
        flag.newFlag(`ex${i}`, `f${i}v${v1}`, `f${i}v${v2}`);
      });
    } else {
      flag.newFace(i, f.map(v => `v${v}`)); // same old flag, if non-n
    }
  }

  if (!foundAny) {
    console.log(`No ${n}-fold components were found.`);
  }

  const newpoly = flag.topoly();
  newpoly.name = `n${n === 0 ? '' : n}${poly.name}`;
  return newpoly;
};

// extrudeN
// ------------------------------------------------------------------------------------------
// for compatibility with older operator spec
export const extrudeN = function (poly, n) {
  const newpoly = insetN(poly, n, 0.0, 0.3);
  newpoly.name = `x${n === 0 ? '' : n}${poly.name}`;
  return newpoly;
};

// loft
// ------------------------------------------------------------------------------------------
export const loft = function (poly, n, alpha) {
  const newpoly = insetN(poly, n, alpha, 0.0);
  newpoly.name = `l${n === 0 ? '' : n}${poly.name}`;
  return newpoly;
};


// Lace
// ----------------------------------------------------------------------------------------------
export const lace = function (poly) {
  console.log(`Taking lace of ${poly.name}...`);
  const flag = new Polyflag();
  const centers = poly.centers();

  // For each face f in the original poly
  for (let i = 0; i < poly.faces.length; i++) {
    const centroid = centers[i];
    forEachCorner(poly.faces[i], (v1, v2, v3) => {
      const t12 = midName(v1, v2);
      const t23 = midName(v2, v3);
      const in12 = `inner_${i}_${t12}`;
      const in23 = `inner_${i}_${t23}`;

      // for each face-corner, we make one new point:
      const midpt = midpoint(poly.vertices[v1], poly.vertices[v2]);
      flag.newV(in12, midpoint(midpt, centroid));
      // and add the old corner-vertex
      flag.newV(`v${v2}`, poly.vertices[v2]);

      // triangle at vertex for each vertex in original face
      flag.newFace(`f${i}_${v2}`, [in12, `v${v2}`, in23]);

      // triangle along original edge
      flag.newFace(`f${i}_${t12}`, [in12, `v${v1}`, `v${v2}`]);

      // inner rotated face of same vertex-number as original
      flag.newFlag(`f_in_${i}`, in12, in23);
    });
  }

  const newpoly = flag.topoly();
  newpoly.name = `L${poly.name}`;
  return newpoly;
};


// Join - Lace
// ----------------------------------------------------------------------------------------------
export const joinlace = function (poly) {
  console.log(`Taking joinlace of ${poly.name}...`);
  const flag = new Polyflag();
  const centers = poly.centers();

  // For each face f in the original poly
  for (let i = 0; i < poly.faces.length; i++) {
    const centroid = centers[i];
    forEachCorner(poly.faces[i], (v1, v2, v3) => {
      const t12 = midName(v1, v2);
      const t23 = midName(v2, v3);
      const in12 = `inner_${i}_${t12}`;
      const in23 = `inner_${i}_${t23}`;

      // for each face-corner, we make one new point:
      const midpt = midpoint(poly.vertices[v1], poly.vertices[v2]);
      flag.newV(in12, midpoint(midpt, centroid));
      // and add the old corner-vertex
      flag.newV(`v${v2}`, poly.vertices[v2]);

      // triangle at vertex for each vertex in original face
      flag.newFace(`f${i}_${v2}`, [in12, `v${v2}`, in23]);

      // half of quad joining adjacent faces across original edge (v1, v2)
      flag.newChain(`fj${t12}`, [`v${v2}`, in12, `v${v1}`]);

      // inner rotated face of same vertex-number as original
      flag.newFlag(`f_in_${i}`, in12, in23);
    });
  }

  const newpoly = flag.topoly();
  newpoly.name = `I${poly.name}`;
  return newpoly;
};


// Stake
// ----------------------------------------------------------------------------------------------
export const stake = function (poly) {
  console.log(`Taking stake of ${poly.name}...`);
  const flag = new Polyflag();
  const centers = poly.centers();

  // For each face f in the original poly
  for (let i = 0; i < poly.faces.length; i++) {
    const centroid = centers[i];
    // face centroid is a new point
    flag.newV(`center${i}`, centroid);

    forEachCorner(poly.faces[i], (v1, v2, v3) => {
      const t12 = midName(v1, v2);
      const t23 = midName(v2, v3);
      const in12 = `inner_${i}_${t12}`;
      const in23 = `inner_${i}_${t23}`;

      // for each face-corner, we make one new point:
      const midpt = midpoint(poly.vertices[v1], poly.vertices[v2]);
      flag.newV(in12, midpoint(midpt, centroid));
      // and add the old corner-vertex
      flag.newV(`v${v2}`, poly.vertices[v2]);

      // quad around corner v2 connected to face center
      flag.newFace(`f${i}_${v2}`, [in12, `v${v2}`, in23, `center${i}`]);

      // triangle along original edge
      flag.newFace(`f${i}_${t12}`, [in12, `v${v1}`, `v${v2}`]);
    });
  }

  const newpoly = flag.topoly();
  newpoly.name = `K${poly.name}`;
  return newpoly;
};


// Perspectiva 1
// ------------------------------------------------------------------------------------------
// an operation reverse-engineered from Perspectiva Corporum Regularium
export const perspectiva1 = function (poly) {
  console.log(`Taking stella of ${poly.name}...`);

  const centers = poly.centers(); // calculate face centers
  const flag = new Polyflag();

  for (let i = 0; i < poly.vertices.length; i++) {
    // each old vertex is a new vertex
    flag.newV(`v${i}`, poly.vertices[i]);
  }

  // iterate over triplets of faces v1,v2,v3
  for (let i = 0; i < poly.faces.length; i++) {
    forEachCorner(poly.faces[i], (v1, v2, v3) => {
      const v12 = `v${v1}~v${v2}`; // names for "oriented" midpoints
      const v21 = `v${v2}~v${v1}`;
      const v23 = `v${v2}~v${v3}`;

      // on each Nface, N new points inset from edge midpoints towards center = "stellated" points
      flag.newV(v12, midpoint(midpoint(poly.vertices[v1], poly.vertices[v2]), centers[i]));

      // inset Nface made of new, stellated points
      flag.newFlag(`in${i}`, v12, v23);

      // new tri face constituting the remainder of the stellated Nface
      flag.newFace(`f${i}v${v2}`, [v23, v12, `v${v2}`]);

      // one of the two new triangles replacing old edge between v1->v2
      flag.newFace(`f${v12}`, [`v${v1}`, v21, v12]);
    });
  }

  const newpoly = flag.topoly();
  newpoly.name = `P${poly.name}`;
  return newpoly;
};


// Triangulate - general polygon triangularization algorithm - defined in separate file
// ------------------------------------------------------------------------------------------
export { triangulate } from './triangulate.js';


// Hollow (skeletonize)
// ------------------------------------------------------------------------------------------
export const hollow = function (poly, inset_dist = 0.5, thickness = 0.2) {
  if (inset_dist === undefined) { inset_dist = 0.5; }
  if (thickness === undefined) { thickness = 0.2; }

  console.log(`Hollowing ${poly.name}...`);

  const dualnormals = dual(poly).normals();
  const normals = poly.normals();
  const centers = poly.centers();

  const flag = new Polyflag();
  for (let i = 0; i < poly.vertices.length; i++) {
    // each old vertex is a new vertex
    const p = poly.vertices[i];
    flag.newV(`v${i}`, p);
    flag.newV(`downv${i}`, add(p, mult(-1 * thickness, dualnormals[i])));
  }

  // new inset vertex for every vert in face
  for (let i = 0; i < poly.faces.length; i++) {
    const f = poly.faces[i];
    const downOffset = mult(-1 * thickness, normals[i]);
    for (let v of f) {
      const tw = tween(poly.vertices[v], centers[i], inset_dist);
      flag.newV(`fin${i}v${v}`, tw);
      flag.newV(`findown${i}v${v}`, add(tw, downOffset));
    }

    forEachEdge(f, (v1, v2) => {
      // top rim quad
      flag.newFace(`${i}v${v1}`, [
        `v${v1}`,
        `v${v2}`,
        `fin${i}v${v2}`,
        `fin${i}v${v1}`
      ]);
      // inner wall quad
      flag.newFace(`sides${i}v${v1}`, [
        `fin${i}v${v1}`,
        `fin${i}v${v2}`,
        `findown${i}v${v2}`,
        `findown${i}v${v1}`
      ]);
      // bottom rim quad
      flag.newFace(`bottom${i}v${v1}`, [
        `downv${v2}`,
        `downv${v1}`,
        `findown${i}v${v1}`,
        `findown${i}v${v2}`
      ]);
    });
  }

  const newpoly = flag.topoly();
  newpoly.name = `H${poly.name}`;
  return newpoly;
};


// ===================================================================================================
// Goldberg-Coxeter Operators  (in progress...)
// ===================================================================================================

// Triangular Subdivision Operator
// ----------------------------------------------------------------------------------------------
// Limited version of the Goldberg-Coxeter u_n operator for triangular meshes.
// Uses exact O(F * n^2) topological indexing across shared vertices and edges.
export const trisub = function (poly, n = 2) {
  console.log(`Taking trisub of ${poly.name}...`);
  if (!n) { n = 2; }

  // No-Op for non-triangular meshes.
  for (let fn = 0; fn < poly.faces.length; fn++) {
    if (poly.faces[fn].length !== 3) {
      return poly;
    }
  }

  const uniqVs = [];
  const vertKeyToIdx = new Map();

  // Returns the unique vertex index for barycentric grid point (i, j) on face fn,
  // creating the 3D vertex only on first visit.
  const getOrCreateVertex = (fn, i1, i2, i3, v1, v21, v31, i, j) => {
    const k = n - i - j;
    let key;
    if (i === 0 && j === 0) {
      key = `v_${i1}`;
    } else if (i === n && j === 0) {
      key = `v_${i2}`;
    } else if (i === 0 && j === n) {
      key = `v_${i3}`;
    } else if (j === 0) {
      // On edge (i1, i2) at step i from i1
      key = i1 < i2 ? `e_${i1}_${i2}_${i}` : `e_${i2}_${i1}_${n - i}`;
    } else if (i === 0) {
      // On edge (i1, i3) at step j from i1
      key = i1 < i3 ? `e_${i1}_${i3}_${j}` : `e_${i3}_${i1}_${n - j}`;
    } else if (k === 0) {
      // On edge (i2, i3) at step j from i2
      key = i2 < i3 ? `e_${i2}_${i3}_${j}` : `e_${i3}_${i2}_${n - j}`;
    } else {
      // Strictly interior to face fn
      key = `f_${fn}_${i}_${j}`;
    }

    let idx = vertKeyToIdx.get(key);
    if (idx === undefined) {
      idx = uniqVs.length;
      vertKeyToIdx.set(key, idx);
      uniqVs.push(add(add(v1, mult(i / n, v21)), mult(j / n, v31)));
    }
    return idx;
  };

  const faces = [];
  for (let fn = 0; fn < poly.faces.length; fn++) {
    const [i1, i2, i3] = poly.faces[fn];
    const v1 = poly.vertices[i1];
    const v2 = poly.vertices[i2];
    const v3 = poly.vertices[i3];
    const v21 = sub(v2, v1);
    const v31 = sub(v3, v1);

    // Build local (i, j) lookup table for face fn
    const grid = new Array(n + 1);
    for (let i = 0; i <= n; i++) {
      grid[i] = new Array(n - i + 1);
      for (let j = 0; j + i <= n; j++) {
        grid[i][j] = getOrCreateVertex(fn, i1, i2, i3, v1, v21, v31, i, j);
      }
    }

    for (let i = 0; i < n; i++) {
      for (let j = 0; j + i < n; j++) {
        faces.push([grid[i][j], grid[i + 1][j], grid[i][j + 1]]);
      }
    }
    for (let i = 1; i < n; i++) {
      for (let j = 0; j + i < n; j++) {
        faces.push([grid[i][j], grid[i][j + 1], grid[i - 1][j + 1]]);
      }
    }
  }

  return new Polyhedron(uniqVs, faces, `u${n}${poly.name}`);
};
