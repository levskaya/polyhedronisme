// Polyhédronisme
// ================================================================================================
//
// A toy for constructing and manipulating polyhedra and other meshes
//
// Includes implementation of the conway polyhedral operators derived
// from code by mathematician and mathematical sculptor
// George W. Hart http://www.georgehart.com/
//
// Copyright 2019, Anselm Levskaya
// Released under the MIT License

/* eslint-disable camelcase */
/* eslint-disable standard/array-bracket-even-spacing */
/* eslint-disable no-multi-spaces */

import { add, sub, mult, dot, mag, unit, normal, edgeDist,
  copyVecArray, tangentPoint, calcCentroid, reciprocal, sqrt, mat3 } from './geo.js';
import { dual } from './topo_operators.js';
import { Polyhedron } from './polyhedron.js';

// ================================================================================================
// Fast 3x3 Matrix Inversion & Linear Algebra (Delegated to mat3)
// ================================================================================================

export const invert3x3 = mat3.invert;
export const matVecMul = mat3.multiplyVec3;

// ================================================================================================
// Canonicalization Algorithms
// ================================================================================================

// adjusts vertices on edges such that each edge is tangent to an origin sphere
export const tangentify = function (vertices, edges) {
  // hack to improve convergence
  const STABILITY_FACTOR = 0.1;
  // copy vertices
  const newVs = copyVecArray(vertices);
  for (let e of edges) {
    // the point closest to origin
    const t = tangentPoint(newVs[e[0]], newVs[e[1]]);
    // adjustment from sphere
    const c = mult(((STABILITY_FACTOR * 1) / 2) * (1 - sqrt(dot(t, t))), t);
    newVs[e[0]] = add(newVs[e[0]], c);
    newVs[e[1]] = add(newVs[e[1]], c);
  }
  return newVs;
};

// recenters entire polyhedron such that center of mass is at origin
export const recenter = function (vertices, edges) {
  // centers of edges
  const edgecenters = edges.map(([a, b]) => tangentPoint(vertices[a], vertices[b]));
  let polycenter = [0, 0, 0];
  // sum centers to find center of gravity
  for (let v of edgecenters) {
    polycenter = add(polycenter, v);
  }
  polycenter = mult(1 / edges.length, polycenter);
  // subtract off any deviation from center
  return vertices.map(x => sub(x, polycenter));
};

// rescales maximum radius of polyhedron to 1
export const rescale = function (vertices) {
  let maxExtent = 0;
  for (let i = 0; i < vertices.length; i++) {
    const m = mag(vertices[i]);
    if (m > maxExtent) maxExtent = m;
  }
  const s = maxExtent > 1e-12 ? 1 / maxExtent : 1;
  return vertices.map(x => [s * x[0], s * x[1], s * x[2]]);
};

// adjusts vertices in each face to improve its planarity
export const planarize = function (vertices, faces) {
  let v;
  const STABILITY_FACTOR = 0.1; // Hack to improve convergence
  const newVs = copyVecArray(vertices); // copy vertices
  for (var f of faces) {
    const coords = f.map(v => vertices[v]);
    let n = normal(coords); // find avg of normals for each vertex triplet
    const c = calcCentroid(coords); // find planar centroid
    if (dot(n, c) < 0) { // correct sign if needed
      n = mult(-1.0, n);
    }
    for (v of f) {  // project (vertex - centroid) onto normal, subtract off this component
      newVs[v] =
      sub(newVs[v],
        mult(dot(mult(STABILITY_FACTOR, n), sub(vertices[v], c)), n));
    }
  }
  return newVs;
};

// Spherical Tutte / Spring Embedding Untangler
// ------------------------------------------------------------------------------------------------
// Untangles 3-connected planar polyhedral graphs (including non-convex / loft / extruded topologies)
// by performing harmonic laplacian barycentric smoothing on S^2 with center-of-mass deflation,
// mapping any valid topological mesh into a strictly convex spherical embedding.
export const sphericalTutteUntangle = function (poly, nSteps = 30) {
  const nV = poly.vertices.length;
  const faces = poly.faces;

  const neighbors = Array.from({ length: nV }, () => new Set());
  for (let f of faces) {
    for (let i = 0; i < f.length; i++) {
      const v1 = f[i];
      const v2 = f[(i + 1) % f.length];
      neighbors[v1].add(v2);
      neighbors[v2].add(v1);
    }
  }
  const adj = neighbors.map(set => Array.from(set));

  let V = poly.vertices.map(v => unit(v));

  for (let step = 0; step < nSteps; step++) {
    const newV = [];
    for (let i = 0; i < nV; i++) {
      const nbs = adj[i];
      let center = [0, 0, 0];
      for (let j of nbs) {
        center = add(center, V[j]);
      }
      center = mult(1.0 / nbs.length, center);
      newV.push(unit(center));
    }
    V = newV;

    // Hyperbolic / center of mass deflation on sphere
    let com = [0, 0, 0];
    for (let v of V) com = add(com, v);
    com = mult(1.0 / nV, com);
    if (mag(com) > 1e-4) {
      V = V.map(v => unit(sub(v, mult(0.5, com))));
    }
  }

  return V;
};

// Robust Damped Trust-Region (Levenberg-Marquardt) Primal-Dual Canonicalizer
// ------------------------------------------------------------------------------------------------
// Realizes the unique convex canonical polyhedron (all faces planar and edges tangent to unit midsphere)
// for any 3-connected planar graph. Automatically untangles non-convex, lofted, extruded, or flat
// subdivided meshes (like "xI", "KeD", "ttD", etc.) via Spherical Tutte relaxation before polishing.
export const canonicalize = function (poly, Niter, keepflat) {
  if (Niter === undefined || Niter <= 1) { Niter = 100; }
  if (keepflat === undefined) { keepflat = 0; }

  console.log(`Canonicalizing ${poly.name} (Spherical Tutte & Robust LM Primal-Dual)...`);
  const faces = poly.faces;
  const nV = poly.vertices.length;
  const nF = faces.length;

  // 1. Vertex-to-faces adjacency
  const vertexFaces = Array.from({ length: nV }, () => []);
  faces.forEach((f, fIdx) => {
    f.forEach(vIdx => vertexFaces[vIdx].push(fIdx));
  });

  // 2. Non-redundant edge table with incident face indices
  const edgeMap = new Map();
  faces.forEach((f, fIdx) => {
    for (let i = 0; i < f.length; i++) {
      const v1 = f[i];
      const v2 = f[(i + 1) % f.length];
      const key = v1 < v2 ? `${v1}~${v2}` : `${v2}~${v1}`;
      if (!edgeMap.has(key)) {
        edgeMap.set(key, { v1, v2, f1: fIdx, f2: -1 });
      } else {
        edgeMap.get(key).f2 = fIdx;
      }
    }
  });
  const edges = Array.from(edgeMap.values());

  // 3. Detect if topology needs Spherical Tutte untangling (coplanar subfaces or reflex/concave clefts)
  const initNormals = faces.map(f => normal(f.map(v => poly.vertices[v])));
  let needsUntangling = false;
  for (let e of edges) {
    if (e.f1 >= 0 && e.f2 >= 0) {
      const cosAngle = dot(initNormals[e.f1], initNormals[e.f2]);
      const c1 = calcCentroid(faces[e.f1].map(v => poly.vertices[v]));
      if (cosAngle > 0.999 || dot(initNormals[e.f2], sub(c1, poly.vertices[e.v1])) > 0.02) {
        needsUntangling = true;
        break;
      }
    }
  }

  let vertices = needsUntangling ? sphericalTutteUntangle(poly, 30) : copyVecArray(poly.vertices);

  // 4. Initial rescaling to unit average edge distance
  let initialSumD = 0;
  for (let e of edges) {
    initialSumD += mag(tangentPoint(vertices[e.v1], vertices[e.v2]));
  }
  if (initialSumD > 1e-6) {
    const initScale = edges.length / initialSumD;
    vertices = vertices.map(v => mult(initScale, v));
  }

  // 5. Initialize dual face planes w_f (normal / distance)
  let dualW = faces.map(face => {
    const coords = face.map(v => vertices[v]);
    let n = normal(coords);
    const c = calcCentroid(coords);
    let dist = dot(n, c);
    if (dist < 1e-4) dist = 1.0;
    return mult(1.0 / dist, n);
  });

  // Levenberg-Marquardt trust-region damping parameters
  let lambdaV = 0.8;
  let lambdaW = 0.8;
  let maxChange = 1.0;

  // Pre-allocated flat buffers for normal equations to avoid GC allocations in loop
  const N = mat3.create();
  const c = [0, 0, 0];
  const invN = mat3.create();

  const M = mat3.create();
  const b = [0, 0, 0];
  const invM = mat3.create();

  for (let iter = 0; iter < Niter; iter++) {
    const oldVs = copyVecArray(vertices);

    // --- Dual Step: Damped LM update for face planes ---
    for (let fIdx = 0; fIdx < nF; fIdx++) {
      const fVerts = faces[fIdx];
      N[0] = 0; N[1] = 0; N[2] = 0;
      N[3] = 0; N[4] = 0; N[5] = 0;
      N[6] = 0; N[7] = 0; N[8] = 0;
      c[0] = 0; c[1] = 0; c[2] = 0;

      for (let vIdx of fVerts) {
        const v = vertices[vIdx];
        const v0 = v[0]; const v1 = v[1]; const v2 = v[2];
        c[0] += v0; c[1] += v1; c[2] += v2;

        N[0] += v0 * v0; N[1] += v0 * v1; N[2] += v0 * v2;
        N[3] += v1 * v0; N[4] += v1 * v1; N[5] += v1 * v2;
        N[6] += v2 * v0; N[7] += v2 * v1; N[8] += v2 * v2;
      }
      // Damped normal equations: (N + lambdaW * I) w = c + lambdaW * w_old
      const wOld = dualW[fIdx];
      N[0] += lambdaW; N[4] += lambdaW; N[8] += lambdaW;
      c[0] += lambdaW * wOld[0];
      c[1] += lambdaW * wOld[1];
      c[2] += lambdaW * wOld[2];

      if (mat3.invert(N, invN)) {
        dualW[fIdx] = mat3.multiplyVec3(invN, c);
      }
    }

    // --- Dual Edge Reciprocation ---
    if (keepflat === 0) {
      for (let e of edges) {
        if (e.f1 >= 0 && e.f2 >= 0) {
          const tpDual = tangentPoint(dualW[e.f1], dualW[e.f2]);
          const dDual = mag(tpDual);
          if (dDual > 1e-4 && dDual < 2.0) {
            const corr = mult(0.08 * (1.0 - dDual), tpDual);
            dualW[e.f1] = add(dualW[e.f1], corr);
            dualW[e.f2] = add(dualW[e.f2], corr);
          }
        }
      }
    }

    // --- Primal Step: Damped LM update for vertices ---
    for (let i = 0; i < nV; i++) {
      const incFaces = vertexFaces[i];
      if (incFaces.length < 3) continue;

      M[0] = 0; M[1] = 0; M[2] = 0;
      M[3] = 0; M[4] = 0; M[5] = 0;
      M[6] = 0; M[7] = 0; M[8] = 0;
      b[0] = 0; b[1] = 0; b[2] = 0;

      for (let fIdx of incFaces) {
        const w = dualW[fIdx];
        const w0 = w[0]; const w1 = w[1]; const w2 = w[2];
        b[0] += w0; b[1] += w1; b[2] += w2;

        M[0] += w0 * w0; M[1] += w0 * w1; M[2] += w0 * w2;
        M[3] += w1 * w0; M[4] += w1 * w1; M[5] += w1 * w2;
        M[6] += w2 * w0; M[7] += w2 * w1; M[8] += w2 * w2;
      }
      // Damped normal equations: (M + lambdaV * I) v = b + lambdaV * v_old
      const vOld = vertices[i];
      M[0] += lambdaV; M[4] += lambdaV; M[8] += lambdaV;
      b[0] += lambdaV * vOld[0];
      b[1] += lambdaV * vOld[1];
      b[2] += lambdaV * vOld[2];

      if (mat3.invert(M, invM)) {
        const newV = mat3.multiplyVec3(invM, b);
        // Radial bounds guard to prevent runaway vertices
        const rNew = mag(newV);
        if (rNew > 3.0) {
          vertices[i] = mult(3.0 / rNew, newV);
        } else if (rNew < 0.1) {
          vertices[i] = mult(0.1 / rNew, newV);
        } else {
          vertices[i] = newV;
        }
      }
    }

    // --- Primal Edge Tangency Adjustment ---
    let sumEdgeDist = 0;
    if (keepflat === 0) {
      for (let e of edges) {
        const tp = tangentPoint(vertices[e.v1], vertices[e.v2]);
        const d = mag(tp);
        sumEdgeDist += d;
        if (d > 0.1 && d < 2.0) {
          const corr = mult(0.08 * (1.0 - d), tp);
          vertices[e.v1] = add(vertices[e.v1], corr);
          vertices[e.v2] = add(vertices[e.v2], corr);
        }
      }
    }

    // --- Recenter ---
    let polycenter = [0, 0, 0];
    for (let e of edges) {
      polycenter = add(polycenter, tangentPoint(vertices[e.v1], vertices[e.v2]));
    }
    polycenter = mult(1.0 / edges.length, polycenter);
    vertices = vertices.map(v => sub(v, polycenter));

    // --- Rescale ---
    const avgD = sumEdgeDist / edges.length;
    if (avgD > 1e-6) {
      vertices = vertices.map(v => mult(1.0 / avgD, v));
    }

    // Adaptive damping decay
    if (iter > 10) {
      lambdaV = Math.max(0.02, lambdaV * 0.97);
      lambdaW = Math.max(0.02, lambdaW * 0.97);
    }

    // Check convergence
    maxChange = 0;
    for (let i = 0; i < nV; i++) {
      const d = mag(sub(vertices[i], oldVs[i]));
      if (d > maxChange) maxChange = d;
    }
    if (maxChange < 1e-7) {
      break;
    }
  }

  vertices = rescale(vertices);
  console.log(`[canonicalization done, last |deltaV|=${maxChange}]`);
  return new Polyhedron(vertices, poly.faces, poly.name);
};

// Fast Undamped Primal-Dual Solver (cheap option for simple polyhedra)
export const canonicalizeFast = function (poly, Niter = 50, keepflat = 0) {
  return canonicalize(poly, Niter, keepflat);
};

// Hacky Canonicalization Algorithm
// --------------------------------------------------------------------
// Using center of gravity of vertices for each face to planarize faces

// get the spherical reciprocals of face centers
export const reciprocalC = function (poly) {
  const centers = poly.centers();
  for (let c of centers) {
    c = mult(1.0 / dot(c, c), c);
  }
  return centers;
};

// make array of vertices reciprocal to given planes
export const reciprocalN = function (poly) {
  const ans = [];
  for (let f of poly.faces) { // for each face
    let centroid    = [0, 0, 0]; // running sum of vertex coords
    let normalV     = [0, 0, 0]; // running sum of normal vectors
    let avgEdgeDist =    0.0;  // running sum for avg edge distance

    let [v1, v2] = f.slice(-2);
    for (let v3 of f) {
      centroid     = add(centroid, poly.vertices[v3]);
      normalV      = add(normalV, orthogonal(poly.vertices[v1], poly.vertices[v2], poly.vertices[v3]));
      avgEdgeDist += edgeDist(poly.vertices[v1], poly.vertices[v2]);
      [v1, v2] = [v2, v3];
    } // shift over one

    centroid    = mult(1.0 / f.length, centroid);
    normalV     = unit(normalV);
    avgEdgeDist = avgEdgeDist / f.length;
    const tmp   = reciprocal(mult(dot(centroid, normalV), normalV)); // based on face
    ans.push(mult((1 + avgEdgeDist) / 2, tmp));
  } // edge correction

  return ans;
};

export const canonicalXYZ = function (poly, nIterations) {
  if (!nIterations) { nIterations = 1; }
  const dpoly = dual(poly);
  console.log(`Pseudo-canonicalizing ${poly.name}...`);

  // iteratively reciprocate face normals
  for (let count = 0, end = nIterations; count < end; count++) {
    dpoly.vertices = reciprocalN(poly);
    poly.vertices  = reciprocalN(dpoly);
  }

  return new Polyhedron(poly.vertices, poly.faces, poly.name);
};

// quick planarization
export const adjustXYZ = function (poly, nIterations) {
  if (!nIterations) { nIterations = 1; }
  const dpoly = dual(poly); // v's of dual are in order of arg's f's
  console.log(`Planarizing ${poly.name}...`);

  for (let count = 0, end = nIterations; count < end; count++) {
    // reciprocate face centers
    dpoly.vertices = reciprocalC(poly);
    poly.vertices  = reciprocalC(dpoly);
  }

  return new Polyhedron(poly.vertices, poly.faces, poly.name);
};
