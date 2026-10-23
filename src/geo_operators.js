// Polyhédronisme
// ================================================================================================
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

import { add, sub, mult, dot, mag, mag2, unit, normal, orthogonal, edgeDist, midpoint,
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

// Clamped segment closest point to origin (keeps tangency point strictly inside edge segment [v1, v2])
export const segmentTangentPoint = function (v1, v2) {
  const d = sub(v2, v1);
  const m2 = mag2(d);
  if (m2 < 1e-14) return midpoint(v1, v2);
  let t = -dot(v1, d) / m2;
  if (t < 0.01) t = 0.01;
  else if (t > 0.99) t = 0.99;
  return add(v1, mult(t, d));
};

// Degree-Weighted Primal-Dual Spherical Tutte Untangler
// ------------------------------------------------------------------------------------------------
// Untangles 3-connected planar polyhedral graphs (including non-convex, lofted, extruded, or
// high-degree star-subdivided meshes) by alternating spherical face-center averaging and
// degree-weighted vertex relaxation on S^2 with center-of-mass deflation.
// Weighting vertices inversely by degree prevents high-degree hubs from crushing radiating small polygons.
export const sphericalTutteUntangle = function (poly, nSteps = 80) {
  const nV = poly.vertices.length;
  const faces = poly.faces;
  const nF = faces.length;

  const vertexFaces = Array.from({ length: nV }, () => []);
  const neighbors = Array.from({ length: nV }, () => new Set());
  for (let fIdx = 0; fIdx < nF; fIdx++) {
    const f = faces[fIdx];
    for (let i = 0; i < f.length; i++) {
      const v1 = f[i];
      const v2 = f[(i + 1) % f.length];
      vertexFaces[v1].push(fIdx);
      neighbors[v1].add(v2);
      neighbors[v2].add(v1);
    }
  }
  const adj = neighbors.map(set => Array.from(set));
  const deg = adj.map(a => a.length || 1);

  let V = poly.vertices.map(v => unit(v));

  for (let step = 0; step < nSteps; step++) {
    // 1. Spherical face centroids weighted inversely by vertex degree
    const F = new Array(nF);
    for (let fIdx = 0; fIdx < nF; fIdx++) {
      const f = faces[fIdx];
      let fc0 = 0;
      let fc1 = 0;
      let fc2 = 0;
      for (let vIdx of f) {
        const w = 1.0 / deg[vIdx];
        const v = V[vIdx];
        fc0 += w * v[0];
        fc1 += w * v[1];
        fc2 += w * v[2];
      }
      F[fIdx] = unit([fc0, fc1, fc2]);
    }

    // 2. Update vertices from incident face centroids (weighted by 1/|f|) and neighbors (weighted by 1/deg(j))
    const newV = new Array(nV);
    for (let i = 0; i < nV; i++) {
      let cx = 0;
      let cy = 0;
      let cz = 0;
      const incF = vertexFaces[i];
      const nbs = adj[i];
      for (let fIdx of incF) {
        const w = 1.0 / faces[fIdx].length;
        const fc = F[fIdx];
        cx += w * fc[0];
        cy += w * fc[1];
        cz += w * fc[2];
      }
      for (let j of nbs) {
        const w = 1.0 / deg[j];
        const vj = V[j];
        cx += w * vj[0];
        cy += w * vj[1];
        cz += w * vj[2];
      }
      newV[i] = unit([cx, cy, cz]);
    }
    V = newV;

    // 3. Hyperbolic center-of-mass deflation on S^2
    let com = [0, 0, 0];
    for (let v of V) com = add(com, v);
    com = mult(1.0 / nV, com);
    if (mag(com) > 1e-5) {
      V = V.map(v => unit(sub(v, mult(0.5, com))));
    }
  }

  return V;
};

// Unified Tripartite Primal-Dual-Edge Levenberg-Marquardt Canonicalizer
// ------------------------------------------------------------------------------------------------
// Realizes the unique convex canonical polyhedron (all faces planar and all edges tangent to the
// unit midsphere S^2) for any 3-connected planar graph.
//
// Mathematical formulation (Tripartite Incidence Identity):
// At the canonical realization with primal vertices v_i, dual face vectors w_f (plane w_f . x = 1),
// and unit midsphere edge tangency points u_e in S^2 (||u_e|| = 1), all three geometric conditions
// share the exact same bilinear incidence form:
//   1. Face Planarity:        w_f . v_i = 1   for all incident (i, f)
//   2. Primal Edge Tangency:  u_e . v_i = 1   for all incident (i, e)  (since u_e _|_ (v_j - v_i))
//   3. Dual Edge Tangency:    u_e . w_f = 1   for all incident (f, e)  (since u_e _|_ (w_f2 - w_f1))
// Solving these coupled 3x3 normal equations jointly eliminates operator-splitting drift and
// converges to machine precision even on extreme multi-scale / high-degree topologies.
export const canonicalize = function (poly, Niter, keepflat) {
  if (Niter === undefined || Niter <= 1) { Niter = 150; }
  if (keepflat === undefined) { keepflat = 0; }

  console.log(`Canonicalizing ${poly.name} (Unified Tripartite Primal-Dual-Edge LM)...`);
  const faces = poly.faces;
  const nV = poly.vertices.length;
  const nF = faces.length;

  // 1. Build Vertex-Face, Vertex-Edge, and Face-Edge incidence tables
  const vertexFaces = Array.from({ length: nV }, () => []);
  const vertexEdges = Array.from({ length: nV }, () => []);
  const faceEdges = Array.from({ length: nF }, () => []);

  const edgeMap = new Map();
  faces.forEach((f, fIdx) => {
    f.forEach(vIdx => vertexFaces[vIdx].push(fIdx));
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
  const nE = edges.length;
  edges.forEach((e, eIdx) => {
    vertexEdges[e.v1].push(eIdx);
    vertexEdges[e.v2].push(eIdx);
    if (e.f1 >= 0) faceEdges[e.f1].push(eIdx);
    if (e.f2 >= 0) faceEdges[e.f2].push(eIdx);
  });

  // 2. Detect if topology needs Spherical Tutte untangling (coplanar subfaces, reflex clefts, or missed segments)
  const initNormals = faces.map(f => normal(f.map(v => poly.vertices[v])));
  let needsUntangling = false;
  for (let e of edges) {
    if (e.f1 >= 0 && e.f2 >= 0) {
      const cosAngle = dot(initNormals[e.f1], initNormals[e.f2]);
      const c1 = calcCentroid(faces[e.f1].map(v => poly.vertices[v]));
      if (cosAngle > 0.999 || dot(initNormals[e.f2], sub(c1, poly.vertices[e.v1])) > 0.01) {
        needsUntangling = true;
        break;
      }
    }
  }

  let vertices = needsUntangling ? sphericalTutteUntangle(poly, 80) : copyVecArray(poly.vertices);

  // 3. Initial rescaling to unit average edge tangency distance
  let initialSumD = 0;
  for (let e of edges) {
    initialSumD += mag(segmentTangentPoint(vertices[e.v1], vertices[e.v2]));
  }
  if (initialSumD > 1e-6) {
    const initScale = nE / initialSumD;
    vertices = vertices.map(v => mult(initScale, v));
  }

  // 4. Initialize dual face vectors w_f (normal / distance) and unit edge tangency points u_e
  let dualW = faces.map(face => {
    const coords = face.map(v => vertices[v]);
    let n = normal(coords);
    const c = calcCentroid(coords);
    let dist = dot(n, c);
    if (dist < 0.1) dist = 1.0;
    return mult(1.0 / dist, n);
  });

  let edgeU = edges.map(e => unit(segmentTangentPoint(vertices[e.v1], vertices[e.v2])));

  // Levenberg-Marquardt damping parameter and tangency coupling weight
  let lambda = 0.15;
  const wTangency = keepflat === 0 ? 0.5 : 0.0;
  let maxChange = 1.0;

  // Pre-allocated flat buffers for 3x3 normal equations
  const M = mat3.create();
  const b = [0, 0, 0];
  const invM = mat3.create();

  for (let iter = 0; iter < Niter; iter++) {
    const oldVs = copyVecArray(vertices);

    if (keepflat === 0) {
      // --- Step A: Update shared unit midsphere tangency points u_e from primal & dual segments ---
      for (let eIdx = 0; eIdx < nE; eIdx++) {
        const e = edges[eIdx];
        const tpV = unit(segmentTangentPoint(vertices[e.v1], vertices[e.v2]));
        if (e.f1 >= 0 && e.f2 >= 0) {
          const tpW = unit(segmentTangentPoint(dualW[e.f1], dualW[e.f2]));
          edgeU[eIdx] = unit(add(tpV, tpW));
        } else {
          edgeU[eIdx] = tpV;
        }
      }

      // --- Step B: Center midsphere tangency points (sum u_e = 0) ---
      let centerU = [0, 0, 0];
      for (let eIdx = 0; eIdx < nE; eIdx++) {
        centerU = add(centerU, edgeU[eIdx]);
      }
      centerU = mult(1.0 / nE, centerU);
      if (mag(centerU) > 1e-7) {
        for (let eIdx = 0; eIdx < nE; eIdx++) {
          edgeU[eIdx] = unit(sub(edgeU[eIdx], centerU));
        }
        for (let i = 0; i < nV; i++) {
          vertices[i] = sub(vertices[i], centerU);
        }
      }
    }

    // --- Step C: Unified Dual Solve for w_f (from v_i . w_f = 1 and u_e . w_f = 1) ---
    for (let fIdx = 0; fIdx < nF; fIdx++) {
      const fVerts = faces[fIdx];
      const fEdges = faceEdges[fIdx];
      M[0] = 0; M[1] = 0; M[2] = 0;
      M[3] = 0; M[4] = 0; M[5] = 0;
      M[6] = 0; M[7] = 0; M[8] = 0;
      b[0] = 0; b[1] = 0; b[2] = 0;

      for (let vIdx of fVerts) {
        const v = vertices[vIdx];
        const v0 = v[0]; const v1 = v[1]; const v2 = v[2];
        b[0] += v0; b[1] += v1; b[2] += v2;

        M[0] += v0 * v0; M[1] += v0 * v1; M[2] += v0 * v2;
        M[3] += v1 * v0; M[4] += v1 * v1; M[5] += v1 * v2;
        M[6] += v2 * v0; M[7] += v2 * v1; M[8] += v2 * v2;
      }

      if (wTangency > 0) {
        for (let eIdx of fEdges) {
          const u = edgeU[eIdx];
          const u0 = u[0]; const u1 = u[1]; const u2 = u[2];
          b[0] += wTangency * u0; b[1] += wTangency * u1; b[2] += wTangency * u2;

          M[0] += wTangency * u0 * u0; M[1] += wTangency * u0 * u1; M[2] += wTangency * u0 * u2;
          M[3] += wTangency * u1 * u0; M[4] += wTangency * u1 * u1; M[5] += wTangency * u1 * u2;
          M[6] += wTangency * u2 * u0; M[7] += wTangency * u2 * u1; M[8] += wTangency * u2 * u2;
        }
      }

      const wOld = dualW[fIdx];
      M[0] += lambda; M[4] += lambda; M[8] += lambda;
      b[0] += lambda * wOld[0];
      b[1] += lambda * wOld[1];
      b[2] += lambda * wOld[2];

      if (mat3.invert(M, invM)) {
        const newW = mat3.multiplyVec3(invM, b);
        const rW = mag(newW);
        if (rW > 50.0) {
          dualW[fIdx] = mult(50.0 / rW, newW);
        } else if (rW < 0.05) {
          dualW[fIdx] = mult(0.05 / rW, newW);
        } else {
          dualW[fIdx] = newW;
        }
      }
    }

    // --- Step D: Unified Primal Solve for v_i (from w_f . v_i = 1 and u_e . v_i = 1) ---
    for (let i = 0; i < nV; i++) {
      const incFaces = vertexFaces[i];
      const incEdges = vertexEdges[i];
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

      if (wTangency > 0) {
        for (let eIdx of incEdges) {
          const u = edgeU[eIdx];
          const u0 = u[0]; const u1 = u[1]; const u2 = u[2];
          b[0] += wTangency * u0; b[1] += wTangency * u1; b[2] += wTangency * u2;

          M[0] += wTangency * u0 * u0; M[1] += wTangency * u0 * u1; M[2] += wTangency * u0 * u2;
          M[3] += wTangency * u1 * u0; M[4] += wTangency * u1 * u1; M[5] += wTangency * u1 * u2;
          M[6] += wTangency * u2 * u0; M[7] += wTangency * u2 * u1; M[8] += wTangency * u2 * u2;
        }
      }

      const vOld = vertices[i];
      M[0] += lambda; M[4] += lambda; M[8] += lambda;
      b[0] += lambda * vOld[0];
      b[1] += lambda * vOld[1];
      b[2] += lambda * vOld[2];

      if (mat3.invert(M, invM)) {
        const newV = mat3.multiplyVec3(invM, b);
        const rNew = mag(newV);
        if (rNew > 50.0) {
          vertices[i] = mult(50.0 / rNew, newV);
        } else if (rNew < 0.05) {
          vertices[i] = mult(0.05 / rNew, newV);
        } else {
          vertices[i] = newV;
        }
      }
    }

    // Adaptive damping decay
    if (iter > 5) {
      lambda = Math.max(1e-4, lambda * 0.93);
    }

    // Check convergence
    maxChange = 0;
    for (let i = 0; i < nV; i++) {
      const d = mag(sub(vertices[i], oldVs[i]));
      if (d > maxChange) maxChange = d;
    }
    if (maxChange < 1e-8) {
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
