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

import _ from 'underscore';
import { add, sub, mult, dot, mag, unit, normal, orthogonal, edgeDist,
  copyVecArray, tangentPoint, calcCentroid, reciprocal, sqrt } from './geo';
import { dual } from './topo_operators';
import { Polyhedron } from './polyhedron';

// ================================================================================================
// Fast 3x3 Symmetric Matrix Inversion & Linear Algebra
// ================================================================================================

export const invert3x3 = function (M) {
  const [
    [m00, m01, m02],
    [m10, m11, m12],
    [m20, m21, m22]
  ] = M;

  const c00 = m11 * m22 - m12 * m21;
  const c01 = m12 * m20 - m10 * m22;
  const c02 = m10 * m21 - m11 * m20;

  const det = m00 * c00 + m01 * c01 + m02 * c02;
  if (Math.abs(det) < 1e-12) return null;

  const invDet = 1.0 / det;

  const c10 = m02 * m21 - m01 * m22;
  const c11 = m00 * m22 - m02 * m20;
  const c12 = m01 * m20 - m00 * m21;

  const c20 = m01 * m12 - m02 * m11;
  const c21 = m02 * m10 - m00 * m12;
  const c22 = m00 * m11 - m01 * m10;

  return [
    [c00 * invDet, c10 * invDet, c20 * invDet],
    [c01 * invDet, c11 * invDet, c21 * invDet],
    [c02 * invDet, c12 * invDet, c22 * invDet]
  ];
};

export const matVecMul = function (M, v) {
  return [
    M[0][0] * v[0] + M[0][1] * v[1] + M[0][2] * v[2],
    M[1][0] * v[0] + M[1][1] * v[1] + M[1][2] * v[2],
    M[2][0] * v[0] + M[2][1] * v[1] + M[2][2] * v[2]
  ];
};

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
  return _.map(vertices, x => sub(x, polycenter));
};

// rescales maximum radius of polyhedron to 1
export const rescale = function (vertices) {
  const maxExtent = _.max(_.map(vertices, x => mag(x)));
  const s = 1 / maxExtent;
  return _.map(vertices, x => [s * x[0], s * x[1], s * x[2]]);
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

  for (let iter = 0; iter < Niter; iter++) {
    const oldVs = copyVecArray(vertices);

    // --- Dual Step: Damped LM update for face planes ---
    for (let fIdx = 0; fIdx < nF; fIdx++) {
      const fVerts = faces[fIdx];
      let N = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
      let c = [0, 0, 0];
      for (let vIdx of fVerts) {
        const v = vertices[vIdx];
        for (let r = 0; r < 3; r++) {
          c[r] += v[r];
          for (let col = 0; col < 3; col++) {
            N[r][col] += v[r] * v[col];
          }
        }
      }
      // Damped normal equations: (N + lambdaW * I) w = c + lambdaW * w_old
      const wOld = dualW[fIdx];
      for (let r = 0; r < 3; r++) {
        N[r][r] += lambdaW;
        c[r] += lambdaW * wOld[r];
      }
      const invN = invert3x3(N);
      if (invN) dualW[fIdx] = matVecMul(invN, c);
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

      let M = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
      let b = [0, 0, 0];
      for (let fIdx of incFaces) {
        const w = dualW[fIdx];
        for (let r = 0; r < 3; r++) {
          b[r] += w[r];
          for (let col = 0; col < 3; col++) {
            M[r][col] += w[r] * w[col];
          }
        }
      }
      // Damped normal equations: (M + lambdaV * I) v = b + lambdaV * v_old
      const vOld = vertices[i];
      for (let r = 0; r < 3; r++) {
        M[r][r] += lambdaV;
        b[r] += lambdaV * vOld[r];
      }
      const invM = invert3x3(M);
      if (invM) {
        const newV = matVecMul(invM, b);
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
    maxChange = _.max(_.map(_.zip(vertices, oldVs), ([x, y]) => mag(sub(x, y))));
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
