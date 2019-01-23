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

import { _ } from 'underscore';
import { add, sub, mult, dot, mag, unit, normal, orthogonal, edgeDist,
  copyVecArray, tangentPoint, calcCentroid, reciprocal, sqrt } from './geo';
import { dual } from './topo_operators';
import { Polyhedron } from './polyhedron';

// ================================================================================================
// Canonicalization Algorithms
// ================================================================================================

// Slow Canonicalization Algorithm
// ----------------------------------------------------------------
//
// This algorithm has some convergence problems, what really needs to be done is to
// sum the three forcing factors together as a conherent force and to use a half-decent
// integrator to make sure that it converges well as opposed to the current hack of
// ad-hoc stability multipliers.  Ideally one would implement a conjugate gradient
// descent or similar pretty thing.
//
// Only try to use this on convex polyhedra that have a chance of being canonicalized,
// otherwise it will probably blow up the geometry.  A much trickier / smarter seed-symmetry
// based geometrical regularizer should be used for fancier/weirder polyhedra.

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
  // const polycenter = [0, 0, 0];
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
    const coords = f.map(v => vertices[v])
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

// combines above three constraint adjustments in iterative cycle
export const canonicalize = function (poly, Niter, keepflat) {
  if (Niter === undefined) { Niter = 1; }
  if (keepflat === undefined) { keepflat = 0; }

  console.log(`Canonicalizing ${poly.name}...`);
  const faces = poly.faces;
  const edges = poly.edges();
  let newVs = poly.vertices;
  let maxChange = 1.0; // convergence tracker
  for (let i = 0; i <= Niter; i++) {
    const oldVs = copyVecArray(newVs); // copy vertices
    if (keepflat === 0) {
      newVs = tangentify(newVs, edges);
    }
    newVs = recenter(newVs, edges);
    newVs = planarize(newVs, faces);
    maxChange = _.max(_.map(_.zip(newVs, oldVs), ([x, y]) => mag(sub(x, y))));
    if (maxChange < 1e-8) {
      break;
    }
  }
  // one should now rescale, but not rescaling here makes for very interesting numerical
  // instabilities that make interesting mutants on multiple applications...
  newVs = rescale(newVs)
  console.log(`[canonicalization done, last |deltaV|=${maxChange}]`);
  const newpoly = new Polyhedron(newVs, poly.faces, poly.name);
  console.log('canonicalize', newpoly);
  return newpoly;
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
