// Polyhédronisme
// ===================================================================================================
//
// A toy for constructing and manipulating polyhedra and other meshes
//
// Copyright 2019, Anselm Levskaya
// Released under the MIT License

/* eslint-disable-camelcase */
/* eslint-disable standard/array-bracket-even-spacing */
/* eslint-disable no-multi-spaces */

import * as polyhedron from './polyhedron';
import * as topo from './topo_operators';
import * as geo from './geo_operators';
// Load peg.js parser for operator-chain-on-base-polyhedra recipes
var opParser = require('./original.pegjs');

// Parser Routines
// ===================================================================================================

const primitiveMap = {
  'T': polyhedron.tetrahedron,
  'O': polyhedron.octahedron,
  'C': polyhedron.cube,
  'I': polyhedron.icosahedron,
  'D': polyhedron.dodecahedron,
  'P': polyhedron.prism,      // takes integer arg
  'A': polyhedron.antiprism,  // takes integer arg
  'Y': polyhedron.pyramid,    // takes integer arg
  'J': polyhedron.johnson,    // takes integer arg
  'U': polyhedron.cupola,     // takes integer arg
  'V': polyhedron.anticupola,  // takes integer arg
  'Q': polyhedron.Qgrid  // takes integer arg
};

const opMap = {
  'd': topo.dual,
  'a': topo.ambo,
  'k': topo.kisN,
  'g': topo.gyro,
  'p': topo.propellor,
  'r': topo.reflect,
  'c': topo.chamfer,
  'w': topo.whirl,
  'n': topo.insetN, // -->needle
  'x': topo.extrudeN,
  'l': topo.loft,
  'P': topo.perspectiva1,
  'q': topo.quinto,
  'u': topo.trisub,
  // 'O': topo.quadsub,
  // z --> topo.zip
  'L': topo.lace,
  'I': topo.joinlace, // change symbol...
  'K': topo.stake,
  'H': topo.hollow,
  // join kis-kis ?
  'Z': topo.triangulate,
  'C': geo.canonicalize,
  'A': geo.adjustXYZ
};
// unclaimed: yihfzv

// List of basic equivalences, easier to replace before parsing.
const specReplacements = [
  [/e/g, 'aa'],   // e --> aa   (abbr. for explode)
  [/b/g, 'ta'],   // b --> ta   (abbr. for bevel)
  [/o/g, 'jj'],   // o --> jj   (abbr. for ortho)
  [/m/g, 'kj'],   // m --> kj   (abbr. for meta)
  [/t(\d*)/g, 'dk$1d'],  // t(n) --> dk(n)d  (dual operations)
  [/j/g, 'dad'],  // j --> dad  (dual operations) # Why not j --> da ?
  [/s/g, 'dgd'],  // s --> dgd  (dual operations) # Why not s --> dg ?
  [/dd/g, ''],    // dd --> null  (order 2)
  [/ad/g, 'a'],   // ad --> a   (a_ = ad_)
  [/gd/g, 'g'],   // gd --> g   (g_ = gd_)
  [/aO$/g, 'aC'],  // aO --> aC  (for uniqueness)
  [/aI$/g, 'aD'],  // aI --> aD  (for uniqueness)
  [/gO$/g, 'gC'],  // gO --> gC  (for uniqueness)
  [/gI$/g, 'gD']];  // gI --> gD  (for uniqueness)

// Execute rewrite rules.
const transformPass = function (notation) {
  let expanded = notation;
  for (let [orig, equiv] of specReplacements) {
    expanded = expanded.replace(orig, equiv);
  }
  console.log(`${notation} executed as ${expanded}`);
  return expanded;
};

// Applies func fn to array args i.e. f, [1,2,3] -> f(1,2,3)
const dispatch = function (fn, args) {
  return fn.apply(this, args || []);
};

// Create polyhedron from notation.
export const generatePoly = function (spec) {
  const transformedSpec = transformPass(spec);
  const opList = opParser.parse(transformedSpec).reverse();

  let op = opList.shift();
  const baseFunc = primitiveMap[op['op']];
  const baseArgs = op['args'];
  let poly = dispatch(baseFunc, baseArgs);

  for (op of opList) {
    const opFunc = opMap[op['op']];
    const opArgs = [poly].concat(op['args']);
    poly = dispatch(opFunc, opArgs);
  }

  // Recenter polyhedra at origin and rescale.
  poly.vertices = geo.recenter(poly.vertices, poly.edges());
  poly.vertices = geo.rescale(poly.vertices);

  // Color the faces of the polyhedra for display.
  poly = polyhedron.paintPolyhedron(poly);

  // Return the poly object.
  return poly;
};
