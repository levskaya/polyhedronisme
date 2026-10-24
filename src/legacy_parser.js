// Polyhédronisme
// ===================================================================================================
//
// A toy for constructing and manipulating polyhedra and other meshes
//
// Copyright 2026, Anselm Levskaya
// Released under the MIT License
//
// Legacy (v1) Parser & Recipe Rewriter for backward compatibility with old ?recipe= URLs.

/* eslint-disable camelcase */
/* eslint-disable standard/array-bracket-even-spacing */
/* eslint-disable no-multi-spaces */

import * as polyhedron from './polyhedron.js';
import * as topo from './topo_operators.js';
import * as geo from './geo_operators.js';

// Load peg.js parser for legacy v1 operator-chain-on-base-polyhedra recipes
let legacyOpParser = null;

export const setLegacyOpParser = function (parser) {
  legacyOpParser = parser;
};

// Legacy v1 -> v2 Recipe Rewriter
// ===================================================================================================

const legacyToV2OpSwaps = {
  'l': 'P', // old stellation 'l' -> 'P' (perspectiva1); 'l' is now loft in v2
  'z': 'Z', // old triangulate 'z' -> 'Z'
  'h': 'H'  // old hollow 'h' -> 'H'
};

// Format parsed opspec list back into recipe notation string
export const formatOps = function (opList) {
  return opList.map(({ op, args }) => {
    if (!args || args.length === 0) {
      return op;
    }
    if (args.length === 1 && args[0] >= 0 && Number.isInteger(args[0])) {
      return `${op}${args[0]}`;
    }
    return `${op}(${args.join(',')})`;
  }).join('');
};

// Rewrites a legacy v1 recipe string (from ?recipe=...) into modern v2 recipe syntax:
// - Old broken 'KN' canonicalizer -> 'C' (drops N, uses 1 round of modern canonicalizer)
// - Old slow 'CN' canonicalizer   -> 'C' (drops N, uses 1 round of modern canonicalizer)
// - Old 'l' (stellaN)             -> 'P' (perspectiva1)
// - Old 'z' (triangulate)         -> 'Z'
// - Old 'h(n,inset,depth)' / 'h'  -> 'H(inset,depth)' / 'H'
export const rewriteLegacyRecipe = function (spec) {
  let ops;
  try {
    ops = legacyOpParser.parse(spec);
  } catch (e) {
    console.warn(`Failed to parse legacy recipe "${spec}" via grammar, falling back to regex rewrite:`, e);
    return spec
      .replace(/h\(\d+,([^)]+)\)/g, 'H($1)')
      .replace(/h/g, 'H')
      .replace(/l/g, 'P')
      .replace(/z/g, 'Z')
      .replace(/[KC]\d+(?=.+)/g, 'C');
  }

  if (!ops || ops.length === 0) {
    return spec;
  }

  // Last element is the base polyhedron seed (e.g. T, C, O, D, I, P5, A4, Y5, J37)
  const baseOp = ops[ops.length - 1];
  const prefixOps = ops.slice(0, -1).map(({ op, args }) => {
    // 1. Old 'KN' and 'CN' geometric canonicalizers -> single round of modern 'C'
    if (op === 'K' || op === 'C') {
      return { op: 'C', args: [] };
    }
    // 2. Old 'h(0,inset,depth)' hollow operator -> 'H(inset,depth)' (drop leading n arg)
    if (op === 'h') {
      if (args.length >= 3) {
        return { op: 'H', args: args.slice(1) };
      }
      return { op: 'H', args };
    }
    // 3. Old 'u(m,0)' -> 'uM' shorthand while preserving general 'u(m,n)' and 'uM'
    if (op === 'u' && args.length === 2 && args[1] === 0 && args[0] > 0) {
      return { op: 'u', args: [args[0]] };
    }
    // 4. Single-letter renamed operators ('l' -> 'P', 'z' -> 'Z')
    if (op in legacyToV2OpSwaps) {
      return { op: legacyToV2OpSwaps[op], args };
    }
    return { op, args };
  });

  const rewritten = formatOps([...prefixOps, baseOp]);
  console.log(`Legacy recipe "${spec}" rewritten to v2 "${rewritten}"`);
  return rewritten;
};


// Standalone Legacy (v1) Parser & Executor
// ===================================================================================================

const legacyPrimitiveMap = {
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
  'V': polyhedron.anticupola  // takes integer arg
};

// Wrapper for legacy hollow(poly, n, inset_dist, thickness) signature
const legacyHollow = function (poly, n, inset_dist, thickness) {
  if (thickness !== undefined) {
    return topo.hollow(poly, inset_dist, thickness);
  }
  return topo.hollow(poly, inset_dist);
};

// Wrapper for legacy C/K canonicalization using modern canonicalize (1 round)
const legacyCanonicalize = function (poly) {
  return geo.canonicalize(poly);
};

const legacyOpMap = {
  'd': topo.dual,
  'a': topo.ambo,
  'k': topo.kisN,
  'g': topo.gyro,
  'p': topo.propellor,
  'r': topo.reflect,
  'h': legacyHollow,
  'c': topo.chamfer,
  'w': topo.whirl,
  'n': topo.insetN,
  'x': topo.extrudeN,
  'l': topo.perspectiva1, // legacy stellaN
  'z': topo.triangulate,
  'K': legacyCanonicalize,
  'C': legacyCanonicalize,
  'A': geo.adjustXYZ,
  'u': topo.trisub,
  'q': topo.quinto
};

const legacySpecReplacements = [
  [/e/g, 'aa'],   // e --> aa   (abbr. for explode)
  [/b/g, 'ta'],   // b --> ta   (abbr. for bevel)
  [/o/g, 'jj'],   // o --> jj   (abbr. for ortho)
  [/m/g, 'kj'],   // m --> kj   (abbr. for meta)
  [/t(\d*)/g, 'dk$1d'],  // t(n) --> dk(n)d  (dual operations)
  [/j/g, 'da'],   // j --> da   (dual-ambo, since ad = a)
  [/s/g, 'dg'],   // s --> dg   (dual-gyro, since gd = g)
  [/dd/g, ''],    // dd --> null  (order 2)
  [/ad/g, 'a'],   // ad --> a   (a_ = ad_)
  [/gd/g, 'g'],   // gd --> g   (g_ = gd_)
  [/aO$/g, 'aC'], // aO --> aC  (for uniqueness)
  [/aI$/g, 'aD'], // aI --> aD  (for uniqueness)
  [/gO$/g, 'gC'], // gO --> gC  (for uniqueness)
  [/gI$/g, 'gD']  // gI --> gD  (for uniqueness)
];

const legacyTransformPass = function (notation) {
  let expanded = notation;
  for (let [orig, equiv] of legacySpecReplacements) {
    expanded = expanded.replace(orig, equiv);
  }
  console.log(`Legacy ${notation} executed as ${expanded}`);
  return expanded;
};

const dispatch = function (fn, args) {
  return fn.apply(this, args || []);
};

export const generatePolyLegacy = function (spec) {
  const transformedSpec = legacyTransformPass(spec);
  const opList = legacyOpParser.parse(transformedSpec).reverse();

  let op = opList.shift();
  const baseFunc = legacyPrimitiveMap[op['op']];
  const baseArgs = op['args'];
  let poly = dispatch(baseFunc, baseArgs);

  for (op of opList) {
    const opFunc = legacyOpMap[op['op']];
    const opArgs = [poly].concat(op['args']);
    poly = dispatch(opFunc, opArgs);
  }

  poly.vertices = geo.recenter(poly.vertices, poly.edges());
  poly.vertices = geo.rescale(poly.vertices);
  poly = polyhedron.paintPolyhedron(poly);

  return poly;
};
