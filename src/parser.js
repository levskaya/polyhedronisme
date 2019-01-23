// Polyhédronisme
// ===================================================================================================
//
// A toy for constructing and manipulating polyhedra and other meshes
//
// Copyright 2019, Anselm Levskaya
// Released under the MIT License

/* eslint-disable camelcase */
/* eslint-disable standard/array-bracket-even-spacing */
/* eslint-disable no-multi-spaces */

import * as polyhedron from './polyhedron';
import * as topo from './topo_operators';
import * as geo from './geo_operators';


// old PEG uses "eval" which can't be used in strict-mode, so have to load it in index.html as a static
// import { PEG } from './libs/peg-0.6.2.js';
// require('./libs/peg-0.6.2.js');

// Parser Routines
// ===================================================================================================
// fairly straightforward Parser Expression Grammar spec for simple
// operator-chain-on-base-polyhedra recipes
const PEG_parser_spec = `\
/* series of opspecs */
start  = opspec+

/* opspec one of:
 A  - single letter
 A3 - single letter and float
 B(5,4.3,3) - function call format w. float args
*/
opspec =
   let:opcode args:opargs {return {"op":let,"args":args};}
/ let:opcode float:float     {return {"op":let,"args":[float]};}
/ let:opcode                     {return {"op":let,"args":[]};}

/*
parentheses surrounding comma-delimited list of floats i.e.
( 1 , 3.2, 4 ) or (1) or (2,3)
*/
opargs = "("
           num:( float:float ","? {return float} )+
         ")" {return num;}

/* just a letter */
opcode = op:[a-zA-Z] {return op;}

/* standard numerical types */
int   = digits:[0-9-]+   { return parseInt(digits.join(""), 10);  }
float = digits:[0-9.-]+  { return parseFloat(digits.join(""), 10); }\
`;

// old library - PEG object loaded statically as global
const op_parser = PEG.buildParser(PEG_parser_spec);
// const op_parser = PEG.generate(PEG_parser_spec);

// applies func fn to array args i.e. f, [1,2,3] -> f(1,2,3)
const dispatch = function (fn, args) {
  return fn.apply(this, args || []);
};

const basemap = {
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

const opmap = {
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

// list of basic equivalences, easier to replace before parsing
const specreplacements = [
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

const getOps = function (notation) {
  let expanded = notation;
  for (let [orig, equiv] of specreplacements) {
    expanded = expanded.replace(orig, equiv);
  }
  console.log(`${notation} executed as ${expanded}`);
  return expanded;
};


// create polyhedron from notation
export const newgeneratePoly = function (notation) {
  const ops_spec = getOps(notation);
  const oplist = op_parser.parse(ops_spec).reverse();

  let op = oplist.shift();
  const basefunc = basemap[op['op']];
  const baseargs = op['args'];
  let poly = dispatch(basefunc, baseargs);

  for (op of oplist) {
    const opfunc = opmap[op['op']];
    const opargs = [poly].concat(op['args']);
    poly = dispatch(opfunc, opargs);
  }

  // Recenter polyhedra at origin (rarely needed)
  poly.vertices = geo.recenter(poly.vertices, poly.edges());
  poly.vertices = geo.rescale(poly.vertices);

  // Color the faces of the polyhedra for display
  poly = polyhedron.paintPolyhedron(poly);

  // return the poly object
  return poly;
};
