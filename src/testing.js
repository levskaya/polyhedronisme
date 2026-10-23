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

// import { generatePoly } from './parser';

// Testing Functions
// ===================================================================================================

// report on face topology
export const topolog = function (poly) {
  let str = '';
  for (let f of poly.faces) {
    str += `${f.length}: `
    for (let v of f) {
      str += `${v}->`;
    }
    str += '\n';
  }
  console.log(str);
};
