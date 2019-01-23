// Polyhédronisme
// ===================================================================================================
//
// A toy for constructing and manipulating polyhedra and other meshes
//
// Includes implementation of the conway polyhedral operators derived
// from code by mathematician and mathematical sculptor
// George W. Hart http://www.georgehart.com/
//
// Copyright 2019, Anselm Levskaya
// Released under the MIT License

// import { newgeneratePoly } from './parser';

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

// test basic cross of all ops against all seeds
/*
export const testrig = function () {
  const seeds = ['T', 'O', 'C', 'I', 'D', 'P3', 'P4', 'A4', 'A5', 'Y3', 'Y4'];
  const ops = ['k', 'a', 'g', 'p', 'd', 'r', 'e', 'b', 'o', 'm', 't', 'j',
    's', 'p', 'c', 'w', 'l', 'n', 'x', 'Z', 'H'];
  console.log('===== Test Basic Ops =====');
  for (let op of ops) {
    console.log(`Operator ${op}`);
    for (let seed of seeds) {
      console.log(op + seed + ':', newgeneratePoly(op + seed));
    }
  }
  console.log('===== Done Testing Basic Ops =====');
};
*/
