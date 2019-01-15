/*
 * decaffeinate suggestions:
 * DS102: Remove unnecessary code created because of implicit returns
 * Full docs: https://github.com/decaffeinate/decaffeinate/blob/master/docs/suggestions.md
 */
// Polyhédronisme
//===================================================================================================
//
// A toy for constructing and manipulating polyhedra and other meshes
//
// Includes implementation of the conway polyhedral operators derived
// from code by mathematician and mathematical sculptor
// George W. Hart http://www.georgehart.com/
//
// Copyright 2019, Anselm Levskaya
// Released under the MIT License


// Testing Functions
//===================================================================================================

//report on face topology
const topolog = function(poly) {
  let str="";
  for (let f of poly.face) {
    for (let v of f) {
      str+=`${v}->`;
    }
    str+="\n";
  }
  return console.log(str);
};

const testrig = function() {
  const seeds=["T","O","C","I","D","P3","P4","A4","A5","Y3","Y4"];
  const ops = ["k","a","g","p","d","n","x","*"];
  console.log("===== Test Basic Ops =====");
  for (let op of ops) {
    console.log(`Operator ${op}`);
    for (let seed of seeds) {
      console.log(op+seed+":", generatePoly(op+seed));
    }
  }
  return console.log("===== Done Testing Basic Ops =====");
};
