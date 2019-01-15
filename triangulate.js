// Polyhédronisme
//===================================================================================================
//
// A toy for constructing and manipulating polyhedra and other meshes
//
// Copyright 2019, Anselm Levskaya
// Released under the MIT License
//
/*
 * decaffeinate suggestions:
 * DS102: Remove unnecessary code created because of implicit returns
 * Full docs: https://github.com/decaffeinate/decaffeinate/blob/master/docs/suggestions.md
 */

// Polyhedra Functions
//===================================================================================================
//
// Set of routines for transforming N-face meshes into triangular meshes, necessary for exporting
// STL or VRML for 3D Printing.
//

// Ear-based triangulation of 2d faces, takes array of 2d coords in the face ordering
// Returns indices of the new diagonal lines to cut.
//
// assumes planarity of course, so this isn't the ideal algo for making aesthetically pleasing
// "flattening" choices in distorted polyhedral planes.
//
const getDiagonals = function(verts){
  let v0, v2;
  const limiter = 999;
  const diagonals = [];
  const ear = [];
  let facelen = verts.length;

  const XOR = (x, y) => (x || y) && !(x && y);
  const Area2     = (Va,Vb,Vc)   => ((Vb[0]-Va[0])*(Vc[1]-Va[1])) - ((Vc[0]-Va[0])*(Vb[1]-Va[1]));
  const Left      = (Va, Vb, Vc) => Area2(Va, Vb, Vc) > 0;
  const LeftOn    = (Va, Vb, Vc) => Area2(Va, Vb, Vc) >= 0;
  const Collinear = (Va, Vb, Vc) => Area2(Va, Vb, Vc) === 0;

  const Between   = function(Va, Vb, Vc) {
    if (Collinear(Va, Vb, Vc)) { return false; }
    if (Va[0] !== Vb[0]) {
      return ((Va[0] <= Vc[0]) && (Vc[0] <= Vb[0])) || ((Va[0] >= Vc[0]) && (Vc[0] >= Vb[0]));
    } else {
      return ((Va[1] <= Vc[1]) && (Vc[1] <= Vb[1])) || ((Va[1] >= Vc[1]) && (Vc[1] >= Vb[1]));
    }
  };

  const IntersectProp = function(Va, Vb, Vc, Vd) {
    if (Collinear(Va, Vb, Vc) || Collinear(Va, Vb, Vd) || Collinear(Vc, Vd, Va) || Collinear(Vc, Vd, Vb)) { return false; }
    return XOR(Left(Va, Vb, Vc), Left(Va, Vb, Vd)) && XOR(Left(Vc, Vd, Va), Left(Vc, Vd, Vb));
  };

  const Intersect = function(Va, Vb, Vc, Vd) {
    if (IntersectProp(Va, Vb, Vc, Vd)) {
      return true;
    } else {
      if (Between(Va, Vb, Vc) || Between(Va, Vb, Vd) || Between(Vc, Vd, Va) || Between(Vc, Vd, Vb)) {
        return true;
      } else {
        return false;
      }
    }
  };

  const InCone = function(a, b) {
    const a1 = (a+1+facelen)%facelen;
    const a0 = ((a-1)+facelen)%facelen;
    if (LeftOn(verts[a], verts[a1], verts[a0])) {
      return (Left(verts[a], verts[b], verts[a0]) && Left(verts[b], verts[a], verts[a1]));
    }
    return !(LeftOn(verts[a], verts[b], verts[a1]) && LeftOn(verts[b], verts[a], verts[a0]));
  };

  const Diagonalie = function(a, b) {
    let c = 0;
    while (true) {
      const c1 = (c+1+facelen)%facelen;
      if ((c !== a) && (c1 !== a) && (c !== b) && (c1 !== b) && IntersectProp(verts[a], verts[b], verts[c], verts[c1])) {
        return false;
      }
      c  = (c+1+facelen)%facelen;
      if (c === 0) { break; }
    }
    return true;
  };

  const Diagonal = (a, b) => InCone(a, b) && InCone(b, a) && Diagonalie(a, b);

  let v1 = 0;
  while (true) {
    v2 = (v1+1+facelen)%facelen;//v1.next
    v0 = ((v1-1)+facelen)%facelen;//v1.prev
    ear[v1] = Diagonal(v0, v2);
    v1 = (v1+1+facelen)%facelen;
    if (v1 === 0) { break; }
  }

  let origIdx = __range__(0, facelen-1, true);
  let n = facelen;//verts.length
  let z = limiter;
  let head = 0; //??
  while ((z > 0) && (n > 3)) {
    z -= 1;
    v2 = head;
    let y = limiter;
    while (true) {
      y -= 1;
      let broke = false;
      if (ear[v2]) {
        let v3 = (v2+1+facelen)%facelen;//v2.next
        let v4 = (v3+1+facelen)%facelen;//v3.next
        v1 = ((v2-1)+facelen)%facelen;//v2.prev
        v0 = ((v1-1)+facelen)%facelen;//v1.prev
        diagonals.push([ origIdx[v1], origIdx[v3] ]);
        ear[v1] = Diagonal(v0, v3);
        ear[v3] = Diagonal(v1, v4);
        //v1.next = v3
        //v3.prev = v1
        verts   = verts.slice(0, +v2 + 1 || undefined).concat(verts.slice(v3));
        origIdx = origIdx.slice(0, +v2 + 1 || undefined).concat(origIdx.slice(v3));
        if (v0>v2) { v0 -= 1; }
        if (v1>v2) { v1 -= 1; }
        if (v3>v2) { v3 -= 1; }
        if (v4>v2) { v4 -= 1; }
        facelen--;
        head = v3;
        n--;
        broke = true;
      }
      v2 = (v2+1+facelen)%facelen;//v2.next
      if ((y <= 0) || !!broke || (v2 === head)) { break; }
    }
  }

  //return diagonals
  return diagonals;
};

// equates triplets of numbers if they can be rotated into identity
const triEq = function(tri1,tri2){
    if (((tri1[0] === tri2[0]) && (tri1[1] === tri2[1]) && (tri1[2] === tri2[2]))
    ||  ((tri1[0] === tri2[1]) && (tri1[1] === tri2[2]) && (tri1[2] === tri2[0]))
    ||  ((tri1[0] === tri2[2]) && (tri1[1] === tri2[0]) && (tri1[2] === tri2[1]))) {
      return true;
    } else {
      return false;
    }
  };

// god-awful but working hack to turn diagonals into triangles
// switch to an edge-matching algo, it would be 10x simpler
const diagsToTris = function(f,diags){
  let d;
  const edges = [];
  const redges = [];
  // get edges from faces as assoc arrays
  for (let [v1,v2] of (__range__(0, f.length-1, true).map((i) => [i,(i+1)%f.length]))) {
    edges[v1]  = [v2];
    redges[v2] = [v1];
  }
  for (d of diags) { // push the diagonals into the assoc arrays in both directions!
    edges[d[0]].push(d[1]);
    edges[d[1]].push(d[0]);
    redges[d[0]].push(d[1]);
    redges[d[1]].push(d[0]);
  }
  const tris=[];
  for (d of diags) {  //orig N-face, N-2 triangles from the N-3 diagonals
    var e1, e2;
    for (e1 of edges[d[1]]) { // edge after diag
      for (e2 of redges[d[0]]) { // edge before diag
        if (e1 === e2) { // if they meet we have a triangle!
          tris.push([d[0],d[1],e1]);
        }
      }
    }
    for (e1 of edges[d[0]]) { // same as above for other dir along diagonal
      for (e2 of redges[d[1]]) {
        if (e1 === e2) {
          tris.push([d[1],d[0],e1]);
        }
      }
    }
  }
  // unfortunately the above duplicates triangles, so filter out repeats
  const uniques = [tris.pop()];
  for (let tri of tris) {
    let already_present = false;
    for (let extant_tri of uniques) {
      if (triEq(tri, extant_tri)) {
        already_present=true;
        break;
      }
    }
    if (!already_present) { uniques.push(tri); }
  }

  return uniques;
};

// driver routine, projects 3d face to 2d, get diagonals then triangles,
// then builds new polyhedron out of them, preserving original face colors
const triangulate = function(poly, colors){
  colors = colors || false;
  console.log(`Triangulating faces of ${poly.name}...`);

  const newpoly = new polyhedron();
  newpoly.xyz = clone(poly.xyz);
  newpoly.face_class = [ ];
  // iterate over triplets of faces v1,v2,v3
  for (let i = 0; i < poly.face.length; i++) {
    const f = poly.face[i];
    if (f.length > 3) {
      const TwoDface = project2dface(f.map((v) => poly.xyz[v]));
      const diags = getDiagonals(TwoDface);
      const tris  = diagsToTris(f,diags);
      for (let j = 0; j < tris.length; j++) {
        const tri = tris[j];
        newpoly.face.push([ f[tri[0]], f[tri[1]], f[tri[2]] ]);
        if (colors) { newpoly.face_class.push(poly.face_class[i]); }
      }
    } else {
      newpoly.face.push([ f[0], f[1], f[2] ]);
      if (colors) { newpoly.face_class.push(poly.face_class[i]); }
    }
  }
  newpoly.name = poly.name; // don't change the name for export
  return newpoly;
};

function __range__(left, right, inclusive) {
  let range = [];
  let ascending = left < right;
  let end = !inclusive ? right : ascending ? right + 1 : right - 1;
  for (let i = left; ascending ? i < end : i > end; ascending ? i++ : i--) {
    range.push(i);
  }
  return range;
}