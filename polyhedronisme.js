/*
 * decaffeinate suggestions:
 * DS101: Remove unnecessary use of Array.from
 * DS102: Remove unnecessary code created because of implicit returns
 * DS202: Simplify dynamic range loops
 * DS207: Consider shorter variations of null checks
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


// Math / Vector / Matrix Functions
//===================================================================================================

// Math is primal, people.
// import math functions to local namespace
const { random, round, floor, sqrt, sin, cos, tan, asin, acos, atan, abs,
        PI, LN10, pow, log
      } = Math;
const log10 = x=> log(x)/LN10;

//returns string w. nsigs digits ignoring magnitude
const sigfigs = function(N, nsigs){
  const normed = pow(10, log10(N)-floor(log10(N)));
  return `${round(normed*(nsigs-1))}`;
};

// for python-style enumerated for-in loops
//  - should use "for [i,x] in AR then do (i,x)->" idiom instead
//  - !! actually even easier:  "for val,idx in array" works!
//enumerate = (ar) -> [i,ar[i]] for i in [0...ar.length]

// general recursive deep-copy function
var clone = function(obj) {
  if ((obj == null) || (typeof obj !== 'object')) {
    return obj;
  }
  const newInstance = new obj.constructor();
  for (let key in obj) {
    newInstance[key] = clone(obj[key]);
  }
  return newInstance;
};

// often useful
const randomchoice = function(array){
  const n = floor(random()*array.length);
  return array[n];
};

// 3d scalar multiplication
const mult = (c, vec) => [c*vec[0],c*vec[1],c*vec[2]];

// 3d element-wise multiply
const _mult = (vec1, vec2) => [vec1[0]*vec2[0],vec1[1]*vec2[1],vec1[2]*vec2[2]];

// 3d vector addition
const add = (vec1, vec2) => [vec1[0]+vec2[0],vec1[1]+vec2[1],vec1[2]+vec2[2]];

// 3d vector subtraction
const sub = (vec1, vec2) => [vec1[0]-vec2[0],vec1[1]-vec2[1],vec1[2]-vec2[2]];

// 3d dot product
const dot = (vec1, vec2) => (vec1[0]*vec2[0]) + (vec1[1]*vec2[1]) + (vec1[2]*vec2[2]);

// 3d cross product d1 x d2
const cross = (d1, d2) => [ (d1[1]*d2[2])-(d1[2]*d2[1]), (d1[2]*d2[0])-(d1[0]*d2[2]),  (d1[0]*d2[1])-(d1[1]*d2[0]) ];

// vector norm
const mag = vec => sqrt(dot(vec,vec));

// vector magnitude squared
const mag2 = vec => dot(vec,vec);

// makes vector unit length
const unit = vec => mult( 1/sqrt(mag2(vec)), vec);

// midpoint between vec1, vec2
const midpoint = (vec1, vec2) => mult(1/2.0,add(vec1,vec2));

// parametric segment between vec1, vec2 w. parameter t ranging from 0 to 1
const tween = (vec1,vec2,t) => [ ((1-t)*vec1[0]) + (t*vec2[0]), ((1-t)*vec1[1]) + (t*vec2[1]), ((1-t)*vec1[2]) + (t*vec2[2]) ];

// uses above to go one-third of the way along vec1->vec2 line
const oneThird = (vec1, vec2) => tween(vec1, vec2, 1/3.0);

// reflect 3vec in unit sphere, spherical reciprocal
const reciprocal = vec => mult( 1.0/mag2(vec), vec);

// point where line v1...v2 tangent to an origin sphere
const tangentPoint= function(v1,v2) {
  const d = sub(v2, v1);
  return sub(v1, mult(dot(d,v1)/mag2(d),d));
};

// distance of line v1...v2 to origin
const edgeDist = (v1,v2) => sqrt(mag2(tangentPoint(v1, v2)));

// square of distance from point v3 to line segment v1...v2
// http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
const linePointDist2 = function(v1,v2,v3) {
  let result;
  const d21 = sub(v2, v1);
  const d13 = sub(v1, v3);
  const m2 = mag2(d21);
  const t = -dot(d13, d21)/m2;
  if (t <= 0) {
    return mag2(d13);
  } else if (t >= 1) {
    result = mag2(sub(v2, v3));
  }

  result = mag2(cross(d21, d13))/m2;
  return result;
};
  
// find vector orthogonal to plane of 3 pts
// -- do the below algos assume this be normalized or not?
const orthogonal = function(v1,v2,v3) {
  // adjacent edge vectors
  const d1 = sub(v2, v1);
  const d2 = sub(v3, v2);
  // cross product
  return cross(d1, d2);
};

// find first element common to 3 sets by brute force search
const intersect = function(set1, set2, set3) {
  for (let s1 of set1) {
    for (let s2 of set2) {
      if (s1 === s2) {
        for (let s3 of set3) {
          if (s1 === s3) {
            return s1;
          }
        }
      }
    }
  }
  return null; // oh noes!
};

// calculate centroid of array of vertices
const calcCentroid = function(xyzs) {
    let centroidV = [0,0,0]; // running sum of vertex coords
    for (let v of xyzs) {
      centroidV = add(centroidV, v);
    }
    return mult(1 / xyzs.length, centroidV );
  };

// calculate average normal vector for array of vertices
const normal = function(xyzs) {
    let normalV = [0,0,0]; // running sum of normal vectors
    let [v1,v2] = Array.from(xyzs.slice(-2));
    for (let v3 of xyzs) {
      normalV = add(normalV, orthogonal(v1, v2, v3));
      [v1,v2] = Array.from([v2,v3]);
    } // shift over one
    return unit(normalV);
  };

// calculates area planar face by summing over subtriangle areas
//  _Assumes_ Convexity!
const convexarea = function(xyzs) {
    let area = 0.0;
    let [v1,v2] = Array.from(xyzs.slice(0, 2));
    for (let v3 of xyzs.slice(2)) {
      //area of sub-triangle
      area += mag( cross(sub(v2,v1), sub(v3,v1)) );
      v2 = v3;
    } // shift over one
    return area;
  };

//returns array of ~3sigfig angle
const faceSignature = function(xyzs) {
    let x;
    const cross_array = [];
    let [v1,v2] = Array.from(xyzs.slice(0, 2));
    for (let v3 of xyzs.slice(2)) {
      //area of sub-triangle
      cross_array.push(mag( cross(sub(v2,v1), sub(v3,v1)) ));
      v2 = v3;
    } // shift over one

    cross_array.sort((a,b)=> a-b); //sort for uniqueness

    let sig=""; // turn it into a string
    for (x of cross_array) { sig+=sigfigs(x,2); }
    // hack to make reflected faces share the same signature
    for (x of cross_array.reverse()) { sig+=sigfigs(x,2); }
    return sig;
  };

// projects 3d polyhedral face to 2d polygon
// for triangulation and face display
const project2dface = function(verts){
  let tmpverts = clone(verts);
  const v0=verts[0];
  tmpverts = _.map(tmpverts, x=> x-v0);

  const n = normal(verts);
  const c = unit(calcCentroid(verts));
  const p = cross(n,c);

  return tmpverts.map((v) => [dot(n,v),dot(p,v)]);
};

// copies array of arrays by value (deep copy)
const copyVecArray = function(vecArray){
  const newVecArray = new Array(vecArray.length);
  for (let i = 0, end = vecArray.length; i < end; i++) {
    newVecArray[i] = vecArray[i].slice(0);
  }
  return newVecArray;
};

// 3d matrix vector multiply
const mv3 = (mat,vec) =>
  //Ghetto custom def of matrix-vector mult
  //example matrix: [[a,b,c],[d,e,f],[g,h,i]]
  [(mat[0][0]*vec[0])+(mat[0][1]*vec[1])+(mat[0][2]*vec[2]),
   (mat[1][0]*vec[0])+(mat[1][1]*vec[1])+(mat[1][2]*vec[2]),
   (mat[2][0]*vec[0])+(mat[2][1]*vec[1])+(mat[2][2]*vec[2])]
;

// 3d matrix matrix multiply
const mm3 = (A,B) =>
  [[(A[0][0]*B[0][0])+(A[0][1]*B[1][0])+(A[0][2]*B[2][0]),
   (A[0][0]*B[0][1])+(A[0][1]*B[1][1])+(A[0][2]*B[2][1]),
   (A[0][0]*B[0][2])+(A[0][1]*B[1][2])+(A[0][2]*B[2][2])],
  [(A[1][0]*B[0][0])+(A[1][1]*B[1][0])+(A[1][2]*B[2][0]),
   (A[1][0]*B[0][1])+(A[1][1]*B[1][1])+(A[1][2]*B[2][1]),
   (A[1][0]*B[0][2])+(A[1][1]*B[1][2])+(A[1][2]*B[2][2])],
  [(A[2][0]*B[0][0])+(A[2][1]*B[1][0])+(A[2][2]*B[2][0]),
   (A[2][0]*B[0][1])+(A[2][1]*B[1][1])+(A[2][2]*B[2][1]),
   (A[2][0]*B[0][2])+(A[2][1]*B[1][2])+(A[2][2]*B[2][2])]]
;

const eye3 = [[1,0,0],[0,1,0],[0,0,1]];

// Rotation Matrix
// Totally ghetto, not at all in agreement with euler angles!
// use quaternions instead
const rotm = function(phi,theta,psi){
    const xy_mat = [
      [cos(phi), -1.0*sin(phi),  0.0],
      [sin(phi),      cos(phi),  0.0],
      [0.0,                0.0,  1.0]];
    const yz_mat = [
      [cos(theta), 0, -1.0*sin(theta)],
      [         0, 1,               0],
      [sin(theta), 0,      cos(theta)]];
    const xz_mat = [
      [1.0,        0,             0],
      [  0, cos(psi), -1.0*sin(psi)],
      [  0, sin(psi),      cos(psi)]];

    return mm3(xz_mat, mm3(yz_mat,xy_mat));
  };


// Rotation Matrix defined by rotation about (unit) axis [x,y,z] for angle radians
const vec_rotm = function(angle, x, y, z) {
  let m;
  angle /= 2;
  const sinA = sin(angle);
  const cosA = cos(angle);
  const sinA2 = sinA*sinA;
  const length = mag([x,y,z]);
  if (length === 0) {
    [x,y,z] = Array.from([0,0,1]);
  }
  if (length !== 1) {
    [x,y,z] = Array.from(unit([x,y,z]));
  }

  //console.log "vec_rotm args",angle,x,y,z,"vars",sinA,cosA

  if ((x === 1) && (y === 0) && (z === 0)) {
      m=[[1,            0,           0],
         [0,    1-(2*sinA2), 2*sinA*cosA],
         [0, -2*sinA*cosA,   1-(2*sinA2)]];
  } else if ((x === 0) && (y === 1) && (z === 0)) {
      m=[[  1-(2*sinA2), 0,  -2*sinA*cosA],
         [          0, 1,             0],
         [2*sinA*cosA, 0,     1-(2*sinA2)]];
  } else if ((x === 0) && (y === 0) && (z === 1)) {
      m=[[   1-(2*sinA2), 2*sinA*cosA, 0],
         [-2*sinA*cosA,   1-(2*sinA2), 0],
         [           0,           0, 1]];
  } else {
      const x2 = x*x;
      const y2 = y*y;
      const z2 = z*z;
      m=[[1-(2*(y2+z2)*sinA2),         2*((x*y*sinA2)+(z*sinA*cosA)), 2*((x*z*sinA2)-(y*sinA*cosA))],
         [2*((y*x*sinA2)-(z*sinA*cosA)),         1-(2*(z2+x2)*sinA2), 2*((y*z*sinA2)+(x*sinA*cosA))],
         [2*((z*x*sinA2)+(y*sinA*cosA)), 2*((z*y*sinA2)-(x*sinA*cosA)),         1-(2*(x2+y2)*sinA2)]];
    }

  //console.log "vec_rotm m", m[0],m[1],m[2],m[3],m[4],m[5],m[6],m[7],m[8]
  //return matrix
  return m;
};

// Perspective Transform
// assumes world's been rotated appropriately such that Z is depth
// scales perspective such that inside depth regions min_real_depth <--> max_real_depth
// perspective lengths vary no more than:   desired_ratio
// with target dimension of roughly length: desired_length
const perspT = function(vec3, max_real_depth, min_real_depth, desired_ratio, desired_length) {
  const z0          = ((max_real_depth * desired_ratio) - min_real_depth)/(1-desired_ratio);
  const scalefactor =  (desired_length * desired_ratio)/(1-desired_ratio);

  // projected [X, Y]
  return [(scalefactor*vec3[0])/(vec3[2]+z0), (scalefactor*vec3[1])/(vec3[2]+z0)];
};

// Inverses perspective transform by projecting plane onto a unit sphere at origin
const invperspT = 
function(x, y, dx, dy, max_real_depth, min_real_depth, desired_ratio, desired_length) {
  const z0 = ((max_real_depth * desired_ratio) - min_real_depth)/(1-desired_ratio);
  const s =  (desired_length * desired_ratio)/(1-desired_ratio);
  const xp = x-dx;
  const yp = y-dy;
  const s2 = s*s;
  const z02 = z0*z0;
  const xp2 = xp*xp;
  const yp2 = yp*yp;

  const xsphere = ((2*s*xp*z0) + sqrt((4*s2*xp2*z02) + (4*xp2*(s2+xp2+yp2)*(1-z02)) ) )/(2.0*(s2+xp2+yp2));
  const ysphere = (((s*yp*z0)/(s2+xp2+yp2)) + ((yp*sqrt((4*s2*z02) + (4*(s2+xp2+yp2)*(1-z02))))/(2.0*(s2+xp2+yp2))));
  const zsphere = sqrt(1-(xsphere*xsphere)-(ysphere*ysphere));

  //console.log  "invperspT", xsphere, ysphere, zsphere, mag([xsphere, ysphere, zsphere])
  return [xsphere, ysphere, zsphere];
};

// Returns rotation matrix that takes vec1 to vec2
const getVec2VecRotM = function(vec1, vec2){
  const axis    = cross(vec1, vec2);
  const angle   = acos(dot(vec1, vec2));

  //console.log  "getVec2VecRotM", angle, axis[0],axis[1],axis[2]
  return vec_rotm(-1*angle, axis[0], axis[1], axis[2]);
};
/*
 * decaffeinate suggestions:
 * DS101: Remove unnecessary use of Array.from
 * DS102: Remove unnecessary code created because of implicit returns
 * DS202: Simplify dynamic range loops
 * DS205: Consider reworking code to avoid use of IIFEs
 * Full docs: https://github.com/decaffeinate/decaffeinate/blob/master/docs/suggestions.md
 */
// Polyhédronisme
//===================================================================================================
//
// A toy for constructing and manipulating polyhedra and other meshes
//
// Copyright 2019, Anselm Levskaya
// Released under the MIT License
//

// Polyhedra Functions
//===================================================================================================
//
// Topology stored as set of "faces."  Each face is list of n vertex indices
// corresponding to one n-sided face.  Vertices listed clockwise as seen from outside.

// Generate an array of edges [v1,v2] for the face.
const faceToEdges = function(face) {
  const edges = [];
  let [v1] = Array.from(face.slice(-1));
  for (let v2 of face) {
    edges.push([v1,v2]);
    v1 = v2;
  }
  return edges;
};

const vertColors = function(poly) {
  const vertcolors=[];
  for (let i = 0; i < poly.face.length; i++) {
    const f = poly.face[i];
    for (let v of f) {
      vertcolors[v] = poly.face_class[i];
    }
  }
  return vertcolors;
};

// Polyhedra Coloring Functions
//===================================================================================================

const rwb_palette  = ["#ff7777","#dddddd","#889999","#fff0e5","#aa3333","#ff0000","#ffffff","#aaaaaa"];

// converts #xxxxxx / #xxx format into list of [r,g,b] floats
const hextofloats = function(hexstr){
  let rgb;
  if (hexstr[0] === "#") {
    hexstr = hexstr.slice(1);
  }
  if (hexstr.length === 3) {
    rgb = hexstr.split('').map(       c=> parseInt(c+c, 16)/255);
  } else {
    rgb = hexstr.match(/.{2}/g).map(  c=> parseInt(c, 16)/255);
  }
  return rgb;
};

const PALETTE = rwb_palette; //GLOBAL
const palette = function(n) {
  if (n < PALETTE.length) {
    return hextofloats(PALETTE[n]);
  } else {
    return hextofloats(PALETTE[PALETTE.length-1]);
  }
};

const paintPolyhedron = function(poly) {
  // Color the faces of the polyhedra for display
  let v;
  poly.face_class = [];
  const colormemory={};

  //memoized color assignment to faces of similar areas
  const colorassign = function(ar, colormemory) {
    const hash = round(100*ar);
    if (hash in colormemory) {
      return colormemory[hash];
    } else {
      const fclr = _.toArray(colormemory).length; //palette _.toArray(colormemory).length
      colormemory[hash] = fclr;
      return fclr;
    }
  };

  for (var f of poly.face) {
    var clr, face_verts;
    if (COLOR_METHOD === "area") {
      // color by face area (quick proxy for different kinds of faces) convexarea
      face_verts = ((() => {
        const result = [];
        for (v of f) {           result.push(poly.xyz[v]);
        }
        return result;
      })());
      clr = colorassign(convexarea(face_verts), colormemory);
    } else if (COLOR_METHOD === "signature") {
      face_verts = ((() => {
        const result1 = [];
        for (v of f) {           result1.push(poly.xyz[v]);
        }
        return result1;
      })());
      clr = colorassign(faceSignature(face_verts), colormemory);
    } else {
      // color by face-sidedness
      clr = f.length-3;
    }

    poly.face_class.push(clr);
  }
  console.log(_.toArray(colormemory).length+" face classes");
  return poly;
};

// z sorts faces of poly
// -------------------------------------------------------------------------
const sortfaces = function(poly) {
  //smallestZ = (x) -> _.sortBy(x,(a,b)->a[2]-b[2])[0]
  //closests = (smallestZ(poly.xyz[v] for v in f) for f in poly.face)
  let idx;
  const centroids  = poly.centers();
  const normals    = poly.normals();
  const ray_origin = [0,0, ((persp_z_max * persp_ratio) - persp_z_min)/(1-persp_ratio)];
  //console.log ray_origin

  // sort by binary-space partition: are you on same side as view-origin or not?
  // !!! there is something wrong with this. even triangulated surfaces have artifacts.
  const planesort = (a,b)=>
    //console.log dot(sub(ray_origin,a[0]),a[1]), dot(sub(b[0],a[0]),a[1])
    -dot(sub(ray_origin,a[0]),a[1])*dot(sub(b[0],a[0]),a[1])
  ;

  // sort by centroid z-depth: not correct but more stable heuristic w. weird non-planar "polygons"
  const zcentroidsort = (a,b)=> a[0][2]-b[0][2];

  const zsortIndex = _.zip(centroids, normals, __range__(0, poly.face.length, false))
    //.sort(planesort)
    .sort(zcentroidsort)
    .map(x=> x[2]);

  // sort all face-associated properties
  poly.face = ((() => {
    const result = [];
    for (idx of zsortIndex) {       result.push(poly.face[idx]);
    }
    return result;
  })());
  return poly.face_class = ((() => {
    const result1 = [];
    for (idx of zsortIndex) {       result1.push(poly.face_class[idx]);
    }
    return result1;
  })());
};


class polyhedron {
  constructor(verts,faces,name) {      // constructor of initially null polyhedron
    this.face = faces || new Array();   // array of faces.          face.length = # faces
    this.xyz  = verts || new Array();   // array of vertex coords.  xyz.length = # of vertices
    this.name = name  || "null polyhedron";
  }

  data() {   // informative string
    const nEdges = (this.face.length + this.xyz.length) - 2; // E = V + F - 2
    return `${this.face.length} faces, ${nEdges} edges, ${this.xyz.length} vertices`;
  }

  moreData() {
    return `min. edge length ${this.minEdgeLength().toPrecision(2)}; min. face radius ${this.minFaceRadius().toPrecision(2)}`;
  }
    
  edges() {
    let e;
    const finalset={};
    const uniqedges=[];
    const alledges = _.map(this.face, faceToEdges);
    for (let edgeset of alledges) {
      for (e of edgeset) {
        var a, b;
        if (e[0] < e[1]) {
          [a,b] = Array.from(e);
        } else {
          [b,a] = Array.from(e);
        }
        finalset[a+'~'+b] = e;
      }
    }
    for (let hash in finalset) {
      e = finalset[hash];
      uniqedges.push(e);
    }

    //return edges
    return uniqedges;
  }

  minEdgeLength() {
    let min2 = Number.MAX_VALUE;
    // Compute minimum edge length
    for (let e of this.edges()) {
      const d2 = mag2(sub(this.xyz[e[0]], this.xyz[e[1]])); // square of edge length
      if (d2 < min2) {
        min2 = d2;
      }
    }
    return sqrt(min2); // This is normalized if rescaling has happened.
  }
    
  minFaceRadius() {
    let min2 = Number.MAX_VALUE;
    const nFaces = this.face.length;
    const centers = this.centers();
    for (let f = 0, end = nFaces, asc = 0 <= end; asc ? f < end : f > end; asc ? f++ : f--) {
      const c = centers[f];
      for (let e of faceToEdges(this.face[f])) {
        // Check distance from center to each edge.
        const de2 = linePointDist2(this.xyz[e[0]], this.xyz[e[1]], c);
        if (de2 < min2) {
          min2 = de2;
        }
      }
    }

    return sqrt(min2);
  }
      
  centers() {
    // get array of face centers
    const centers_array = [];
    for (let f of this.face) {
      let fcenter = [0,0,0];
      for (let v of f) { //avg vertex coords
        fcenter = add(fcenter, this.xyz[v]);
      } // add
      centers_array.push(mult(1.0/f.length, fcenter));
    } // div by n
    // return face-ordered array of centroids
    return centers_array;
  }

  normals() {
  // get array of face normals
    const normals_array = [];
    for (let f of this.face) {
      normals_array.push(normal(f.map((v) => this.xyz[v])));
    }
    return normals_array;
  }

  // Export / Formatting Routines --------------------------------------------------

  // produces vanilla OBJ files for import into 3d apps
  toOBJ() {
    let f;
    let v;
    let objstr="#Produced by polyHédronisme http://levskaya.github.com/polyhedronisme\n";
    objstr+=`group ${this.name}\n`;
    objstr+="#vertices\n";
    for (v of this.xyz) {
      objstr += `v ${v[0]} ${v[1]} ${v[2]}\n`;
    }

    objstr += "#normal vector defs \n";
    for (f of this.face) {
      const norm = normal((() => {
        const result = [];
        for (v of f) {           result.push(this.xyz[v]);
        }
        return result;
      })());
      objstr += `vn ${norm[0]} ${norm[1]} ${norm[2]}\n`;
    }

    objstr += "#face defs \n";
    for (let i = 0; i < this.face.length; i++) {
      f = this.face[i];
      objstr += "f ";
      for (v of f) {
        objstr += `${v+1}//${i+1} `;
      }
      objstr += "\n";
    }

    return objstr;
  }

  toX3D() {
    let v;
    const SCALE_FACTOR = .01; //ShapeWays uses 1unit = 1meter, so reduce to 1cm scale
    // opening cruft
    let x3dstr=`\
<?xml version="1.0" encoding ="UTF-8"?>
<X3D profile="Interchange" version="3.0">
<head>
<component name="Rendering" level="3"/>
<meta name="generator" content="Polyhedronisme"/>
<meta name="version" content="0.1.0"/>
</head>
<Scene>
<Shape>
<IndexedFaceSet normalPerVertex="false" coordIndex="\
`;
    // face indices
    for (let f of this.face) {
      for (v of f) {
        x3dstr+=`${v} `;
      }
      x3dstr+='-1\n';
    }
    x3dstr+='">\n';

    // per-face Color
    x3dstr+='<Color color="';
    for (let cl of vertColors(this)) {//@face_class
      const clr=palette(cl);
      x3dstr+=`${clr[0]} ${clr[1]} ${clr[2]} `;
    }
    x3dstr+='"/>';

    // re-scaled xyz coordinates
    x3dstr+='<Coordinate point="';
    for (v of this.xyz) {
      x3dstr+=`${v[0]*SCALE_FACTOR} ${v[1]*SCALE_FACTOR} ${v[2]*SCALE_FACTOR} `;
    }
    x3dstr+='"/>\n';

      // end cruft
    x3dstr+=`\
</IndexedFaceSet>
</Shape>
</Scene>
</X3D>`;

    return x3dstr;
  }

  toVRML() {
    let v;
    const SCALE_FACTOR = .01; //ShapeWays uses 1unit = 1meter, so reduce to 1cm scale
    // opening cruft
    let x3dstr=`\
#VRML V2.0 utf8
#Generated by Polyhedronisme
NavigationInfo {
	type [ "EXAMINE", "ANY" ]
}
Transform {
  scale 1 1 1
  translation 0 0 0
  children
  [
    Shape
    {
      geometry IndexedFaceSet
      {
        creaseAngle .5
        solid FALSE
        coord Coordinate
        {
          point
          [\
`;
    // re-scaled xyz coordinates
    for (v of this.xyz) {
      x3dstr+=`${v[0]*SCALE_FACTOR} ${v[1]*SCALE_FACTOR} ${v[2]*SCALE_FACTOR},`;
    }
    x3dstr=x3dstr.slice(0, +-2 + 1 || undefined);
    x3dstr+=`\
    ]
}
color Color
{
  color
  [\
`;
    // per-face Color
    for (let cl of this.face_class) {
      const clr=palette(cl);
      x3dstr+=`${clr[0]} ${clr[1]} ${clr[2]} ,`;
    }
    x3dstr=x3dstr.slice(0, +-2 + 1 || undefined);
    x3dstr+=`\
  ]
}
colorPerVertex FALSE
coordIndex
[\
`;
    // face indices
    for (let f of this.face) {
      for (v of f) {
        x3dstr+=`${v}, `;
      }
      x3dstr+='-1,';
    }
    x3dstr=x3dstr.slice(0, +-2 + 1 || undefined);
    x3dstr+=`\
          ]
      }
      appearance Appearance
      {
        material Material
        {
	       ambientIntensity 0.2
	       diffuseColor 0.9 0.9 0.9
	       specularColor .1 .1 .1
	       shininess .5
        }
      }
    }
  ]
}\
`;
    return x3dstr;
  }
}

//===================================================================================================
// Primitive Polyhedra Seeds
//===================================================================================================

const tetrahedron = function() {
  const poly = new polyhedron();
  poly.name = "T";
  poly.face = [ [0,1,2], [0,2,3], [0,3,1], [1,3,2] ];
  poly.xyz  = [ [1.0,1.0,1.0], [1.0,-1.0,-1.0], [-1.0,1.0,-1.0], [-1.0,-1.0,1.0] ];
  return poly;
};

const octahedron = function() {
  const poly = new polyhedron();
  poly.name = "O";
  poly.face = [ [0,1,2], [0,2,3], [0,3,4], [0,4,1], [1,4,5], [1,5,2], [2,5,3], [3,5,4] ];
  poly.xyz  = [ [0,0,1.414], [1.414,0,0], [0,1.414,0], [-1.414,0,0], [0,-1.414,0], [0,0,-1.414] ];
  return poly;
};

const cube = function() {
  const poly = new polyhedron();
  poly.name = "C";
  poly.face = [ [3,0,1,2], [3,4,5,0], [0,5,6,1], [1,6,7,2], [2,7,4,3], [5,4,7,6] ];
  poly.xyz  = [ [0.707,0.707,0.707], [-0.707,0.707,0.707], [-0.707,-0.707,0.707], [0.707,-0.707,0.707],
                [0.707,-0.707,-0.707], [0.707,0.707,-0.707], [-0.707,0.707,-0.707], [-0.707,-0.707,-0.707] ];
  return poly;
};

const icosahedron = function() {
  const poly = new polyhedron();
  poly.name = "I";
  poly.face = [ [0,1,2], [0,2,3], [0,3,4], [0,4,5],
    [0,5,1], [1,5,7], [1,7,6], [1,6,2],
    [2,6,8], [2,8,3], [3,8,9], [3,9,4],
    [4,9,10], [4,10,5], [5,10,7], [6,7,11],
    [6,11,8], [7,10,11], [8,11,9], [9,11,10] ];

  poly.xyz = [ [0,0,1.176], [1.051,0,0.526],
    [0.324,1.0,0.525], [-0.851,0.618,0.526],
    [-0.851,-0.618,0.526], [0.325,-1.0,0.526],
    [0.851,0.618,-0.526], [0.851,-0.618,-0.526],
    [-0.325,1.0,-0.526], [-1.051,0,-0.526],
    [-0.325,-1.0,-0.526], [0,0,-1.176] ];
  return poly;
};

const dodecahedron = function() {
   const poly = new polyhedron();
   poly.name = "D";
   poly.face = [ [0,1,4,7,2], [0,2,6,9,3], [0,3,8,5,1],
      [1,5,11,10,4], [2,7,13,12,6], [3,9,15,14,8],
      [4,10,16,13,7], [5,8,14,17,11], [6,12,18,15,9],
      [10,11,17,19,16], [12,13,16,19,18], [14,15,18,19,17] ];
   poly.xyz = [ [0,0,1.07047], [0.713644,0,0.797878],
      [-0.356822,0.618,0.797878], [-0.356822,-0.618,0.797878],
      [0.797878,0.618034,0.356822], [0.797878,-0.618,0.356822],
      [-0.934172,0.381966,0.356822], [0.136294,1.0,0.356822],
      [0.136294,-1.0,0.356822], [-0.934172,-0.381966,0.356822],
      [0.934172,0.381966,-0.356822], [0.934172,-0.381966,-0.356822],
      [-0.797878,0.618,-0.356822], [-0.136294,1.0,-0.356822],
      [-0.136294,-1.0,-0.356822], [-0.797878,-0.618034,-0.356822],
      [0.356822,0.618,-0.797878], [0.356822,-0.618,-0.797878],
      [-0.713644,0,-0.797878], [0,0,-1.07047] ];
   return poly;
 };

const prism = function(n) {
  let i;
  let asc, end;
  let asc1, end1;
  let asc2, end2;
  const theta = (2*PI)/n; // pie angle
  const h = Math.sin(theta/2); // half-edge
  let poly = new polyhedron();
  poly.name = `P${n}`;

  for (i = 0, end = n, asc = 0 <= end; asc ? i < end : i > end; asc ? i++ : i--) { // vertex #'s 0 to n-1 around one face
    poly.xyz.push([-cos(i*theta), -sin(i*theta),  -h]);
  }
  for (i = 0, end1 = n, asc1 = 0 <= end1; asc1 ? i < end1 : i > end1; asc1 ? i++ : i--) { // vertex #'s n to 2n-1 around other
    poly.xyz.push([-cos(i*theta), -sin(i*theta), h]);
  }

  poly.face.push(__range__(n-1, 0, true));   //top
  poly.face.push(__range__(n, 2*n, false)); //bottom
  for (i = 0, end2 = n, asc2 = 0 <= end2; asc2 ? i < end2 : i > end2; asc2 ? i++ : i--) { //n square sides
    poly.face.push([i, (i+1)%n, ((i+1)%n)+n, i+n]);
  }

  poly = adjustXYZ(poly,1);
  return poly;
};

const antiprism = function(n) {
  let i;
  let asc, end;
  let asc1, end1;
  let asc2, end2;
  const theta = (2*PI)/n; // pie angle
  let h = sqrt(1-(4/((4+(2*cos(theta/2)))-(2*cos(theta)))));
  let r = sqrt(1-(h*h));
  const f = sqrt((h*h) + pow(r*cos(theta/2),2) );
  // correction so edge midpoints (not vertices) on unit sphere
  r = -r/f;
  h = -h/f;
  let poly = new polyhedron();
  poly.name = `A${n}`;

  for (i = 0, end = n, asc = 0 <= end; asc ? i < end : i > end; asc ? i++ : i--) { // vertex #'s 0...n-1 around one face
    poly.xyz.push([r * cos(i*theta), r * sin(i*theta), h]);
  }
  for (i = 0, end1 = n, asc1 = 0 <= end1; asc1 ? i < end1 : i > end1; asc1 ? i++ : i--) { // vertex #'s n...2n-1 around other
    poly.xyz.push([r * cos((i+0.5)*theta), r * sin((i+0.5)*theta), -h]);
  }

  poly.face.push(__range__(n-1, 0, true));   //top
  poly.face.push(__range__(n, (2*n)-1, true)); //bottom
  for (i = 0, end2 = n-1, asc2 = 0 <= end2; asc2 ? i <= end2 : i >= end2; asc2 ? i++ : i--) { //2n triangular sides
    poly.face.push([i, (i+1)%n, i+n]);
    poly.face.push([i, i+n, ((((n+i)-1)%n)+n)]);
  }

  poly = adjustXYZ(poly,1);
  return poly;
};

const pyramid = function(n) {
  let i;
  let asc, end;
  let asc1, end1;
  const theta = (2*PI)/n; // pie angle
  const height = 1;
  let poly = new polyhedron();
  poly.name = `Y${n}`;

  for (i = 0, end = n, asc = 0 <= end; asc ? i < end : i > end; asc ? i++ : i--) { // vertex #'s 0...n-1 around one face
    poly.xyz.push([-cos(i*theta), -sin(i*theta), -0.2]);
  }
  poly.xyz.push([0,0, height]); // apex

  poly.face.push(__range__(n-1, 0, true)); // base
  for (i = 0, end1 = n, asc1 = 0 <= end1; asc1 ? i < end1 : i > end1; asc1 ? i++ : i--) { // n triangular sides
    poly.face.push([i, (i+1)%n, n]);
  }

  poly = canonicalXYZ(poly,3);
  return poly;
};


function __range__(left, right, inclusive) {
  let range = [];
  let ascending = left < right;
  let end = !inclusive ? right : ascending ? right + 1 : right - 1;
  for (let i = left; ascending ? i < end : i > end; ascending ? i++ : i--) {
    range.push(i);
  }
  return range;
}/*
 * decaffeinate suggestions:
 * DS101: Remove unnecessary use of Array.from
 * DS102: Remove unnecessary code created because of implicit returns
 * DS202: Simplify dynamic range loops
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


//===================================================================================================
// Polyhedron Flagset Construct
//
// A Flag is an associative triple of a face index and two adjacent vertex indices,
// listed in geometric clockwise order (staring into the normal)
//
// Face_i -> V_i -> V_j
//
// They are a useful abstraction for defining topological transformations of the polyhedral mesh, as
// one can refer to vertices and faces that don't yet exist or haven't been traversed yet in the
// transformation code.
//
// A flag is similar in concept to a directed halfedge in halfedge data structures.
//
class polyflag {
  constructor() {
    this.flags= new Object(); // flags[face][vertex] = next vertex of flag; symbolic triples
    this.verts= new Object(); // XYZ coordinates
    this.xyzs = new Object(); // [symbolic names] holds vertex index
  }

  // Add a new vertex named "name" with coordinates "xyz".
  newV(name, xyz) {
    if (this.verts[name] === undefined) {
      this.verts[name] = 0;
      return this.xyzs[name] = xyz;
    }
  }

  newFlag(facename, v1, v2) {
    if (this.flags[facename] === undefined) {
      this.flags[facename] = {};
    }
    return this.flags[facename][v1] = v2;
  }

  topoly() {
    let i, v;
    const poly = new polyhedron();

    let ctr=0; // first number the vertices
    for (i in this.verts) {
      v = this.verts[i];
      poly.xyz[ctr]=this.xyzs[i]; // store in array
      this.verts[i] = ctr;
      ctr++;
    }

    ctr=0;
    for (i in this.flags) {
      var v0;
      const f = this.flags[i];
      poly.face[ctr] = []; // new face
      // grab _any_ vertex as starting point
      for (let j in f) {
        v = f[j];
        v0 = v;
        break;
      }  // need just one
      // build face out of all the edge relations in the flag assoc array
      v = v0; // v moves around face
      poly.face[ctr].push(this.verts[v]); //record index
      v = this.flags[i][v]; // goto next vertex
      let faceCTR=0;
      while (v !== v0) { // loop until back to start
        poly.face[ctr].push(this.verts[v]);
        v = this.flags[i][v];
        faceCTR++;
        if (faceCTR>1000) { // necessary to prevent browser hangs on badly formed flagsets!
          console.log("Bad flag spec, have a neverending face:", i, this.flags[i]);
          break;
        }
      }
      ctr++;
    }

    poly.name = "unknown polyhedron";
    return poly;
  }
}


//===================================================================================================
// Polyhedron Operators
//===================================================================================================
//          for each vertex of new polyhedron:
//              call newV(Vname, xyz) with a symbolic name and coordinates
//          for each flag of new polyhedron:
//              call newFlag(Fname, Vname1, Vname2) with a symbolic name for the new face
//              and the symbolic name for two vertices forming an oriented edge
//          ORIENTATION -must- be dealt with properly to make a manifold (correct) mesh.
//          Specifically, no edge v1->v2 can ever be crossed in the -same direction- by
//          two different faces
//
//          call topoly() to assemble flags into polyhedron structure by following the orbits
//          of the vertex mapping stored in the flagset for each new face
//
//          set name as appropriate

// Kis(N)
// ------------------------------------------------------------------------------------------
// Kis (abbreviated from triakis) transforms an N-sided face into an N-pyramid rooted at the
// same base vertices.
// only kis n-sided faces, but n==0 means kiss all.
//
const kisN = function(poly, n, apexdist){
  let i;
  if (!n) { n = 0; }
  if (!apexdist) { apexdist = 0.1; }
  console.log(`Taking kis of ${n===0 ? "" : n}-sided faces of ${poly.name}...`);

  const flag = new polyflag();
  for (i = 0; i < poly.xyz.length; i++) {
    // each old vertex is a new vertex
    const p = poly.xyz[i];
    flag.newV(`v${i}`, p);
  }

  const normals = poly.normals();
  const centers = poly.centers();
  let foundAny = false;                 // alert if don't find any
  for (i = 0; i < poly.face.length; i++) {
    const f = poly.face[i];
    let v1 = `v${f[f.length-1]}`;
    for (let v of f) {
      const v2 = `v${v}`;
      if ((f.length === n) || (n === 0)) {
        foundAny = true;
        const apex = `apex${i}`;
        const fname = `${i}${v1}`;
        flag.newV(apex, add(centers[i],mult(apexdist,normals[i]))); // new vertices in centers of n-sided face
        flag.newFlag(fname,   v1,   v2); // the old edge of original face
        flag.newFlag(fname,   v2, apex); // up to apex of pyramid
        flag.newFlag(fname, apex,   v1); // and back down again
      } else {
        flag.newFlag(`${i}`, v1, v2);  // same old flag, if non-n
      }
      v1=v2;
    }
  }  // current becomes previous

  if (!foundAny) {
    console.log(`No ${n}-fold components were found.`);
  }

  const newpoly = flag.topoly();
  newpoly.name = `k${n === 0 ? "" : n}${poly.name}`;
  //newpoly.xyz = adjustXYZ(newpoly, 3)
  //newpoly.xyz = canonicalXYZ(newpoly, 3)  # this tends to make results look like shit
  return newpoly;
};


// Ambo
// ------------------------------------------------------------------------------------------
// The best way to think of the ambo operator is as a topological "tween" between a polyhedron
// and its dual polyhedron.  Thus the ambo of a dual polyhedron is the same as the ambo of the
// original. Also called "Rectify".
//
const ambo = function(poly){
  console.log(`Taking ambo of ${poly.name}...`);

  // helper func to insure unique names of midpoints
  const midName = function(v1, v2) { if (v1<v2) { return v1+"_"+v2; } else { return v2+"_"+v1; } };

  const flag = new polyflag();

  // For each face f in the original poly
  for (let i = 0; i < poly.face.length; i++) {
    const f = poly.face[i];
    let [v1, v2] = Array.from(f.slice(-2));
    for (let v3 of f) {
      if (v1 < v2) { // vertices are the midpoints of all edges of original poly
        flag.newV(midName(v1,v2), midpoint(poly.xyz[v1], poly.xyz[v2]));
      }
      // two new flags:
      // One whose face corresponds to the original f:
      flag.newFlag(`orig${i}`,  midName(v1,v2), midName(v2,v3));
      // Another flag whose face  corresponds to (the truncated) v2:
      flag.newFlag(`dual${v2}`, midName(v2,v3), midName(v1,v2));
      // shift over one
      [v1, v2] = Array.from([v2, v3]);
    }
  }

  const newpoly = flag.topoly();
  newpoly.name = `a${poly.name}`;
  //newpoly.xyz = adjustXYZ(newpoly, 2)
  return newpoly;
};


// Gyro
// ----------------------------------------------------------------------------------------------
// This is the dual operator to "snub", i.e dual*Gyro = Snub.  It is a bit easier to implement
// this way.
//
// Snub creates at each vertex a new face, expands and twists it, and adds two new triangles to
// replace each edge.
//
const gyro = function(poly){
  let f, i, v;
  console.log(`Taking gyro of ${poly.name}...`);

  const flag = new polyflag();

  for (i = 0; i < poly.xyz.length; i++) {
    v = poly.xyz[i];
    flag.newV(`v${i}`, unit(v));
  }  // each old vertex is a new vertex

  const centers = poly.centers(); // new vertices in center of each face
  for (i = 0; i < poly.face.length; i++) {
    f = poly.face[i];
    flag.newV(`center${i}`, unit(centers[i]));
  }

  for (i = 0; i < poly.face.length; i++) {
    f = poly.face[i];
    let [v1, v2] = Array.from(f.slice(-2));
    for (let j = 0; j < f.length; j++) {
      v = f[j];
      const v3 = v;
      flag.newV(v1+"~"+v2, oneThird(poly.xyz[v1],poly.xyz[v2]));  // new v in face
      const fname = i+"f"+v1;
      flag.newFlag(fname, `center${i}`,      v1+"~"+v2); // five new flags
      flag.newFlag(fname, v1+"~"+v2,  v2+"~"+v1);
      flag.newFlag(fname, v2+"~"+v1,  `v${v2}`);
      flag.newFlag(fname, `v${v2}`,     v2+"~"+v3);
      flag.newFlag(fname, v2+"~"+v3,  `center${i}`);
      [v1, v2] = Array.from([v2, v3]);
    }
  }                       // shift over one

  const newpoly = flag.topoly();
  newpoly.name = `g${poly.name}`;
  //newpoly.xyz = adjustXYZ(newpoly, 3)
  return newpoly;
};


// Propellor
// ------------------------------------------------------------------------------------------
// builds a new 'skew face' by making new points along edges, 1/3rd the distance from v1->v2,
// then connecting these into a new inset face.  This breaks rotational symmetry about the
// faces, whirling them into gyres
//
const propellor = function(poly) {
  let i, v;
  console.log(`Taking propellor of ${poly.name}...`);

  const flag = new polyflag();

  for (i = 0; i < poly.xyz.length; i++) {
    v = poly.xyz[i];
    flag.newV(`v${i}`, unit(v));
  }  // each old vertex is a new vertex

  for (i = 0; i < poly.face.length; i++) {
    const f = poly.face[i];
    let [v1, v2] = Array.from(f.slice(-2));
    for (v of f) {
      const v3 = `${v}`;
      flag.newV(v1+"~"+v2, oneThird(poly.xyz[v1], poly.xyz[v2]));  // new v in face, 1/3rd along edge
      const fname = `${i}f${v2}`;
      flag.newFlag(`v${i}`, v1+"~"+v2,  v2+"~"+v3); // five new flags
      flag.newFlag(fname,   v1+"~"+v2,  v2+"~"+v1);
      flag.newFlag(fname,   v2+"~"+v1,     `v${v2}`);
      flag.newFlag(fname,      `v${v2}`,  v2+"~"+v3);
      flag.newFlag(fname,   v2+"~"+v3,  v1+"~"+v2);
      [v1, v2] = Array.from([v2, v3]);
    }
  }                       // shift over one

  const newpoly = flag.topoly();
  newpoly.name = `p${poly.name}`;
  //newpoly.xyz  = adjustXYZ(newpoly, 3)
  return newpoly;
};


// Reflection
// ------------------------------------------------------------------------------------------
// geometric reflection through origin
const reflect = function(poly) {
  let i;
  let asc, end;
  let asc1, end1;
  console.log(`Taking reflection of ${poly.name}...`);
  for (i = 0, end = poly.xyz.length-1, asc = 0 <= end; asc ? i <= end : i >= end; asc ? i++ : i--) {
     poly.xyz[i] = mult(-1, poly.xyz[i]);
  }         // reflect each point through origin
  for (i = 0, end1 = poly.face.length-1, asc1 = 0 <= end1; asc1 ? i <= end1 : i >= end1; asc1 ? i++ : i--) {
     poly.face[i] = poly.face[i].reverse();
  }       // repair clockwise-ness of faces!
  poly.name = `r${poly.name}`;
  return poly;
};


// Dual
// ------------------------------------------------------------------------------------------------
// The dual of a polyhedron is another mesh wherein:
// - every face in the original becomes a vertex in the dual
// - every vertex in the original becomes a face in the dual
//
// So N_faces, N_vertices = N_dualfaces, N_dualvertices
//
// The new vertex coordinates are convenient to set to the original face centroids.
//
const dual = function(poly) {
  let f, i, v1, v2;
  let asc, end;
  let asc1, end1;
  console.log(`Taking dual of ${poly.name}...`);

  const flag = new polyflag();

  const face = []; // make table of face as fn of edge
  for (i = 0, end = poly.xyz.length-1, asc = 0 <= end; asc ? i <= end : i >= end; asc ? i++ : i--) {
    face[i] = {};
  } // create empty associative table

  for (i = 0; i < poly.face.length; i++) {
    f = poly.face[i];
    v1 = f[f.length-1]; //previous vertex
    for (v2 of f) {
      // THIS ASSUMES that no 2 faces that share an edge share it in the same orientation!
      // which of course never happens for proper manifold meshes, so get your meshes right.
      face[v1][`v${v2}`] = `${i}`;
      v1=v2;
    }
  } // current becomes previous

  const centers = poly.centers();
  for (i = 0, end1 = poly.face.length-1, asc1 = 0 <= end1; asc1 ? i <= end1 : i >= end1; asc1 ? i++ : i--) {
    flag.newV(`${i}`,centers[i]);
  }

  for (i = 0; i < poly.face.length; i++) {
    f = poly.face[i];
    v1 = f[f.length-1]; //previous vertex
    for (v2 of f) {
      flag.newFlag(v1, face[v2][`v${v1}`], `${i}`);
      v1=v2;
    }
  } // current becomes previous

  const dpoly = flag.topoly(); // build topological dual from flags

  // match F index ordering to V index ordering on dual
  const sortF = [];
  for (f of dpoly.face) {
    const k = intersect(poly.face[f[0]],poly.face[f[1]],poly.face[f[2]]);
    sortF[k] = f;
  }
  dpoly.face = sortF;

  if (poly.name[0] !== "d") {
    dpoly.name = `d${poly.name}`;
  } else {
    dpoly.name = poly.name.slice(1);
  }

  return dpoly;
};


// Chamfer
// ----------------------------------------------------------------------------------------
// A truncation along a polyhedron's edges.
// Chamfering or edge-truncation is similar to expansion, moving faces apart and outward,
// but also maintains the original vertices. Adds a new hexagonal face in place of each
// original edge.
// A polyhedron with e edges will have a chamfered form containing 2e new vertices,
// 3e new edges, and e new hexagonal faces. -- Wikipedia
// See also http://dmccooey.com/polyhedra/Chamfer.html
//
// The dist parameter could control how deeply to chamfer.
// But I'm not sure about implementing that yet.
//
// Q: what is the dual operation of chamfering? I.e.
// if cX = dxdX, and xX = dcdX, what operation is x?

// We could "almost" do this in terms of already-implemented operations:
// cC = t4daC = t4jC, cO = t3daO, cD = t5daD, cI = t3daI
// But it doesn't work for cases like T.

const chamfer = function(poly, dist) {
  console.log(`Taking chamfer of ${poly.name}...`);

  if (!dist) { dist = 0.5; }

  const flag = new polyflag();

  const normals = poly.normals();

  // For each face f in the original poly
  for (let i = 0; i < poly.face.length; i++) {
    const f = poly.face[i];
    let v1 = f[f.length-1];
    let v1new = i + "_" + v1;

    for (let v2 of f) {
      // TODO: figure out what distances will give us a planar hex face.
      // Move each old vertex further from the origin.
      flag.newV(v2, mult(1.0 + dist, poly.xyz[v2]));
      // Add a new vertex, moved parallel to normal.
      const v2new = i + "_" + v2;
      flag.newV(v2new, add(poly.xyz[v2], mult(dist*1.5, normals[i])));
      // Four new flags:
      // One whose face corresponds to the original face:
      flag.newFlag(`orig${i}`, v1new, v2new);
      // And three for the edges of the new hexagon:
      const facename = (v1<v2 ? `hex${v1}_${v2}` : `hex${v2}_${v1}`);
      flag.newFlag(facename, v2, v2new);
      flag.newFlag(facename, v2new, v1new);
      flag.newFlag(facename, v1new, v1);
      v1 = v2;
      v1new = v2new;
    }
  }

  const newpoly = flag.topoly();
  newpoly.name = `c${poly.name}`;
  return newpoly;
};


// Whirl
// ----------------------------------------------------------------------------------------------
// Gyro followed by truncation of vertices centered on original faces.
// This create 2 new hexagons for every original edge.
// (https://en.wikipedia.org/wiki/Conway_polyhedron_notation#Operations_on_polyhedra)
//
// Possible extension: take a parameter n that means only whirl n-sided faces.
// If we do that, the flags marked #* below will need to have their other sides
// filled in one way or another, depending on whether the adjacent face is
// whirled or not.

const whirl = function(poly, n) {
  let i, v;
  console.log(`Taking whirl of ${poly.name}...`);
  if (!n) { n = 0; }
  
  const flag = new polyflag();

  for (i = 0; i < poly.xyz.length; i++) {
    v = poly.xyz[i];
    flag.newV(`v${i}`, unit(v));
  }  // each old vertex is a new vertex

  const centers = poly.centers(); // new vertices around center of each face
  //for f,i in poly.face
  //  # Whirl: use "center"+i+"~"+v1
  //  flag.newV "center"+i+"~"+v1, unit(centers[i])

  for (i = 0; i < poly.face.length; i++) {
    const f = poly.face[i];
    let [v1, v2] = Array.from(f.slice(-2));
    for (let j = 0; j < f.length; j++) {
      v = f[j];
      const v3 = v;
      // New vertex along edge
      const v1_2 = oneThird(poly.xyz[v1],poly.xyz[v2]);
      flag.newV(v1+"~"+v2, v1_2);
      // New vertices near center of face
      const cv1name = `center${i}~${v1}`;
      const cv2name = `center${i}~${v2}`;
      flag.newV(cv1name, unit(oneThird(centers[i], v1_2))); 
      const fname = i+"f"+v1;
      // New hexagon for each original edge
      flag.newFlag(fname, cv1name,      v1+"~"+v2);
      flag.newFlag(fname, v1+"~"+v2,  v2+"~"+v1); //*
      flag.newFlag(fname, v2+"~"+v1,  `v${v2}`);    //*
      flag.newFlag(fname, `v${v2}`,     v2+"~"+v3); //*
      flag.newFlag(fname, v2+"~"+v3,  cv2name);
      flag.newFlag(fname, cv2name, cv1name);
      // New face in center of each old face      
      flag.newFlag(`c${i}`, cv1name, cv2name);
      
      [v1, v2] = Array.from([v2, v3]);
    }
  }                       // shift over one

  const newpoly = flag.topoly();
  newpoly.name = `w${poly.name}`;
  //newpoly.xyz = adjustXYZ(newpoly, 3)
  return newpoly;
};
  

// insetN
// ------------------------------------------------------------------------------------------
const insetN = function(poly, n, inset_dist, popout_dist){
  let f, i, v;
  if (!n) { n = 0; }
  if (!inset_dist) { inset_dist = 0.5; }
  if (!popout_dist) { popout_dist = -0.2; }

  console.log(`Taking inset of ${n===0 ? "" : n}-sided faces of ${poly.name}...`);

  const flag = new polyflag();
  for (i = 0; i < poly.xyz.length; i++) {
    // each old vertex is a new vertex
    const p = poly.xyz[i];
    flag.newV(`v${i}`, p);
  }

  const normals = poly.normals();
  const centers = poly.centers();
  for (i = 0; i < poly.face.length; i++) { //new inset vertex for every vert in face
    f = poly.face[i];
    if ((f.length === n) || (n === 0)) {
      for (v of f) {
        flag.newV(`f${i}v${v}`, add(tween(poly.xyz[v],centers[i],inset_dist),mult(popout_dist,normals[i])));
      }
    }
  }

  let foundAny = false;                 // alert if don't find any
  for (i = 0; i < poly.face.length; i++) {
    f = poly.face[i];
    let v1 = `v${f[f.length-1]}`;
    for (v of f) {
      const v2 = `v${v}`;
      if ((f.length === n) || (n === 0)) {
        foundAny = true;
        const fname = i + v1;
        flag.newFlag(fname,      v1,       v2);
        flag.newFlag(fname,      v2,       `f${i}${v2}`);
        flag.newFlag(fname, `f${i}${v2}`,  `f${i}${v1}`);
        flag.newFlag(fname, `f${i}${v1}`,  v1);
        //new inset, extruded face
        flag.newFlag(`ex${i}`, `f${i}${v1}`,  `f${i}${v2}`);
      } else {
        flag.newFlag(i, v1, v2);  // same old flag, if non-n
      }
      v1=v2;
    }
  }  // current becomes previous

  if (!foundAny) {
    console.log(`No ${n}-fold components were found.`);
  }

  const newpoly = flag.topoly();
  newpoly.name = `n${n === 0 ? "" : n}${poly.name}`;
  //newpoly.xyz = adjustXYZ(newpoly, 3)
  //newpoly.xyz = canonicalXYZ(newpoly, 3)  # this tends to make results look like shit
  return newpoly;
};

// ExtrudeN
// ------------------------------------------------------------------------------------------
const extrudeN = function(poly, n){
  let f, i, v;
  if (!n) { n = 0; }
  console.log(`Taking extrusion of ${n===0 ? "" : n}-sided faces of ${poly.name}...`);

  const flag = new polyflag();
  for (i = 0; i < poly.xyz.length; i++) {
    // each old vertex is a new vertex
    const p = poly.xyz[i];
    flag.newV(`v${i}`, p);
  }

  const normals = poly.normals();
  const centers = poly.centers();
  for (i = 0; i < poly.face.length; i++) { //new inset vertex for every vert in face
    f = poly.face[i];
    if ((f.length === n) || (n === 0)) {
      for (v of f) {
        //flag.newV "f"+i+"v"+v, add(midpoint(poly.xyz[v],centers[i]),mult(-0.2,normals[i]))
        flag.newV(`f${i}v${v}`, add(poly.xyz[v], mult(0.3,normals[i])));
      }
    }
  }

  let foundAny = false;                 // alert if don't find any
  for (i = 0; i < poly.face.length; i++) {
    f = poly.face[i];
    let v1 = `v${f[f.length-1]}`;
    for (v of f) {
      const v2 = `v${v}`;
      if ((f.length === n) || (n === 0)) {
        foundAny = true;
        //fname = i+v1
        flag.newFlag(i+v1,       v1,       v2);
        flag.newFlag(i+v1,       v2, `f${i}${v2}`);
        flag.newFlag(i+v1, `f${i}${v2}`, `f${i}${v1}`);
        flag.newFlag(i+v1, `f${i}${v1}`,       v1);
        //new inset, extruded face
        flag.newFlag(`ex${i}`, `f${i}${v1}`,  `f${i}${v2}`);
      } else {
        flag.newFlag(i, v1, v2);  // same old flag, if non-n
      }
      v1=v2;
    }
  }  // current becomes previous

  if (!foundAny) {
    console.log(`No ${n}-fold components were found.`);
  }

  const newpoly = flag.topoly();
  newpoly.name = `x${n === 0 ? "" : n}${poly.name}`;
  //console.log newpoly
  //newpoly.xyz = adjustXYZ(newpoly, 3)
  //newpoly.xyz = canonicalXYZ(newpoly, 3)  # this tends to make results look like shit
  return newpoly;
};

// hollow / skeletonize
// ------------------------------------------------------------------------------------------
const hollow = function(poly, n, inset_dist, thickness){
  let f, i, v;
  if (!n) { n = 0; }
  if (!inset_dist) { inset_dist = 0.5; }
  if (!thickness) { thickness = 0.2; }

  console.log(`Skeletonizing ${n===0 ? "" : n}-sided faces of ${poly.name}...`);

  const dualnormals = dual(poly).normals();
  const normals = poly.normals();
  const centers = poly.centers();

  const flag = new polyflag();
  for (i = 0; i < poly.xyz.length; i++) {
    // each old vertex is a new vertex
    const p = poly.xyz[i];
    flag.newV(`v${i}`, p);
    flag.newV(`downv${i}`,  add(p,mult(-1*thickness,dualnormals[i])));
  }

  for (i = 0; i < poly.face.length; i++) { //new inset vertex for every vert in face
    //if f.length is n or n is 0
    f = poly.face[i];
    for (v of f) {
      flag.newV(`fin${i}v${v}`, tween(poly.xyz[v],centers[i],inset_dist));
      flag.newV(`findown${i}v${v}`, add(tween(poly.xyz[v],centers[i],inset_dist),mult(-1*thickness,normals[i])));
    }
  }

  //foundAny = false                 # alert if don't find any
  for (i = 0; i < poly.face.length; i++) {
    f = poly.face[i];
    let v1 = `v${f[f.length-1]}`;
    for (v of f) {
      const v2 = `v${v}`;
      //if f.length is n or n is 0
      const foundAny = true;
      let fname = i + v1;
      flag.newFlag(fname,      v1,       v2);
      flag.newFlag(fname,      v2,       `fin${i}${v2}`);
      flag.newFlag(fname, `fin${i}${v2}`,  `fin${i}${v1}`);
      flag.newFlag(fname, `fin${i}${v1}`,  v1);

      fname = `sides${i}${v1}`;
      flag.newFlag(fname, `fin${i}${v1}`,     `fin${i}${v2}`);
      flag.newFlag(fname, `fin${i}${v2}`,     `findown${i}${v2}`);
      flag.newFlag(fname, `findown${i}${v2}`, `findown${i}${v1}`);
      flag.newFlag(fname, `findown${i}${v1}`, `fin${i}${v1}`);

      fname = `bottom${i}${v1}`;
      flag.newFlag(fname,  `down${v2}`,      `down${v1}`);
      flag.newFlag(fname,  `down${v1}`,      `findown${i}${v1}`);
      flag.newFlag(fname,  `findown${i}${v1}`, `findown${i}${v2}`);
      flag.newFlag(fname,  `findown${i}${v2}`, `down${v2}`);

      //new inset, extruded face
      //flag.newFlag "ex"+i, "fin"+i+v1,  "fin"+i+v2
      //else
      //  flag.newFlag i, v1, v2  # same old flag, if non-n
      v1=v2;
    }
  }  // current becomes previous

  //if not foundAny
  //  console.log "No #{n}-fold components were found."

  const newpoly = flag.topoly();
  //newpoly.name = "h" + (if n is 0 then "" else n) + poly.name
  newpoly.name = `h${poly.name}`;
  //newpoly.xyz = adjustXYZ(newpoly, 3)
  //newpoly.xyz = canonicalXYZ(newpoly, 3)  # this tends to make results look like shit
  return newpoly;
};


// StellaN
// ------------------------------------------------------------------------------------------
const stellaN = function(poly){
  let i;
  console.log(`Taking stella of ${poly.name}...`);

  const centers = poly.centers();  // calculate face centers

  const flag = new polyflag();
  for (i = 0; i < poly.xyz.length; i++) {
    const p = poly.xyz[i];
    flag.newV(`v${i}`, p);
  }      // each old vertex is a new vertex

  // iterate over triplets of faces v1,v2,v3
  for (i = 0; i < poly.face.length; i++) {
    const f = poly.face[i];
    let v1 = `v${f[f.length-2]}`;
    let v2 = `v${f[f.length-1]}`;
    let vert1 = poly.xyz[f[f.length-2]];
    let vert2 = poly.xyz[f[f.length-1]];
    for (let v of f) {
      const v3 = `v${v}`;
      const vert3 = poly.xyz[v];
      const v12=v1+"~"+v2; // names for "oriented" midpoints
      const v21=v2+"~"+v1;
      const v23=v2+"~"+v3;

      // on each Nface, N new points inset from edge midpoints towards center = "stellated" points
      flag.newV(v12, midpoint( midpoint(vert1,vert2), centers[i] ));

      // inset Nface made of new, stellated points
      flag.newFlag(`in${i}`,      v12,       v23);

      // new tri face constituting the remainder of the stellated Nface
      flag.newFlag(`f${i}${v2}`,      v23,      v12);
      flag.newFlag(`f${i}${v2}`,       v12,      v2);
      flag.newFlag(`f${i}${v2}`,      v2,      v23);

      // one of the two new triangles replacing old edge between v1->v2
      flag.newFlag(`f${v12}`,     v1,        v21);
      flag.newFlag(`f${v12}`,     v21,       v12);
      flag.newFlag(`f${v12}`,      v12,       v1);

      [v1,v2]=Array.from([v2,v3]);  // current becomes previous
      [vert1,vert2]=Array.from([vert2,vert3]);
    }
  }

  const newpoly = flag.topoly();
  newpoly.name = `l${poly.name}`;
  //newpoly.xyz = adjustXYZ(newpoly, 3)
  //newpoly.xyz = canonicalXYZ(newpoly, 3)  # this tends to make results look like shit
  return newpoly;
};
/*
 * decaffeinate suggestions:
 * DS102: Remove unnecessary code created because of implicit returns
 * DS201: Simplify complex destructure assignments
 * DS202: Simplify dynamic range loops
 * DS205: Consider reworking code to avoid use of IIFEs
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


//===================================================================================================
// Canonicalization Algorithms
//===================================================================================================

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
const tangentify = function(xyzs, edges) {
  // hack to improve convergence
  const STABILITY_FACTOR = 0.1; 
  // copy vertices
  const newVs = copyVecArray(xyzs);
  for (let e of edges) {
    // the point closest to origin
    const t = tangentPoint( newVs[e[0]], newVs[e[1]] );
    // adjustment from sphere
    const c = mult(((STABILITY_FACTOR*1)/2)*(1-sqrt(dot(t,t))), t);
    newVs[e[0]] = add(newVs[e[0]], c);
    newVs[e[1]] = add(newVs[e[1]], c);
  }
  return newVs;
};

// recenters entire polyhedron such that center of mass is at origin
const recenter = function(xyzs, edges) {
  const edgecenters = ((() => {
    const result = [];
    for (let [a,b] of edges) {
      result.push(tangentPoint(xyzs[a], xyzs[b]));
    }
    return result;
  })()); //centers of edges
  let polycenter = [0,0,0];
  // sum centers to find center of gravity
  for (let v of edgecenters) { 
    polycenter = add(polycenter, v);
  }
  polycenter = mult(1/edges.length, polycenter);
  // subtract off any deviation from center
  return _.map(xyzs, x=> sub(x, polycenter));
};

// rescales maximum radius of polyhedron to 1
const rescale = function(xyzs) {
  const polycenter = [0,0,0];
  const maxExtent = _.max(_.map(xyzs, x=>mag(x)));
  const s = 1 / maxExtent;
  return _.map(xyzs, x=>[s*x[0],s*x[1],s*x[2]]);
};

// adjusts vertices in each face to improve its planarity
const planarize = function(xyzs, faces) {
  let v;
  const STABILITY_FACTOR = 0.1; // Hack to improve convergence
  const newVs = copyVecArray(xyzs); // copy vertices
  for (var f of faces) {
    const coords = ((() => {
      const result = [];
      for (v of f) {         
        result.push(xyzs[v]);
      }
      return result;
    })());
    let n = normal(coords); // find avg of normals for each vertex triplet
    const c = calcCentroid(coords); // find planar centroid
    if (dot(n,c) < 0) { // correct sign if needed
      n = mult(-1.0,n);
    }
    for (v of f) {  // project (vertex - centroid) onto normal, subtract off this component
      newVs[v] = add(newVs[v], mult(dot(mult(STABILITY_FACTOR, n), sub(c, xyzs[v])), n) );
    }
  }
  return newVs;
};

// combines above three constraint adjustments in iterative cycle
const canonicalize = function(poly, Niter) {
  if (!Niter) { Niter = 1; }
  console.log(`Canonicalizing ${poly.name}...`);
  const faces = poly.face;
  const edges = poly.edges();
  let newVs = poly.xyz;
  let maxChange=1.0; // convergence tracker
  for (let i = 0, end = Niter, asc = 0 <= end; asc ? i <= end : i >= end; asc ? i++ : i--) {
    const oldVs = copyVecArray(newVs); //copy vertices
    newVs = tangentify(newVs, edges);
    newVs = recenter(newVs, edges);
    newVs = planarize(newVs, faces);
    maxChange = _.max(_.map(_.zip(newVs,oldVs), 
    function(...args){ 
      const [x,y] = args[0]; 
      return mag(sub(x,y));  
    }));
    if (maxChange < 1e-8) {
      break;
    }
  }
  // one should now rescale, but not rescaling here makes for very interesting numerical
  // instabilities that make interesting mutants on multiple applications...
  // more experience will tell what to do
  //newVs = rescale(newVs)
  console.log(`[canonicalization done, last |deltaV|=${maxChange}]`);
  const newpoly = new polyhedron(newVs, poly.face, poly.name);
  console.log("canonicalize" , newpoly);
  return newpoly;
};

// Hacky Canonicalization Algorithm
// --------------------------------------------------------------------
// Using center of gravity of vertices for each face to planarize faces

// get the spherical reciprocals of face centers
const reciprocalC = function(poly) {
  const centers = poly.centers();
  for (let c of centers) {
    c = mult(1.0/dot(c,c), c);
  }
  return centers;
};

// make array of vertices reciprocal to given planes
const reciprocalN = function(poly) {
  const ans = [];
  for (let f of poly.face) { //for each face
    let centroid    = [0,0,0]; // running sum of vertex coords
    let normalV     = [0,0,0]; // running sum of normal vectors
    let avgEdgeDist =    0.0;  // running sum for avg edge distance

    let [v1, v2] = f.slice(-2);
    for (let v3 of f) {
      centroid     = add(centroid, poly.xyz[v3]);
      normalV      = add(normalV, orthogonal(poly.xyz[v1], poly.xyz[v2], poly.xyz[v3]));
      avgEdgeDist += edgeDist(poly.xyz[v1], poly.xyz[v2]);
      [v1, v2] = [v2, v3];
    } // shift over one

    centroid    = mult(1.0/f.length, centroid);
    normalV     = unit(normalV);
    avgEdgeDist = avgEdgeDist / f.length;
    const tmp   = reciprocal(mult(dot(centroid, normalV), normalV)); // based on face
    ans.push(mult((1 + avgEdgeDist) / 2, tmp));
  } // edge correction

  return ans;
};

const canonicalXYZ = function(poly, nIterations) {
  if (!nIterations) { nIterations = 1; }
  const dpoly = dual(poly);
  console.log(`Pseudo-canonicalizing ${poly.name}...`);

  // iteratively reciprocate face normals
  for (let count = 0, end = nIterations; count < end; count++) {
    dpoly.xyz = reciprocalN(poly);
    poly.xyz  = reciprocalN(dpoly);
  }

  return new polyhedron(poly.xyz, poly.face, poly.name);
};


// quick planarization
const adjustXYZ = function(poly, nIterations) {
  if (!nIterations) { nIterations = 1; }
  const dpoly = dual(poly); // v's of dual are in order of arg's f's
  console.log(`Planarizing ${poly.name}...`);

  for (let count = 0, end = nIterations; count < end; count++) {
    // reciprocate face centers
    dpoly.xyz = reciprocalC(poly);
    poly.xyz  = reciprocalC(dpoly);
  }

  return new polyhedron(poly.xyz, poly.face, poly.name);
};


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
// Copyright 2019, Anselm Levskaya
// Released under the MIT License
//

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
}/*
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
// Copyright 2019, Anselm Levskaya
// Released under the MIT License

// Parser Routines
//===================================================================================================

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
const op_parser = PEG.buildParser(PEG_parser_spec);

//applies func fn to array args i.e. f, [1,2,3] -> f(1,2,3)
const dispatch = function(fn, args) { return fn.apply(this, args || []); };

const basemap = {
  "T": tetrahedron,
  "O": octahedron,
  "C": cube,
  "I": icosahedron,
  "D": dodecahedron,
  "P": prism,     //takes integer arg
  "A": antiprism, //takes integer arg
  "Y": pyramid   //takes integer arg
  };

const opmap = {
  "d": dual,
  "k": kisN,
  "a": ambo,
  "g": gyro,
  "p": propellor,
  "r": reflect,
  "h": hollow,
  "c": chamfer,
  "w": whirl,
  "n": insetN,
  "x": extrudeN,
  "l": stellaN,
  "z": triangulate,
  "K": canonicalXYZ,
  "C": canonicalize,
  "A": adjustXYZ
  };

//list of basic equivalences, easier to replace before parsing
const specreplacements = [
  [/e/g, "aa"],   // e --> aa   (abbr. for explode)
  [/b/g, "ta"],   // b --> ta   (abbr. for bevel)
  [/o/g, "jj"],   // o --> jj   (abbr. for ortho)
  [/m/g, "kj"],   // m --> kj   (abbr. for meta)
  [/t(\d*)/g, "dk$1d"],  // t(n) --> dk(n)d  (dual operations)
  [/j/g, "dad"],  // j --> dad  (dual operations) # Why not j --> da ?
  [/s/g, "dgd"],  // s --> dgd  (dual operations) # Why not s --> dg ?
  [/dd/g, ""],    // dd --> null  (order 2)
  [/ad/g, "a"],   // ad --> a   (a_ = ad_)
  [/gd/g, "g"],   // gd --> g   (g_ = gd_)
  [/aO/g, "aC"],  // aO --> aC  (for uniqueness)
  [/aI/g, "aD"],  // aI --> aD  (for uniqueness)
  [/gO/g, "gC"],  // gO --> gC  (for uniqueness)
  [/gI/g, "gD"]];  // gI --> gD  (for uniqueness)

const getOps = function(notation) {
  let expanded = notation;
  for (let [orig,equiv] of specreplacements) {
    expanded = expanded.replace(orig,equiv);
  }
  console.log(`${notation} executed as ${expanded}`);
  return expanded;
};

// create polyhedron from notation
const newgeneratePoly = function(notation) {
  //poly = new polyhedron()

  const ops_spec = getOps(notation);
  const oplist = op_parser.parse(ops_spec).reverse();

  let op = oplist.shift();
  const basefunc = basemap[op["op"]];
  const baseargs = op["args"];
  let poly     = dispatch(basefunc, baseargs);

  //console.log "base", poly

  for (op of oplist) {
    const opfunc = opmap[op["op"]];
    const opargs = [poly].concat(op["args"]);
    //console.log opargs
    poly   = dispatch(opfunc, opargs);
  }

  //console.log "final", poly

  // Recenter polyhedra at origin (rarely needed)
  poly.xyz = recenter(poly.xyz, poly.edges());
  poly.xyz = rescale(poly.xyz);

  // Color the faces of the polyhedra for display
  poly = paintPolyhedron(poly);

  // return the poly object
  return poly;
};
/*
 * decaffeinate suggestions:
 * DS102: Remove unnecessary code created because of implicit returns
 * DS205: Consider reworking code to avoid use of IIFEs
 * Full docs: https://github.com/decaffeinate/decaffeinate/blob/master/docs/suggestions.md
 */
// Polyhédronisme
//===================================================================================================
//
// A toy for constructing and manipulating polyhedra and other meshes
//
// Copyright 2019, Anselm Levskaya
// Released under the MIT License


// GLOBALS
//===================================================================================================
let ctx={}; // for global access to canvas context
const CANVAS_WIDTH  = 500; //canvas dims
const CANVAS_HEIGHT = 400; //canvas dims
let globPolys={}; // constructed polyhedras

let globRotM = clone(eye3);
let globlastRotM = clone(eye3);
//globtheta = 0 # rotation and projective mapping parameters
//globphi   = 0
let perspective_scale = 800;
const persp_z_max = 5;
const persp_z_min = 0;
const persp_ratio = 0.8;
const _2d_x_offset = CANVAS_WIDTH/2; //300
const _2d_y_offset = CANVAS_HEIGHT/2; //140

const globtime = new Date(); // for animation

const BG_CLEAR = true; // clear background or colored?
const BG_COLOR = "rgba(255,255,255,1.0)"; // background color
const COLOR_METHOD = "signature"; //"area"
let PaintMode = "fillstroke";
const ctx_linewidth = 0.5; // for outline of faces

// Mouse Event Variables
let MOUSEDOWN=false;
let LastMouseX=0;
let LastMouseY=0;
let LastSphVec=[1,0,0]; //for 3d trackball

// random grabbag of polyhedra
const DEFAULT_RECIPES = [
  "C2dakD","oC20kkkT","kn4C40A0dA4","opD",
  "lT","lK5oC","knD","dn6x4K5bT","oox4P7",
  "n18n18n9n9n9soxY9"];

// File-saving objects used to export txt/canvas-png
// const saveText = function(text, filename) {
//   const BB = window.BlobBuilder || window.WebKitBlobBuilder || window.MozBlobBuilder;
//   const bb = new BB();
//   bb.append(text);
//   return saveAs(bb.getBlob(`text/plain;charset=${document.characterSet}`), filename);
// };
const saveText = function(text, filename) {
  const blb = new Blob([text], 
    {type: `text/plain;charset=${document.characterSet}`});
  return saveAs(blb, filename);
}

// parses URL string for polyhedron recipe, for bookmarking
// should use #! href format instead
const parseurl = function() {
  let e;
  const urlParams = {};
  const a = /\+/g;  // Regex for replacing addition symbol with a space
  const r = /([^&=]+)=?([^&]*)/g;
  const d = s => decodeURIComponent(s.replace(a, " "));
  const q = window.location.search.substring(1);

  while ((e=r.exec(q))) {
    urlParams[d(e[1])] = d(e[2]);
  }
  return urlParams;
};

// update the shareable link URL with the current recipe and palette
const setlink = function() {
  const specs = $("#spec").val().split(/\s+/g).slice(0, 2);
  // strip any existing parameters
  let link = location.protocol + '//' + location.host + location.pathname;
  link += `?recipe=${encodeURIComponent(specs[0])}`;
  if (PALETTE !== rwb_palette) {
    link += `&palette=${encodeURIComponent(PALETTE.reduce((x,y)=> x+" "+y))}`;
  }
  return $("#link").attr("href", link);
};


// Drawing Functions
//==================================================================================================

// init canvas element
// -------------------------------------------------------------------------------
const init = function() {
  const canvas = $('#poly');
  canvas.width(CANVAS_WIDTH);
  canvas.height(CANVAS_HEIGHT);

  ctx = canvas[0].getContext("2d");
  ctx.lineWidth = ctx_linewidth;

  if (BG_CLEAR) {
    ctx.clearRect(0, 0, CANVAS_WIDTH, CANVAS_HEIGHT);
  } else {
    ctx.clearRect(0, 0, CANVAS_WIDTH, CANVAS_HEIGHT);
    ctx.fillStyle = BG_COLOR;
    ctx.fillRect(0, 0, CANVAS_WIDTH, CANVAS_HEIGHT);
  }
    
  const exp = $('#expandcollapse');
  return exp.click(function() {
    if (/minus/.test(exp.attr('src'))) {  // Contains 'minus'
      $('#morestats').hide();
      return exp.attr('src', 'media/plus.png');
    } else {
      $('#morestats').show();      
      return exp.attr('src', 'media/minus.png');
    }
  });
};

// clear canvas
// -----------------------------------------------------------------------------------
const clear = function() {
  if (BG_CLEAR) {
    return ctx.clearRect(0, 0, CANVAS_WIDTH, CANVAS_HEIGHT);
  } else {
    ctx.clearRect(0, 0, CANVAS_WIDTH, CANVAS_HEIGHT);
    ctx.fillStyle = BG_COLOR;
    return ctx.fillRect(0, 0, CANVAS_WIDTH, CANVAS_HEIGHT);
  }
};


// main drawing routine for polyhedra
//===================================================================================================
const drawpoly = function(poly, tvec) {
  let v;
  if (!tvec) { tvec = [3, 3, 3]; }


  // rotate poly in 3d
  const oldxyz = _.map(poly.xyz, x=> x);
  poly.xyz = _.map(poly.xyz, x=> mv3(globRotM,x));

  // z sort faces
  sortfaces(poly);

  for (let fno = 0; fno < poly.face.length; fno++) {
    var face = poly.face[fno];
    ctx.beginPath();
    // move to first vertex of face
    const v0 = face[face.length-1];
    let [x,y] = perspT(add(tvec,poly.xyz[v0]), persp_z_max,persp_z_min,persp_ratio,perspective_scale);
    ctx.moveTo(x+_2d_x_offset, y+_2d_y_offset);
    // loop around face, defining polygon
    for (v of face) {
      [x,y] = perspT(add(tvec,poly.xyz[v]),persp_z_max,persp_z_min,persp_ratio,perspective_scale);
      ctx.lineTo(x+_2d_x_offset, y+_2d_y_offset);
    }

    // use pre-computed colors
    let clr = palette(poly.face_class[fno]);

    // shade based on simple cosine illumination factor
    const face_verts = ((() => {
      const result = [];
      for (v of face) {
        result.push(poly.xyz[v]);
      }
      return result;
    })());
    const illum = dot(normal(face_verts), unit([1, -1, 0]));
    clr = mult((((illum / 2.0) + 0.5) * 0.7) + 0.3, clr);

    if ((PaintMode === "fill") || (PaintMode === "fillstroke")) {
      ctx.fillStyle = `rgba(${round(clr[0]*255)}, ${round(clr[1]*255)}, ${round(clr[2]*255)}, ${1.0})`;
      ctx.fill();
    // make cartoon stroke (=black) / realistic stroke an option (=below)
      ctx.strokeStyle = `rgba(${round(clr[0]*255)}, ${round(clr[1]*255)}, ${round(clr[2]*255)}, ${1.0})`;
      ctx.stroke();
    }
    if (PaintMode === "fillstroke") {
      ctx.fillStyle = `rgba(${round(clr[0]*255)}, ${round(clr[1]*255)}, ${round(clr[2]*255)}, ${1.0})`;
      ctx.fill();
      ctx.strokeStyle = "rgba(0,0,0, .3)";  // light lines, less cartoony, more render-y
      ctx.stroke();
    }
    if (PaintMode === "stroke") {
      ctx.strokeStyle = "rgba(0,0,0, .8)";
      ctx.stroke();
    }
  }

  // reset coords, for setting absolute rotation, as poly is passed by ref
  return poly.xyz = oldxyz;
};


// draw polyhedra just once
// -----------------------------------------------------------------------------------
const drawShape = function() {
  clear();
  return globPolys.map((p, i) =>
    drawpoly(p,[0+(3*i),0,3]));
};

// update V E F stats on page
// -----------------------------------------------------------------------------------
const updateStats = () =>
  (() => {
    const result = [];
    for (let i = 0; i < globPolys.length; i++) {
      const p = globPolys[i];
      $("#basicstats").text(p.data());
      result.push($("#morestats").text(p.moreData()));
    }
    return result;
  })()
;

// loop for animation
// -----------------------------------------------------------------------------------
var animateShape = function() {
  clear();
  const globtheta=((2*Math.PI)/180.0)*globtime.getSeconds()*0.1;
  for (let i = 0; i < globPolys.length; i++) {
    const p = globPolys[i];
    drawpoly(p,[0+(3*i),0,3]);
  }
  return setTimeout(animateShape, 100);
};


// Initialization and Basic UI
//===================================================================================================

$( function() { //wait for page to load

  // incorrectly added by decaffeinate I believe, shadows these globals:
  // let PALETTE, specs;

  init(); //init canvas

  const urlParams = parseurl(); //see if recipe is spec'd in URL
  if ("recipe" in urlParams) {
    specs=[urlParams["recipe"]];
    $("#spec").val(specs);
  } else {
    specs=[randomchoice(DEFAULT_RECIPES)];
    $("#spec").val(specs);
    setlink();
  }

  // set initial palette spec
  if ("palette" in urlParams) {
    PALETTE = urlParams["palette"].split(/\s+/g);
    setlink();
  }
  $("#palette").val( PALETTE.reduce((x,y)=> x+" "+y) );

  // construct the polyhedra from spec
  globPolys = _.map(specs, x=> newgeneratePoly(x));
  updateStats();

  // draw it
  drawShape();


  // Event Handlers
  // ----------------------------------------------------

  // when spec changes in input, parse and draw new polyhedra
  $("#spec").change(function(e) {
    specs = $("#spec").val().split(/\s+/g).slice(0, 2); //only allow one recipe for now
    globPolys = _.map(specs, x=> newgeneratePoly(x));
    updateStats();
    //animateShape()
    setlink();
    return drawShape();
  });

  // when palette changes in input, redraw polyhedra
  $("#palette").change(function(e) {
    PALETTE = $(this).val().split(/\s+/g);
    setlink();
    return drawShape();
  });

  // Basic manipulation: rotation and scaling of geometry
  // ----------------------------------------------------

  // mousewheel changes scale of drawing
  $("#poly").mousewheel( function(e,delta, deltaX, deltaY){
    e.preventDefault();
    perspective_scale*=(10+delta)/10;
    return drawShape();
  });

  // Implement standard trackball routines
  // ---------------------------------------
  $("#poly").mousedown( function(e){
    e.preventDefault();
    MOUSEDOWN=true;
    LastMouseX=e.clientX-$(this).offset().left; //relative mouse coords
    LastMouseY=e.clientY-($(this).offset().top-$(window).scrollTop());
    //calculate inverse projection of point to sphere
    const tmpvec=invperspT(LastMouseX,LastMouseY,_2d_x_offset,_2d_y_offset,persp_z_max,persp_z_min,persp_ratio,perspective_scale);
    if ((tmpvec[0]*tmpvec[1]*tmpvec[2]*0) === 0) {  //quick NaN check
      LastSphVec=tmpvec;
    }
    return globlastRotM = clone(globRotM); //copy last transform state
  });
  $("#poly").mouseup( function(e){
    e.preventDefault();
    return MOUSEDOWN=false;
  });
  $("#poly").mouseleave( function(e){
    e.preventDefault();
    return MOUSEDOWN=false;
  });
  $("#poly").mousemove( function(e){
    e.preventDefault();
    if (MOUSEDOWN) {
      const MouseX=e.clientX-$(this).offset().left;
      const MouseY=e.clientY-($(this).offset().top-$(window).scrollTop());
      const SphVec=invperspT(MouseX,MouseY,_2d_x_offset,_2d_y_offset,persp_z_max,persp_z_min,persp_ratio,perspective_scale);

      // quick NaN check
      if (((SphVec[0]*SphVec[1]*SphVec[2]*0) === 0) && ((LastSphVec[0]*LastSphVec[1]*LastSphVec[2]*0) === 0)) {
        globRotM = mm3(getVec2VecRotM(LastSphVec,SphVec),globlastRotM);
      }

      return drawShape();
    }
  });

  // State control via some buttons
  // ---------------------------------------

  $("#strokeonly").click(function(e) {
    PaintMode = "stroke";
    return drawShape();
  });
  
  $("#fillonly").click(function(e) {
    PaintMode = "fill";
    return drawShape();
  });
  
  $("#fillandstroke").click(function(e) {
    PaintMode = "fillstroke";
    return drawShape();
  });

  $("#siderot").click(function(e) {
    globRotM = vec_rotm(PI/2,0,1,0);
    return drawShape();
  });
  
  $("#toprot").click(function(e) {
    globRotM = vec_rotm(PI/2,1,0,0);
    return drawShape();
  });

  $("#frontrot").click(function(e) {
    globRotM = rotm(0,0,0);
    return drawShape();
  });

  // Export Options
  // ---------------------------------------
  $("#pngsavebutton").click(function(e){
    const canvas=$("#poly")[0];
    //this works, but is janky
    //window.location = canvas.toDataURL("image/png")
    const spec = $("#spec").val().split(/\s+/g)[0];
    const filename = `polyhedronisme-${spec.replace(/\([^\)]+\)/g, "")}.png`;
    //blobUtil.canvasToBlob(canvas, 'image/png').then(
    //  blob=>saveAs(blob, filename))
    return canvas.toBlobHD(blob => saveAs(blob, filename));
    //return canvas.toBlob( blob=> saveAs(blob, filename));
  });

  $("#objsavebutton").click(function(e){
    const objtxt = globPolys[0].toOBJ();
    const spec = $("#spec").val().split(/\s+/g)[0];
    const filename = `polyhedronisme-${spec.replace(/\([^\)]+\)/g, "")}.obj`;
    return saveText(objtxt,filename);
  });

  return $("#x3dsavebutton").click(function(e){
    const triangulated = triangulate(globPolys[0],true); //triangulate to preserve face_colors for 3d printing
    const x3dtxt = triangulated.toVRML();
    const spec = $("#spec").val().split(/\s+/g)[0];
    //filename = "polyhedronisme-"+spec+".x3d"
    const filename = `polyhedronisme-${spec.replace(/\([^\)]+\)/g, "")}.wrl`;
    return saveText(x3dtxt,filename);
  });
});
