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

// import math functions to local namespace
const { random, round, floor, sqrt, 
        sin, cos, tan, asin, acos, atan, 
        abs, pow, log,
        PI, LN10
      } = Math;
const log10 = x=> log(x)/LN10;

//returns string w. nsigs digits ignoring magnitude
const sigfigs = function(N, nsigs){
  const normed = pow(10, log10(N)-floor(log10(N)));
  return `${round(normed*(nsigs-1))}`;
};

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
const mult = (c, vec) => 
  [c*vec[0], c*vec[1], c*vec[2]];

// 3d element-wise multiply
const _mult = (vec1, vec2) => 
  [vec1[0]*vec2[0], vec1[1]*vec2[1], vec1[2]*vec2[2]];

// 3d vector addition
const add = (vec1, vec2) => 
  [vec1[0]+vec2[0], vec1[1]+vec2[1], vec1[2]+vec2[2]];

// 3d vector subtraction
const sub = (vec1, vec2) => 
  [vec1[0]-vec2[0], vec1[1]-vec2[1], vec1[2]-vec2[2]];

// 3d dot product
const dot = (vec1, vec2) => 
  (vec1[0]*vec2[0]) + (vec1[1]*vec2[1]) + (vec1[2]*vec2[2]);

// 3d cross product d1 x d2
const cross = (d1, d2) => 
  [(d1[1]*d2[2]) - (d1[2]*d2[1]), 
   (d1[2]*d2[0]) - (d1[0]*d2[2]),  
   (d1[0]*d2[1]) - (d1[1]*d2[0]) ];

// vector norm
const mag = vec => sqrt(dot(vec, vec));

// vector magnitude squared
const mag2 = vec => dot(vec, vec);

// makes vector unit length
const unit = vec => mult(1 / sqrt(mag2(vec)), vec);

// midpoint between vec1, vec2
const midpoint = (vec1, vec2) => mult(1/2.0, add(vec1, vec2));

// parametric segment between vec1, vec2 w. parameter t ranging from 0 to 1
const tween = (vec1, vec2, t) => 
  [((1-t)*vec1[0]) + (t*vec2[0]), 
   ((1-t)*vec1[1]) + (t*vec2[1]), 
   ((1-t)*vec1[2]) + (t*vec2[2])];

// uses above to go one-third of the way along vec1->vec2 line
const oneThird = (vec1, vec2) => tween(vec1, vec2, 1/3.0);

// reflect 3vec in unit sphere, spherical reciprocal
const reciprocal = vec => mult(1.0 / mag2(vec), vec);

// point where line v1...v2 tangent to an origin sphere
const tangentPoint= function(v1, v2) {
  const d = sub(v2, v1);
  return sub(v1, mult(dot(d, v1)/mag2(d), d));
};

// distance of line v1...v2 to origin
const edgeDist = (v1, v2) => sqrt(mag2(tangentPoint(v1, v2)));

// square of distance from point v3 to line segment v1...v2
// http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
// calculates min distance from 
// point v3 to finite line segment between v1 and v2
const linePointDist2 = function(v1, v2, v3) {
  let result;
  const d21 = sub(v2, v1);
  const d13 = sub(v1, v3);
  const d23 = sub(v2, v3);
  const m2 = mag2(d21);
  const t = -dot(d13, d21)/m2;
  if (t <= 0) { 
    // closest to point beyond v1, clip to |v3-v1|^2
    result = mag2(d13);
  } else if (t >= 1) { 
    // closest to point beyond v2, clip to |v3-v2|^2
    result = mag2(d23);
  } else {
    // closest in-between v1, v2
    result = mag2(cross(d21, d13))/m2;
  }
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
  return null; // empty intersection
};

// calculate centroid of array of vertices
const calcCentroid = function(xyzs) {
  // running sum of vertex coords
  let centroidV = [0,0,0];
    for (let v of xyzs) {
      centroidV = add(centroidV, v);
    }
    return mult(1 / xyzs.length, centroidV );
  };

// calculate average normal vector for array of vertices
const normal = function(xyzs) {
  // running sum of normal vectors
  let normalV = [0,0,0]; 
    let [v1,v2] = xyzs.slice(-2);
    for (let v3 of xyzs) {
      normalV = add(normalV, orthogonal(v1, v2, v3));
      [v1,v2] = [v2,v3];
    } // shift over one
    return unit(normalV);
  };

// calculates area planar face by summing over subtriangle areas
//  _Assumes_ Convexity!
const convexarea = function(xyzs) {
    let area = 0.0;
    let [v1,v2] = xyzs.slice(0, 2);
    for (let v3 of xyzs.slice(2)) {
      //area of sub-triangle
      area += mag( cross(sub(v2, v1), sub(v3, v1)) );
      v2 = v3;
    } // shift over one
    return area;
  };

//returns array of ~3sigfig angle
const faceSignature = function(xyzs) {
    let x;
    const cross_array = [];
    let [v1,v2] = xyzs.slice(0, 2);
    for (let v3 of xyzs.slice(2)) {
      //area of sub-triangle
      cross_array.push(mag( cross(sub(v2, v1), sub(v3, v1)) ));
      v2 = v3;
    } // shift over one

    cross_array.sort((a,b)=>a-b); //sort for uniqueness

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
    [x,y,z] = [0,0,1];
  }
  if (length !== 1) {
    [x,y,z] = unit([x,y,z]);
  }

  if ((x === 1) && (y === 0) && (z === 0)) {
      m = [[1,              0,           0],
          [0,    1-(2*sinA2), 2*sinA*cosA],
          [0,   -2*sinA*cosA, 1-(2*sinA2)]];
  } else if ((x === 0) && (y === 1) && (z === 0)) {
      m = [[1-(2*sinA2), 0,  -2*sinA*cosA],
          [          0, 1,             0],
          [2*sinA*cosA, 0,   1-(2*sinA2)]];
  } else if ((x === 0) && (y === 0) && (z === 1)) {
      m = [[   1-(2*sinA2),   2*sinA*cosA, 0],
          [  -2*sinA*cosA,   1-(2*sinA2), 0],
          [             0,             0, 1]];
  } else {
      const x2 = x*x;
      const y2 = y*y;
      const z2 = z*z;
      m = 
        [[1-(2*(y2+z2)*sinA2), 2*((x*y*sinA2)+(z*sinA*cosA)), 2*((x*z*sinA2)-(y*sinA*cosA))],
        [2*((y*x*sinA2)-(z*sinA*cosA)), 1-(2*(z2+x2)*sinA2), 2*((y*z*sinA2)+(x*sinA*cosA))],
        [2*((z*x*sinA2)+(y*sinA*cosA)), 2*((z*y*sinA2)-(x*sinA*cosA)), 1-(2*(x2+y2)*sinA2)]];
    }

  return m;
};

// Perspective Transform
// assumes world's been rotated appropriately such that Z is depth
// scales perspective such that inside depth regions min_real_depth <--> max_real_depth
// perspective lengths vary no more than:   desired_ratio
// with target dimension of roughly length: desired_length
const perspT = function(vec3, max_real_depth, min_real_depth, 
                        desired_ratio, desired_length) {
  const z0 = 
    ((max_real_depth * desired_ratio) - min_real_depth) / (1-desired_ratio);
  const scalefactor = 
    (desired_length * desired_ratio) / (1-desired_ratio);
  // projected [X, Y]
  return [(scalefactor*vec3[0])/(vec3[2]+z0), (scalefactor*vec3[1])/(vec3[2]+z0)];
};

// Inverses perspective transform by projecting plane onto a unit sphere at origin
const invperspT = 
  function(x, y, dx, dy, max_real_depth, min_real_depth, 
           desired_ratio, desired_length) {
  const z0 = 
    ((max_real_depth * desired_ratio) - min_real_depth)/(1-desired_ratio);
  const s = (desired_length * desired_ratio)/(1-desired_ratio);
  const xp = x-dx;
  const yp = y-dy;
  const s2 = s*s;
  const z02 = z0*z0;
  const xp2 = xp*xp;
  const yp2 = yp*yp;

  const xsphere = ((2*s*xp*z0) 
                    + sqrt((4*s2*xp2*z02) 
                    + (4*xp2*(s2+xp2+yp2)*(1-z02))))/(2.0*(s2+xp2+yp2));
  const ysphere = (((s*yp*z0)/(s2+xp2+yp2)) 
                   + ((yp*sqrt((4*s2*z02) 
                   + (4*(s2+xp2+yp2)*(1-z02))))/(2.0*(s2+xp2+yp2))));
  const zsphere = sqrt(1 - (xsphere*xsphere) - (ysphere*ysphere));

  return [xsphere, ysphere, zsphere];
};

// Returns rotation matrix that takes vec1 to vec2
const getVec2VecRotM = function(vec1, vec2){
  const axis  = cross(vec1, vec2);
  const angle = acos(dot(vec1, vec2));
  return vec_rotm(-1*angle, axis[0], axis[1], axis[2]);
};
