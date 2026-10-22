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

/* eslint-disable camelcase */
/* eslint-disable standard/array-bracket-even-spacing */
/* eslint-disable no-multi-spaces */

import * as vec2 from './math/vec2.js';
import * as vec3 from './math/vec3.js';
import * as mat3 from './math/mat3.js';
import * as quat from './math/quat.js';

export { vec2, vec3, mat3, quat };

// Math Functions & Constants
// ===================================================================================================
const { random, round, floor, sqrt,
  sin, cos, tan, asin, acos, atan, PI,
  abs, pow, log, LN10 } = Math;

export { random, round, floor, sqrt,
  sin, cos, tan, asin, acos, atan, PI,
  abs, pow, log, LN10 };

export const log10 = Math.log10;
export const sign = Math.sign;

export function __range__ (left, right, inclusive) {
  const range = [];
  const ascending = left < right;
  const end = !inclusive ? right : ascending ? right + 1 : right - 1;
  for (let i = left; ascending ? i < end : i > end; ascending ? i++ : i--) {
    range.push(i);
  }
  return range;
}

// returns string w. nsigs digits ignoring magnitude
export const sigfigs = function (N, nsigs) {
  if (N === 0) return '0';
  const mantissa = N / pow(10, floor(log10(abs(N))));
  const truncated_mantissa = round(mantissa * pow(10, (nsigs - 1)));
  return `${truncated_mantissa}`;
};

// general recursive deep-copy function
export const clone = function (obj) {
  if ((obj == null) || (typeof obj !== 'object')) {
    return obj;
  }
  if (Array.isArray(obj)) {
    return obj.map(clone);
  }
  const newInstance = new obj.constructor();
  for (let key of Object.keys(obj)) {
    newInstance[key] = clone(obj[key]);
  }
  return newInstance;
};

// often useful
export const randomchoice = function (array) {
  const n = floor(random() * array.length);
  return array[n];
};

// 3D Vector Math Primitives
// ===================================================================================================
export const mult = (c, vec) => vec3.scale(vec, c);
export const add = vec3.add;
export const sub = vec3.sub;
export const _mult = vec3.multiply;
export const dot = vec3.dot;
export const cross = vec3.cross;
export const mag = vec3.mag;
export const mag2 = vec3.magSq;
export const unit = vec3.normalize;
export const midpoint = vec3.midpoint;
export const tween = vec3.lerp;
export const oneThird = vec3.oneThird;
export const reciprocal = vec3.reciprocal;
export const tangentPoint = vec3.tangentPoint;
export const edgeDist = vec3.edgeDist;
export const linePointDist2 = vec3.linePointDistSq;
export const orthogonal = vec3.orthogonal;
export const calcCentroid = vec3.calcCentroid;
export const normal = vec3.normal;
export const planararea = vec3.planarArea;

// 2D Vector Math Primitives
// ===================================================================================================
export const mult2d = (c, vec) => vec2.scale(vec, c);
export const add2d = vec2.add;
export const sub2d = vec2.sub;
export const dot2d = vec2.dot;
export const cross2d = vec2.cross;
export const mag2d = vec2.mag;
export const mag22d = vec2.magSq;
export const unit2d = vec2.normalize;
export const midpoint2d = vec2.midpoint;
export const tween2d = vec2.lerp;
export const calcCentroid2d = vec2.calcCentroid;
export const planararea2d = vec2.planarArea;

// find first element common to 3 sets by brute force search
export const intersect = function (set1, set2, set3) {
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

// congruence signature for assigning same colors to congruent faces
export const faceSignature = function (vertices, sensitivity) {
  let x;
  const cross_array = [];
  let [v1, v2] = vertices.slice(-2);
  for (let v3 of vertices) {
    // accumulate inner angles
    cross_array.push(mag(cross(sub(v1, v2), sub(v3, v2))));
    [v1, v2] = [v2, v3];
  }
  // sort angles to create unique sequence
  cross_array.sort((a, b) => a - b);

  // render sorted angles as quantized digit strings
  // this is the congruence signature
  let sig = '';
  for (x of cross_array) { sig += sigfigs(x, sensitivity); }
  // hack to make reflected faces share the same signature
  for (x of cross_array.slice().reverse()) { sig += sigfigs(x, sensitivity); }

  return sig;
};

// projects 3d polyhedral face to 2d polygon for triangulation and face display
export const project2dface = function (verts) {
  if (verts.length < 3) return verts.map(() => [0, 0]);
  const v0 = verts[0];
  const n = normal(verts);

  // Form an orthonormal in-plane basis (u, v)
  let u = sub(verts[1], v0);
  if (mag2(u) < 1e-12) {
    u = sub(verts[2], v0);
  }
  u = unit(u);
  const v = cross(n, u);

  return verts.map((pt) => {
    const d = sub(pt, v0);
    return [dot(u, d), dot(v, d)];
  });
};

// copies array of arrays by value (deep copy)
export const copyVecArray = function (vecArray) {
  const newVecArray = new Array(vecArray.length);
  for (let i = 0, end = vecArray.length; i < end; i++) {
    newVecArray[i] = vecArray[i].slice(0);
  }
  return newVecArray;
};

// 3D Matrix & Linear Algebra
// ===================================================================================================
export const mv3 = (mat, vec) => mat3.multiplyVec3(mat, vec);
export const mm3 = (A, B) => mat3.multiply(A, B);
export const eye3 = [1, 0, 0, 0, 1, 0, 0, 0, 1];

export const rotm = (phi, theta, psi) => mat3.fromEuler(phi, theta, psi);
export const vec_rotm = (angle, x, y, z) => mat3.fromAxisAngle([x, y, z], angle);

// Returns rotation matrix taking vec1 to vec2
export const getVec2VecRotM = function (vec1, vec2) {
  const q = quat.fromUnitVectors(unit(vec1), unit(vec2));
  return mat3.fromQuat(q);
};
if (typeof window !== 'undefined') {
  window.getVec2VecRotM = getVec2VecRotM;
}

// Perspective Transform
// ===================================================================================================
export const persp_z_max = 5;
export const persp_z_min = 0;
export const persp_ratio = 0.8;
export const perspective_scale = 1200;

// assumes world's been rotated appropriately such that Z is depth
export const perspT = function (vec3, max_real_depth = persp_z_max, min_real_depth = persp_z_min, desired_ratio = persp_ratio, desired_length = perspective_scale) {
  const z0 = ((max_real_depth * desired_ratio) - min_real_depth) / (1 - desired_ratio);
  const scalefactor = (desired_length * desired_ratio) / (1 - desired_ratio);
  // projected [X, Y]
  return [(scalefactor * vec3[0]) / (vec3[2] + z0), (scalefactor * vec3[1]) / (vec3[2] + z0)];
};

// z-depth to use where [screen_x, screen_y, magic_z_depth] corresponds to correct world-ray from origin
export const magic_depth = function (max_real_depth, min_real_depth, desired_ratio, desired_length) {
  const z0 = ((max_real_depth * desired_ratio) - min_real_depth) / (1 - desired_ratio);
  const scalefactor = (desired_length * desired_ratio) / (1 - desired_ratio);
  return scalefactor - z0;
};

// Inverses perspective transform by projecting plane onto a unit sphere at origin
export const invperspT = function (x, y, dx, dy, max_real_depth, min_real_depth, desired_ratio, desired_length) {
  const z0 = ((max_real_depth * desired_ratio) - min_real_depth) / (1 - desired_ratio);
  const s = (desired_length * desired_ratio) / (1 - desired_ratio);
  const xp = x - dx;
  const yp = y - dy;
  const s2 = s * s;
  const z02 = z0 * z0;
  const xp2 = xp * xp;
  const yp2 = yp * yp;
  const denom = 2.0 * (s2 + xp2 + yp2);

  const xsphere = ((2 * s * xp * z0) +
                  sqrt((4 * s2 * xp2 * z02) +
                  (4 * xp2 * (s2 + xp2 + yp2) * (1 - z02)))) / denom;
  const ysphere = (((s * yp * z0) / (s2 + xp2 + yp2)) +
                 ((yp * sqrt((4 * s2 * z02) +
                 (4 * (s2 + xp2 + yp2) * (1 - z02)))) / denom));
  const zsphere = sqrt(Math.max(0, 1 - (xsphere * xsphere) - (ysphere * ysphere)));

  return [xsphere, ysphere, zsphere];
};

// 2D point inclusion in 2D poly by winding number
// see http://geomalgorithms.com/a03-_inclusion.html
export const pointInPolygon = function (point, vs) {
  const [ x, y ] = point;
  let wn = 0;
  const isLeft = function (P0, P1, P2) {
    return ((P1[0] - P0[0]) * (P2[1] - P0[1])) - ((P2[0] - P0[0]) * (P1[1] - P0[1]));
  };
  for (var i = 0, j = vs.length - 1; i < vs.length; j = i++) {
    const [ xi, yi ] = vs[i];
    const [ xj, yj ] = vs[j];
    if (yj <= y) {
      if (yi > y) {
        if (isLeft([xj, yj], [xi, yi], [x, y]) > 0) {
          wn++;
        }
      }
    } else {
      if (yi <= y) {
        if (isLeft([xj, yj], [xi, yi], [x, y]) < 0) {
          wn--;
        }
      }
    }
  }
  return wn !== 0;
};
