// Polyhédronisme - Flat 3x3 Matrix Math (Row-Major)
// ===================================================================================================
// Layout:
// [ 0, 1, 2,
//   3, 4, 5,
//   6, 7, 8 ]

export const create = () => [
  1, 0, 0,
  0, 1, 0,
  0, 0, 1
];

export const identity = (out = create()) => {
  out[0] = 1; out[1] = 0; out[2] = 0;
  out[3] = 0; out[4] = 1; out[5] = 0;
  out[6] = 0; out[7] = 0; out[8] = 1;
  return out;
};

export const clone = (a) => [
  a[0], a[1], a[2],
  a[3], a[4], a[5],
  a[6], a[7], a[8]
];

export const copy = (out, a) => {
  out[0] = a[0]; out[1] = a[1]; out[2] = a[2];
  out[3] = a[3]; out[4] = a[4]; out[5] = a[5];
  out[6] = a[6]; out[7] = a[7]; out[8] = a[8];
  return out;
};

// Convert legacy nested matrix [[a,b,c],[d,e,f],[g,h,i]] to flat 9 array
export const fromNested = (nested, out = create()) => {
  out[0] = nested[0][0]; out[1] = nested[0][1]; out[2] = nested[0][2];
  out[3] = nested[1][0]; out[4] = nested[1][1]; out[5] = nested[1][2];
  out[6] = nested[2][0]; out[7] = nested[2][1]; out[8] = nested[2][2];
  return out;
};

// Convert flat 9 array to legacy nested matrix [[a,b,c],[d,e,f],[g,h,i]]
export const toNested = (flat) => [
  [flat[0], flat[1], flat[2]],
  [flat[3], flat[4], flat[5]],
  [flat[6], flat[7], flat[8]]
];

// Matrix multiplication: out = A * B
export const multiply = (A, B, out = create()) => {
  // Support both flat arrays and legacy nested arrays transparently
  const a00 = A[0] !== undefined && typeof A[0] === 'number' ? A[0] : A[0][0];
  const a01 = A[1] !== undefined && typeof A[1] === 'number' ? A[1] : A[0][1];
  const a02 = A[2] !== undefined && typeof A[2] === 'number' ? A[2] : A[0][2];
  const a10 = A[3] !== undefined && typeof A[3] === 'number' ? A[3] : A[1][0];
  const a11 = A[4] !== undefined && typeof A[4] === 'number' ? A[4] : A[1][1];
  const a12 = A[5] !== undefined && typeof A[5] === 'number' ? A[5] : A[1][2];
  const a20 = A[6] !== undefined && typeof A[6] === 'number' ? A[6] : A[2][0];
  const a21 = A[7] !== undefined && typeof A[7] === 'number' ? A[7] : A[2][1];
  const a22 = A[8] !== undefined && typeof A[8] === 'number' ? A[8] : A[2][2];

  const b00 = B[0] !== undefined && typeof B[0] === 'number' ? B[0] : B[0][0];
  const b01 = B[1] !== undefined && typeof B[1] === 'number' ? B[1] : B[0][1];
  const b02 = B[2] !== undefined && typeof B[2] === 'number' ? B[2] : B[0][2];
  const b10 = B[3] !== undefined && typeof B[3] === 'number' ? B[3] : B[1][0];
  const b11 = B[4] !== undefined && typeof B[4] === 'number' ? B[4] : B[1][1];
  const b12 = B[5] !== undefined && typeof B[5] === 'number' ? B[5] : B[1][2];
  const b20 = B[6] !== undefined && typeof B[6] === 'number' ? B[6] : B[2][0];
  const b21 = B[7] !== undefined && typeof B[7] === 'number' ? B[7] : B[2][1];
  const b22 = B[8] !== undefined && typeof B[8] === 'number' ? B[8] : B[2][2];

  out[0] = (a00 * b00) + (a01 * b10) + (a02 * b20);
  out[1] = (a00 * b01) + (a01 * b11) + (a02 * b21);
  out[2] = (a00 * b02) + (a01 * b12) + (a02 * b22);

  out[3] = (a10 * b00) + (a11 * b10) + (a12 * b20);
  out[4] = (a10 * b01) + (a11 * b11) + (a12 * b21);
  out[5] = (a10 * b02) + (a11 * b12) + (a12 * b22);

  out[6] = (a20 * b00) + (a21 * b10) + (a22 * b20);
  out[7] = (a20 * b01) + (a21 * b11) + (a22 * b21);
  out[8] = (a20 * b02) + (a21 * b12) + (a22 * b22);
  return out;
};

// Matrix-vector multiply: out = M * v
export const multiplyVec3 = (M, v, out = [0, 0, 0]) => {
  const m00 = M[0] !== undefined && typeof M[0] === 'number' ? M[0] : M[0][0];
  const m01 = M[1] !== undefined && typeof M[1] === 'number' ? M[1] : M[0][1];
  const m02 = M[2] !== undefined && typeof M[2] === 'number' ? M[2] : M[0][2];
  const m10 = M[3] !== undefined && typeof M[3] === 'number' ? M[3] : M[1][0];
  const m11 = M[4] !== undefined && typeof M[4] === 'number' ? M[4] : M[1][1];
  const m12 = M[5] !== undefined && typeof M[5] === 'number' ? M[5] : M[1][2];
  const m20 = M[6] !== undefined && typeof M[6] === 'number' ? M[6] : M[2][0];
  const m21 = M[7] !== undefined && typeof M[7] === 'number' ? M[7] : M[2][1];
  const m22 = M[8] !== undefined && typeof M[8] === 'number' ? M[8] : M[2][2];

  const x = v[0];
  const y = v[1];
  const z = v[2];

  out[0] = (m00 * x) + (m01 * y) + (m02 * z);
  out[1] = (m10 * x) + (m11 * y) + (m12 * z);
  out[2] = (m20 * x) + (m21 * y) + (m22 * z);
  return out;
};

// Transpose: out = M^T
export const transpose = (M, out = create()) => {
  const m00 = M[0]; const m01 = M[1]; const m02 = M[2];
  const m10 = M[3]; const m11 = M[4]; const m12 = M[5];
  const m20 = M[6]; const m21 = M[7]; const m22 = M[8];

  out[0] = m00; out[1] = m10; out[2] = m20;
  out[3] = m01; out[4] = m11; out[5] = m21;
  out[6] = m02; out[7] = m12; out[8] = m22;
  return out;
};

// Fast 3x3 Matrix Inversion via cofactor expansion. Returns null if singular.
export const invert = (M, out = create()) => {
  const m00 = M[0] !== undefined && typeof M[0] === 'number' ? M[0] : M[0][0];
  const m01 = M[1] !== undefined && typeof M[1] === 'number' ? M[1] : M[0][1];
  const m02 = M[2] !== undefined && typeof M[2] === 'number' ? M[2] : M[0][2];
  const m10 = M[3] !== undefined && typeof M[3] === 'number' ? M[3] : M[1][0];
  const m11 = M[4] !== undefined && typeof M[4] === 'number' ? M[4] : M[1][1];
  const m12 = M[5] !== undefined && typeof M[5] === 'number' ? M[5] : M[1][2];
  const m20 = M[6] !== undefined && typeof M[6] === 'number' ? M[6] : M[2][0];
  const m21 = M[7] !== undefined && typeof M[7] === 'number' ? M[7] : M[2][1];
  const m22 = M[8] !== undefined && typeof M[8] === 'number' ? M[8] : M[2][2];

  const c00 = (m11 * m22) - (m12 * m21);
  const c01 = (m12 * m20) - (m10 * m22);
  const c02 = (m10 * m21) - (m11 * m20);

  const det = (m00 * c00) + (m01 * c01) + (m02 * c02);
  if (Math.abs(det) < 1e-12) return null;

  const invDet = 1.0 / det;

  const c10 = (m02 * m21) - (m01 * m22);
  const c11 = (m00 * m22) - (m02 * m20);
  const c12 = (m01 * m20) - (m00 * m21);

  const c20 = (m01 * m12) - (m02 * m11);
  const c21 = (m02 * m10) - (m00 * m12);
  const c22 = (m00 * m11) - (m01 * m10);

  out[0] = c00 * invDet; out[1] = c10 * invDet; out[2] = c20 * invDet;
  out[3] = c01 * invDet; out[4] = c11 * invDet; out[5] = c21 * invDet;
  out[6] = c02 * invDet; out[7] = c12 * invDet; out[8] = c22 * invDet;
  return out;
};

// Rotation matrix from axis and angle (radians)
export const fromAxisAngle = (axis, angle, out = create()) => {
  let [x, y, z] = axis;
  const len = Math.hypot(x, y, z);
  if (len < 1e-12) {
    return identity(out);
  }
  x /= len;
  y /= len;
  z /= len;

  const s = Math.sin(angle);
  const c = Math.cos(angle);
  const t = 1 - c;

  out[0] = (x * x * t) + c;
  out[1] = (x * y * t) - (z * s);
  out[2] = (x * z * t) + (y * s);

  out[3] = (y * x * t) + (z * s);
  out[4] = (y * y * t) + c;
  out[5] = (y * z * t) - (x * s);

  out[6] = (z * x * t) - (y * s);
  out[7] = (z * y * t) + (x * s);
  out[8] = (z * z * t) + c;
  return out;
};

// Rotation matrix from unit quaternion [x, y, z, w]
export const fromQuat = (q, out = create()) => {
  const x = q[0];
  const y = q[1];
  const z = q[2];
  const w = q[3];

  const x2 = x + x;
  const y2 = y + y;
  const z2 = z + z;

  const xx = x * x2;
  const xy = x * y2;
  const xz = x * z2;
  const yy = y * y2;
  const yz = y * z2;
  const zz = z * z2;
  const wx = w * x2;
  const wy = w * y2;
  const wz = w * z2;

  out[0] = 1 - (yy + zz);
  out[1] = xy - wz;
  out[2] = xz + wy;

  out[3] = xy + wz;
  out[4] = 1 - (xx + zz);
  out[5] = yz - wx;

  out[6] = xz - wy;
  out[7] = yz + wx;
  out[8] = 1 - (xx + yy);
  return out;
};

// Euler rotation matrix matching legacy rotm(phi, theta, psi)
export const fromEuler = (phi, theta, psi, out = create()) => {
  const cPhi = Math.cos(phi);
  const sPhi = Math.sin(phi);
  const cThe = Math.cos(theta);
  const sThe = Math.sin(theta);
  const cPsi = Math.cos(psi);
  const sPsi = Math.sin(psi);

  // XY rotation
  const xy = [
    cPhi, -sPhi, 0,
    sPhi, cPhi, 0,
    0, 0, 1
  ];
  // YZ rotation
  const yz = [
    cThe, 0, -sThe,
    0, 1, 0,
    sThe, 0, cThe
  ];
  // XZ rotation
  const xz = [
    1, 0, 0,
    0, cPsi, -sPsi,
    0, sPsi, cPsi
  ];

  const tmp = create();
  multiply(yz, xy, tmp);
  multiply(xz, tmp, out);
  return out;
};
