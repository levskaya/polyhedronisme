// Polyhédronisme - 3D Vector Math
// ===================================================================================================

export const create = (x = 0, y = 0, z = 0) => [x, y, z];

export const set = (out, x, y, z) => {
  out[0] = x;
  out[1] = y;
  out[2] = z;
  return out;
};

export const copy = (out, a) => {
  out[0] = a[0];
  out[1] = a[1];
  out[2] = a[2];
  return out;
};

export const clone = (a) => [a[0], a[1], a[2]];

export const add = (a, b, out = [0, 0, 0]) => {
  out[0] = a[0] + b[0];
  out[1] = a[1] + b[1];
  out[2] = a[2] + b[2];
  return out;
};

export const sub = (a, b, out = [0, 0, 0]) => {
  out[0] = a[0] - b[0];
  out[1] = a[1] - b[1];
  out[2] = a[2] - b[2];
  return out;
};

export const scale = (a, s, out = [0, 0, 0]) => {
  out[0] = a[0] * s;
  out[1] = a[1] * s;
  out[2] = a[2] * s;
  return out;
};

// Element-wise multiply
export const multiply = (a, b, out = [0, 0, 0]) => {
  out[0] = a[0] * b[0];
  out[1] = a[1] * b[1];
  out[2] = a[2] * b[2];
  return out;
};

export const dot = (a, b) => (a[0] * b[0]) + (a[1] * b[1]) + (a[2] * b[2]);

export const cross = (a, b, out = [0, 0, 0]) => {
  const ax = a[0];
  const ay = a[1];
  const az = a[2];
  const bx = b[0];
  const by = b[1];
  const bz = b[2];
  out[0] = (ay * bz) - (az * by);
  out[1] = (az * bx) - (ax * bz);
  out[2] = (ax * by) - (ay * bx);
  return out;
};

export const magSq = (a) => (a[0] * a[0]) + (a[1] * a[1]) + (a[2] * a[2]);

export const mag = (a) => Math.hypot(a[0], a[1], a[2]);

export const normalize = (a, out = [0, 0, 0]) => {
  const m = Math.hypot(a[0], a[1], a[2]);
  if (m > 1e-12) {
    const inv = 1 / m;
    out[0] = a[0] * inv;
    out[1] = a[1] * inv;
    out[2] = a[2] * inv;
  } else {
    out[0] = 0;
    out[1] = 0;
    out[2] = 0;
  }
  return out;
};

export const midpoint = (a, b, out = [0, 0, 0]) => {
  out[0] = (a[0] + b[0]) * 0.5;
  out[1] = (a[1] + b[1]) * 0.5;
  out[2] = (a[2] + b[2]) * 0.5;
  return out;
};

export const lerp = (a, b, t, out = [0, 0, 0]) => {
  out[0] = ((1 - t) * a[0]) + (t * b[0]);
  out[1] = ((1 - t) * a[1]) + (t * b[1]);
  out[2] = ((1 - t) * a[2]) + (t * b[2]);
  return out;
};

export const oneThird = (a, b, out = [0, 0, 0]) => lerp(a, b, 1 / 3, out);

// Spherical reciprocal: reflect in unit sphere
export const reciprocal = (a, out = [0, 0, 0]) => {
  const m2 = magSq(a);
  const s = m2 > 1e-12 ? 1 / m2 : 0;
  return scale(a, s, out);
};

// Point on line segment (v1...v2) closest to the origin
export const tangentPoint = (v1, v2, out = [0, 0, 0]) => {
  const dx = v2[0] - v1[0];
  const dy = v2[1] - v1[1];
  const dz = v2[2] - v1[2];
  const m2 = (dx * dx) + (dy * dy) + (dz * dz);
  if (m2 < 1e-14) {
    out[0] = v1[0];
    out[1] = v1[1];
    out[2] = v1[2];
    return out;
  }
  const dDotV1 = (dx * v1[0]) + (dy * v1[1]) + (dz * v1[2]);
  const factor = dDotV1 / m2;
  out[0] = v1[0] - (factor * dx);
  out[1] = v1[1] - (factor * dy);
  out[2] = v1[2] - (factor * dz);
  return out;
};

// Distance of line v1...v2 to origin
export const edgeDist = (v1, v2) => {
  const tp = tangentPoint(v1, v2);
  return mag(tp);
};

// Minimum squared distance from point v3 to finite line segment between v1 and v2
export const linePointDistSq = (v1, v2, v3) => {
  const d21x = v2[0] - v1[0];
  const d21y = v2[1] - v1[1];
  const d21z = v2[2] - v1[2];

  const d13x = v1[0] - v3[0];
  const d13y = v1[1] - v3[1];
  const d13z = v1[2] - v3[2];

  const m2 = (d21x * d21x) + (d21y * d21y) + (d21z * d21z);
  if (m2 < 1e-14) {
    return (d13x * d13x) + (d13y * d13y) + (d13z * d13z);
  }

  const t = -((d13x * d21x) + (d13y * d21y) + (d13z * d21z)) / m2;
  if (t <= 0) {
    return (d13x * d13x) + (d13y * d13y) + (d13z * d13z);
  } else if (t >= 1) {
    const d23x = v2[0] - v3[0];
    const d23y = v2[1] - v3[1];
    const d23z = v2[2] - v3[2];
    return (d23x * d23x) + (d23y * d23y) + (d23z * d23z);
  } else {
    // cross(d21, d13)
    const cx = (d21y * d13z) - (d21z * d13y);
    const cy = (d21z * d13x) - (d21x * d13z);
    const cz = (d21x * d13y) - (d21y * d13x);
    return ((cx * cx) + (cy * cy) + (cz * cz)) / m2;
  }
};

// Vector orthogonal to plane of 3 points: (v2 - v1) x (v3 - v2)
export const orthogonal = (v1, v2, v3, out = [0, 0, 0]) => {
  const d1x = v2[0] - v1[0];
  const d1y = v2[1] - v1[1];
  const d1z = v2[2] - v1[2];

  const d2x = v3[0] - v2[0];
  const d2y = v3[1] - v2[1];
  const d2z = v3[2] - v2[2];

  out[0] = (d1y * d2z) - (d1z * d2y);
  out[1] = (d1z * d2x) - (d1x * d2z);
  out[2] = (d1x * d2y) - (d1y * d2x);
  return out;
};

// Centroid of array of 3D vertices
export const calcCentroid = (vertices, out = [0, 0, 0]) => {
  if (vertices.length === 0) {
    out[0] = 0;
    out[1] = 0;
    out[2] = 0;
    return out;
  }
  let x = 0;
  let y = 0;
  let z = 0;
  for (let i = 0; i < vertices.length; i++) {
    x += vertices[i][0];
    y += vertices[i][1];
    z += vertices[i][2];
  }
  const inv = 1 / vertices.length;
  out[0] = x * inv;
  out[1] = y * inv;
  out[2] = z * inv;
  return out;
};

// Newell's method / sum of edge cross products for average polygon normal
export const normal = (vertices, out = [0, 0, 0]) => {
  if (vertices.length < 3) {
    out[0] = 0;
    out[1] = 0;
    out[2] = 1;
    return out;
  }
  let nx = 0;
  let ny = 0;
  let nz = 0;
  let [v1, v2] = vertices.slice(-2);
  for (let i = 0; i < vertices.length; i++) {
    const v3 = vertices[i];
    const d1x = v2[0] - v1[0];
    const d1y = v2[1] - v1[1];
    const d1z = v2[2] - v1[2];

    const d2x = v3[0] - v2[0];
    const d2y = v3[1] - v2[1];
    const d2z = v3[2] - v2[2];

    nx += (d1y * d2z) - (d1z * d2y);
    ny += (d1z * d2x) - (d1x * d2z);
    nz += (d1x * d2y) - (d1y * d2x);

    v1 = v2;
    v2 = v3;
  }
  const m = Math.hypot(nx, ny, nz);
  if (m > 1e-12) {
    const inv = 1 / m;
    out[0] = nx * inv;
    out[1] = ny * inv;
    out[2] = nz * inv;
  } else {
    out[0] = 0;
    out[1] = 0;
    out[2] = 1;
  }
  return out;
};

// Planar area of polygon in 3D
export const planarArea = (vertices) => {
  if (vertices.length < 3) return 0;
  let vx = 0;
  let vy = 0;
  let vz = 0;
  let [v1, v2] = vertices.slice(-2);
  for (let i = 0; i < vertices.length; i++) {
    const v3 = vertices[i];
    vx += (v1[1] * v2[2]) - (v1[2] * v2[1]);
    vy += (v1[2] * v2[0]) - (v1[0] * v2[2]);
    vz += (v1[0] * v2[1]) - (v1[1] * v2[0]);
    v1 = v2;
    v2 = v3;
  }
  const n = normal(vertices);
  return Math.abs(((n[0] * vx) + (n[1] * vy) + (n[2] * vz)) * 0.5);
};
