// Polyhédronisme - 2D Vector Math
// ===================================================================================================

export const create = (x = 0, y = 0) => [x, y];

export const set = (out, x, y) => {
  out[0] = x;
  out[1] = y;
  return out;
};

export const copy = (out, a) => {
  out[0] = a[0];
  out[1] = a[1];
  return out;
};

export const clone = (a) => [a[0], a[1]];

export const add = (a, b, out = [0, 0]) => {
  out[0] = a[0] + b[0];
  out[1] = a[1] + b[1];
  return out;
};

export const sub = (a, b, out = [0, 0]) => {
  out[0] = a[0] - b[0];
  out[1] = a[1] - b[1];
  return out;
};

export const scale = (a, s, out = [0, 0]) => {
  out[0] = a[0] * s;
  out[1] = a[1] * s;
  return out;
};

export const multiply = (a, b, out = [0, 0]) => {
  out[0] = a[0] * b[0];
  out[1] = a[1] * b[1];
  return out;
};

export const dot = (a, b) => (a[0] * b[0]) + (a[1] * b[1]);

// 2D scalar cross product (determinant)
export const cross = (a, b) => (a[0] * b[1]) - (a[1] * b[0]);

export const magSq = (a) => (a[0] * a[0]) + (a[1] * a[1]);

export const mag = (a) => Math.hypot(a[0], a[1]);

export const normalize = (a, out = [0, 0]) => {
  const m = Math.hypot(a[0], a[1]);
  if (m > 1e-12) {
    const inv = 1 / m;
    out[0] = a[0] * inv;
    out[1] = a[1] * inv;
  } else {
    out[0] = 0;
    out[1] = 0;
  }
  return out;
};

export const midpoint = (a, b, out = [0, 0]) => {
  out[0] = (a[0] + b[0]) * 0.5;
  out[1] = (a[1] + b[1]) * 0.5;
  return out;
};

export const lerp = (a, b, t, out = [0, 0]) => {
  out[0] = ((1 - t) * a[0]) + (t * b[0]);
  out[1] = ((1 - t) * a[1]) + (t * b[1]);
  return out;
};

export const calcCentroid = (vertices, out = [0, 0]) => {
  if (vertices.length === 0) {
    out[0] = 0;
    out[1] = 0;
    return out;
  }
  let x = 0;
  let y = 0;
  for (let i = 0; i < vertices.length; i++) {
    x += vertices[i][0];
    y += vertices[i][1];
  }
  const inv = 1 / vertices.length;
  out[0] = x * inv;
  out[1] = y * inv;
  return out;
};

export const planarArea = (vertices) => {
  if (vertices.length < 3) return 0;
  let vsum = 0;
  let [v1, v2] = vertices.slice(-2);
  for (let i = 0; i < vertices.length; i++) {
    const v3 = vertices[i];
    vsum += cross(v1, v2);
    v1 = v2;
    v2 = v3;
  }
  return Math.abs(vsum * 0.5);
};
