// Polyhédronisme - Quaternion Math [x, y, z, w]
// ===================================================================================================

export const create = () => [0, 0, 0, 1];

export const identity = (out = create()) => {
  out[0] = 0;
  out[1] = 0;
  out[2] = 0;
  out[3] = 1;
  return out;
};

export const copy = (out, a) => {
  out[0] = a[0];
  out[1] = a[1];
  out[2] = a[2];
  out[3] = a[3];
  return out;
};

export const clone = (a) => [a[0], a[1], a[2], a[3]];

export const set = (out, x, y, z, w) => {
  out[0] = x;
  out[1] = y;
  out[2] = z;
  out[3] = w;
  return out;
};

export const setAxisAngle = (axis, rad, out = create()) => {
  let [x, y, z] = axis;
  const len = Math.hypot(x, y, z);
  if (len < 1e-12) {
    return identity(out);
  }
  const half = rad * 0.5;
  const s = Math.sin(half) / len;
  out[0] = x * s;
  out[1] = y * s;
  out[2] = z * s;
  out[3] = Math.cos(half);
  return out;
};

// Quaternion multiplication: out = a * b
export const multiply = (a, b, out = create()) => {
  const ax = a[0]; const ay = a[1]; const az = a[2]; const aw = a[3];
  const bx = b[0]; const by = b[1]; const bz = b[2]; const bw = b[3];

  out[0] = (aw * bx) + (ax * bw) + (ay * bz) - (az * by);
  out[1] = (aw * by) - (ax * bz) + (ay * bw) + (az * bx);
  out[2] = (aw * bz) + (ax * by) - (ay * bx) + (az * bw);
  out[3] = (aw * bw) - (ax * bx) - (ay * by) - (az * bz);
  return out;
};

export const normalize = (a, out = create()) => {
  const len = Math.hypot(a[0], a[1], a[2], a[3]);
  if (len > 1e-12) {
    const inv = 1 / len;
    out[0] = a[0] * inv;
    out[1] = a[1] * inv;
    out[2] = a[2] * inv;
    out[3] = a[3] * inv;
  } else {
    identity(out);
  }
  return out;
};

// Computes the rotation quaternion that takes unit vector u to unit vector v
export const fromUnitVectors = (u, v, out = create()) => {
  const ux = u[0]; const uy = u[1]; const uz = u[2];
  const vx = v[0]; const vy = v[1]; const vz = v[2];

  const dot = (ux * vx) + (uy * vy) + (uz * vz);
  const r = dot + 1;

  if (r < 1e-7) {
    // Vectors are collinear and opposite: 180 degree rotation around an orthogonal axis
    let orthX, orthY, orthZ;
    if (Math.abs(ux) > Math.abs(uz)) {
      orthX = -uy; orthY = ux; orthZ = 0;
    } else {
      orthX = 0; orthY = -uz; orthZ = uy;
    }
    return setAxisAngle([orthX, orthY, orthZ], Math.PI, out);
  }

  // Cross product
  const cx = (uy * vz) - (uz * vy);
  const cy = (uz * vx) - (ux * vz);
  const cz = (ux * vy) - (uy * vx);

  out[0] = cx;
  out[1] = cy;
  out[2] = cz;
  out[3] = r;
  return normalize(out, out);
};

// Fast vector rotation by unit quaternion: v' = v + 2*r*(u x v) + 2*(u x (u x v))
export const rotateVec3 = (q, v, out = [0, 0, 0]) => {
  const qx = q[0]; const qy = q[1]; const qz = q[2]; const qw = q[3];
  const vx = v[0]; const vy = v[1]; const vz = v[2];

  // t = 2 * cross(q.xyz, v)
  const tx = 2 * ((qy * vz) - (qz * vy));
  const ty = 2 * ((qz * vx) - (qx * vz));
  const tz = 2 * ((qx * vy) - (qy * vx));

  // cross(q.xyz, t)
  const cx = (qy * tz) - (qz * ty);
  const cy = (qz * tx) - (qx * tz);
  const cz = (qx * ty) - (qy * tx);

  out[0] = vx + (qw * tx) + cx;
  out[1] = vy + (qw * ty) + cy;
  out[2] = vz + (qw * tz) + cz;
  return out;
};

// Spherical linear interpolation between quaternions a and b
export const slerp = (a, b, t, out = create()) => {
  let bx = b[0]; let by = b[1]; let bz = b[2]; let bw = b[3];
  let cosHalfTheta = (a[0] * bx) + (a[1] * by) + (a[2] * bz) + (a[3] * bw);

  if (cosHalfTheta < 0) {
    cosHalfTheta = -cosHalfTheta;
    bx = -bx; by = -by; bz = -bz; bw = -bw;
  }

  if (cosHalfTheta >= 1.0) {
    out[0] = a[0]; out[1] = a[1]; out[2] = a[2]; out[3] = a[3];
    return out;
  }

  const sqrSinHalfTheta = 1.0 - (cosHalfTheta * cosHalfTheta);
  if (sqrSinHalfTheta <= 1e-12) {
    // Linear interpolation when very close
    out[0] = a[0] + (t * (bx - a[0]));
    out[1] = a[1] + (t * (by - a[1]));
    out[2] = a[2] + (t * (bz - a[2]));
    out[3] = a[3] + (t * (bw - a[3]));
    return normalize(out, out);
  }

  const sinHalfTheta = Math.sqrt(sqrSinHalfTheta);
  const halfTheta = Math.atan2(sinHalfTheta, cosHalfTheta);
  const ratioA = Math.sin((1 - t) * halfTheta) / sinHalfTheta;
  const ratioB = Math.sin(t * halfTheta) / sinHalfTheta;

  out[0] = (a[0] * ratioA) + (bx * ratioB);
  out[1] = (a[1] * ratioA) + (by * ratioB);
  out[2] = (a[2] * ratioA) + (bz * ratioB);
  out[3] = (a[3] * ratioA) + (bw * ratioB);
  return out;
};
