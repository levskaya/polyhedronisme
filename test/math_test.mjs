// Test suite for modernized math module and polyhedral operations
import * as assert from 'node:assert';
import * as vec2 from '../src/math/vec2.js';
import * as vec3 from '../src/math/vec3.js';
import * as mat3 from '../src/math/mat3.js';
import * as quat from '../src/math/quat.js';
import * as geo from '../src/geo.js';

console.log('=== Running Math Subsystem Tests ===\n');

// 1. vec2 tests
console.log('Testing vec2...');
{
  const a = [1, 2];
  const b = [3, 4];
  assert.deepStrictEqual(vec2.add(a, b), [4, 6]);

  const out = [0, 0];
  vec2.add(a, b, out);
  assert.deepStrictEqual(out, [4, 6]);

  assert.deepStrictEqual(vec2.sub(a, b), [-2, -2]);
  assert.deepStrictEqual(vec2.scale(a, 2.5), [2.5, 5]);
  assert.strictEqual(vec2.dot(a, b), 11);
  assert.strictEqual(vec2.cross(a, b), -2); // 1*4 - 2*3 = -2
  assert.strictEqual(vec2.magSq(a), 5);
  assert.strictEqual(vec2.mag([3, 4]), 5);
  const n2 = vec2.normalize([3, 4]);
  assert.strictEqual(Math.abs(n2[0] - 0.6) < 1e-12, true);
  assert.strictEqual(Math.abs(n2[1] - 0.8) < 1e-12, true);
  assert.deepStrictEqual(vec2.midpoint(a, b), [2, 3]);
  assert.deepStrictEqual(vec2.lerp(a, b, 0.5), [2, 3]);

  const poly2d = [[0, 0], [4, 0], [4, 3], [0, 3]];
  assert.deepStrictEqual(vec2.calcCentroid(poly2d), [2, 1.5]);
  assert.strictEqual(vec2.planarArea(poly2d), 12);
}
console.log('✓ vec2 passed');

// 2. vec3 tests
console.log('Testing vec3...');
{
  const a = [1, 2, 3];
  const b = [4, 5, 6];
  assert.deepStrictEqual(vec3.add(a, b), [5, 7, 9]);

  const out = [0, 0, 0];
  vec3.add(a, b, out);
  assert.deepStrictEqual(out, [5, 7, 9]);

  assert.deepStrictEqual(vec3.sub(b, a), [3, 3, 3]);
  assert.deepStrictEqual(vec3.scale(a, 2), [2, 4, 6]);
  assert.deepStrictEqual(vec3.multiply(a, b), [4, 10, 18]);
  assert.strictEqual(vec3.dot(a, b), 32);

  // Cross product
  assert.deepStrictEqual(vec3.cross([1, 0, 0], [0, 1, 0]), [0, 0, 1]);
  assert.deepStrictEqual(vec3.cross([0, 1, 0], [0, 0, 1]), [1, 0, 0]);
  assert.deepStrictEqual(vec3.cross([0, 0, 1], [1, 0, 0]), [0, 1, 0]);

  assert.strictEqual(vec3.magSq([2, 3, 6]), 49);
  assert.strictEqual(vec3.mag([2, 3, 6]), 7);
  const n3 = vec3.normalize([2, 3, 6]);
  assert.strictEqual(Math.abs(n3[0] - (2 / 7)) < 1e-12, true);
  assert.strictEqual(Math.abs(n3[1] - (3 / 7)) < 1e-12, true);
  assert.strictEqual(Math.abs(n3[2] - (6 / 7)) < 1e-12, true);
  assert.deepStrictEqual(vec3.midpoint([0, 0, 0], [2, 4, 6]), [1, 2, 3]);

  // Reciprocal
  const r = vec3.reciprocal([2, 0, 0]);
  assert.deepStrictEqual(r, [0.5, 0, 0]);

  // Tangent point
  const tp = vec3.tangentPoint([1, 1, 0], [1, -1, 0]);
  assert.strictEqual(Math.abs(tp[0] - 1) < 1e-12, true);
  assert.strictEqual(Math.abs(tp[1]) < 1e-12, true);
  assert.strictEqual(Math.abs(tp[2]) < 1e-12, true);

  // Newell's normal on a flat square in XY plane
  const square = [[0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]];
  const n = vec3.normal(square);
  assert.strictEqual(Math.abs(n[0]) < 1e-12, true);
  assert.strictEqual(Math.abs(n[1]) < 1e-12, true);
  assert.strictEqual(Math.abs(n[2] - 1) < 1e-12, true);

  // Planar area
  assert.strictEqual(Math.abs(vec3.planarArea(square) - 1.0) < 1e-12, true);
}
console.log('✓ vec3 passed');

// 3. mat3 tests
console.log('Testing mat3...');
{
  const I = mat3.create();
  assert.deepStrictEqual(I, [1, 0, 0, 0, 1, 0, 0, 0, 1]);

  const v = [3, 4, 5];
  assert.deepStrictEqual(mat3.multiplyVec3(I, v), [3, 4, 5]);

  // Test matrix multiplication
  const A = [
    1, 2, 3,
    0, 1, 4,
    5, 6, 0
  ];
  const invA = mat3.invert(A);
  assert.notStrictEqual(invA, null);

  const product = mat3.multiply(A, invA);
  // product should be identity matrix
  for (let r = 0; r < 3; r++) {
    for (let c = 0; c < 3; c++) {
      const expected = r === c ? 1 : 0;
      const actual = product[(r * 3) + c];
      assert.strictEqual(Math.abs(actual - expected) < 1e-10, true, `M * M^-1 failed at (${r},${c}): ${actual}`);
    }
  }

  // Rotation around Z by 90 degrees (pi/2)
  const rotZ = mat3.fromAxisAngle([0, 0, 1], Math.PI / 2);
  const rotV = mat3.multiplyVec3(rotZ, [1, 0, 0]);
  assert.strictEqual(Math.abs(rotV[0]) < 1e-12, true);
  assert.strictEqual(Math.abs(rotV[1] - 1) < 1e-12, true);
  assert.strictEqual(Math.abs(rotV[2]) < 1e-12, true);
}
console.log('✓ mat3 passed');

// 4. quat tests
console.log('Testing quat...');
{
  const q = quat.create();
  assert.deepStrictEqual(q, [0, 0, 0, 1]);

  // Rotate [1, 0, 0] by 90 deg around Z
  const qRotZ = quat.setAxisAngle([0, 0, 1], Math.PI / 2);
  const rotV = quat.rotateVec3(qRotZ, [1, 0, 0]);
  assert.strictEqual(Math.abs(rotV[0]) < 1e-12, true);
  assert.strictEqual(Math.abs(rotV[1] - 1) < 1e-12, true);
  assert.strictEqual(Math.abs(rotV[2]) < 1e-12, true);

  // fromUnitVectors: taking [1, 0, 0] to [0, 1, 0]
  const qArc = quat.fromUnitVectors([1, 0, 0], [0, 1, 0]);
  const arcV = quat.rotateVec3(qArc, [1, 0, 0]);
  assert.strictEqual(Math.abs(arcV[0]) < 1e-12, true);
  assert.strictEqual(Math.abs(arcV[1] - 1) < 1e-12, true);
  assert.strictEqual(Math.abs(arcV[2]) < 1e-12, true);

  // Check mat3.fromQuat correspondence
  const mFromQ = mat3.fromQuat(qRotZ);
  const vViaM = mat3.multiplyVec3(mFromQ, [1, 0, 0]);
  assert.strictEqual(Math.abs(vViaM[0]) < 1e-12, true);
  assert.strictEqual(Math.abs(vViaM[1] - 1) < 1e-12, true);
  assert.strictEqual(Math.abs(vViaM[2]) < 1e-12, true);
}
console.log('✓ quat passed');

// 5. geo.js backward compatibility & project2dface bug fix
console.log('Testing geo.js and project2dface...');
{
  assert.strictEqual(typeof geo.add, 'function');
  assert.strictEqual(typeof geo.sub, 'function');
  assert.strictEqual(typeof geo.mult, 'function');
  assert.strictEqual(typeof geo.dot, 'function');
  assert.strictEqual(typeof geo.cross, 'function');
  assert.strictEqual(typeof geo.unit, 'function');
  assert.strictEqual(typeof geo.tangentPoint, 'function');
  assert.strictEqual(typeof geo.calcCentroid, 'function');
  assert.strictEqual(typeof geo.normal, 'function');
  assert.strictEqual(typeof geo.planararea, 'function');

  // Verify project2dface does NOT produce NaN and gives valid 2D coordinates
  const face3d = [[1, 0, 0], [1, 2, 0], [1, 2, 2], [1, 0, 2]]; // square of side 2 in x=1 plane
  const face2d = geo.project2dface(face3d);
  assert.strictEqual(face2d.length, 4);
  for (let pt of face2d) {
    assert.strictEqual(isNaN(pt[0]), false, `project2dface produced NaN for x: ${pt}`);
    assert.strictEqual(isNaN(pt[1]), false, `project2dface produced NaN for y: ${pt}`);
  }
  // Area in 2D should equal area in 3D (4.0)
  const area2d = vec2.planarArea(face2d);
  assert.strictEqual(Math.abs(area2d - 4.0) < 1e-10, true, `Expected area 4.0, got ${area2d}`);
}
console.log('✓ geo.js and project2dface passed');

// 6. Polyhedra & Solvers Tests
console.log('Testing Polyhedron primitives and canonicalization...');
{
  const { cube, tetrahedron, octahedron } = await import('../src/polyhedron.js');
  const { dual, kisN, ambo } = await import('../src/topo_operators.js');
  const { canonicalize, sphericalTutteUntangle } = await import('../src/geo_operators.js');

  const c = cube();
  assert.strictEqual(c.vertices.length, 8);
  assert.strictEqual(c.faces.length, 6);
  assert.strictEqual(c.edges().length, 12);
  // Euler characteristic: V - E + F = 2
  assert.strictEqual(c.vertices.length - c.edges().length + c.faces.length, 2);

  // Dual of Cube should be Octahedron (6 vertices, 8 faces, 12 edges)
  const dCube = dual(c);
  assert.strictEqual(dCube.vertices.length, 6);
  assert.strictEqual(dCube.faces.length, 8);
  assert.strictEqual(dCube.edges().length, 12);
  assert.strictEqual(dCube.vertices.length - dCube.edges().length + dCube.faces.length, 2);

  // Kis of Octahedron
  const kOct = kisN(octahedron(), 0, 0.1);
  assert.strictEqual(kOct.faces.length, 24); // 8 * 3 triangular faces
  assert.strictEqual(kOct.vertices.length - kOct.edges().length + kOct.faces.length, 2);

  // Canonicalize Cube
  const canonC = canonicalize(c, 20);
  assert.strictEqual(canonC.vertices.length, 8);
  assert.strictEqual(canonC.faces.length, 6);
  // Verify all vertices are valid numbers
  for (let v of canonC.vertices) {
    assert.strictEqual(isNaN(v[0]), false);
    assert.strictEqual(isNaN(v[1]), false);
    assert.strictEqual(isNaN(v[2]), false);
  }

  // Spherical Tutte untangler
  const untangled = sphericalTutteUntangle(c, 10);
  assert.strictEqual(untangled.length, 8);
  for (let v of untangled) {
    const r = vec3.mag(v);
    assert.strictEqual(Math.abs(r - 1.0) < 1e-6, true, `Expected unit radius, got ${r}`);
  }
  // Dodecahedron & Icosahedron
  const { dodecahedron, icosahedron } = await import('../src/polyhedron.js');
  const { gyro } = await import('../src/topo_operators.js');

  const d = dodecahedron();
  assert.strictEqual(d.faces.length, 12);
  assert.strictEqual(d.vertices.length, 20);
  assert.strictEqual(d.edges().length, 30);
  assert.strictEqual(d.vertices.length - d.edges().length + d.faces.length, 2);

  const ic = icosahedron();
  assert.strictEqual(ic.faces.length, 20);
  assert.strictEqual(ic.vertices.length, 12);
  assert.strictEqual(ic.edges().length, 30);
  assert.strictEqual(ic.vertices.length - ic.edges().length + ic.faces.length, 2);

  // Ambo of Cube (Cuboctahedron: 14 faces, 12 vertices, 24 edges)
  const aCube = ambo(c);
  assert.strictEqual(aCube.faces.length, 14);
  assert.strictEqual(aCube.vertices.length, 12);
  assert.strictEqual(aCube.edges().length, 24);
  assert.strictEqual(aCube.vertices.length - aCube.edges().length + aCube.faces.length, 2);

  // Gyro of Tetrahedron (3 * 4 = 12 pentagonal faces)
  const gTet = gyro(tetrahedron());
  assert.strictEqual(gTet.faces.length, 12);
  assert.strictEqual(gTet.vertices.length - gTet.edges().length + gTet.faces.length, 2);

  // Canonicalize Dodecahedron and check midsphere edge tangency
  const canonD = canonicalize(d, 30);
  assert.strictEqual(canonD.faces.length, 12);
  assert.strictEqual(canonD.vertices.length, 20);
  for (let e of canonD.edges()) {
    const tp = geo.tangentPoint(canonD.vertices[e[0]], canonD.vertices[e[1]]);
    const dist = geo.mag(tp);
    // Average edge distance should be centered around unit sphere
    assert.strictEqual(dist > 0.5 && dist < 1.5, true, `Edge distance out of range: ${dist}`);
  }
}
console.log('✓ Polyhedron primitives, topological operators, and canonicalizer passed');

// 7. Legacy Parser & Recipe Rewriter Tests
console.log('Testing legacy v1 parser and recipe rewriter...');
{
  const { createRequire } = await import('node:module');
  const req = createRequire(import.meta.url);
  const fs = req('node:fs');
  const path = req('node:path');
  const pegjs = req('pegjs');
  const grammar = fs.readFileSync(path.join(path.dirname(new URL(import.meta.url).pathname), '../src/legacy.pegjs'), 'utf8');
  const { rewriteLegacyRecipe, generatePolyLegacy, setLegacyOpParser } = await import('../src/legacy_parser.js');
  setLegacyOpParser(pegjs.generate(grammar));

  // Verify KN and CN rewrites (and preservation of base C)
  assert.strictEqual(rewriteLegacyRecipe('C2dakD'), 'CdakD');
  assert.strictEqual(rewriteLegacyRecipe('K50dakD'), 'CdakD');
  assert.strictEqual(rewriteLegacyRecipe('C200tC'), 'CtC');
  assert.strictEqual(rewriteLegacyRecipe('K10C'), 'CC');
  assert.strictEqual(rewriteLegacyRecipe('oC20kkkT'), 'oCkkkT');
  assert.strictEqual(rewriteLegacyRecipe('kn4C40A0dA4'), 'kn4CA0dA4');

  // Verify renamed operator rewrites: l -> P, z -> Z, h(0,inset,depth) -> H(inset,depth), h -> H
  assert.strictEqual(rewriteLegacyRecipe('lC'), 'PC');
  assert.strictEqual(rewriteLegacyRecipe('zD'), 'ZD');
  assert.strictEqual(rewriteLegacyRecipe('h(0,0.4,0.1)C'), 'H(0.4,0.1)C');
  assert.strictEqual(rewriteLegacyRecipe('hT'), 'HT');
  assert.strictEqual(rewriteLegacyRecipe('C100zh(0,0.5,0.2)lD'), 'CZH(0.5,0.2)PD');

  // Verify standalone legacy polyhedron generation
  const p = generatePolyLegacy('C2dakD');
  assert.strictEqual(p.faces.length > 0, true);
  assert.strictEqual(p.vertices.length > 0, true);
}
console.log('✓ Legacy v1 parser and recipe rewriter passed');

console.log('\nAll math and geometry tests passed successfully!');
