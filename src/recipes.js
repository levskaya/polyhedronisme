// Polyhédronisme
// ===================================================================================================
//
// Curated Polyhedron Recipe Gallery & Smart Procedural Recipe Generator
//
// Copyright 2026, Anselm Levskaya
// Released under the MIT License

import { random, randomchoice } from './geo.js';
import { generatePoly } from './parser.js';

// Categorized gallery of classic and modern polyhedral recipes
export const PRESET_RECIPE_GROUPS = [
  {
    label: 'Classic Conway & Catalan',
    recipes: [
      { spec: 'CdakD', name: 'CdakD — Canonical Ambo-Kis Dodecahedron' },
      { spec: 'CsD', name: 'CsD — Canonical Snub Dodecahedron' },
      { spec: 'CgD', name: 'CgD — Pentagonal Hexecontahedron' },
      { spec: 'CbD', name: 'CbD — Great Rhombicosidodecahedron' },
      { spec: 'CwdD', name: 'CwdD — Canonical Whirled Dual Dodecahedron' },
      { spec: 'Ctk5dI', name: 'Ctk5dI — Truncated Pentakis Dodecahedron' },
      { spec: 'aobD', name: 'aobD — Ambo-Ortho-Bevel Dodecahedron' },
      { spec: 'oCkkkT', name: 'oCkkkT — Ortho-Canonical Triple-Kis Tetrahedron' }
    ]
  },
  {
    label: 'Sculptural & Extruded',
    recipes: [
      { spec: 'n(5,0.4,0.35)n(6,0.45,-0.25)tI', name: 'n(5,…)n(6,…)tI — Towers & Dimples Buckyball' },
      { spec: 'k(5,-0.35)n(6,0.4,0.25)tI', name: 'k(5,-0.35)n(6,…)tI — Cratered Hex-Rim Fullerene' },
      { spec: 'n(6,0.35,0.45)cD', name: 'n(6,0.35,0.45)cD — Hex-Turret Chamfered Dodecahedron' },
      { spec: 'k(5,-0.4)D', name: 'k(5,-0.4)D — Concave Star-Dimpled Dodecahedron' },
      { spec: 'n(5,0.35,0.45)n(5,0.4,-0.3)D', name: 'n(5,…)n(5,…)D — Stepped Pentagonal Ziggurat' },
      { spec: 'n(8,0.45,-0.38)CbC', name: 'n(8,0.45,-0.38)CbC — Octagonal-Well Bevel Cube' },
      { spec: 'knD', name: 'knD — Kis-Inset Dodecahedron' }
    ]
  },
  {
    label: 'Skeletal Frames (Hollow)',
    recipes: [
      { spec: 'H(0.55,0.12)CbD', name: 'H(0.55,0.12)CbD — Leonardo Frame (Bevel D)' },
      { spec: 'H(0.5,0.15)CcI', name: 'H(0.5,0.15)CcI — Hollow Chamfered Icosahedron' },
      { spec: 'H(0.58,0.14)aD', name: 'H(0.58,0.14)aD — Open-Frame Icosidodecahedron' },
      { spec: 'H(0.48,0.16)CsC', name: 'H(0.48,0.16)CsC — Chiral Snub-Cube Cage' },
      { spec: 'H(0.5,0.12)qD', name: 'H(0.5,0.12)qD — Quinto Stained-Glass Lattice' },
      { spec: 'H(0.52,0.14)CtI', name: 'H(0.52,0.14)CtI — Skeletal Buckyball Cage' },
      { spec: 'H(0.45,0.15)wC', name: 'H(0.45,0.15)wC — Whirled Cube Lattice' }
    ]
  },
  {
    label: 'Goldberg, Fullerenes & Geodesics',
    recipes: [
      { spec: 'Cdu3I', name: 'Cdu3I — Goldberg GP(3,0) Fullerene (92 faces)' },
      { spec: 'Cdu4I', name: 'Cdu4I — Goldberg GP(4,0) Fullerene (162 faces)' },
      { spec: 'Cdu5T', name: 'Cdu5T — Tetrahedral Goldberg Sphere' },
      { spec: 'Cdu4O', name: 'Cdu4O — Octahedral Goldberg Sphere' },
      { spec: 'k6Cdu3I', name: 'k6Cdu3I — Triangulated-Hex Geodesic Dome' },
      { spec: 'n(6,0.45,-0.15)Cdu3I', name: 'n(6,0.45,-0.15)Cdu3I — Dimpled Viral Capsid' },
      { spec: 'H(0.55,0.1)Cdu3I', name: 'H(0.55,0.1)Cdu3I — Geodesic Fullerene Cage' }
    ]
  },
  {
    label: 'Chiral Swirls & Gyres',
    recipes: [
      { spec: 'n5wD', name: 'n5wD — Whirled Dodecahedron w/ Inset Eyes' },
      { spec: 'k5gD', name: 'k5gD — Pentagonal-Starburst Gyro Dodecahedron' },
      { spec: 'pwI', name: 'pwI — Double-Twist Propellor-Whirl Icosahedron' },
      { spec: 'opD', name: 'opD — Ortho-Propellor Dodecahedron' },
      { spec: 'n5pD', name: 'n5pD — Recessed-Pentagon Propellor Dodecahedron' },
      { spec: 'CpdA6', name: 'CpdA6 — Canonical Propellor Hexagonal Bipyramid' },
      { spec: 'waA5', name: 'waA5 — Whirled Ambo Pentagonal Antiprism' }
    ]
  },
  {
    label: 'Renaissance & Filigree',
    recipes: [
      { spec: 'PPD', name: 'PPD — Jamnitzer Double-Perspectiva Dodecahedron' },
      { spec: 'PT', name: 'PT — Jamnitzer Perspectiva Tetrahedron' },
      { spec: 'PqD', name: 'PqD — Perspectiva Quinto-Dodecahedron Star' },
      { spec: 'qaxI', name: 'qaxI — Quinto-Ambo-Extrude Icosahedron' },
      { spec: 'CKdD', name: 'CKdD — Canonical Stake-Dual Dodecahedron' },
      { spec: 'LbaC', name: 'LbaC — Laced Bevel-Ambo Cube Filigree' },
      { spec: 'IqI', name: 'IqI — Woven Join-Lace Quinto Icosahedron' },
      { spec: 'n5qD', name: 'n5qD — Popped-Pentagon Quinto Dodecahedron' },
      { spec: 'qA40PA40oC', name: 'qA40PA40oC — Quinto-Perspectiva Ortho-Cube' }
    ]
  },
  {
    label: 'Exotic Seeds (Johnson, Cupolas, Prisms)',
    recipes: [
      { spec: 'Ck3J91', name: 'Ck3J91 — Canonical Kis-3 Bilunabirotunda (J91)' },
      { spec: 'qqJ37', name: 'qqJ37 — Double-Quinto Gyrate Rhombicuboctahedron' },
      { spec: 'n5qJ72', name: 'n5qJ72 — Inset-Quinto Gyrate Rhombicosidodecahedron' },
      { spec: 'CwJ84', name: 'CwJ84 — Canonical Whirled Snub Disphenoid (J84)' },
      { spec: 'n4qJ26', name: 'n4qJ26 — Inset-Quinto Gyrobifastigium (J26)' },
      { spec: 'k5n10U5', name: 'k5n10U5 — Crowned Pentagonal Cupola (U5)' },
      { spec: 'kn4CA0dA4', name: 'kn4CA0dA4 — Kis-Inset Canonical Dual Antiprism' },
      { spec: 'oox4P7', name: 'oox4P7 — Double-Ortho Extruded Heptagonal Prism' }
    ]
  }
];

// Flat list of all curated recipes
export const DEFAULT_RECIPES = PRESET_RECIPE_GROUPS.flatMap(g => g.recipes.map(r => r.spec));

// Curated subset of Johnson solids with high symmetry or compelling face mixtures
const INTERESTING_JOHNSON_IDS = [
  12, 13, 14, 15, 16, 17, 18, 22, 26, 27, 28, 29, 30, 31,
  34, 35, 37, 44, 45, 48, 51, 58, 63, 65, 66, 67, 68, 72,
  73, 74, 75, 84, 85, 86, 88, 89, 90, 91, 92
];

// Inspect distinct face sizes and vertex degrees of a polyhedron
const analyzeTopology = function (poly) {
  const faceSizeCounts = new Map();
  const vertDegCounts = new Array(poly.vertices.length).fill(0);

  for (const f of poly.faces) {
    const len = f.length;
    faceSizeCounts.set(len, (faceSizeCounts.get(len) || 0) + 1);
    for (const v of f) {
      vertDegCounts[v]++;
    }
  }

  const degSet = new Set(vertDegCounts.filter(d => d >= 3));
  const faceSizes = Array.from(faceSizeCounts.keys()).sort((a, b) => a - b);
  const vertDegrees = Array.from(degSet).sort((a, b) => a - b);
  const allTriangles = faceSizes.length === 1 && faceSizes[0] === 3;

  return { faceSizes, faceSizeCounts, vertDegrees, allTriangles };
};

// Pick a random base seed polyhedron string
const randomSeedSpec = function () {
  const roll = random();
  if (roll < 0.52) {
    return randomchoice(['T', 'C', 'O', 'I', 'D', 'C', 'I', 'D']);
  } else if (roll < 0.72) {
    const family = randomchoice(['P', 'A', 'Y']);
    const n = family === 'Y'
      ? randomchoice([3, 4, 5, 6])
      : randomchoice([3, 4, 5, 6, 7, 8]);
    return `${family}${n}`;
  } else if (roll < 0.82) {
    const family = randomchoice(['U', 'V']);
    const n = randomchoice([3, 4, 5]);
    return `${family}${n}`;
  } else {
    const jId = randomchoice(INTERESTING_JOHNSON_IDS);
    return `J${jId}`;
  }
};

// Smart procedural recipe generator:
// Chains a random seed + 1-3 compatible topological operators (respecting face-count budgets
// and actual face-sidedness / vertex degrees), optional canonicalization, and optional 3D
// sculptural finishers (parameterized inset/extrude, negative kis dimples, hollow frames, etc.).
export const generateRandomRecipe = function (avoidSpec = '') {
  for (let attempt = 0; attempt < 12; attempt++) {
    try {
      const seed = randomSeedSpec();
      let spec = seed;
      let poly = generatePoly(spec);
      let topo = analyzeTopology(poly);

      // Path A: Goldberg / Subdivision path on all-triangle seeds (~15% chance when triangular)
      if (topo.allTriangles && poly.faces.length <= 20 && random() < 0.25) {
        const subN = poly.faces.length <= 8 ? randomchoice([3, 4, 5]) : randomchoice([2, 3, 4]);
        spec = `Cdu${subN}${spec}`;
        poly = generatePoly(spec);
        topo = analyzeTopology(poly);

        // Optional finisher on Goldberg polyhedron (e.g. hex inset, kis6, or hollow)
        const finRoll = random();
        if (finRoll < 0.35 && topo.faceSizes.includes(6)) {
          const inset = randomchoice([0.35, 0.4, 0.45]);
          const depth = randomchoice([-0.2, -0.15, 0.2, 0.25]);
          spec = `n(6,${inset},${depth})${spec}`;
        } else if (finRoll < 0.55 && topo.faceSizes.includes(6)) {
          spec = `k6${spec}`;
        } else if (finRoll < 0.72 && poly.faces.length <= 120) {
          const inset = randomchoice([0.5, 0.55]);
          const thick = randomchoice([0.1, 0.12]);
          spec = `H(${inset},${thick})${spec}`;
        }

        if (spec !== avoidSpec) {
          const finalPoly = generatePoly(spec);
          if (finalPoly.faces.length >= 12 && finalPoly.faces.length <= 900) {
            return spec;
          }
        }
        continue;
      }

      // Path B: Conway / Topological operator chain
      // Keep chains shorter on lower-symmetry seeds (P, A, Y, U, V, J) so face classes stay crisp
      const isPlatonic = ['T', 'C', 'O', 'I', 'D'].includes(seed);
      const targetSteps = isPlatonic
        ? (random() < 0.3 ? 1 : (random() < 0.8 ? 2 : 3))
        : (random() < 0.65 ? 1 : 2);
      let lastOp = '';

      for (let step = 0; step < targetSteps; step++) {
        const nF = poly.faces.length;
        if (nF > 150) break;

        const candidates = [];

        // Light multipliers (~1x-2.5x)
        candidates.push('a', 'j', 't');
        if (lastOp !== 'd') candidates.push('d');

        // Selective kis / truncate based on actual face sizes & vertex degrees
        if (topo.faceSizes.length >= 2) {
          const targetF = randomchoice(topo.faceSizes);
          candidates.push(`k${targetF}`, `k${targetF}`);
        } else {
          candidates.push('k');
        }
        if (topo.vertDegrees.length >= 2) {
          const targetD = randomchoice(topo.vertDegrees);
          candidates.push(`t${targetD}`);
        }

        // Medium multipliers (~3x-4.5x) if face count is moderate
        if (nF <= 80) {
          candidates.push('c', 'p', 'L', 'I', 'K', 'e', 'o');
        }

        // Heavy multipliers (~5x-7x) only when face count is small
        if (nF <= 42) {
          candidates.push('g', 's', 'q', 'b', 'm', 'w');
        }

        const chosenOp = randomchoice(candidates);
        spec = `${chosenOp}${spec}`;
        lastOp = chosenOp;
        poly = generatePoly(spec);
        topo = analyzeTopology(poly);
      }

      // Optional Canonicalization before 3D sculptural relief
      const nF = poly.faces.length;
      let canonicalized = false;
      if (nF >= 10 && nF <= 220 && random() < 0.55) {
        const cOp = nF <= 90 ? 'C' : 'C25';
        spec = `${cOp}${spec}`;
        poly = generatePoly(spec);
        topo = analyzeTopology(poly);
        canonicalized = true;
      }

      // Optional Sculptural Finisher (or mandatory if only 1 light op was applied and not canonicalized)
      const FinisherProb = (!canonicalized && targetSteps === 1) ? 0.85 : 0.55;
      if (random() < FinisherProb && poly.faces.length <= 260) {
        const finRoll = random();
        const nonTriSizes = topo.faceSizes.filter(s => s >= 4);
        const pickSize = nonTriSizes.length > 0 ? randomchoice(nonTriSizes) : randomchoice(topo.faceSizes);

        if (finRoll < 0.38 && pickSize) {
          // Parameterized inset/extrude on a valid face size
          const inset = randomchoice([0.35, 0.4, 0.45, 0.5]);
          const depth = randomchoice([-0.35, -0.25, -0.18, 0.22, 0.32, 0.42]);
          spec = `n(${pickSize},${inset},${depth})${spec}`;

          // If there's a second distinct polygon size, occasionally add a contrasting extrusion/dimple!
          const otherSizes = nonTriSizes.filter(s => s !== pickSize);
          if (otherSizes.length > 0 && poly.faces.length <= 120 && random() < 0.45) {
            const pickSize2 = randomchoice(otherSizes);
            const inset2 = randomchoice([0.35, 0.4, 0.45]);
            const depth2 = depth > 0 ? randomchoice([-0.3, -0.2]) : randomchoice([0.25, 0.35]);
            spec = `n(${pickSize2},${inset2},${depth2})${spec}`;
          }
        } else if (finRoll < 0.56 && pickSize) {
          // Parameterized kis (either outward spikes or inward concave dimples)
          const depth = randomchoice([-0.4, -0.3, -0.22, 0.25, 0.35]);
          spec = `k(${pickSize},${depth})${spec}`;
        } else if (finRoll < 0.72 && pickSize) {
          // Classic concise extrudeN or insetN shorthand
          const op = randomchoice(['n', 'x']);
          spec = `${op}${pickSize}${spec}`;
        } else if (finRoll < 0.88 && poly.faces.length <= 130) {
          // Hollow skeletal frame
          const inset = randomchoice([0.45, 0.5, 0.55, 0.6]);
          const thick = randomchoice([0.1, 0.12, 0.15]);
          spec = `H(${inset},${thick})${spec}`;
        } else if (poly.faces.length <= 75) {
          // Jamnitzer Perspectiva stellation
          spec = `P${spec}`;
        }
      }

      if (spec === avoidSpec) continue;

      const finalPoly = generatePoly(spec);
      const validVertices = finalPoly.vertices.every(
        v => Number.isFinite(v[0]) && Number.isFinite(v[1]) && Number.isFinite(v[2])
      );
      const nClasses = finalPoly.num_face_classes || 1;
      if (
        validVertices &&
        finalPoly.faces.length >= 12 &&
        finalPoly.faces.length <= 950 &&
        nClasses >= 2 &&
        nClasses <= 14
      ) {
        return spec;
      }
    } catch (e) {
      // Retry if any random combination hit a degenerate edge case
    }
  }

  // Fallback to a random curated recipe
  return randomchoice(DEFAULT_RECIPES.filter(r => r !== avoidSpec));
};
