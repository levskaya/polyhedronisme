// Polyhédronisme
// ===================================================================================================
//
// Color Palettes, Perceptual Color Math (OKLCH), and Harmonic Random Generators
//
// Copyright 2026, Anselm Levskaya
// Released under the MIT License

import { random, PI, cos, sin } from './geo.js';

// Default landing palette (Classic Red-White-Black/Slate on warm paper)
export const rwb_palette = [
  '#ff7777', '#dddddd', '#889999', '#fff0e5',
  '#aa3333', '#ff0000', '#ffffff', '#aaaaaa'
];

// Curated preset palettes tailored for 3D cosine-shaded polyhedra on #f0f0e5 paper
export const PRESET_PALETTES = [
  {
    id: 'rwb',
    name: 'Classic RWB',
    colors: [...rwb_palette]
  },
  {
    id: 'rwb-refined',
    name: 'Warm Editorial',
    colors: ['#e86f68', '#d9d7ce', '#7c9094', '#f7efe4', '#9e3333', '#c99b52', '#525c60', '#b5b2a8']
  },
  {
    id: 'ukiyoe',
    name: 'Ukiyo-e Woodblock',
    colors: ['#d15b47', '#e8dfd1', '#3b5d6e', '#e3b464', '#8fa382', '#6e6a64', '#c78c78', '#253842']
  },
  {
    id: 'bauhaus',
    name: 'Bauhaus',
    colors: ['#d14938', '#ede6d6', '#2d5a7b', '#e5a93c', '#9e9a90', '#4a4843', '#7a9a95', '#c97b4b']
  },
  {
    id: 'celadon',
    name: 'Celadon & Clay',
    colors: ['#8ab09c', '#efe5d5', '#c86d51', '#68787d', '#d9a05b', '#8c3f2b', '#b8c4b8', '#4d585a']
  },
  {
    id: 'delft',
    name: 'Delft Ceramic',
    colors: ['#2b5782', '#f5f2e9', '#678fb8', '#c9d9e8', '#8c9ba5', '#19344f', '#d9b478', '#46698c']
  },
  {
    id: 'mineral',
    name: 'Mineral & Stone',
    colors: ['#c29b78', '#e6e2d8', '#6e7f80', '#a36652', '#8f9982', '#534f4a', '#d6c5b3', '#7a6354']
  },
  {
    id: 'origami',
    name: 'Origami Pastel',
    colors: ['#e8998d', '#eef0e8', '#89b0ae', '#f2d096', '#b6a6ca', '#a3b899', '#d6c2b2', '#7c8483']
  }
];

// Normalize any 3-digit or 6-digit hex string to lowercase #rrggbb
export const normalizeHex = function (hexstr, fallback = '#cccccc') {
  if (typeof hexstr !== 'string') return fallback;
  let s = hexstr.trim();
  if (s.startsWith('#')) s = s.slice(1);
  if (/^[0-9a-fA-F]{3}$/.test(s)) {
    return '#' + s.split('').map(c => c + c).join('').toLowerCase();
  }
  if (/^[0-9a-fA-F]{6}$/.test(s)) {
    return '#' + s.toLowerCase();
  }
  return fallback;
};

// Converts #xxxxxx / #xxx format into list of [r,g,b] floats in [0, 1]
export const hextofloats = function (hexstr) {
  const norm = normalizeHex(hexstr, '#cccccc').slice(1);
  return [
    parseInt(norm.slice(0, 2), 16) / 255,
    parseInt(norm.slice(2, 4), 16) / 255,
    parseInt(norm.slice(4, 6), 16) / 255
  ];
};

// Clamp float to [0, 1]
const clamp01 = x => Math.max(0, Math.min(1, x));

// Converts [r,g,b] floats in [0, 1] to #xxxxxx form
export const floatstohex = function (rgb) {
  const r_hex = Math.round(clamp01(rgb[0]) * 255).toString(16).padStart(2, '0');
  const g_hex = Math.round(clamp01(rgb[1]) * 255).toString(16).padStart(2, '0');
  const b_hex = Math.round(clamp01(rgb[2]) * 255).toString(16).padStart(2, '0');
  return '#' + r_hex + g_hex + b_hex;
};

// Converts [h,s,l] float args in [0, 1] to [r,g,b] list
export const hsl2rgb = function (h, s, l) {
  let r, g, b;
  h = ((h % 1) + 1) % 1;
  s = clamp01(s);
  l = clamp01(l);
  if (s === 0) {
    r = g = b = l;
  } else {
    const hue2rgb = function (p, q, t) {
      if (t < 0) t += 1;
      if (t > 1) t -= 1;
      if (t < 1 / 6) return p + (q - p) * 6 * t;
      if (t < 1 / 2) return q;
      if (t < 2 / 3) return p + (q - p) * (2 / 3 - t) * 6;
      return p;
    };
    const q = l < 0.5 ? l * (1 + s) : l + s - l * s;
    const p = 2 * l - q;
    r = hue2rgb(p, q, h + 1 / 3);
    g = hue2rgb(p, q, h);
    b = hue2rgb(p, q, h - 1 / 3);
  }
  return [r, g, b];
};

// Perceptually uniform OKLCH -> sRGB [r, g, b] conversion
// L in [0, 1], C in [0, ~0.35], hDeg in [0, 360)
const linearToSrgb = c => (c <= 0.0031308 ? 12.92 * c : 1.055 * Math.pow(c, 1 / 2.4) - 0.055);

const oklchToLinearRgb = function (L, C, hRad) {
  const a = C * cos(hRad);
  const b = C * sin(hRad);

  const l_ = L + 0.3963377774 * a + 0.2158037573 * b;
  const m_ = L - 0.1055613458 * a - 0.0638541728 * b;
  const s_ = L - 0.0894841775 * a - 1.2914855480 * b;

  const l = l_ * l_ * l_;
  const m = m_ * m_ * m_;
  const s = s_ * s_ * s_;

  return [
    +4.0767416621 * l - 3.3077115913 * m + 0.2309699292 * s,
    -1.2684380046 * l + 2.6097574011 * m - 0.3413193965 * s,
    -0.0041960863 * l - 0.7034186147 * m + 1.7076147010 * s
  ];
};

export const oklchToHex = function (L, C, hDeg) {
  const hRad = (((hDeg % 360) + 360) % 360) * (PI / 180);
  const targetL = clamp01(L);
  // Binary-search chroma reduction if outside sRGB gamut so hue & lightness stay true
  let cLow = 0;
  let cHigh = Math.max(0, C);
  let rgbLin = oklchToLinearRgb(targetL, cHigh, hRad);
  const inGamut = rgb =>
    rgb[0] >= -0.001 && rgb[0] <= 1.001 &&
    rgb[1] >= -0.001 && rgb[1] <= 1.001 &&
    rgb[2] >= -0.001 && rgb[2] <= 1.001;

  if (!inGamut(rgbLin)) {
    for (let i = 0; i < 8; i++) {
      const cMid = 0.5 * (cLow + cHigh);
      const test = oklchToLinearRgb(targetL, cMid, hRad);
      if (inGamut(test)) {
        cLow = cMid;
        rgbLin = test;
      } else {
        cHigh = cMid;
      }
    }
  }

  return floatstohex([
    clamp01(linearToSrgb(clamp01(rgbLin[0]))),
    clamp01(linearToSrgb(clamp01(rgbLin[1]))),
    clamp01(linearToSrgb(clamp01(rgbLin[2])))
  ]);
};

// Generate a tasteful, color-theory-driven palette of `count` hex colors
export const rndcolors = function (count = 8) {
  const baseHue = random() * 360;
  const strategy = Math.floor(random() * 4);
  const specs = [];

  if (strategy === 0) {
    // 1. Editorial Accent + Neutrals (in the spirit of Classic RWB)
    const compHue = (baseHue + 165 + random() * 30) % 360;
    specs.push(
      [0.68, 0.14, baseHue],                      // soft primary accent
      [0.88, 0.015, baseHue],                     // light warm/cool grey
      [0.62, 0.035, compHue],                     // muted slate complement
      [0.95, 0.025, (baseHue + 40) % 360],        // warm ivory / paper tint
      [0.48, 0.13, baseHue],                      // deep primary shade
      [0.74, 0.11, (baseHue + 35) % 360],         // secondary warm/analogous accent
      [0.54, 0.08, compHue],                      // deeper complement
      [0.76, 0.02, compHue]                       // mid neutral pewter
    );
  } else if (strategy === 1) {
    // 2. Split-Complementary Ceramic
    const h2 = (baseHue + 150) % 360;
    const h3 = (baseHue + 210) % 360;
    specs.push(
      [0.64, 0.13, baseHue],
      [0.92, 0.02, baseHue],
      [0.60, 0.10, h2],
      [0.78, 0.09, h3],
      [0.48, 0.11, baseHue],
      [0.72, 0.03, h2],
      [0.52, 0.09, h3],
      [0.84, 0.05, (baseHue + 30) % 360]
    );
  } else if (strategy === 2) {
    // 3. Triadic Mineral (muted saturation, distinct lightness steps)
    const dir = random() < 0.5 ? 1 : -1;
    const h2 = (baseHue + dir * 120 + 360) % 360;
    const h3 = (baseHue + dir * 240 + 360) % 360;
    specs.push(
      [0.65, 0.12, baseHue],
      [0.91, 0.025, (baseHue + 20) % 360],
      [0.58, 0.10, h2],
      [0.76, 0.11, h3],
      [0.49, 0.09, baseHue],
      [0.74, 0.07, h2],
      [0.55, 0.03, h3],
      [0.83, 0.06, baseHue]
    );
  } else {
    // 4. Analogous Woodblock / Golden-Angle Harmony
    const step = 32 + random() * 20;
    const lightnessRing = [0.66, 0.90, 0.56, 0.78, 0.48, 0.72, 0.60, 0.84];
    const chromaRing = [0.12, 0.02, 0.10, 0.11, 0.11, 0.04, 0.09, 0.06];
    for (let i = 0; i < 8; i++) {
      specs.push([
        lightnessRing[i],
        chromaRing[i],
        (baseHue + i * step) % 360
      ]);
    }
  }

  const out = [];
  for (let i = 0; i < count; i++) {
    const [L, C, H] = specs[i % specs.length];
    const jitterL = clamp01(L + (random() - 0.5) * 0.04);
    const jitterC = Math.max(0.01, C + (random() - 0.5) * 0.02);
    out.push(oklchToHex(jitterL, jitterC, H));
  }
  return out;
};

// Permute palette colors so face-class assignments change visibly
export const shufflePalette = function (pal) {
  if (!Array.isArray(pal) || pal.length <= 1) return [...(pal || rwb_palette)];
  const next = [...pal];
  // Fisher-Yates shuffle, ensuring at least the first few slots change
  for (let attempt = 0; attempt < 5; attempt++) {
    for (let i = next.length - 1; i > 0; i--) {
      const j = Math.floor(random() * (i + 1));
      const tmp = next[i];
      next[i] = next[j];
      next[j] = tmp;
    }
    const changed = next.some((c, idx) => c !== pal[idx]);
    if (changed) break;
  }
  return next;
};
