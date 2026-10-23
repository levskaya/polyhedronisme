// Polyhédronisme
// ===================================================================================================
//
// A toy for constructing and manipulating polyhedra and other meshes
//
// Copyright 2026, Anselm Levskaya
// Released under the MIT License

/* eslint-disable camelcase */

import { add, clone, eye3, mv3, mm3, getVec2VecRotM, perspT, invperspT,
  vec_rotm, rotm, mult, dot, unit, normal, round, randomchoice, PI, quat, mat3 } from './geo.js';
import { rwb_palette, sortFaces, palette, rndcolors } from './polyhedron.js';
import { generatePoly } from './parser.js';
import { rewriteLegacyRecipe, setLegacyOpParser } from './legacy_parser.js';
import legacyParser from './legacy.pegjs';
import { triangulate } from './topo_operators.js';
import { toOBJ, toVRML, toSVG } from './exporters.js';
import { sortFacesCarefully } from './depth.js';

// docs
import frontHtml from './components/front.md';

// global css
import './assets/polyhedronisme.css';

// Initialize legacy recipe parser
setLegacyOpParser(legacyParser);

// GLOBALS
// ===================================================================================================
export let ctx = {}; // for global access to canvas context
export let globPolys = {}; // constructed polyhedras

const CANVAS_WIDTH = 720; // canvas dims
const CANVAS_HEIGHT = 500; // canvas dims
const PNG_EXPORT_SCALE = 4; // scale factor for high-res PNG exports
let globQuat = quat.create();
let globLastQuat = quat.create();
let globRotM = mat3.create();
const updateRotM = function () {
  mat3.fromQuat(globQuat, globRotM);
};
export let perspective_scale = 1200;
export const persp_z_max = 5;
export const persp_z_min = 0;
export const persp_ratio = 0.8;
const _2d_x_offset = CANVAS_WIDTH / 2;
const _2d_y_offset = CANVAS_HEIGHT / 2;

export let SORTING_ALGO = 'FAST'; // or 'CAREFUL';
// window.set_algo = x => (SORTING_ALGO = x);

// Back-face culling cutoff disabled for now (set to -Infinity so all faces pass through)
export const VISIBILITY_ZNORMAL_CUTOFF = -Infinity; // originally -0.1
export var PALETTE = rwb_palette;
const BG_CLEAR = true; // clear background or colored?
const BG_COLOR = 'rgba(255,255,255,1.0)'; // background color
export let COLOR_METHOD = 'signature'; // "area", "edges"
export let COLOR_SENSITIVITY = 2; // color sensitivity to variation
// in congruence signature or planar area
const ctx_linewidth = 0.5; // for outline of faces
let PaintMode = 'fillstroke';

// mouse event variables
let MOUSEDOWN = false;
let LastMouseX = 0;
let LastMouseY = 0;
// state variable for 3d trackball
let LastSphVec = [1, 0, 0];

// random grabbag of polyhedra
const DEFAULT_RECIPES = [
  'CdakD', 'oCkkkT', 'kn4CA0dA4', 'opD',
  'PT', 'qA40PA40oC', 'knD', 'dn6x4bT', 'oox4P7',
  'qqJ37', 'aobD', 'qaxI'];

// File-saving helpers used to export txt/svg/obj/wrl/canvas-png
const downloadBlob = function (blob, filename) {
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  URL.revokeObjectURL(url);
};

const saveText = function (text, filename) {
  const blb = new Blob([text], { type: 'text/plain;charset=utf-8' });
  downloadBlob(blb, filename);
};

// parses URL string for polyhedron recipe, for bookmarking
const parseURL = function () {
  const urlParams = {};
  const params = new URLSearchParams(window.location.search);
  for (const [key, value] of params.entries()) {
    urlParams[key] = value;
  }
  return urlParams;
};

// update the shareable link URL with the current recipe and palette
const setLink = function () {
  const specEl = document.getElementById('spec');
  const linkEl = document.getElementById('link');
  const specs = specEl.value.split(/\s+/g).slice(0, 2);
  // strip any existing parameters
  let link = location.protocol + '//' + location.host + location.pathname;
  link += `?recipe_v2=${encodeURIComponent(specs[0])}`;
  if (PALETTE !== rwb_palette) {
    link += `&palette=${encodeURIComponent(PALETTE.reduce((x, y) => x + ' ' + y))}`;
  }
  linkEl.setAttribute('href', link);
};

// Drawing Functions
// ==================================================================================================

// init canvas element
// -------------------------------------------------------------------------------
const init = function () {
  const canvas = document.getElementById('poly');
  canvas.width = CANVAS_WIDTH;
  canvas.height = CANVAS_HEIGHT;

  ctx = canvas.getContext('2d');
  ctx.lineWidth = ctx_linewidth;

  if (BG_CLEAR) {
    ctx.clearRect(0, 0, CANVAS_WIDTH, CANVAS_HEIGHT);
  } else {
    ctx.clearRect(0, 0, CANVAS_WIDTH, CANVAS_HEIGHT);
    ctx.fillStyle = BG_COLOR;
    ctx.fillRect(0, 0, CANVAS_WIDTH, CANVAS_HEIGHT);
  }
};

// clear canvas
// -----------------------------------------------------------------------------------
const clear = function () {
  if (BG_CLEAR) {
    ctx.clearRect(0, 0, CANVAS_WIDTH, CANVAS_HEIGHT);
  } else {
    ctx.clearRect(0, 0, CANVAS_WIDTH, CANVAS_HEIGHT);
    ctx.fillStyle = BG_COLOR;
    ctx.fillRect(0, 0, CANVAS_WIDTH, CANVAS_HEIGHT);
  }
};

// main drawing routine for polyhedra
// ===================================================================================================
export const drawpoly = function (poly, tvec) {
  let v;
  if (!tvec) { tvec = [0, 0, 0]; }

  // rotate poly in 3d
  const oldxyz = poly.vertices.map(x => x);
  poly.vertices = poly.vertices.map(x => add(mv3(globRotM, x), tvec));

  try {
    // z sort faces
    sortFaces(poly, SORTING_ALGO);

    for (let fno = 0; fno < poly.faces.length; fno++) {
      // Back-face culling disabled for now (was leading to artifacts on certain polyhedra/orientations)
      // if (!poly.visibility[fno]) { continue; }
      var face = poly.faces[fno];
      ctx.beginPath();
      // move to first vertex of face
      const v0 = face[face.length - 1];
      let [x, y] = perspT(poly.vertices[v0], persp_z_max, persp_z_min, persp_ratio, perspective_scale);
      ctx.moveTo(x + _2d_x_offset, y + _2d_y_offset);
      // loop around face, defining polygon
      for (v of face) {
        [x, y] = perspT(poly.vertices[v], persp_z_max, persp_z_min, persp_ratio, perspective_scale);
        ctx.lineTo(x + _2d_x_offset, y + _2d_y_offset);
      }

      // use pre-computed colors
      let clr = palette(poly.face_classes[fno], PALETTE);

      // shade based on simple cosine illumination factor
      const face_verts = face.map((v) => poly.vertices[v]);
      // TODO: these magic illumination parameters should be global constants or parameters
      const illum = dot(normal(face_verts), unit([1, -1, 0]));
      clr = mult((((illum / 2.0) + 0.5) * 0.7) + 0.3, clr);

      if ((PaintMode === 'fill') || (PaintMode === 'fillstroke')) {
        ctx.fillStyle =
        `rgba(${round(clr[0] * 255)}, ${round(clr[1] * 255)}, ${round(clr[2] * 255)}, ${1.0})`;
        ctx.fill();
        // make cartoon stroke (=black) / realistic stroke an option (=below)
        ctx.strokeStyle =
        `rgba(${round(clr[0] * 255)}, ${round(clr[1] * 255)}, ${round(clr[2] * 255)}, ${1.0})`;
        ctx.stroke();
      }
      if (PaintMode === 'fillstroke') {
        ctx.fillStyle =
        `rgba(${round(clr[0] * 255)}, ${round(clr[1] * 255)}, ${round(clr[2] * 255)}, ${1.0})`;
        ctx.fill();
        ctx.strokeStyle = 'rgba(0,0,0,0.3)'; // light lines, less cartoony, more render-y
        ctx.stroke();
      }
      if (PaintMode === 'stroke') {
        ctx.strokeStyle = 'rgba(0,0,0,0.8)';
        ctx.stroke();
      }
    }

    ctx.save();

    ctx.fillStyle = 'red';
    if (window._pts) {
      for (let i = 0; i < window._pts.length; i++) {
        let [ x, y ] = window._pts[i][0];
        ctx.beginPath();
        ctx.arc(x + _2d_x_offset, y + _2d_y_offset, 3, 0, Math.PI * 2);
        ctx.fill();
        let d = Math.round(window._pts[i][3] * 1000) / 1000;
        ctx.fillText(`${window._pts[i][1]},${window._pts[i][2]}:${d}`, x + _2d_x_offset + 5, y + _2d_y_offset);
      }
    }
    ctx.lineWidth = 2.0;
    ctx.strokeStyle = 'black';
    ctx.fillStyle = 'rgba(255,0,0,0.5)';
    if (window._traces) {
      for (let i = 0; i < window._traces.length; i++) {
        let trace = window._traces[i];
        let [ x, y ] = trace[trace.length - 1];
        ctx.beginPath();
        ctx.moveTo(x + _2d_x_offset, y + _2d_y_offset);
        for (let j = 0; j < window._traces[i].length; j++) {
          // loop around face, defining polygon
          let [ x, y ] = trace[j];
          ctx.lineTo(x + _2d_x_offset, y + _2d_y_offset);
        }
      }
      ctx.stroke();
      ctx.fill();
    }
    // labelFaces(poly);
    ctx.restore();

  } catch (e) {
    console.log(e);
  } finally {
    // reset coords, for setting absolute rotation, as poly is passed by ref
    poly.vertices = oldxyz;
  }
};

const labelFaces = function (poly, show_hidden = false) {
  // Label all visible faces
  ctx.font = '10px Arial';
  ctx.fillStyle = 'black';
  let centers = poly.centers();
  let normals = poly.normals();
  for (let fno = 0; fno < poly.faces.length; fno++) {
    let [x, y] = perspT(centers[fno], persp_z_max, persp_z_min, persp_ratio, perspective_scale);
    if (normals[fno][2] < 0 && show_hidden) {
      ctx.fillStyle = 'gray';
      ctx.fillText(`${fno}`, x + _2d_x_offset, y + _2d_y_offset);
    } else {
      ctx.fillStyle = 'black';
      ctx.fillText(`${fno}`, x + _2d_x_offset, y + _2d_y_offset);
    }
  }
  // window._fvs = poly.faces.map(f => f.map(v => poly.vertices[v]));
};

// draw polyhedra just once
// -----------------------------------------------------------------------------------
const drawShape = function () {
  clear();
  globPolys.map((p, i) => drawpoly(p, [0 + (3 * i), 0, 3]));
};

// update V E F stats on page
// -----------------------------------------------------------------------------------
const updateStats = function () {
  const basicStatsEl = document.getElementById('basicstats');
  const moreStatsEl = document.getElementById('morestats');
  for (let i = 0; i < globPolys.length; i++) {
    const p = globPolys[i];
    basicStatsEl.innerHTML = p.data();
    moreStatsEl.innerHTML = p.moreData();
  }
};

// Initialization and Basic UI
// ===================================================================================================

const setupApp = function () {
  // Render markdown documentation
  const docsEl = document.getElementById('docs');
  if (docsEl) {
    docsEl.innerHTML = frontHtml;
  }

  init(); // init canvas

  const specEl = document.getElementById('spec');
  const paletteEl = document.getElementById('palette');
  const canvasEl = document.getElementById('poly');

  let specs;
  const urlParams = parseURL(); // see if recipe is spec'd in URL
  if ('recipe_v2' in urlParams) {
    specs = [urlParams['recipe_v2']];
    specEl.value = specs[0];
    setLink();
  } else if ('recipe' in urlParams) {
    specs = urlParams['recipe'].trim().split(/\s+/g).slice(0, 2).map(rewriteLegacyRecipe);
    specEl.value = specs.join(' ');
    setLink();
  } else {
    specs = [randomchoice(DEFAULT_RECIPES)];
    specEl.value = specs[0];
    setLink();
  }

  // set initial palette spec
  if ('palette' in urlParams) {
    PALETTE = urlParams['palette'].trim().split(/\s+/g);
    setLink();
  }
  paletteEl.value = PALETTE.reduce((x, y) => x + ' ' + y);

  // construct the polyhedra from spec
  globPolys = specs.map(x => generatePoly(x));
  updateStats();

  // draw it
  drawShape();

  // Event Handlers
  // ----------------------------------------------------

  // when spec changes in input, parse and draw new polyhedra
  specEl.addEventListener('change', function () {
    // only allow one recipe for now
    let specs = specEl.value.trim().split(/\s+/g).slice(0, 2);
    globPolys = specs.map(x => generatePoly(x));
    updateStats();
    setLink();
    drawShape();
  });

  // when palette changes in input, redraw polyhedra
  paletteEl.addEventListener('change', function () {
    PALETTE = paletteEl.value.trim().split(/\s+/g);
    setLink();
    drawShape();
  });

  // randomize palette
  document.getElementById('rndcolors').addEventListener('click', function () {
    let newpalette = rndcolors();
    PALETTE = newpalette;
    paletteEl.value = newpalette.join(' ');
    setLink();
    drawShape();
  });

  // Basic manipulation: rotation and scaling of geometry
  // ----------------------------------------------------

  // mousewheel changes scale of drawing
  canvasEl.addEventListener('wheel', function (e) {
    e.preventDefault();
    const delta = -Math.sign(e.deltaY);
    perspective_scale *= (10 + delta) / 10;
    drawShape();
  }, { passive: false });

  // Implement standard trackball routines
  // ---------------------------------------
  canvasEl.addEventListener('mousedown', function (e) {
    e.preventDefault();
    MOUSEDOWN = true;
    const rect = canvasEl.getBoundingClientRect();
    LastMouseX = e.clientX - rect.left;
    LastMouseY = e.clientY - rect.top;
    // calculate inverse projection of point to sphere
    const tmpvec = invperspT(LastMouseX, LastMouseY,
      _2d_x_offset, _2d_y_offset,
      persp_z_max, persp_z_min,
      persp_ratio, perspective_scale);
    // quick NaN check
    if ((tmpvec[0] * tmpvec[1] * tmpvec[2] * 0) === 0) {
      LastSphVec = tmpvec;
    }
    // copy last transform state
    quat.copy(globLastQuat, globQuat);
  });

  canvasEl.addEventListener('mouseup', function (e) {
    e.preventDefault();
    MOUSEDOWN = false;
  });

  canvasEl.addEventListener('mouseleave', function (e) {
    e.preventDefault();
    MOUSEDOWN = false;
  });

  canvasEl.addEventListener('mousemove', function (e) {
    e.preventDefault();
    if (MOUSEDOWN) {
      const rect = canvasEl.getBoundingClientRect();
      const MouseX = e.clientX - rect.left;
      const MouseY = e.clientY - rect.top;
      const SphVec = invperspT(MouseX, MouseY,
        _2d_x_offset, _2d_y_offset,
        persp_z_max, persp_z_min,
        persp_ratio, perspective_scale);

      // quick NaN check
      if (((SphVec[0] * SphVec[1] * SphVec[2] * 0) === 0) &&
           ((LastSphVec[0] * LastSphVec[1] * LastSphVec[2] * 0) === 0)) {
        const deltaQ = quat.fromUnitVectors(unit(LastSphVec), unit(SphVec));
        quat.multiply(deltaQ, globLastQuat, globQuat);
        quat.normalize(globQuat, globQuat);
        updateRotM();
      }
      drawShape();
    }
  });

  // State control via some buttons
  // ---------------------------------------
  document.getElementById('strokeonly').addEventListener('click', function () {
    PaintMode = 'stroke';
    drawShape();
  });

  document.getElementById('fillonly').addEventListener('click', function () {
    PaintMode = 'fill';
    drawShape();
  });

  document.getElementById('fillandstroke').addEventListener('click', function () {
    PaintMode = 'fillstroke';
    drawShape();
  });

  document.getElementById('siderot').addEventListener('click', function () {
    quat.setAxisAngle([0, 1, 0], PI / 2, globQuat);
    updateRotM();
    drawShape();
  });

  document.getElementById('toprot').addEventListener('click', function () {
    quat.setAxisAngle([1, 0, 0], PI / 2, globQuat);
    updateRotM();
    drawShape();
  });

  document.getElementById('frontrot').addEventListener('click', function () {
    quat.identity(globQuat);
    updateRotM();
    drawShape();
  });

  // Export Options
  // ---------------------------------------
  document.getElementById('pngsavebutton').addEventListener('click', function () {
    const spec = specEl.value.split(/\s+/g)[0];
    const filename = `polyhedronisme-${spec.replace(/\([^)]+\)/g, '')}.png`;
    const exportCanvas = document.createElement('canvas');
    exportCanvas.width = CANVAS_WIDTH * PNG_EXPORT_SCALE;
    exportCanvas.height = CANVAS_HEIGHT * PNG_EXPORT_SCALE;
    const exportCtx = exportCanvas.getContext('2d');
    exportCtx.scale(PNG_EXPORT_SCALE, PNG_EXPORT_SCALE);
    exportCtx.lineWidth = ctx_linewidth;

    const prevCtx = ctx;
    try {
      ctx = exportCtx;
      drawShape();
    } finally {
      ctx = prevCtx;
    }

    exportCanvas.toBlob(blob => {
      if (blob) downloadBlob(blob, filename);
    }, 'image/png');
  });

  document.getElementById('svgsavebutton').addEventListener('click', function () {
    const svgtxt = toSVG(globPolys[0], {
      rotM: globRotM,
      palette: PALETTE,
      paintMode: PaintMode,
      cullHidden: true,
      perspective_scale: perspective_scale,
      persp_z_max: persp_z_max,
      persp_z_min: persp_z_min,
      persp_ratio: persp_ratio
    });
    const spec = specEl.value.split(/\s+/g)[0];
    const filename = `polyhedronisme-${spec.replace(/\([^)]+\)/g, '')}.svg`;
    saveText(svgtxt, filename);
  });

  document.getElementById('objsavebutton').addEventListener('click', function () {
    const objtxt = toOBJ(globPolys[0]);
    const spec = specEl.value.split(/\s+/g)[0];
    const filename = `polyhedronisme-${spec.replace(/\([^)]+\)/g, '')}.obj`;
    saveText(objtxt, filename);
  });

  document.getElementById('x3dsavebutton').addEventListener('click', function () {
    const triangulated = triangulate(globPolys[0], true); // triangulate to preserve face_colors for 3d printing
    const x3dtxt = toVRML(triangulated, PALETTE);
    const spec = specEl.value.split(/\s+/g)[0];
    const filename = `polyhedronisme-${spec.replace(/\([^)]+\)/g, '')}.wrl`;
    saveText(x3dtxt, filename);
  });
};

if (document.readyState === 'loading') {
  document.addEventListener('DOMContentLoaded', setupApp);
} else {
  setupApp();
}
