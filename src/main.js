// Polyhédronisme
// ===================================================================================================
//
// A toy for constructing and manipulating polyhedra and other meshes
//
// Copyright 2019, Anselm Levskaya
// Released under the MIT License

/* eslint-disable camelcase */

// polyfills for saving png/obj/vrml
import 'file-saver';
import 'blobjs';
import 'canvas-toBlob';

import Vue from 'vue'
import { _ } from 'underscore';

import { add, clone, eye3, mv3, mm3, getVec2VecRotM, perspT, invperspT,
  vec_rotm, rotm, mult, dot, unit, normal, round, randomchoice, PI } from './geo';
import { rwb_palette, sortFaces, palette, rndcolors } from './polyhedron';
import { generatePoly } from './parser';
import { triangulate } from './topo_operators';
import { toOBJ, toVRML } from './exporters';
// import * as testing from './testing'
import { sortFacesCarefully } from './depth';

// docs
import frontText from './components/front.md'

// global css
require('./assets/polyhedronisme.css');

// debug
// window.testrig = testing.testrig;
Vue.config.productionTip = false;

// eslint-disable-next-line no-new
new Vue({
  el: '#docs',
  components: { frontText },
  template: '<front-text/>'
})

// GLOBALS
// ===================================================================================================
export let ctx = {}; // for global access to canvas context
export let globPolys = {}; // constructed polyhedras

const CANVAS_WIDTH = 720; // canvas dims
const CANVAS_HEIGHT = 500; // canvas dims
let globRotM = clone(eye3);
let globLastRotM = clone(eye3);
export let perspective_scale = 1200;
export const persp_z_max = 5;
export const persp_z_min = 0;
export const persp_ratio = 0.8;
const _2d_x_offset = CANVAS_WIDTH / 2;
const _2d_y_offset = CANVAS_HEIGHT / 2;

export let SORTING_ALGO = 'FAST'; // or 'CAREFUL';
window.set_algo = x => (SORTING_ALGO = x);

export const VISIBILITY_ZNORMAL_CUTOFF = -0.1;
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
  'C2dakD', 'oC20kkkT', 'kn4C40A0dA4', 'opD',
  'PT', 'qA40PA40oC', 'knD', 'dn6x4bT', 'oox4P7',
  'qqJ37', 'aobD', 'qaxI'];

// File-saving objects used to export txt/canvas-png
const saveText = function (text, filename) {
  const blb = new Blob([text],
    {type: `text/plain;charset=${document.characterSet}`});
  /* eslint-disable no-undef */
  saveAs(blb, filename);
}

// parses URL string for polyhedron recipe, for bookmarking
// should use #! href format instead
const parseURL = function () {
  let e;
  const urlParams = {};
  const a = /\+/g; // Regex for replacing addition symbol with a space
  const r = /([^&=]+)=?([^&]*)/g;
  const d = s => decodeURIComponent(s.replace(a, ' '));
  const q = window.location.search.substring(1);

  while ((e = r.exec(q))) {
    urlParams[d(e[1])] = d(e[2]);
  }
  return urlParams;
};

// update the shareable link URL with the current recipe and palette
const setLink = function () {
  const specs = $('#spec').val().split(/\s+/g).slice(0, 2);
  // strip any existing parameters
  let link = location.protocol + '//' + location.host + location.pathname;
  link += `?recipe=${encodeURIComponent(specs[0])}`;
  if (PALETTE !== rwb_palette) {
    link += `&palette=${encodeURIComponent(PALETTE.reduce((x, y) => x + ' ' + y))}`;
  }
  $('#link').attr('href', link);
};

// Drawing Functions
// ==================================================================================================

// init canvas element
// -------------------------------------------------------------------------------
const init = function () {
  const canvas = $('#poly');
  canvas.width(CANVAS_WIDTH);
  canvas.height(CANVAS_HEIGHT);

  ctx = canvas[0].getContext('2d');
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
  const oldxyz = _.map(poly.vertices, x => x);
  poly.vertices = _.map(poly.vertices, x => add(mv3(globRotM, x), tvec));

  try {
    // z sort faces
    sortFaces(poly, SORTING_ALGO);

    for (let fno = 0; fno < poly.faces.length; fno++) {
      if (!poly.visibility[fno]) { continue; }
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
      let clr = palette(poly.face_classes[fno]);

      // shade based on simple cosine illumination factor
      const face_verts = face.map((v) => poly.vertices[v])
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
        let d = Math.round(window._pts[i][3]*1000)/1000;
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
    console.log(e)
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
}

// draw polyhedra just once
// -----------------------------------------------------------------------------------
const drawShape = function () {
  clear();
  globPolys.map((p, i) => drawpoly(p, [0 + (3 * i), 0, 3]));
};

// update V E F stats on page
// -----------------------------------------------------------------------------------
const updateStats = function () {
  for (let i = 0; i < globPolys.length; i++) {
    const p = globPolys[i];
    $('#basicstats').html(p.data());
    $('#morestats').html(p.moreData());
  }
}

// Initialization and Basic UI
// ===================================================================================================

$(function () { // wait for page to load
  init(); // init canvas

  let specs;
  const urlParams = parseURL(); // see if recipe is spec'd in URL
  if ('recipe' in urlParams) {
    specs = [urlParams['recipe']];
    $('#spec').val(specs);
  } else {
    specs = [randomchoice(DEFAULT_RECIPES)];
    $('#spec').val(specs);
    setLink();
  }

  // set initial palette spec
  if ('palette' in urlParams) {
    PALETTE = urlParams['palette'].split(/\s+/g);
    setLink();
  }
  $('#palette').val(PALETTE.reduce((x, y) => x + ' ' + y));

  // construct the polyhedra from spec
  globPolys = _.map(specs, x => generatePoly(x));
  updateStats();

  // draw it
  drawShape();

  // Event Handlers
  // ----------------------------------------------------

  // when spec changes in input, parse and draw new polyhedra
  $('#spec').change(function (e) {
    // only allow one recipe for now
    let specs = $('#spec').val().split(/\s+/g).slice(0, 2);
    globPolys = _.map(specs, x => generatePoly(x));
    updateStats();
    // animateShape()
    setLink();
    drawShape();
  });

  // when palette changes in input, redraw polyhedra
  $('#palette').change(function (e) {
    PALETTE = $(this).val().split(/\s+/g);
    setLink();
    drawShape();
  });

  // randomize palette
  $('#rndcolors').click(function (e) {
    let newpalette = rndcolors();
    PALETTE = newpalette;
    $('#palette').val(newpalette.join(' '));
    setLink();
    drawShape();
  });

  // Basic manipulation: rotation and scaling of geometry
  // ----------------------------------------------------

  // mousewheel changes scale of drawing
  $('#poly').mousewheel(function (e, delta, deltaX, deltaY) {
    e.preventDefault();
    perspective_scale *= (10 + delta) / 10;
    drawShape();
  });

  // Implement standard trackball routines
  // ---------------------------------------
  $('#poly').mousedown(function (e) {
    e.preventDefault();
    MOUSEDOWN = true;
    // relative mouse coords
    LastMouseX = e.clientX - $(this).offset().left;
    LastMouseY = e.clientY - ($(this).offset().top - $(window).scrollTop());
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
    globLastRotM = clone(globRotM);
  });
  $('#poly').mouseup(function (e) {
    e.preventDefault();
    MOUSEDOWN = false;
  });
  $('#poly').mouseleave(function (e) {
    e.preventDefault();
    MOUSEDOWN = false;
  });
  $('#poly').mousemove(function (e) {
    e.preventDefault();
    if (MOUSEDOWN) {
      const MouseX = e.clientX - $(this).offset().left;
      const MouseY = e.clientY - ($(this).offset().top - $(window).scrollTop());
      const SphVec = invperspT(MouseX, MouseY,
        _2d_x_offset, _2d_y_offset,
        persp_z_max, persp_z_min,
        persp_ratio, perspective_scale);

      // quick NaN check
      if (((SphVec[0] * SphVec[1] * SphVec[2] * 0) === 0) &&
           ((LastSphVec[0] * LastSphVec[1] * LastSphVec[2] * 0) === 0)) {
        globRotM = mm3(getVec2VecRotM(LastSphVec, SphVec), globLastRotM);
      }
      drawShape();
    }
  });

  // State control via some buttons
  // ---------------------------------------
  $('#strokeonly').click(function (e) {
    PaintMode = 'stroke';
    drawShape();
  });

  $('#fillonly').click(function (e) {
    PaintMode = 'fill';
    drawShape();
  });

  $('#fillandstroke').click(function (e) {
    PaintMode = 'fillstroke';
    drawShape();
  });

  $('#siderot').click(function (e) {
    globRotM = vec_rotm(PI / 2, 0, 1, 0);
    drawShape();
  });

  $('#toprot').click(function (e) {
    globRotM = vec_rotm(PI / 2, 1, 0, 0);
    drawShape();
  });

  $('#frontrot').click(function (e) {
    globRotM = rotm(0, 0, 0);
    drawShape();
  });

  // Export Options
  // ---------------------------------------
  $('#pngsavebutton').click(function (e) {
    const canvas = $('#poly')[0];
    const spec = $('#spec').val().split(/\s+/g)[0];
    const filename = `polyhedronisme-${spec.replace(/\([^)]+\)/g, '')}.png`;
    canvas.toBlobHD(blob => saveAs(blob, filename));
  });

  $('#objsavebutton').click(function (e) {
    const objtxt = toOBJ(globPolys[0]);
    const spec = $('#spec').val().split(/\s+/g)[0];
    const filename = `polyhedronisme-${spec.replace(/\([^)]+\)/g, '')}.obj`;
    saveText(objtxt, filename);
  });

  $('#x3dsavebutton').click(function (e) {
    const triangulated = triangulate(globPolys[0], true); // triangulate to preserve face_colors for 3d printing
    const x3dtxt = toVRML(triangulated);
    const spec = $('#spec').val().split(/\s+/g)[0];
    const filename = `polyhedronisme-${spec.replace(/\([^)]+\)/g, '')}.wrl`;
    saveText(x3dtxt, filename);
  });
});