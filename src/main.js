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
import { rwb_palette, PRESET_PALETTES, sortFaces, palette, rndcolors,
  shufflePalette, normalizeHex } from './polyhedron.js';
import { PRESET_RECIPE_GROUPS, DEFAULT_RECIPES, generateRandomRecipe } from './recipes.js';
import { generatePoly, setOpParser } from './parser.js';
import defaultOpParser from './original.pegjs';
import { rewriteLegacyRecipe, setLegacyOpParser } from './legacy_parser.js';
import legacyParser from './legacy.pegjs';
import { triangulate } from './topo_operators.js';
import { toOBJ, toVRML, toSVG, toSTL, toGLB } from './exporters.js';
import { sortFacesCarefully } from './depth.js';

// docs
import frontHtml from './components/front.md';

// global css
import './assets/polyhedronisme.css';

// Initialize recipe parsers
setOpParser(defaultOpParser);
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
const isDefaultPalette = function (pal) {
  if (!Array.isArray(pal) || pal.length !== rwb_palette.length) return false;
  return pal.every((c, i) => normalizeHex(c) === normalizeHex(rwb_palette[i]));
};

const setLink = function () {
  const specEl = document.getElementById('spec');
  const linkEl = document.getElementById('link');
  const specs = specEl.value.split(/\s+/g).slice(0, 2);
  // strip any existing parameters
  let link = location.protocol + '//' + location.host + location.pathname;
  link += `?recipe_v2=${encodeURIComponent(specs[0])}`;
  if (!isDefaultPalette(PALETTE)) {
    link += `&palette=${encodeURIComponent(PALETTE.join(' '))}`;
  }
  linkEl.setAttribute('href', link);
};

// Drawing Functions
// ==================================================================================================

// init canvas element
// -------------------------------------------------------------------------------
const init = function () {
  const canvas = document.getElementById('poly');
  const dpr = window.devicePixelRatio || 1;
  canvas.width = CANVAS_WIDTH * dpr;
  canvas.height = CANVAS_HEIGHT * dpr;
  canvas.style.width = `${CANVAS_WIDTH}px`;
  canvas.style.height = `${CANVAS_HEIGHT}px`;

  ctx = canvas.getContext('2d');
  ctx.scale(dpr, dpr);
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
  for (let i = 0; i < globPolys.length; i++) {
    const p = globPolys[i];
    basicStatsEl.innerHTML = p.data();
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
  const recipePresetEl = document.getElementById('recipe-preset');
  const examplePolyBtn = document.getElementById('examplepoly');
  const rndPolyBtn = document.getElementById('rndpoly');
  const paletteEl = document.getElementById('palette');
  const presetEl = document.getElementById('palette-preset');
  const swatchesEl = document.getElementById('palette-swatches');
  const addColorBtn = document.getElementById('add-color');
  const removeColorBtn = document.getElementById('remove-color');
  const shuffleBtn = document.getElementById('shufflecolors');
  const rndBtn = document.getElementById('rndcolors');
  const toggleHexBtn = document.getElementById('togglehex');
  const hexRowEl = document.getElementById('palette-hex-row');
  const canvasEl = document.getElementById('poly');

  // Populate curated recipe preset selector options
  if (recipePresetEl) {
    for (const group of PRESET_RECIPE_GROUPS) {
      const optgroup = document.createElement('optgroup');
      optgroup.label = group.label;
      for (const r of group.recipes) {
        const opt = document.createElement('option');
        opt.value = r.spec;
        opt.textContent = r.name;
        optgroup.appendChild(opt);
      }
      recipePresetEl.appendChild(optgroup);
    }
    const customRecipeOpt = document.createElement('option');
    customRecipeOpt.value = 'custom';
    customRecipeOpt.textContent = 'Custom / Procedural';
    customRecipeOpt.disabled = true;
    customRecipeOpt.hidden = true;
    recipePresetEl.appendChild(customRecipeOpt);
  }

  // Sync recipe preset dropdown to match current spec
  const syncRecipePresetSelect = function (currentSpec) {
    if (!recipePresetEl) return;
    const trimmed = (currentSpec || '').trim();
    recipePresetEl.value = DEFAULT_RECIPES.includes(trimmed) ? trimmed : 'custom';
  };

  // Populate preset selector options
  if (presetEl) {
    for (const preset of PRESET_PALETTES) {
      const opt = document.createElement('option');
      opt.value = preset.id;
      opt.textContent = preset.name;
      presetEl.appendChild(opt);
    }
    const customOpt = document.createElement('option');
    customOpt.value = 'custom';
    customOpt.textContent = 'Custom';
    customOpt.disabled = true;
    customOpt.hidden = true;
    presetEl.appendChild(customOpt);
  }

  // Compute which palette indices are actively used by the current polyhedron(s)
  const getUsedPaletteIndices = function () {
    const used = new Set();
    if (Array.isArray(globPolys) && PALETTE.length > 0) {
      for (const poly of globPolys) {
        if (poly && Array.isArray(poly.face_classes)) {
          for (const fc of poly.face_classes) {
            used.add(((fc % PALETTE.length) + PALETTE.length) % PALETTE.length);
          }
        }
      }
    }
    return used;
  };

  // Sync preset dropdown to match current PALETTE
  const syncPresetSelect = function () {
    if (!presetEl) return;
    const normCurrent = PALETTE.map(c => normalizeHex(c));
    const matched = PRESET_PALETTES.find(p =>
      p.colors.length === normCurrent.length &&
      p.colors.every((c, i) => normalizeHex(c) === normCurrent[i])
    );
    presetEl.value = matched ? matched.id : 'custom';
  };

  // Sync raw hex input box and +/- button states
  const syncPaletteAuxUI = function () {
    if (paletteEl && document.activeElement !== paletteEl) {
      paletteEl.value = PALETTE.join(' ');
    }
    syncPresetSelect();
    if (removeColorBtn) removeColorBtn.disabled = PALETTE.length <= 2;
    if (addColorBtn) addColorBtn.disabled = PALETTE.length >= 16;
  };

  // Update visual styling of existing swatch nodes without rebuilding DOM
  const refreshSwatchVisuals = function () {
    if (!swatchesEl) return;
    const used = getUsedPaletteIndices();
    const items = swatchesEl.querySelectorAll('.swatch-item');
    items.forEach((item, idx) => {
      const hex = normalizeHex(PALETTE[idx]);
      const isUsed = used.size === 0 || used.has(idx);
      item.style.backgroundColor = hex;
      item.classList.toggle('swatch-unused', !isUsed);
      item.title = isUsed
        ? `Color #${idx + 1} (${hex}) — click to change, drag to reorder`
        : `Color #${idx + 1} (${hex}, unused on this shape) — drag left to use, or click to edit`;
      const input = item.querySelector('.swatch-input');
      if (input && input.value !== hex) {
        input.value = hex;
      }
    });
  };

  // Render interactive color swatch squares with click-to-pick and drag-to-reorder
  const renderSwatches = function () {
    if (!swatchesEl) return;
    swatchesEl.innerHTML = '';

    PALETTE.forEach((color, initialIdx) => {
      const hex = normalizeHex(color);
      PALETTE[initialIdx] = hex;

      const item = document.createElement('div');
      item.className = 'swatch-item';

      const input = document.createElement('input');
      input.type = 'color';
      input.className = 'swatch-input';
      input.value = hex;

      let suppressNextClick = false;

      // Live color picking
      input.addEventListener('input', function () {
        const currentIdx = Array.from(swatchesEl.children).indexOf(item);
        if (currentIdx === -1) return;
        PALETTE[currentIdx] = normalizeHex(input.value);
        refreshSwatchVisuals();
        syncPaletteAuxUI();
        setLink();
        drawShape();
      });

      // Suppress color picker popup if pointer interaction was a drag-reorder
      input.addEventListener('click', function (e) {
        if (suppressNextClick) {
          e.preventDefault();
          e.stopPropagation();
          suppressNextClick = false;
        }
      }, true);

      // Pointer-based drag reordering across the swatch strip
      item.addEventListener('pointerdown', function (e) {
        if (e.button !== 0) return;
        const startX = e.clientX;
        const startY = e.clientY;
        let didDrag = false;

        const onPointerMove = function (moveEvt) {
          const dist = Math.hypot(moveEvt.clientX - startX, moveEvt.clientY - startY);
          if (!didDrag && dist > 4) {
            didDrag = true;
            suppressNextClick = true;
            item.classList.add('swatch-dragging');
          }
          if (!didDrag) return;

          const children = Array.from(swatchesEl.children);
          const fromIdx = children.indexOf(item);
          if (fromIdx === -1) return;

          for (let toIdx = 0; toIdx < children.length; toIdx++) {
            if (toIdx === fromIdx) continue;
            const rect = children[toIdx].getBoundingClientRect();
            const midX = rect.left + rect.width / 2;
            if ((toIdx < fromIdx && moveEvt.clientX < midX) ||
                (toIdx > fromIdx && moveEvt.clientX > midX)) {
              // Move color in PALETTE and reorder DOM node in-place
              const [movedColor] = PALETTE.splice(fromIdx, 1);
              PALETTE.splice(toIdx, 0, movedColor);
              if (toIdx < fromIdx) {
                swatchesEl.insertBefore(item, children[toIdx]);
              } else {
                swatchesEl.insertBefore(item, children[toIdx].nextSibling);
              }
              refreshSwatchVisuals();
              syncPaletteAuxUI();
              setLink();
              drawShape();
              break;
            }
          }
        };

        const onPointerUp = function () {
          window.removeEventListener('pointermove', onPointerMove);
          window.removeEventListener('pointerup', onPointerUp);
          item.classList.remove('swatch-dragging');
          if (didDrag) {
            setTimeout(() => { suppressNextClick = false; }, 50);
          }
        };

        window.addEventListener('pointermove', onPointerMove);
        window.addEventListener('pointerup', onPointerUp);
      });

      item.appendChild(input);
      swatchesEl.appendChild(item);
    });

    refreshSwatchVisuals();
    syncPaletteAuxUI();
  };

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

  // set initial palette spec (defaulting to a fresh copy of rwb_palette)
  if ('palette' in urlParams) {
    const parsedPal = urlParams['palette'].trim().split(/[\s,]+/g).filter(Boolean).map(c => normalizeHex(c));
    PALETTE = parsedPal.length >= 2 ? parsedPal : [...rwb_palette];
    setLink();
  } else {
    PALETTE = [...rwb_palette];
  }

  // construct the polyhedra from spec
  globPolys = specs.map(x => generatePoly(x));
  syncRecipePresetSelect(specs[0]);
  updateStats();
  renderSwatches();

  // draw it
  drawShape();

  // Apply a recipe spec string, rebuild polyhedron, and sync UI
  const applyRecipeSpec = function (rawSpec) {
    const trimmed = (rawSpec || '').trim();
    if (!trimmed) return;
    specEl.value = trimmed;
    const nextSpecs = trimmed.split(/\s+/g).slice(0, 2);
    globPolys = nextSpecs.map(x => generatePoly(x));
    syncRecipePresetSelect(nextSpecs[0]);
    updateStats();
    refreshSwatchVisuals();
    setLink();
    drawShape();
  };

  // Event Handlers
  // ----------------------------------------------------

  // when spec changes in input, parse and draw new polyhedra
  specEl.addEventListener('change', function () {
    applyRecipeSpec(specEl.value);
  });

  // curated recipe preset selection
  if (recipePresetEl) {
    recipePresetEl.addEventListener('change', function () {
      if (recipePresetEl.value && recipePresetEl.value !== 'custom') {
        applyRecipeSpec(recipePresetEl.value);
      }
    });
  }

  // pick a random curated example
  if (examplePolyBtn) {
    examplePolyBtn.addEventListener('click', function () {
      const current = specEl.value.trim();
      const pool = DEFAULT_RECIPES.filter(r => r !== current);
      const next = randomchoice(pool.length > 0 ? pool : DEFAULT_RECIPES);
      applyRecipeSpec(next);
    });
  }

  // synthesize a procedural random recipe
  if (rndPolyBtn) {
    rndPolyBtn.addEventListener('click', function () {
      const next = generateRandomRecipe(specEl.value.trim());
      applyRecipeSpec(next);
    });
  }

  // preset selection
  if (presetEl) {
    presetEl.addEventListener('change', function () {
      const preset = PRESET_PALETTES.find(p => p.id === presetEl.value);
      if (preset) {
        PALETTE = [...preset.colors];
        renderSwatches();
        setLink();
        drawShape();
      }
    });
  }

  // add/remove color slots
  if (addColorBtn) {
    addColorBtn.addEventListener('click', function () {
      if (PALETTE.length >= 16) return;
      const extra = rndcolors(1)[0];
      PALETTE = [...PALETTE, extra];
      renderSwatches();
      setLink();
      drawShape();
    });
  }

  if (removeColorBtn) {
    removeColorBtn.addEventListener('click', function () {
      if (PALETTE.length <= 2) return;
      PALETTE = PALETTE.slice(0, -1);
      renderSwatches();
      setLink();
      drawShape();
    });
  }

  // shuffle current palette order
  if (shuffleBtn) {
    shuffleBtn.addEventListener('click', function () {
      PALETTE = shufflePalette(PALETTE);
      renderSwatches();
      setLink();
      drawShape();
    });
  }

  // randomize palette using color-theory generator
  if (rndBtn) {
    rndBtn.addEventListener('click', function () {
      const targetLen = Math.max(6, Math.min(12, PALETTE.length || 8));
      PALETTE = rndcolors(targetLen);
      renderSwatches();
      setLink();
      drawShape();
    });
  }

  // toggle raw hex input row
  if (toggleHexBtn && hexRowEl) {
    toggleHexBtn.addEventListener('click', function () {
      const willShow = hexRowEl.hasAttribute('hidden');
      if (willShow) {
        hexRowEl.removeAttribute('hidden');
        toggleHexBtn.classList.add('active');
        paletteEl.value = PALETTE.join(' ');
        paletteEl.focus();
      } else {
        hexRowEl.setAttribute('hidden', '');
        toggleHexBtn.classList.remove('active');
      }
    });
  }

  // when palette changes in raw hex input, sync swatches and redraw polyhedra
  const applyRawHexInput = function () {
    const tokens = paletteEl.value.trim().split(/[\s,]+/g).filter(Boolean);
    const valid = tokens.filter(t => /^#?([0-9a-fA-F]{3}|[0-9a-fA-F]{6})$/.test(t)).map(t => normalizeHex(t));
    if (valid.length >= 1) {
      PALETTE = valid;
      renderSwatches();
      setLink();
      drawShape();
    }
  };

  paletteEl.addEventListener('change', applyRawHexInput);
  paletteEl.addEventListener('input', applyRawHexInput);

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
  document.getElementById('fitscreen').addEventListener('click', function () {
    if (!Array.isArray(globPolys) || globPolys.length === 0) return;
    let maxUnitX = 0;
    let maxUnitY = 0;
    for (let i = 0; i < globPolys.length; i++) {
      const poly = globPolys[i];
      if (!poly || !Array.isArray(poly.vertices)) continue;
      const tvec = [3 * i, 0, 3];
      for (const v of poly.vertices) {
        const tv = add(mv3(globRotM, v), tvec);
        const [ux, uy] = perspT(tv, persp_z_max, persp_z_min, persp_ratio, 1);
        const ax = Math.abs(ux);
        const ay = Math.abs(uy);
        if (ax > maxUnitX) maxUnitX = ax;
        if (ay > maxUnitY) maxUnitY = ay;
      }
    }
    if (maxUnitX > 1e-9 && maxUnitY > 1e-9) {
      const margin = 0.9;
      const scaleX = (_2d_x_offset * margin) / maxUnitX;
      const scaleY = (_2d_y_offset * margin) / maxUnitY;
      perspective_scale = Math.min(scaleX, scaleY);
      drawShape();
    }
  });

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

  document.getElementById('stlsavebutton').addEventListener('click', function () {
    const stlbuf = toSTL(globPolys[0]);
    const spec = specEl.value.split(/\s+/g)[0];
    const filename = `polyhedronisme-${spec.replace(/\([^)]+\)/g, '')}.stl`;
    downloadBlob(new Blob([stlbuf], { type: 'model/stl' }), filename);
  });

  document.getElementById('glbsavebutton').addEventListener('click', function () {
    const glbbuf = toGLB(globPolys[0], PALETTE);
    const spec = specEl.value.split(/\s+/g)[0];
    const filename = `polyhedronisme-${spec.replace(/\([^)]+\)/g, '')}.glb`;
    downloadBlob(new Blob([glbbuf], { type: 'model/gltf-binary' }), filename);
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
