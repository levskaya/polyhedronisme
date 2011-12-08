(function() {
  var BG_CLEAR, BG_COLOR, CANVAS_HEIGHT, CANVAS_WIDTH, COLOR_METHOD, CSGToPoly, DEFAULT_RECIPES, LastMouseX, LastMouseY, LastSphVec, MOUSEDOWN, PALETTE, PI, PaintMode, abs, acos, add, adjustXYZ, ambo, animateShape, antiprism, arrayToCSGVertex, asin, atan, calcCentroid, canonicalXYZ, canonicalize, clear, clone, colorassign, convexarea, copyVecArray, cos, cross, csgIntersect, csgSubtract, csgUnion, ctx, ctx_linewidth, cube, def_palette, diagsToTris, dodecahedron, dot, drawShape, drawpoly, dual, edgeDist, edge_vert_isect, edgesToFace, enumerate, extrudeN, eye3, faceToEdges, faceprint, floor, generatePoly, getDiagonals, getOps, getVec2VecRotM, globPolys, globRotM, globlastRotM, globtime, gyro, hextofloats, icosahedron, init, insetN, intersect, invperspT, joinFaces, kisN, mag, mag2, meshFix, midpoint, mm3, mult, mv3, normal, octahedron, old_centers, old_normals, oneThird, orthogonal, paintPolyhedron, palette, parseurl, perspT, persp_ratio, persp_z_max, persp_z_min, perspective_scale, planarize, polyToCSG, polyflag, polyhedron, pow, prism, project2dface, propellor, pyramid, random, randomchoice, recenter, reciprocal, reciprocalC, reciprocalN, reflect, rescale, rotm, round, rwb_palette, rwbg_palette, saveText, sin, sortfaces, specreplacements, sqrt, stellaN, sub, tan, tangentPoint, tangentify, testrig, tetrahedron, topolog, triEq, triangulate, tween, uniquifyedges, uniquifyverts, unit, uniteFaces, vec_rotm, vertColors, vertIndex, _2d_x_offset, _2d_y_offset, _mult;
  var __indexOf = Array.prototype.indexOf || function(item) {
    for (var i = 0, l = this.length; i < l; i++) {
      if (this[i] === item) return i;
    }
    return -1;
  };
  random = Math.random;
  round = Math.round;
  floor = Math.floor;
  sqrt = Math.sqrt;
  sin = Math.sin;
  cos = Math.cos;
  tan = Math.tan;
  asin = Math.asin;
  acos = Math.acos;
  atan = Math.atan;
  pow = Math.pow;
  abs = Math.abs;
  PI = Math.PI;
  enumerate = function(ar) {
    var i, _ref, _results;
    _results = [];
    for (i = 0, _ref = ar.length - 1; 0 <= _ref ? i <= _ref : i >= _ref; 0 <= _ref ? i++ : i--) {
      _results.push([i, ar[i]]);
    }
    return _results;
  };
  clone = function(obj) {
    var key, newInstance;
    if (!(obj != null) || typeof obj !== 'object') {
      return obj;
    }
    newInstance = new obj.constructor();
    for (key in obj) {
      newInstance[key] = clone(obj[key]);
    }
    return newInstance;
  };
  randomchoice = function(array) {
    var n;
    n = floor(random() * array.length);
    return array[n];
  };
  mult = function(c, vec) {
    return [c * vec[0], c * vec[1], c * vec[2]];
  };
  _mult = function(vec1, vec2) {
    return [vec1[0] * vec2[0], vec1[1] * vec2[1], vec1[2] * vec2[2]];
  };
  add = function(vec1, vec2) {
    return [vec1[0] + vec2[0], vec1[1] + vec2[1], vec1[2] + vec2[2]];
  };
  sub = function(vec1, vec2) {
    return [vec1[0] - vec2[0], vec1[1] - vec2[1], vec1[2] - vec2[2]];
  };
  dot = function(vec1, vec2) {
    return vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];
  };
  cross = function(d1, d2) {
    return [d1[1] * d2[2] - d1[2] * d2[1], d1[2] * d2[0] - d1[0] * d2[2], d1[0] * d2[1] - d1[1] * d2[0]];
  };
  mag = function(vec) {
    return sqrt(dot(vec, vec));
  };
  mag2 = function(vec) {
    return dot(vec, vec);
  };
  unit = function(vec) {
    return mult(1 / sqrt(mag2(vec)), vec);
  };
  midpoint = function(vec1, vec2) {
    return mult(1 / 2.0, add(vec1, vec2));
  };
  tween = function(vec1, vec2, t) {
    return [(1 - t) * vec1[0] + t * vec2[0], (1 - t) * vec1[1] + t * vec2[1], (1 - t) * vec1[2] + t * vec2[2]];
  };
  oneThird = function(vec1, vec2) {
    return tween(vec1, vec2, 1 / 3.0);
  };
  reciprocal = function(vec) {
    return mult(1.0 / mag2(vec), vec);
  };
  tangentPoint = function(v1, v2) {
    var d;
    d = sub(v2, v1);
    return sub(v1, mult(dot(d, v1) / mag2(d), d));
  };
  edgeDist = function(v1, v2) {
    return sqrt(mag2(tangentPoint(v1, v2)));
  };
  orthogonal = function(v1, v2, v3) {
    var d1, d2;
    d1 = sub(v2, v1);
    d2 = sub(v3, v2);
    return cross(d1, d2);
  };
  intersect = function(set1, set2, set3) {
    var s1, s2, s3, _i, _j, _k, _len, _len2, _len3;
    for (_i = 0, _len = set1.length; _i < _len; _i++) {
      s1 = set1[_i];
      for (_j = 0, _len2 = set2.length; _j < _len2; _j++) {
        s2 = set2[_j];
        if (s1 === s2) {
          for (_k = 0, _len3 = set3.length; _k < _len3; _k++) {
            s3 = set3[_k];
            if (s1 === s3) {
              return s1;
            }
          }
        }
      }
    }
    return null;
  };
  calcCentroid = function(xyzs) {
    var centroidV, v, _i, _len;
    centroidV = [0, 0, 0];
    for (_i = 0, _len = xyzs.length; _i < _len; _i++) {
      v = xyzs[_i];
      centroidV = add(centroidV, v);
    }
    return mult(1 / xyzs.length, centroidV);
  };
  normal = function(xyzs) {
    var normalV, v1, v2, v3, _i, _len, _ref, _ref2;
    normalV = [0, 0, 0];
    _ref = xyzs.slice(-2), v1 = _ref[0], v2 = _ref[1];
    for (_i = 0, _len = xyzs.length; _i < _len; _i++) {
      v3 = xyzs[_i];
      normalV = add(normalV, orthogonal(v1, v2, v3));
      _ref2 = [v2, v3], v1 = _ref2[0], v2 = _ref2[1];
    }
    return unit(normalV);
  };
  convexarea = function(xyzs) {
    var area, v1, v2, v3, _i, _len, _ref, _ref2;
    area = 0.0;
    _ref = xyzs.slice(0, 2), v1 = _ref[0], v2 = _ref[1];
    _ref2 = xyzs.slice(2);
    for (_i = 0, _len = _ref2.length; _i < _len; _i++) {
      v3 = _ref2[_i];
      area += mag(cross(sub(v2, v1), sub(v3, v1)));
      v2 = v3;
    }
    return area;
  };
  project2dface = function(verts) {
    var c, n, p, tmpverts, v, v0, _i, _len, _results;
    tmpverts = clone(verts);
    v0 = verts[0];
    tmpverts = _.map(tmpverts, function(x) {
      return x - v0;
    });
    n = normal(verts);
    c = unit(calcCentroid(verts));
    p = cross(n, c);
    _results = [];
    for (_i = 0, _len = tmpverts.length; _i < _len; _i++) {
      v = tmpverts[_i];
      _results.push([dot(n, v), dot(p, v)]);
    }
    return _results;
  };
  copyVecArray = function(vecArray) {
    var i, newVecArray, _ref;
    newVecArray = new Array(vecArray.length);
    for (i = 0, _ref = vecArray.length - 1; 0 <= _ref ? i <= _ref : i >= _ref; 0 <= _ref ? i++ : i--) {
      newVecArray[i] = vecArray[i].slice(0);
    }
    return newVecArray;
  };
  mv3 = function(mat, vec) {
    return [mat[0][0] * vec[0] + mat[0][1] * vec[1] + mat[0][2] * vec[2], mat[1][0] * vec[0] + mat[1][1] * vec[1] + mat[1][2] * vec[2], mat[2][0] * vec[0] + mat[2][1] * vec[1] + mat[2][2] * vec[2]];
  };
  mm3 = function(A, B) {
    return [[A[0][0] * B[0][0] + A[0][1] * B[1][0] + A[0][2] * B[2][0], A[0][0] * B[0][1] + A[0][1] * B[1][1] + A[0][2] * B[2][1], A[0][0] * B[0][2] + A[0][1] * B[1][2] + A[0][2] * B[2][2]], [A[1][0] * B[0][0] + A[1][1] * B[1][0] + A[1][2] * B[2][0], A[1][0] * B[0][1] + A[1][1] * B[1][1] + A[1][2] * B[2][1], A[1][0] * B[0][2] + A[1][1] * B[1][2] + A[1][2] * B[2][2]], [A[2][0] * B[0][0] + A[2][1] * B[1][0] + A[2][2] * B[2][0], A[2][0] * B[0][1] + A[2][1] * B[1][1] + A[2][2] * B[2][1], A[2][0] * B[0][2] + A[2][1] * B[1][2] + A[2][2] * B[2][2]]];
  };
  eye3 = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
  rotm = function(phi, theta, psi) {
    var xy_mat, xz_mat, yz_mat;
    xy_mat = [[cos(phi), -1.0 * sin(phi), 0.0], [sin(phi), cos(phi), 0.0], [0.0, 0.0, 1.0]];
    yz_mat = [[cos(theta), 0, -1.0 * sin(theta)], [0, 1, 0], [sin(theta), 0, cos(theta)]];
    xz_mat = [[1.0, 0, 0], [0, cos(psi), -1.0 * sin(psi)], [0, sin(psi), cos(psi)]];
    return mm3(xz_mat, mm3(yz_mat, xy_mat));
  };
  vec_rotm = function(angle, x, y, z) {
    var cosA, length, m, sinA, sinA2, x2, y2, z2, _ref, _ref2;
    angle /= 2;
    sinA = sin(angle);
    cosA = cos(angle);
    sinA2 = sinA * sinA;
    length = mag([x, y, z]);
    if (length === 0) {
      _ref = [0, 0, 1], x = _ref[0], y = _ref[1], z = _ref[2];
    }
    if (length !== 1) {
      _ref2 = unit([x, y, z]), x = _ref2[0], y = _ref2[1], z = _ref2[2];
    }
    if (x === 1 && y === 0 && z === 0) {
      m = [[1, 0, 0], [0, 1 - 2 * sinA2, 2 * sinA * cosA], [0, -2 * sinA * cosA, 1 - 2 * sinA2]];
    } else if (x === 0 && y === 1 && z === 0) {
      m = [[1 - 2 * sinA2, 0, -2 * sinA * cosA], [0, 1, 0], [2 * sinA * cosA, 0, 1 - 2 * sinA2]];
    } else if (x === 0 && y === 0 && z === 1) {
      m = [[1 - 2 * sinA2, 2 * sinA * cosA, 0], [-2 * sinA * cosA, 1 - 2 * sinA2, 0], [0, 0, 1]];
    } else {
      x2 = x * x;
      y2 = y * y;
      z2 = z * z;
      m = [[1 - 2 * (y2 + z2) * sinA2, 2 * (x * y * sinA2 + z * sinA * cosA), 2 * (x * z * sinA2 - y * sinA * cosA)], [2 * (y * x * sinA2 - z * sinA * cosA), 1 - 2 * (z2 + x2) * sinA2, 2 * (y * z * sinA2 + x * sinA * cosA)], [2 * (z * x * sinA2 + y * sinA * cosA), 2 * (z * y * sinA2 - x * sinA * cosA), 1 - 2 * (x2 + y2) * sinA2]];
    }
    return m;
  };
  perspT = function(vec3, max_real_depth, min_real_depth, desired_ratio, desired_length) {
    var scalefactor, z0;
    z0 = (max_real_depth * desired_ratio - min_real_depth) / (1 - desired_ratio);
    scalefactor = desired_length * desired_ratio / (1 - desired_ratio);
    return [scalefactor * vec3[0] / (vec3[2] + z0), scalefactor * vec3[1] / (vec3[2] + z0)];
  };
  invperspT = function(x, y, dx, dy, max_real_depth, min_real_depth, desired_ratio, desired_length) {
    var s, s2, xp, xp2, xsphere, yp, yp2, ysphere, z0, z02, zsphere;
    z0 = (max_real_depth * desired_ratio - min_real_depth) / (1 - desired_ratio);
    s = desired_length * desired_ratio / (1 - desired_ratio);
    xp = x - dx;
    yp = y - dy;
    s2 = s * s;
    z02 = z0 * z0;
    xp2 = xp * xp;
    yp2 = yp * yp;
    xsphere = (2 * s * xp * z0 + sqrt(4 * s2 * xp2 * z02 + 4 * xp2 * (s2 + xp2 + yp2) * (1 - z02))) / (2.0 * (s2 + xp2 + yp2));
    ysphere = (s * yp * z0) / (s2 + xp2 + yp2) + (yp * sqrt(4 * s2 * z02 + 4 * (s2 + xp2 + yp2) * (1 - z02))) / (2.0 * (s2 + xp2 + yp2));
    zsphere = sqrt(1 - xsphere * xsphere - ysphere * ysphere);
    return [xsphere, ysphere, zsphere];
  };
  getVec2VecRotM = function(vec1, vec2) {
    var angle, axis;
    axis = cross(vec1, vec2);
    angle = acos(dot(vec1, vec2));
    return vec_rotm(-1 * angle, axis[0], axis[1], axis[2]);
  };
  faceToEdges = function(face) {
    var edges, v1, v2, _i, _len;
    edges = [];
    v1 = face.slice(-1)[0];
    for (_i = 0, _len = face.length; _i < _len; _i++) {
      v2 = face[_i];
      edges.push([v1, v2]);
      v1 = v2;
    }
    return edges;
  };
  vertColors = function(poly) {
    var f, i, v, vertcolors, _i, _j, _len, _len2, _ref, _ref2;
    vertcolors = [];
    _ref = enumerate(poly.face);
    for (_i = 0, _len = _ref.length; _i < _len; _i++) {
      _ref2 = _ref[_i], i = _ref2[0], f = _ref2[1];
      for (_j = 0, _len2 = f.length; _j < _len2; _j++) {
        v = f[_j];
        vertcolors[v] = poly.face_colors[i];
      }
    }
    return vertcolors;
  };
  polyhedron = (function() {
    function polyhedron() {
      this.face = new Array();
      this.xyz = new Array();
      this.name = "null polyhedron";
    }
    polyhedron.prototype.data = function() {
      var nEdges;
      nEdges = this.face.length + this.xyz.length - 2;
      return "(" + this.face.length + " faces, " + nEdges + " edges, " + this.xyz.length + " vertices)";
    };
    polyhedron.prototype.getEdges = function() {
      var a, alledges, b, e, edgeset, finalset, hash, uniqedges, _i, _j, _len, _len2;
      finalset = {};
      uniqedges = [];
      alledges = _.map(this.face, faceToEdges);
      for (_i = 0, _len = alledges.length; _i < _len; _i++) {
        edgeset = alledges[_i];
        for (_j = 0, _len2 = edgeset.length; _j < _len2; _j++) {
          e = edgeset[_j];
          if (e[0] < e[1]) {
            a = e[0], b = e[1];
          } else {
            b = e[0], a = e[1];
          }
          finalset[a + '~' + b] = e;
        }
      }
      for (hash in finalset) {
        e = finalset[hash];
        uniqedges.push(e);
      }
      return uniqedges;
    };
    polyhedron.prototype.centers = function() {
      var centers_array, f, fcenter, v, _i, _j, _len, _len2, _ref;
      centers_array = [];
      _ref = this.face;
      for (_i = 0, _len = _ref.length; _i < _len; _i++) {
        f = _ref[_i];
        fcenter = [0, 0, 0];
        for (_j = 0, _len2 = f.length; _j < _len2; _j++) {
          v = f[_j];
          fcenter = add(fcenter, this.xyz[v]);
        }
        centers_array.push(mult(1.0 / f.length, fcenter));
      }
      return centers_array;
    };
    polyhedron.prototype.normals = function() {
      var f, normals_array, v, _i, _len, _ref;
      normals_array = [];
      _ref = this.face;
      for (_i = 0, _len = _ref.length; _i < _len; _i++) {
        f = _ref[_i];
        normals_array.push(normal((function() {
          var _j, _len2, _results;
          _results = [];
          for (_j = 0, _len2 = f.length; _j < _len2; _j++) {
            v = f[_j];
            _results.push(this.xyz[v]);
          }
          return _results;
        }).call(this)));
      }
      return normals_array;
    };
    polyhedron.prototype.rescale = function(s) {
      this.xyz = rescale(this.xyz, s);
      return this;
    };
    polyhedron.prototype.translate = function(x) {
      this.xyz = _.map(this.xyz, function(v) {
        return [v[0] + x[0], v[1] + x[1], v[2] + x[2]];
      });
      return this;
    };
    polyhedron.prototype.toOBJ = function() {
      var f, i, norm, objstr, v, _i, _j, _k, _l, _len, _len2, _len3, _len4, _ref, _ref2, _ref3, _ref4;
      objstr = "#Produced by polyHédronisme http://levskaya.github.com/polyhedronisme\n";
      objstr += "group " + this.name + "\n";
      objstr += "#vertices\n";
      _ref = this.xyz;
      for (_i = 0, _len = _ref.length; _i < _len; _i++) {
        v = _ref[_i];
        objstr += "v " + v[0] + " " + v[1] + " " + v[2] + "\n";
      }
      objstr += "#normal vector defs \n";
      _ref2 = this.face;
      for (_j = 0, _len2 = _ref2.length; _j < _len2; _j++) {
        f = _ref2[_j];
        norm = normal((function() {
          var _k, _len3, _results;
          _results = [];
          for (_k = 0, _len3 = f.length; _k < _len3; _k++) {
            v = f[_k];
            _results.push(this.xyz[v]);
          }
          return _results;
        }).call(this));
        objstr += "vn " + norm[0] + " " + norm[1] + " " + norm[2] + "\n";
      }
      objstr += "#face defs \n";
      _ref3 = enumerate(this.face);
      for (_k = 0, _len3 = _ref3.length; _k < _len3; _k++) {
        _ref4 = _ref3[_k], i = _ref4[0], f = _ref4[1];
        objstr += "f ";
        for (_l = 0, _len4 = f.length; _l < _len4; _l++) {
          v = f[_l];
          objstr += "" + (v + 1) + "//" + (i + 1) + " ";
        }
        objstr += "\n";
      }
      return objstr;
    };
    polyhedron.prototype.toX3D = function() {
      var SCALE_FACTOR, clr, f, v, x3dstr, _i, _j, _k, _l, _len, _len2, _len3, _len4, _ref, _ref2, _ref3;
      SCALE_FACTOR = .01;
      x3dstr = '<?xml version="1.0" encoding ="UTF-8"?>\n<X3D profile="Interchange" version="3.0">\n<head>\n<component name="Rendering" level="3"/>\n<meta name="generator" content="Polyhedronisme"/>\n<meta name="version" content="0.1.0"/>\n</head>\n<Scene>\n<Shape>\n<IndexedFaceSet normalPerVertex="false" coordIndex="';
      _ref = this.face;
      for (_i = 0, _len = _ref.length; _i < _len; _i++) {
        f = _ref[_i];
        for (_j = 0, _len2 = f.length; _j < _len2; _j++) {
          v = f[_j];
          x3dstr += "" + v + " ";
        }
        x3dstr += '-1\n';
      }
      x3dstr += '">\n';
      x3dstr += '<Color color="';
      _ref2 = vertColors(this);
      for (_k = 0, _len3 = _ref2.length; _k < _len3; _k++) {
        clr = _ref2[_k];
        x3dstr += "" + clr[0] + " " + clr[1] + " " + clr[2] + " ";
      }
      x3dstr += '"/>';
      x3dstr += '<Coordinate point="';
      _ref3 = this.xyz;
      for (_l = 0, _len4 = _ref3.length; _l < _len4; _l++) {
        v = _ref3[_l];
        x3dstr += "" + (v[0] * SCALE_FACTOR) + " " + (v[1] * SCALE_FACTOR) + " " + (v[2] * SCALE_FACTOR) + " ";
      }
      x3dstr += '"/>\n';
      x3dstr += '</IndexedFaceSet>\n</Shape>\n</Scene>\n</X3D>';
      return x3dstr;
    };
    polyhedron.prototype.toVRML = function() {
      var SCALE_FACTOR, clr, f, v, x3dstr, _i, _j, _k, _l, _len, _len2, _len3, _len4, _ref, _ref2, _ref3;
      SCALE_FACTOR = .01;
      x3dstr = '#VRML V2.0 utf8\n#Generated by Polyhedronisme\nNavigationInfo {\n	type [ "EXAMINE", "ANY" ]\n}\nTransform {\n  scale 1 1 1\n  translation 0 0 0\n  children\n  [\n    Shape\n    {\n      geometry IndexedFaceSet\n      {\n        creaseAngle .5\n        solid FALSE\n        coord Coordinate\n        {\n          point\n          [';
      _ref = this.xyz;
      for (_i = 0, _len = _ref.length; _i < _len; _i++) {
        v = _ref[_i];
        x3dstr += "" + (v[0] * SCALE_FACTOR) + " " + (v[1] * SCALE_FACTOR) + " " + (v[2] * SCALE_FACTOR) + ",";
      }
      x3dstr = x3dstr.slice(0, -1);
      x3dstr += '    ]\n}\ncolor Color\n{\n  color\n  [';
      _ref2 = this.face_colors;
      for (_j = 0, _len2 = _ref2.length; _j < _len2; _j++) {
        clr = _ref2[_j];
        x3dstr += "" + clr[0] + " " + clr[1] + " " + clr[2] + " ,";
      }
      x3dstr = x3dstr.slice(0, -1);
      x3dstr += '  ]\n}\ncolorPerVertex FALSE\ncoordIndex\n[';
      _ref3 = this.face;
      for (_k = 0, _len3 = _ref3.length; _k < _len3; _k++) {
        f = _ref3[_k];
        for (_l = 0, _len4 = f.length; _l < _len4; _l++) {
          v = f[_l];
          x3dstr += "" + v + ", ";
        }
        x3dstr += '-1,';
      }
      x3dstr = x3dstr.slice(0, -1);
      x3dstr += '          ]\n      }\n      appearance Appearance\n      {\n        material Material\n        {\n	       ambientIntensity 0.2\n	       diffuseColor 0.9 0.9 0.9\n	       specularColor .1 .1 .1\n	       shininess .5\n        }\n      }\n    }\n  ]\n}';
      return x3dstr;
    };
    return polyhedron;
  })();
  tetrahedron = function() {
    var poly;
    poly = new polyhedron();
    poly.name = "T";
    poly.face = [[0, 1, 2], [0, 2, 3], [0, 3, 1], [1, 3, 2]];
    poly.xyz = [[1.0, 1.0, 1.0], [1.0, -1.0, -1.0], [-1.0, 1.0, -1.0], [-1.0, -1.0, 1.0]];
    return poly;
  };
  octahedron = function() {
    var poly;
    poly = new polyhedron();
    poly.name = "O";
    poly.face = [[0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 1], [1, 4, 5], [1, 5, 2], [2, 5, 3], [3, 5, 4]];
    poly.xyz = [[0, 0, 1.414], [1.414, 0, 0], [0, 1.414, 0], [-1.414, 0, 0], [0, -1.414, 0], [0, 0, -1.414]];
    return poly;
  };
  cube = function() {
    var poly;
    poly = new polyhedron();
    poly.name = "C";
    poly.face = [[3, 0, 1, 2], [3, 4, 5, 0], [0, 5, 6, 1], [1, 6, 7, 2], [2, 7, 4, 3], [5, 4, 7, 6]];
    poly.xyz = [[0.707, 0.707, 0.707], [-0.707, 0.707, 0.707], [-0.707, -0.707, 0.707], [0.707, -0.707, 0.707], [0.707, -0.707, -0.707], [0.707, 0.707, -0.707], [-0.707, 0.707, -0.707], [-0.707, -0.707, -0.707]];
    return poly;
  };
  icosahedron = function() {
    var poly;
    poly = new polyhedron();
    poly.name = "I";
    poly.face = [[0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 5], [0, 5, 1], [1, 5, 7], [1, 7, 6], [1, 6, 2], [2, 6, 8], [2, 8, 3], [3, 8, 9], [3, 9, 4], [4, 9, 10], [4, 10, 5], [5, 10, 7], [6, 7, 11], [6, 11, 8], [7, 10, 11], [8, 11, 9], [9, 11, 10]];
    poly.xyz = [[0, 0, 1.176], [1.051, 0, 0.526], [0.324, 1.0, 0.525], [-0.851, 0.618, 0.526], [-0.851, -0.618, 0.526], [0.325, -1.0, 0.526], [0.851, 0.618, -0.526], [0.851, -0.618, -0.526], [-0.325, 1.0, -0.526], [-1.051, 0, -0.526], [-0.325, -1.0, -0.526], [0, 0, -1.176]];
    return poly;
  };
  dodecahedron = function() {
    var poly;
    poly = new polyhedron();
    poly.name = "D";
    poly.face = [[0, 1, 4, 7, 2], [0, 2, 6, 9, 3], [0, 3, 8, 5, 1], [1, 5, 11, 10, 4], [2, 7, 13, 12, 6], [3, 9, 15, 14, 8], [4, 10, 16, 13, 7], [5, 8, 14, 17, 11], [6, 12, 18, 15, 9], [10, 11, 17, 19, 16], [12, 13, 16, 19, 18], [14, 15, 18, 19, 17]];
    poly.xyz = [[0, 0, 1.07047], [0.713644, 0, 0.797878], [-0.356822, 0.618, 0.797878], [-0.356822, -0.618, 0.797878], [0.797878, 0.618034, 0.356822], [0.797878, -0.618, 0.356822], [-0.934172, 0.381966, 0.356822], [0.136294, 1.0, 0.356822], [0.136294, -1.0, 0.356822], [-0.934172, -0.381966, 0.356822], [0.934172, 0.381966, -0.356822], [0.934172, -0.381966, -0.356822], [-0.797878, 0.618, -0.356822], [-0.136294, 1.0, -0.356822], [-0.136294, -1.0, -0.356822], [-0.797878, -0.618034, -0.356822], [0.356822, 0.618, -0.797878], [0.356822, -0.618, -0.797878], [-0.713644, 0, -0.797878], [0, 0, -1.07047]];
    return poly;
  };
  prism = function(n) {
    var h, i, poly, theta, _i, _j, _ref, _ref2, _ref3, _ref4, _ref5, _results, _results2;
    theta = 2 * PI / n;
    h = Math.sin(theta / 2);
    poly = new polyhedron();
    poly.name = "P" + n;
    for (i = 0, _ref = n - 1; 0 <= _ref ? i <= _ref : i >= _ref; 0 <= _ref ? i++ : i--) {
      poly.xyz.push([-cos(i * theta), -sin(i * theta), -h]);
    }
    for (i = 0, _ref2 = n - 1; 0 <= _ref2 ? i <= _ref2 : i >= _ref2; 0 <= _ref2 ? i++ : i--) {
      poly.xyz.push([-cos(i * theta), -sin(i * theta), h]);
    }
    poly.face.push((function() {
      _results = [];
      for (var _i = _ref3 = n - 1; _ref3 <= 0 ? _i <= 0 : _i >= 0; _ref3 <= 0 ? _i++ : _i--){ _results.push(_i); }
      return _results;
    }).apply(this));
    poly.face.push((function() {
      _results2 = [];
      for (var _j = n, _ref4 = 2 * n - 1; n <= _ref4 ? _j <= _ref4 : _j >= _ref4; n <= _ref4 ? _j++ : _j--){ _results2.push(_j); }
      return _results2;
    }).apply(this));
    for (i = 0, _ref5 = n - 1; 0 <= _ref5 ? i <= _ref5 : i >= _ref5; 0 <= _ref5 ? i++ : i--) {
      poly.face.push([i, (i + 1) % n, (i + 1) % n + n, i + n]);
    }
    poly.xyz = adjustXYZ(poly, 1);
    return poly;
  };
  antiprism = function(n) {
    var f, h, i, poly, r, theta, _i, _j, _ref, _ref2, _ref3, _ref4, _ref5, _results, _results2;
    theta = 2 * PI / n;
    h = sqrt(1 - 4 / (4 + 2 * cos(theta / 2) - 2 * cos(theta)));
    r = sqrt(1 - h * h);
    f = sqrt(h * h + pow(r * cos(theta / 2), 2));
    r = -r / f;
    h = -h / f;
    poly = new polyhedron();
    poly.name = "A" + n;
    for (i = 0, _ref = n - 1; 0 <= _ref ? i <= _ref : i >= _ref; 0 <= _ref ? i++ : i--) {
      poly.xyz.push([r * cos(i * theta), r * sin(i * theta), h]);
    }
    for (i = 0, _ref2 = n - 1; 0 <= _ref2 ? i <= _ref2 : i >= _ref2; 0 <= _ref2 ? i++ : i--) {
      poly.xyz.push([r * cos((i + 0.5) * theta), r * sin((i + 0.5) * theta), -h]);
    }
    poly.face.push((function() {
      _results = [];
      for (var _i = _ref3 = n - 1; _ref3 <= 0 ? _i <= 0 : _i >= 0; _ref3 <= 0 ? _i++ : _i--){ _results.push(_i); }
      return _results;
    }).apply(this));
    poly.face.push((function() {
      _results2 = [];
      for (var _j = n, _ref4 = 2 * n - 1; n <= _ref4 ? _j <= _ref4 : _j >= _ref4; n <= _ref4 ? _j++ : _j--){ _results2.push(_j); }
      return _results2;
    }).apply(this));
    for (i = 0, _ref5 = n - 1; 0 <= _ref5 ? i <= _ref5 : i >= _ref5; 0 <= _ref5 ? i++ : i--) {
      poly.face.push([i, (i + 1) % n, i + n]);
      poly.face.push([i, i + n, (n + i - 1) % n + n]);
    }
    poly.xyz = adjustXYZ(poly, 1);
    return poly;
  };
  pyramid = function(n) {
    var height, i, poly, theta, _i, _ref, _ref2, _ref3, _results;
    theta = 2 * PI / n;
    height = 1;
    poly = new polyhedron();
    poly.name = "Y" + n;
    for (i = 0, _ref = n - 1; 0 <= _ref ? i <= _ref : i >= _ref; 0 <= _ref ? i++ : i--) {
      poly.xyz.push([-cos(i * theta), -sin(i * theta), -0.2]);
    }
    poly.xyz.push([0, 0, height]);
    poly.face.push((function() {
      _results = [];
      for (var _i = _ref2 = n - 1; _ref2 <= 0 ? _i <= 0 : _i >= 0; _ref2 <= 0 ? _i++ : _i--){ _results.push(_i); }
      return _results;
    }).apply(this));
    for (i = 0, _ref3 = n - 1; 0 <= _ref3 ? i <= _ref3 : i >= _ref3; 0 <= _ref3 ? i++ : i--) {
      poly.face.push([i, (i + 1) % n, n]);
    }
    poly.xyz = canonicalXYZ(poly, 3);
    return poly;
  };
  polyflag = (function() {
    function polyflag() {
      this.flags = new Object();
      this.verts = new Object();
      this.xyzs = new Object();
    }
    polyflag.prototype.newV = function(name, xyz) {
      if (this.verts[name] === void 0) {
        this.verts[name] = 0;
        return this.xyzs[name] = xyz;
      }
    };
    polyflag.prototype.newFlag = function(face, v1, v2) {
      if (this.flags[face] === void 0) {
        this.flags[face] = {};
      }
      return this.flags[face][v1] = v2;
    };
    polyflag.prototype.topoly = function() {
      var ctr, f, faceCTR, i, j, poly, v, v0, _ref, _ref2;
      poly = new polyhedron();
      ctr = 0;
      _ref = this.verts;
      for (i in _ref) {
        v = _ref[i];
        poly.xyz[ctr] = this.xyzs[i];
        this.verts[i] = ctr;
        ctr++;
      }
      ctr = 0;
      _ref2 = this.flags;
      for (i in _ref2) {
        f = _ref2[i];
        poly.face[ctr] = [];
        for (j in f) {
          v = f[j];
          v0 = v;
          break;
        }
        v = v0;
        poly.face[ctr].push(this.verts[v]);
        v = this.flags[i][v];
        faceCTR = 0;
        while (v !== v0) {
          poly.face[ctr].push(this.verts[v]);
          v = this.flags[i][v];
          faceCTR++;
          if (faceCTR > 1000) {
            console.log("Bad flag spec, have a neverending face:", i, this.flags[i]);
            break;
          }
        }
        ctr++;
      }
      poly.name = "unknown polyhedron";
      return poly;
    };
    return polyflag;
  })();
  kisN = function(poly, n) {
    var apex, centers, f, flag, fname, foundAny, i, newpoly, normals, p, v, v1, v2, _i, _j, _k, _len, _len2, _len3, _ref, _ref2, _ref3, _ref4;
    console.log("Taking kis of " + (n === 0 ? "" : n) + "-sided faces of " + poly.name + "...");
    flag = new polyflag();
    _ref = enumerate(poly.xyz);
    for (_i = 0, _len = _ref.length; _i < _len; _i++) {
      _ref2 = _ref[_i], i = _ref2[0], p = _ref2[1];
      flag.newV("v" + i, p);
    }
    normals = poly.normals();
    centers = poly.centers();
    foundAny = false;
    _ref3 = enumerate(poly.face);
    for (_j = 0, _len2 = _ref3.length; _j < _len2; _j++) {
      _ref4 = _ref3[_j], i = _ref4[0], f = _ref4[1];
      v1 = "v" + f[f.length - 1];
      for (_k = 0, _len3 = f.length; _k < _len3; _k++) {
        v = f[_k];
        v2 = "v" + v;
        if (f.length === n || n === 0) {
          foundAny = true;
          apex = "apex" + i;
          fname = "" + i + v1;
          flag.newV(apex, add(centers[i], mult(0.2, normals[i])));
          flag.newFlag(fname, v1, v2);
          flag.newFlag(fname, v2, apex);
          flag.newFlag(fname, apex, v1);
        } else {
          flag.newFlag("" + i, v1, v2);
        }
        v1 = v2;
      }
    }
    if (!foundAny) {
      console.log("No " + n + "-fold components were found.");
    }
    newpoly = flag.topoly();
    newpoly.name = "k" + (n === 0 ? "" : n) + poly.name;
    return newpoly;
  };
  ambo = function(poly) {
    var f, flag, i, midName, newpoly, v1, v2, v3, _i, _j, _len, _len2, _ref, _ref2, _ref3, _ref4;
    console.log("Taking ambo of " + poly.name + "...");
    midName = function(v1, v2) {
      if (v1 < v2) {
        return v1 + "_" + v2;
      } else {
        return v2 + "_" + v1;
      }
    };
    flag = new polyflag();
    _ref = enumerate(poly.face);
    for (_i = 0, _len = _ref.length; _i < _len; _i++) {
      _ref2 = _ref[_i], i = _ref2[0], f = _ref2[1];
      _ref3 = f.slice(-2), v1 = _ref3[0], v2 = _ref3[1];
      for (_j = 0, _len2 = f.length; _j < _len2; _j++) {
        v3 = f[_j];
        if (v1 < v2) {
          flag.newV(midName(v1, v2), midpoint(poly.xyz[v1], poly.xyz[v2]));
        }
        flag.newFlag("orig" + i, midName(v1, v2), midName(v2, v3));
        flag.newFlag("dual" + v2, midName(v2, v3), midName(v1, v2));
        _ref4 = [v2, v3], v1 = _ref4[0], v2 = _ref4[1];
      }
    }
    newpoly = flag.topoly();
    newpoly.name = "a" + poly.name;
    return newpoly;
  };
  gyro = function(poly) {
    var centers, f, flag, fname, i, j, newpoly, v, v1, v2, v3, _i, _j, _k, _l, _len, _len2, _len3, _len4, _ref, _ref10, _ref2, _ref3, _ref4, _ref5, _ref6, _ref7, _ref8, _ref9;
    console.log("Taking gyro of " + poly.name + "...");
    flag = new polyflag();
    _ref = enumerate(poly.xyz);
    for (_i = 0, _len = _ref.length; _i < _len; _i++) {
      _ref2 = _ref[_i], i = _ref2[0], v = _ref2[1];
      flag.newV("v" + i, unit(v));
    }
    centers = poly.centers();
    _ref3 = enumerate(poly.face);
    for (_j = 0, _len2 = _ref3.length; _j < _len2; _j++) {
      _ref4 = _ref3[_j], i = _ref4[0], f = _ref4[1];
      flag.newV("center" + i, unit(centers[i]));
    }
    _ref5 = enumerate(poly.face);
    for (_k = 0, _len3 = _ref5.length; _k < _len3; _k++) {
      _ref6 = _ref5[_k], i = _ref6[0], f = _ref6[1];
      _ref7 = f.slice(-2), v1 = _ref7[0], v2 = _ref7[1];
      _ref8 = enumerate(f);
      for (_l = 0, _len4 = _ref8.length; _l < _len4; _l++) {
        _ref9 = _ref8[_l], j = _ref9[0], v = _ref9[1];
        v3 = v;
        flag.newV(v1 + "~" + v2, oneThird(poly.xyz[v1], poly.xyz[v2]));
        fname = i + "f" + v1;
        flag.newFlag(fname, "center" + i, v1 + "~" + v2);
        flag.newFlag(fname, v1 + "~" + v2, v2 + "~" + v1);
        flag.newFlag(fname, v2 + "~" + v1, "v" + v2);
        flag.newFlag(fname, "v" + v2, v2 + "~" + v3);
        flag.newFlag(fname, v2 + "~" + v3, "center" + i);
        _ref10 = [v2, v3], v1 = _ref10[0], v2 = _ref10[1];
      }
    }
    newpoly = flag.topoly();
    newpoly.name = "g" + poly.name;
    return newpoly;
  };
  propellor = function(poly) {
    var f, flag, fname, i, newpoly, v, v1, v2, v3, _i, _j, _k, _len, _len2, _len3, _ref, _ref2, _ref3, _ref4, _ref5, _ref6;
    console.log("Taking propellor of " + poly.name + "...");
    flag = new polyflag();
    _ref = enumerate(poly.xyz);
    for (_i = 0, _len = _ref.length; _i < _len; _i++) {
      _ref2 = _ref[_i], i = _ref2[0], v = _ref2[1];
      flag.newV("v" + i, unit(v));
    }
    _ref3 = enumerate(poly.face);
    for (_j = 0, _len2 = _ref3.length; _j < _len2; _j++) {
      _ref4 = _ref3[_j], i = _ref4[0], f = _ref4[1];
      _ref5 = f.slice(-2), v1 = _ref5[0], v2 = _ref5[1];
      for (_k = 0, _len3 = f.length; _k < _len3; _k++) {
        v = f[_k];
        v3 = "" + v;
        flag.newV(v1 + "~" + v2, oneThird(poly.xyz[v1], poly.xyz[v2]));
        fname = "" + i + "f" + v2;
        flag.newFlag("v" + i, v1 + "~" + v2, v2 + "~" + v3);
        flag.newFlag(fname, v1 + "~" + v2, v2 + "~" + v1);
        flag.newFlag(fname, v2 + "~" + v1, "v" + v2);
        flag.newFlag(fname, "v" + v2, v2 + "~" + v3);
        flag.newFlag(fname, v2 + "~" + v3, v1 + "~" + v2);
        _ref6 = [v2, v3], v1 = _ref6[0], v2 = _ref6[1];
      }
    }
    newpoly = flag.topoly();
    newpoly.name = "p" + poly.name;
    return newpoly;
  };
  reflect = function(poly) {
    var i, _ref, _ref2;
    console.log("Taking reflection of " + poly.name + "...");
    for (i = 0, _ref = poly.xyz.length - 1; 0 <= _ref ? i <= _ref : i >= _ref; 0 <= _ref ? i++ : i--) {
      poly.xyz[i] = mult(-1, poly.xyz[i]);
    }
    for (i = 0, _ref2 = poly.face.length - 1; 0 <= _ref2 ? i <= _ref2 : i >= _ref2; 0 <= _ref2 ? i++ : i--) {
      poly.face[i] = poly.face[i].reverse();
    }
    poly.name = "r" + poly.name;
    return poly;
  };
  dual = function(poly) {
    var centers, dpoly, f, face, flag, i, k, sortF, v1, v2, _i, _j, _k, _l, _len, _len2, _len3, _len4, _len5, _m, _ref, _ref2, _ref3, _ref4, _ref5, _ref6, _ref7;
    console.log("Taking dual of " + poly.name + "...");
    flag = new polyflag();
    face = [];
    for (i = 0, _ref = poly.xyz.length - 1; 0 <= _ref ? i <= _ref : i >= _ref; 0 <= _ref ? i++ : i--) {
      face[i] = {};
    }
    _ref2 = enumerate(poly.face);
    for (_i = 0, _len = _ref2.length; _i < _len; _i++) {
      _ref3 = _ref2[_i], i = _ref3[0], f = _ref3[1];
      v1 = f[f.length - 1];
      for (_j = 0, _len2 = f.length; _j < _len2; _j++) {
        v2 = f[_j];
        face[v1]["v" + v2] = "" + i;
        v1 = v2;
      }
    }
    centers = poly.centers();
    for (i = 0, _ref4 = poly.face.length - 1; 0 <= _ref4 ? i <= _ref4 : i >= _ref4; 0 <= _ref4 ? i++ : i--) {
      flag.newV("" + i, centers[i]);
    }
    _ref5 = enumerate(poly.face);
    for (_k = 0, _len3 = _ref5.length; _k < _len3; _k++) {
      _ref6 = _ref5[_k], i = _ref6[0], f = _ref6[1];
      v1 = f[f.length - 1];
      for (_l = 0, _len4 = f.length; _l < _len4; _l++) {
        v2 = f[_l];
        flag.newFlag(v1, face[v2]["v" + v1], "" + i);
        v1 = v2;
      }
    }
    console.log("dual flag:", flag);
    dpoly = flag.topoly();
    sortF = [];
    _ref7 = dpoly.face;
    for (_m = 0, _len5 = _ref7.length; _m < _len5; _m++) {
      f = _ref7[_m];
      k = intersect(poly.face[f[0]], poly.face[f[1]], poly.face[f[2]]);
      sortF[k] = f;
    }
    dpoly.face = sortF;
    if (poly.name[0] !== "d") {
      dpoly.name = "d" + poly.name;
    } else {
      dpoly.name = poly.name.slice(1);
    }
    return dpoly;
  };
  insetN = function(poly, n) {
    var centers, f, flag, fname, foundAny, i, newpoly, normals, p, v, v1, v2, _i, _j, _k, _l, _len, _len2, _len3, _len4, _len5, _m, _ref, _ref2, _ref3, _ref4, _ref5, _ref6;
    console.log("Taking inset of " + (n === 0 ? "" : n) + "-sided faces of " + poly.name + "...");
    flag = new polyflag();
    _ref = enumerate(poly.xyz);
    for (_i = 0, _len = _ref.length; _i < _len; _i++) {
      _ref2 = _ref[_i], i = _ref2[0], p = _ref2[1];
      flag.newV("v" + i, p);
    }
    normals = poly.normals();
    centers = poly.centers();
    _ref3 = enumerate(poly.face);
    for (_j = 0, _len2 = _ref3.length; _j < _len2; _j++) {
      _ref4 = _ref3[_j], i = _ref4[0], f = _ref4[1];
      if (f.length === n || n === 0) {
        for (_k = 0, _len3 = f.length; _k < _len3; _k++) {
          v = f[_k];
          flag.newV("f" + i + "v" + v, add(midpoint(poly.xyz[v], centers[i]), mult(-0.2, normals[i])));
        }
      }
    }
    foundAny = false;
    _ref5 = enumerate(poly.face);
    for (_l = 0, _len4 = _ref5.length; _l < _len4; _l++) {
      _ref6 = _ref5[_l], i = _ref6[0], f = _ref6[1];
      v1 = "v" + f[f.length - 1];
      for (_m = 0, _len5 = f.length; _m < _len5; _m++) {
        v = f[_m];
        v2 = "v" + v;
        if (f.length === n || n === 0) {
          foundAny = true;
          fname = i + v1;
          flag.newFlag(fname, v1, v2);
          flag.newFlag(fname, v2, "f" + i + v2);
          flag.newFlag(fname, "f" + i + v2, "f" + i + v1);
          flag.newFlag(fname, "f" + i + v1, v1);
          flag.newFlag("ex" + i, "f" + i + v1, "f" + i + v2);
        } else {
          flag.newFlag(i, v1, v2);
        }
        v1 = v2;
      }
    }
    if (!foundAny) {
      console.log("No " + n + "-fold components were found.");
    }
    newpoly = flag.topoly();
    newpoly.name = "n" + (n === 0 ? "" : n) + poly.name;
    return newpoly;
  };
  extrudeN = function(poly, n) {
    var centers, f, flag, foundAny, i, newpoly, normals, p, v, v1, v2, _i, _j, _k, _l, _len, _len2, _len3, _len4, _len5, _m, _ref, _ref2, _ref3, _ref4, _ref5, _ref6;
    console.log("Taking extrusion of " + (n === 0 ? "" : n) + "-sided faces of " + poly.name + "...");
    flag = new polyflag();
    _ref = enumerate(poly.xyz);
    for (_i = 0, _len = _ref.length; _i < _len; _i++) {
      _ref2 = _ref[_i], i = _ref2[0], p = _ref2[1];
      flag.newV("v" + i, p);
    }
    normals = poly.normals();
    centers = poly.centers();
    _ref3 = enumerate(poly.face);
    for (_j = 0, _len2 = _ref3.length; _j < _len2; _j++) {
      _ref4 = _ref3[_j], i = _ref4[0], f = _ref4[1];
      if (f.length === n || n === 0) {
        for (_k = 0, _len3 = f.length; _k < _len3; _k++) {
          v = f[_k];
          flag.newV("f" + i + "v" + v, add(poly.xyz[v], mult(0.3, normals[i])));
        }
      }
    }
    foundAny = false;
    _ref5 = enumerate(poly.face);
    for (_l = 0, _len4 = _ref5.length; _l < _len4; _l++) {
      _ref6 = _ref5[_l], i = _ref6[0], f = _ref6[1];
      v1 = "v" + f[f.length - 1];
      for (_m = 0, _len5 = f.length; _m < _len5; _m++) {
        v = f[_m];
        v2 = "v" + v;
        if (f.length === n || n === 0) {
          foundAny = true;
          flag.newFlag(i + v1, v1, v2);
          flag.newFlag(i + v1, v2, "f" + i + v2);
          flag.newFlag(i + v1, "f" + i + v2, "f" + i + v1);
          flag.newFlag(i + v1, "f" + i + v1, v1);
          flag.newFlag("ex" + i, "f" + i + v1, "f" + i + v2);
        } else {
          flag.newFlag(i, v1, v2);
        }
        v1 = v2;
      }
    }
    if (!foundAny) {
      console.log("No " + n + "-fold components were found.");
    }
    newpoly = flag.topoly();
    newpoly.name = "x" + (n === 0 ? "" : n) + poly.name;
    console.log(newpoly);
    return newpoly;
  };
  stellaN = function(poly) {
    var centers, f, flag, i, newpoly, p, v, v1, v12, v2, v21, v23, v3, vert1, vert2, vert3, _i, _j, _k, _len, _len2, _len3, _ref, _ref2, _ref3, _ref4, _ref5, _ref6;
    console.log("Taking stella of " + poly.name + "...");
    centers = poly.centers();
    flag = new polyflag();
    _ref = enumerate(poly.xyz);
    for (_i = 0, _len = _ref.length; _i < _len; _i++) {
      _ref2 = _ref[_i], i = _ref2[0], p = _ref2[1];
      flag.newV("v" + i, p);
    }
    _ref3 = enumerate(poly.face);
    for (_j = 0, _len2 = _ref3.length; _j < _len2; _j++) {
      _ref4 = _ref3[_j], i = _ref4[0], f = _ref4[1];
      v1 = "v" + f[f.length - 2];
      v2 = "v" + f[f.length - 1];
      vert1 = poly.xyz[f[f.length - 2]];
      vert2 = poly.xyz[f[f.length - 1]];
      for (_k = 0, _len3 = f.length; _k < _len3; _k++) {
        v = f[_k];
        v3 = "v" + v;
        vert3 = poly.xyz[v];
        v12 = v1 + "~" + v2;
        v21 = v2 + "~" + v1;
        v23 = v2 + "~" + v3;
        flag.newV(v12, midpoint(midpoint(vert1, vert2), centers[i]));
        flag.newFlag("in" + i, v12, v23);
        flag.newFlag("f" + i + v2, v23, v12);
        flag.newFlag("f" + i + v2, v12, v2);
        flag.newFlag("f" + i + v2, v2, v23);
        flag.newFlag("f" + v12, v1, v21);
        flag.newFlag("f" + v12, v21, v12);
        flag.newFlag("f" + v12, v12, v1);
        _ref5 = [v2, v3], v1 = _ref5[0], v2 = _ref5[1];
        _ref6 = [vert2, vert3], vert1 = _ref6[0], vert2 = _ref6[1];
      }
    }
    newpoly = flag.topoly();
    newpoly.name = "*" + poly.name;
    return newpoly;
  };
  tangentify = function(xyzs, edges) {
    var STABILITY_FACTOR, c, e, newVs, t, _i, _len;
    STABILITY_FACTOR = 0.1;
    newVs = copyVecArray(xyzs);
    for (_i = 0, _len = edges.length; _i < _len; _i++) {
      e = edges[_i];
      t = tangentPoint(newVs[e[0]], newVs[e[1]]);
      c = mult(STABILITY_FACTOR * 1 / 2 * (1 - sqrt(dot(t, t))), t);
      newVs[e[0]] = add(newVs[e[0]], c);
      newVs[e[1]] = add(newVs[e[1]], c);
    }
    return newVs;
  };
  recenter = function(xyzs, edges) {
    var a, b, edgecenters, polycenter, v, _i, _len;
    edgecenters = (function() {
      var _i, _len, _ref, _results;
      _results = [];
      for (_i = 0, _len = edges.length; _i < _len; _i++) {
        _ref = edges[_i], a = _ref[0], b = _ref[1];
        _results.push(tangentPoint(xyzs[a], xyzs[b]));
      }
      return _results;
    })();
    polycenter = [0, 0, 0];
    for (_i = 0, _len = edgecenters.length; _i < _len; _i++) {
      v = edgecenters[_i];
      polycenter = add(polycenter, v);
    }
    polycenter = mult(1 / edges.length, polycenter);
    return _.map(xyzs, function(x) {
      return sub(x, polycenter);
    });
  };
  rescale = function(xyzs, scale) {
    var maxExtent, polycenter, s;
    polycenter = [0, 0, 0];
    maxExtent = _.max(_.map(xyzs, function(x) {
      return mag(x);
    }));
    s = 1 / maxExtent;
    scale || (scale = s);
    return _.map(xyzs, function(x) {
      return [scale * x[0], scale * x[1], scale * x[2]];
    });
  };
  planarize = function(xyzs, faces) {
    var STABILITY_FACTOR, c, coords, f, n, newVs, v, _i, _j, _len, _len2;
    STABILITY_FACTOR = 0.1;
    newVs = copyVecArray(xyzs);
    for (_i = 0, _len = faces.length; _i < _len; _i++) {
      f = faces[_i];
      coords = (function() {
        var _j, _len2, _results;
        _results = [];
        for (_j = 0, _len2 = f.length; _j < _len2; _j++) {
          v = f[_j];
          _results.push(xyzs[v]);
        }
        return _results;
      })();
      n = normal(coords);
      c = calcCentroid(coords);
      if (dot(n, c) < 0) {
        n = mult(-1.0, n);
      }
      for (_j = 0, _len2 = f.length; _j < _len2; _j++) {
        v = f[_j];
        newVs[v] = add(newVs[v], mult(dot(mult(STABILITY_FACTOR, n), sub(c, xyzs[v])), n));
      }
    }
    return newVs;
  };
  canonicalize = function(poly, Niter) {
    var edges, faces, i, maxChange, newVs, oldVs;
    console.log("Canonicalizing " + poly.name + "...");
    faces = poly.face;
    edges = poly.getEdges();
    newVs = poly.xyz;
    maxChange = 1.0;
    for (i = 0; 0 <= Niter ? i <= Niter : i >= Niter; 0 <= Niter ? i++ : i--) {
      oldVs = copyVecArray(newVs);
      newVs = tangentify(newVs, edges);
      newVs = recenter(newVs, edges);
      newVs = planarize(newVs, faces);
      maxChange = _.max(_.map(_.zip(newVs, oldVs), function(_arg) {
        var x, y;
        x = _arg[0], y = _arg[1];
        return mag(sub(x, y));
      }));
      if (maxChange < 1e-8) {
        break;
      }
    }
    console.log("[canonicalization done, last |deltaV|=" + maxChange + "]");
    return newVs;
  };
  reciprocalC = function(poly) {
    var c, centers, _i, _len;
    centers = poly.centers();
    for (_i = 0, _len = centers.length; _i < _len; _i++) {
      c = centers[_i];
      c = mult(1.0 / dot(c, c), c);
    }
    return centers;
  };
  reciprocalN = function(poly) {
    var ans, avgEdgeDist, centroid, f, normalV, tmp, v1, v2, v3, _i, _j, _len, _len2, _ref, _ref2, _ref3;
    ans = [];
    _ref = poly.face;
    for (_i = 0, _len = _ref.length; _i < _len; _i++) {
      f = _ref[_i];
      centroid = [0, 0, 0];
      normalV = [0, 0, 0];
      avgEdgeDist = 0.0;
      _ref2 = f.slice(-2), v1 = _ref2[0], v2 = _ref2[1];
      for (_j = 0, _len2 = f.length; _j < _len2; _j++) {
        v3 = f[_j];
        centroid = add(centroid, poly.xyz[v3]);
        normalV = add(normalV, orthogonal(poly.xyz[v1], poly.xyz[v2], poly.xyz[v3]));
        avgEdgeDist += edgeDist(poly.xyz[v1], poly.xyz[v2]);
        _ref3 = [v2, v3], v1 = _ref3[0], v2 = _ref3[1];
      }
      centroid = mult(1.0 / f.length, centroid);
      normalV = unit(normalV);
      avgEdgeDist = avgEdgeDist / f.length;
      tmp = reciprocal(mult(dot(centroid, normalV), normalV));
      ans.push(mult((1 + avgEdgeDist) / 2, tmp));
    }
    return ans;
  };
  canonicalXYZ = function(poly, nIterations) {
    var count, dpoly, _ref;
    dpoly = dual(poly);
    console.log("Pseudo-canonicalizing " + poly.name + "...");
    for (count = 0, _ref = nIterations - 1; 0 <= _ref ? count <= _ref : count >= _ref; 0 <= _ref ? count++ : count--) {
      dpoly.xyz = reciprocalN(poly);
      poly.xyz = reciprocalN(dpoly);
    }
    return poly.xyz;
  };
  adjustXYZ = function(poly, nIterations) {
    var count, dpoly, _ref;
    dpoly = dual(poly);
    console.log("Planarizing " + poly.name + "...");
    for (count = 0, _ref = nIterations - 1; 0 <= _ref ? count <= _ref : count >= _ref; 0 <= _ref ? count++ : count--) {
      dpoly.xyz = reciprocalC(poly);
      poly.xyz = reciprocalC(dpoly);
    }
    return poly.xyz;
  };
  getDiagonals = function(verts) {
    var Area2, Between, Collinear, Diagonal, Diagonalie, InCone, Intersect, IntersectProp, Left, LeftOn, XOR, broke, diagonals, ear, facelen, head, limiter, n, origIdx, v0, v1, v2, v3, v4, y, z, _i, _ref, _results;
    limiter = 999;
    diagonals = [];
    ear = [];
    facelen = verts.length;
    XOR = function(x, y) {
      return (x || y) && !(x && y);
    };
    Area2 = function(Va, Vb, Vc) {
      return (Vb[0] - Va[0]) * (Vc[1] - Va[1]) - (Vc[0] - Va[0]) * (Vb[1] - Va[1]);
    };
    Left = function(Va, Vb, Vc) {
      return Area2(Va, Vb, Vc) > 0;
    };
    LeftOn = function(Va, Vb, Vc) {
      return Area2(Va, Vb, Vc) >= 0;
    };
    Collinear = function(Va, Vb, Vc) {
      return Area2(Va, Vb, Vc) === 0;
    };
    Between = function(Va, Vb, Vc) {
      if (Collinear(Va, Vb, Vc)) {
        return false;
      }
      if (Va[0] !== Vb[0]) {
        return (Va[0] <= Vc[0]) && (Vc[0] <= Vb[0]) || (Va[0] >= Vc[0]) && (Vc[0] >= Vb[0]);
      } else {
        return (Va[1] <= Vc[1]) && (Vc[1] <= Vb[1]) || (Va[1] >= Vc[1]) && (Vc[1] >= Vb[1]);
      }
    };
    IntersectProp = function(Va, Vb, Vc, Vd) {
      if (Collinear(Va, Vb, Vc) || Collinear(Va, Vb, Vd) || Collinear(Vc, Vd, Va) || Collinear(Vc, Vd, Vb)) {
        return false;
      }
      return XOR(Left(Va, Vb, Vc), Left(Va, Vb, Vd)) && XOR(Left(Vc, Vd, Va), Left(Vc, Vd, Vb));
    };
    Intersect = function(Va, Vb, Vc, Vd) {
      if (IntersectProp(Va, Vb, Vc, Vd)) {
        return true;
      } else {
        if (Between(Va, Vb, Vc) || Between(Va, Vb, Vd) || Between(Vc, Vd, Va) || Between(Vc, Vd, Vb)) {
          return true;
        } else {
          return false;
        }
      }
    };
    InCone = function(a, b) {
      var a0, a1;
      a1 = (a + 1 + facelen) % facelen;
      a0 = (a - 1 + facelen) % facelen;
      if (LeftOn(verts[a], verts[a1], verts[a0])) {
        return Left(verts[a], verts[b], verts[a0]) && Left(verts[b], verts[a], verts[a1]);
      }
      return !(LeftOn(verts[a], verts[b], verts[a1]) && LeftOn(verts[b], verts[a], verts[a0]));
    };
    Diagonalie = function(a, b) {
      var c, c1;
      c = 0;
      while (true) {
        c1 = (c + 1 + facelen) % facelen;
        if ((c !== a) && (c1 !== a) && (c !== b) && (c1 !== b) && IntersectProp(verts[a], verts[b], verts[c], verts[c1])) {
          return false;
        }
        c = (c + 1 + facelen) % facelen;
        if (c === 0) {
          break;
        }
      }
      return true;
    };
    Diagonal = function(a, b) {
      return InCone(a, b) && InCone(b, a) && Diagonalie(a, b);
    };
    v1 = 0;
    while (true) {
      v2 = (v1 + 1 + facelen) % facelen;
      v0 = (v1 - 1 + facelen) % facelen;
      ear[v1] = Diagonal(v0, v2);
      v1 = (v1 + 1 + facelen) % facelen;
      if (v1 === 0) {
        break;
      }
    }
    origIdx = (function() {
      _results = [];
      for (var _i = 0, _ref = facelen - 1; 0 <= _ref ? _i <= _ref : _i >= _ref; 0 <= _ref ? _i++ : _i--){ _results.push(_i); }
      return _results;
    }).apply(this);
    n = facelen;
    z = limiter;
    head = 0;
    while (z > 0 && n > 3) {
      z -= 1;
      v2 = head;
      y = limiter;
      while (true) {
        y -= 1;
        broke = false;
        if (ear[v2]) {
          v3 = (v2 + 1 + facelen) % facelen;
          v4 = (v3 + 1 + facelen) % facelen;
          v1 = (v2 - 1 + facelen) % facelen;
          v0 = (v1 - 1 + facelen) % facelen;
          diagonals.push([origIdx[v1], origIdx[v3]]);
          ear[v1] = Diagonal(v0, v3);
          ear[v3] = Diagonal(v1, v4);
          verts = verts.slice(0, (v2 + 1) || 9e9).concat(verts.slice(v3));
          origIdx = origIdx.slice(0, (v2 + 1) || 9e9).concat(origIdx.slice(v3));
          if (v0 > v2) {
            v0 -= 1;
          }
          if (v1 > v2) {
            v1 -= 1;
          }
          if (v3 > v2) {
            v3 -= 1;
          }
          if (v4 > v2) {
            v4 -= 1;
          }
          facelen--;
          head = v3;
          n--;
          broke = true;
        }
        v2 = (v2 + 1 + facelen) % facelen;
        if (!(y > 0 && !broke && v2 !== head)) {
          break;
        }
      }
    }
    return diagonals;
  };
  triEq = function(tri1, tri2) {
    if (((tri1[0] === tri2[0]) && (tri1[1] === tri2[1]) && (tri1[2] === tri2[2])) || (tri1[0] === tri2[1]) && (tri1[1] === tri2[2]) && (tri1[2] === tri2[0]) || (tri1[0] === tri2[2]) && (tri1[1] === tri2[0]) && (tri1[2] === tri2[1])) {
      return true;
    } else {
      return false;
    }
  };
  diagsToTris = function(f, diags) {
    var already_present, d, e1, e2, edges, extant_tri, i, redges, tri, tris, uniques, v1, v2, _i, _j, _k, _l, _len, _len2, _len3, _len4, _len5, _len6, _len7, _len8, _len9, _m, _n, _o, _p, _q, _ref, _ref2, _ref3, _ref4, _ref5, _ref6;
    edges = [];
    redges = [];
    _ref = (function() {
      var _ref, _results;
      _results = [];
      for (i = 0, _ref = f.length - 1; 0 <= _ref ? i <= _ref : i >= _ref; 0 <= _ref ? i++ : i--) {
        _results.push([i, (i + 1) % f.length]);
      }
      return _results;
    })();
    for (_i = 0, _len = _ref.length; _i < _len; _i++) {
      _ref2 = _ref[_i], v1 = _ref2[0], v2 = _ref2[1];
      edges[v1] = [v2];
      redges[v2] = [v1];
    }
    for (_j = 0, _len2 = diags.length; _j < _len2; _j++) {
      d = diags[_j];
      edges[d[0]].push(d[1]);
      edges[d[1]].push(d[0]);
      redges[d[0]].push(d[1]);
      redges[d[1]].push(d[0]);
    }
    tris = [];
    for (_k = 0, _len3 = diags.length; _k < _len3; _k++) {
      d = diags[_k];
      _ref3 = edges[d[1]];
      for (_l = 0, _len4 = _ref3.length; _l < _len4; _l++) {
        e1 = _ref3[_l];
        _ref4 = redges[d[0]];
        for (_m = 0, _len5 = _ref4.length; _m < _len5; _m++) {
          e2 = _ref4[_m];
          if (e1 === e2) {
            tris.push([d[0], d[1], e1]);
          }
        }
      }
      _ref5 = edges[d[0]];
      for (_n = 0, _len6 = _ref5.length; _n < _len6; _n++) {
        e1 = _ref5[_n];
        _ref6 = redges[d[1]];
        for (_o = 0, _len7 = _ref6.length; _o < _len7; _o++) {
          e2 = _ref6[_o];
          if (e1 === e2) {
            tris.push([d[1], d[0], e1]);
          }
        }
      }
    }
    uniques = [tris.pop()];
    for (_p = 0, _len8 = tris.length; _p < _len8; _p++) {
      tri = tris[_p];
      already_present = false;
      for (_q = 0, _len9 = uniques.length; _q < _len9; _q++) {
        extant_tri = uniques[_q];
        if (triEq(tri, extant_tri)) {
          already_present = true;
          break;
        }
      }
      if (!already_present) {
        uniques.push(tri);
      }
    }
    return uniques;
  };
  triangulate = function(poly, colors) {
    var TwoDface, diags, f, i, j, newpoly, tri, tris, v, _i, _j, _len, _len2, _ref, _ref2, _ref3, _ref4;
    colors = colors || false;
    console.log("Triangulating faces of " + poly.name + "...");
    newpoly = new polyhedron();
    newpoly.xyz = clone(poly.xyz);
    newpoly.face_colors = [];
    _ref = enumerate(poly.face);
    for (_i = 0, _len = _ref.length; _i < _len; _i++) {
      _ref2 = _ref[_i], i = _ref2[0], f = _ref2[1];
      if (f.length > 3) {
        TwoDface = project2dface((function() {
          var _j, _len2, _results;
          _results = [];
          for (_j = 0, _len2 = f.length; _j < _len2; _j++) {
            v = f[_j];
            _results.push(poly.xyz[v]);
          }
          return _results;
        })());
        diags = getDiagonals(TwoDface);
        tris = diagsToTris(f, diags);
        _ref3 = enumerate(tris);
        for (_j = 0, _len2 = _ref3.length; _j < _len2; _j++) {
          _ref4 = _ref3[_j], j = _ref4[0], tri = _ref4[1];
          newpoly.face.push([f[tri[0]], f[tri[1]], f[tri[2]]]);
          if (colors) {
            newpoly.face_colors.push(poly.face_colors[i]);
          }
        }
      } else {
        newpoly.face.push([f[0], f[1], f[2]]);
        if (colors) {
          newpoly.face_colors.push(poly.face_colors[i]);
        }
      }
    }
    newpoly.name = poly.name;
    return newpoly;
  };
  topolog = function(poly) {
    var f, str, v, _i, _j, _len, _len2, _ref;
    str = "";
    _ref = poly.face;
    for (_i = 0, _len = _ref.length; _i < _len; _i++) {
      f = _ref[_i];
      for (_j = 0, _len2 = f.length; _j < _len2; _j++) {
        v = f[_j];
        str += "" + v + "->";
      }
      str += "\n";
    }
    return console.log(str);
  };
  testrig = function() {
    var op, ops, seed, seeds, _i, _j, _len, _len2;
    seeds = ["T", "O", "C", "I", "D", "P3", "P4", "A4", "A5", "Y3", "Y4"];
    ops = ["k", "a", "g", "p", "d", "n", "x", "*"];
    console.log("===== Test Basic Ops =====");
    for (_i = 0, _len = ops.length; _i < _len; _i++) {
      op = ops[_i];
      console.log("Operator " + op);
      for (_j = 0, _len2 = seeds.length; _j < _len2; _j++) {
        seed = seeds[_j];
        console.log(op + seed + ":", generatePoly(op + seed));
      }
    }
    return console.log("===== Done Testing Basic Ops =====");
  };
  arrayToCSGVertex = function(pos, norm) {
    var csgvec;
    pos || (pos = [0, 0, 0]);
    norm || (norm = unit(pos));
    csgvec = new CSG.Vertex(new CSG.Vector(pos[0], pos[1], pos[2]), new CSG.Vector(norm[0], norm[1], norm[2]));
    return csgvec;
  };
  polyToCSG = function(poly) {
    var CSGPolys, csg, f, v, _i, _len, _ref;
    CSGPolys = [];
    _ref = poly.face;
    for (_i = 0, _len = _ref.length; _i < _len; _i++) {
      f = _ref[_i];
      CSGPolys.push(new CSG.Polygon((function() {
        var _j, _len2, _results;
        _results = [];
        for (_j = 0, _len2 = f.length; _j < _len2; _j++) {
          v = f[_j];
          _results.push(arrayToCSGVertex(poly.xyz[v]));
        }
        return _results;
      })(), {}));
    }
    csg = CSG.fromPolygons(CSGPolys);
    return csg;
  };
  vertIndex = function(vertlist, vert) {
    var i, v, _i, _len, _ref, _ref2;
    _ref = enumerate(vertlist);
    for (_i = 0, _len = _ref.length; _i < _len; _i++) {
      _ref2 = _ref[_i], i = _ref2[0], v = _ref2[1];
      if (mag(sub(vert, v)) < 1E-6) {
        return i;
      }
    }
    return -1;
  };
  uniquifyverts = function(vertlist) {
    var already_present, u, uniqverts, v, _i, _j, _len, _len2;
    uniqverts = [];
    uniqverts.push(vertlist[0]);
    for (_i = 0, _len = vertlist.length; _i < _len; _i++) {
      v = vertlist[_i];
      already_present = false;
      for (_j = 0, _len2 = uniqverts.length; _j < _len2; _j++) {
        u = uniqverts[_j];
        if (mag(sub(v, u)) < 1E-6) {
          already_present = true;
        }
      }
      if (!already_present) {
        uniqverts.push(v);
      }
    }
    return uniqverts;
  };
  edge_vert_isect = function(evert1, evert2, testvert) {
    var a, b, dT, dv, _ref;
    dv = sub(evert2, evert1);
    dT = sub(testvert, evert1);
    a = dot(unit(dv), dT);
    b = mag(sub(unit(dv), unit(dT)));
    if (b < 1E-6 && (0 < (_ref = mag(dT)) && _ref < mag(dv))) {
      return true;
    } else {
      return false;
    }
  };
  meshFix = function(poly) {
    var edge0, edge1, f, faces, newf, newfaces, newpoly, testv, v, verts, vno, _i, _j, _k, _len, _len2, _len3, _ref, _ref2;
    faces = poly.face;
    verts = poly.xyz;
    newfaces = [];
    for (_i = 0, _len = faces.length; _i < _len; _i++) {
      f = faces[_i];
      newf = [];
      edge0 = verts[f[f.length - 1]];
      for (_j = 0, _len2 = f.length; _j < _len2; _j++) {
        v = f[_j];
        edge1 = verts[v];
        _ref = enumerate(verts);
        for (_k = 0, _len3 = _ref.length; _k < _len3; _k++) {
          _ref2 = _ref[_k], vno = _ref2[0], testv = _ref2[1];
          if (edge_vert_isect(edge0, edge1, testv)) {
            newf.push(vno);
          }
        }
        newf.push(v);
        edge0 = edge1;
      }
      newfaces.push(newf);
    }
    newpoly = new polyhedron();
    newpoly.face = newfaces;
    newpoly.xyz = clone(verts);
    return newpoly;
  };
  CSGToPoly = function(csgpoly) {
    var csg, faces, poly, v, verts, _i, _j, _k, _len, _len2, _len3, _ref, _ref2, _ref3;
    faces = [];
    verts = [];
    _ref = csgpoly.toPolygons();
    for (_i = 0, _len = _ref.length; _i < _len; _i++) {
      csg = _ref[_i];
      _ref2 = csg.vertices;
      for (_j = 0, _len2 = _ref2.length; _j < _len2; _j++) {
        v = _ref2[_j];
        verts.push([v.pos.x, v.pos.y, v.pos.z]);
      }
    }
    verts = uniquifyverts(verts);
    _ref3 = csgpoly.toPolygons();
    for (_k = 0, _len3 = _ref3.length; _k < _len3; _k++) {
      csg = _ref3[_k];
      faces.push((function() {
        var _l, _len4, _ref4, _results;
        _ref4 = csg.vertices;
        _results = [];
        for (_l = 0, _len4 = _ref4.length; _l < _len4; _l++) {
          v = _ref4[_l];
          _results.push(vertIndex(verts, [v.pos.x, v.pos.y, v.pos.z]));
        }
        return _results;
      })());
    }
    poly = new polyhedron();
    poly.face = faces;
    poly.xyz = verts;
    return meshFix(poly);
  };
  csgUnion = function(polyA, polyB) {
    var csgA, csgB, csgresult;
    csgA = polyToCSG(polyA);
    csgB = polyToCSG(polyB);
    csgresult = csgA.union(csgB);
    return CSGToPoly(csgresult);
  };
  csgIntersect = function(polyA, polyB) {
    var csgA, csgB, csgresult;
    csgA = polyToCSG(polyA);
    csgB = polyToCSG(polyB);
    csgresult = csgA.intersect(csgB);
    return CSGToPoly(csgresult);
  };
  csgSubtract = function(polyA, polyB) {
    var csgA, csgB, csgresult;
    csgA = polyToCSG(polyA);
    csgB = polyToCSG(polyB);
    csgresult = csgA.subtract(csgB);
    return CSGToPoly(csgresult);
  };
  edgesToFace = function(edges) {
    var e, edict, face, itCTR, v, v0, _i, _len;
    face = [];
    edict = {};
    for (_i = 0, _len = edges.length; _i < _len; _i++) {
      e = edges[_i];
      edict[e[0]] = e[1];
    }
    v0 = edges[0][0];
    v = edges[0][1];
    face.push(v);
    itCTR = 0;
    while (v !== v0) {
      v = edict[v];
      face.push(v);
      itCTR++;
      if (itCTR > 1000) {
        console.log("Bad edges to face join, have a neverending face:", edges);
        break;
      }
    }
    return face;
  };
  uniquifyedges = function(edgelist) {
    var already_present, u, uniqverts, v, _i, _j, _len, _len2;
    uniqverts = [];
    uniqverts.push(vertlist[0]);
    for (_i = 0, _len = vertlist.length; _i < _len; _i++) {
      v = vertlist[_i];
      already_present = false;
      for (_j = 0, _len2 = uniqverts.length; _j < _len2; _j++) {
        u = uniqverts[_j];
        if (mag(sub(v, u)) < 1E-6) {
          already_present = true;
        }
      }
      if (!already_present) {
        uniqverts.push(v);
      }
    }
    return uniqverts;
  };
  faceprint = function(face) {
    var str, v, _i, _len;
    str = "";
    for (_i = 0, _len = face.length; _i < _len; _i++) {
      v = face[_i];
      str += "" + v + "->";
    }
    return str.slice(0, (str.length - 3 + 1) || 9e9);
  };
  joinFaces = function(faces) {
    var alledges, e, e2, i, newedges, nuked, _i, _len, _ref, _ref2;
    alledges = [];
    alledges = (_.reduce(_.map(faces, faceToEdges), (function(memo, x) {
      return memo.concat(x);
    }), alledges)).slice(0);
    newedges = [];
    while (alledges.length > 0) {
      e = alledges.pop();
      nuked = false;
      if (alledges.length > 0) {
        _ref = enumerate(alledges);
        for (_i = 0, _len = _ref.length; _i < _len; _i++) {
          _ref2 = _ref[_i], i = _ref2[0], e2 = _ref2[1];
          if (e[0] === e2[1] && e[1] === e2[0]) {
            nuked = true;
            alledges.splice(i, 1);
          }
        }
        if (!nuked) {
          newedges.push(e);
        }
      } else {
        newedges.push(e);
      }
    }
    return edgesToFace(newedges);
  };
  uniteFaces = function(poly) {
    var connected_set, dpoly, e, edict, face_idx, face_set, face_sets, faces_to_join, hash, idx, isNeighbor, j, joinset, joinsets, k, ndict, newfaces, newpoly, normals, normhash, t, v, _i, _j, _k, _l, _len, _len2, _len3, _len4, _len5, _len6, _len7, _m, _n, _o, _p, _ref, _ref2, _ref3, _ref4, _ref5, _results;
    dpoly = dual(poly);
    edict = {};
    _ref = dpoly.getEdges();
    for (_i = 0, _len = _ref.length; _i < _len; _i++) {
      e = _ref[_i];
      edict[e[0]] = [];
    }
    _ref2 = dpoly.getEdges();
    for (_j = 0, _len2 = _ref2.length; _j < _len2; _j++) {
      e = _ref2[_j];
      edict[e[0]].push(e[1]);
      edict[e[1]].push(e[0]);
    }
    isNeighbor = function(Fi, Fi2) {
      if (edict[Fi].indexOf(Fi2) !== -1) {
        return true;
      } else {
        return false;
      }
    };
    normals = poly.normals();
    normhash = function(norm) {
      return round(100 * norm[0]) + "~" + round(100 * norm[1]) + "~" + round(100 * norm[2]);
    };
    face_sets = [];
    ndict = {};
    face_idx = (function() {
      _results = [];
      for (var _k = 0, _ref3 = poly.face.length - 1; 0 <= _ref3 ? _k <= _ref3 : _k >= _ref3; 0 <= _ref3 ? _k++ : _k--){ _results.push(_k); }
      return _results;
    }).apply(this);
    for (_l = 0, _len3 = face_idx.length; _l < _len3; _l++) {
      idx = face_idx[_l];
      hash = normhash(normals[idx]);
      if (ndict[hash]) {
        ndict[hash].push(idx);
      } else {
        ndict[hash] = [idx];
      }
    }
    console.log("ndict", ndict);
    for (k in ndict) {
      v = ndict[k];
      face_sets.push(clone(v));
      console.log(clone(v));
    }
    joinsets = [];
    for (_m = 0, _len4 = face_sets.length; _m < _len4; _m++) {
      face_set = face_sets[_m];
      connected_set = [];
      while (face_set.length > 0) {
        connected_set.push(0);
        for (j = 0, _ref4 = face_set.length - 1; 0 <= _ref4 ? j <= _ref4 : j >= _ref4; 0 <= _ref4 ? j++ : j--) {
          if (__indexOf.call(connected_set, j) >= 0) {
            continue;
          }
          for (_n = 0, _len5 = connected_set.length; _n < _len5; _n++) {
            t = connected_set[_n];
            if (isNeighbor(face_set[j], face_set[t])) {
              connected_set.push(j);
              j = 0;
              break;
            }
          }
        }
        console.log(connected_set.sort((function(a, b) {
          return a - b;
        })).reverse(), clone(face_set));
        joinset = [];
        _ref5 = connected_set.sort((function(a, b) {
          return a - b;
        })).reverse();
        for (_o = 0, _len6 = _ref5.length; _o < _len6; _o++) {
          t = _ref5[_o];
          joinset.push(face_set.splice(t, 1)[0]);
        }
        joinsets.push(joinset);
        connected_set = [];
      }
    }
    for (_p = 0, _len7 = joinsets.length; _p < _len7; _p++) {
      j = joinsets[_p];
      console.log(j);
    }
    faces_to_join = _.map(joinsets, (function(x) {
      return _.map(x, function(i) {
        return poly.face[i];
      });
    }));
    newfaces = _.map(faces_to_join, joinFaces);
    newpoly = new polyhedron();
    newpoly.xyz = clone(poly.xyz);
    newpoly.face = newfaces;
    console.log(newpoly);
    return newpoly;
  };
  ctx = {};
  CANVAS_WIDTH = 500;
  CANVAS_HEIGHT = 400;
  globPolys = {};
  globRotM = clone(eye3);
  globlastRotM = clone(eye3);
  perspective_scale = 800;
  persp_z_max = 5;
  persp_z_min = 0;
  persp_ratio = 0.8;
  _2d_x_offset = CANVAS_WIDTH / 2;
  _2d_y_offset = CANVAS_HEIGHT / 2;
  globtime = new Date();
  BG_CLEAR = true;
  BG_COLOR = "rgba(255,255,255,1.0)";
  COLOR_METHOD = "area";
  PaintMode = "fillstroke";
  ctx_linewidth = 0.5;
  MOUSEDOWN = false;
  LastMouseX = 0;
  LastMouseY = 0;
  LastSphVec = [1, 0, 0];
  DEFAULT_RECIPES = ["C2dakD", "oC20kkkT", "kn4C40A0dA4", "opD", "lT", "lK5oC", "knD", "dn6x4K5bT", "oox4P7", "n18n18n9n9n9soxY9"];
  saveText = function(text, filename) {
    var BB, bb;
    BB = window.BlobBuilder || window.WebKitBlobBuilder || window.MozBlobBuilder;
    bb = new BB();
    bb.append(text);
    return saveAs(bb.getBlob("text/plain;charset=" + document.characterSet), filename);
  };
  old_centers = [];
  old_normals = [];
  def_palette = ["#ff3333", "#33ff33", "#3333ff", "#ffff33", "#ff33ff", "#33ffff", "#dddddd", "#555555", "#dd0000", "#00dd00", "#0000dd"];
  rwb_palette = ["#ff8888", "#dddddd", "#777777", "#aa3333", "#ff0000", "#ffffff", "#aaaaaa"];
  rwbg_palette = ["#ff8888", "#ffeeee", "#88ff88", "#dd7777", "#ff2222", "#22ff22", "#ee4422", "#aaaaaa"];
  hextofloats = function(hexstr) {
    var rgb;
    if (hexstr[0] === "#") {
      hexstr = hexstr.slice(1);
    }
    if (hexstr.length === 3) {
      rgb = hexstr.split('').map(function(c) {
        return parseInt(c + c, 16) / 255;
      });
    } else {
      rgb = hexstr.match(/.{2}/g).map(function(c) {
        return parseInt(c, 16) / 255;
      });
    }
    return rgb;
  };
  PALETTE = rwb_palette;
  palette = function(n) {
    if (n < PALETTE.length) {
      return hextofloats(PALETTE[n]);
    } else {
      return hextofloats(PALETTE[PALETTE.length - 1]);
    }
  };
  colorassign = function(ar, colormemory) {
    var fclr, hash;
    hash = round(100 * ar);
    if (hash in colormemory) {
      return colormemory[hash];
    } else {
      fclr = palette(_.toArray(colormemory).length);
      colormemory[hash] = fclr;
      return fclr;
    }
  };
  paintPolyhedron = function(poly) {
    var clr, colormemory, f, face_verts, v, _i, _len, _ref;
    poly.face_colors = [];
    colormemory = {};
    _ref = poly.face;
    for (_i = 0, _len = _ref.length; _i < _len; _i++) {
      f = _ref[_i];
      if (COLOR_METHOD === "area") {
        face_verts = (function() {
          var _j, _len2, _results;
          _results = [];
          for (_j = 0, _len2 = f.length; _j < _len2; _j++) {
            v = f[_j];
            _results.push(poly.xyz[v]);
          }
          return _results;
        })();
        clr = colorassign(convexarea(face_verts), colormemory);
      } else {
        clr = palette(f.length - 3);
      }
      poly.face_colors.push(clr);
    }
    return poly;
  };
  specreplacements = [[/e/g, "aa"], [/b/g, "ta"], [/o/g, "jj"], [/m/g, "kj"], [/t(\d*)/g, "dk$1d"], [/j/g, "dad"], [/s/g, "dgd"], [/dd/g, ""], [/ad/g, "a"], [/gd/g, "g"], [/aO/g, "aC"], [/aI/g, "aD"], [/gO/g, "gC"], [/gI/g, "gD"]];
  getOps = function(notation) {
    var equiv, expanded, orig, _i, _len, _ref;
    expanded = notation;
    for (_i = 0, _len = specreplacements.length; _i < _len; _i++) {
      _ref = specreplacements[_i], orig = _ref[0], equiv = _ref[1];
      expanded = expanded.replace(orig, equiv);
    }
    console.log("" + notation + " executed as " + expanded);
    return expanded;
  };
  generatePoly = function(notation) {
    var n, ops, poly;
    poly = new polyhedron();
    n = 0;
    ops = getOps(notation);
    if (ops.search(/([0-9]+)$/) !== -1) {
      n = 1 * RegExp.lastParen;
      ops = ops.slice(0, -RegExp.lastParen.length);
    }
    switch (ops.slice(-1)) {
      case "T":
        poly = tetrahedron();
        break;
      case "O":
        poly = octahedron();
        break;
      case "C":
        poly = cube();
        break;
      case "I":
        poly = icosahedron();
        break;
      case "D":
        poly = dodecahedron();
        break;
      case "P":
        poly = prism(n);
        break;
      case "A":
        poly = antiprism(n);
        break;
      case "Y":
        poly = pyramid(n);
        break;
      default:
        return;
    }
    while (ops !== "") {
      n = 0;
      if (ops.search(/([0-9]+)$/) !== -1) {
        n = 1 * RegExp.lastParen;
        ops = ops.slice(0, -RegExp.lastParen.length);
      }
      switch (ops.slice(-1)) {
        case "d":
          poly = dual(poly);
          break;
        case "k":
          poly = kisN(poly, n);
          break;
        case "a":
          poly = ambo(poly);
          break;
        case "g":
          poly = gyro(poly);
          break;
        case "p":
          poly = propellor(poly);
          break;
        case "r":
          poly = reflect(poly);
          break;
        case "K":
          poly.xyz = canonicalXYZ(poly, n === 0 ? 1 : n);
          break;
        case "C":
          poly.xyz = canonicalize(poly, n === 0 ? 1 : n);
          break;
        case "A":
          poly.xyz = adjustXYZ(poly, n === 0 ? 1 : n);
          break;
        case "n":
          poly = insetN(poly, n);
          break;
        case "x":
          poly = extrudeN(poly, n);
          break;
        case "l":
          poly = stellaN(poly, n);
          break;
        case "z":
          poly = triangulate(poly, false);
      }
      ops = ops.slice(0, -1);
    }
    poly.xyz = recenter(poly.xyz, poly.getEdges());
    poly.xyz = rescale(poly.xyz);
    poly = paintPolyhedron(poly);
    return poly;
  };
  parseurl = function() {
    var a, d, e, q, r, urlParams;
    urlParams = {};
    a = /\+/g;
    r = /([^&=]+)=?([^&]*)/g;
    d = function(s) {
      return decodeURIComponent(s.replace(a, " "));
    };
    q = window.location.search.substring(1);
    while (e = r.exec(q)) {
      urlParams[d(e[1])] = d(e[2]);
    }
    return urlParams;
  };
  init = function() {
    var canvas;
    canvas = $('#poly');
    canvas.width(CANVAS_WIDTH);
    canvas.height(CANVAS_HEIGHT);
    ctx = canvas[0].getContext("2d");
    ctx.lineWidth = ctx_linewidth;
    if (BG_CLEAR) {
      return ctx.clearRect(0, 0, CANVAS_WIDTH, CANVAS_HEIGHT);
    } else {
      ctx.clearRect(0, 0, CANVAS_WIDTH, CANVAS_HEIGHT);
      ctx.fillStyle = BG_COLOR;
      return ctx.fillRect(0, 0, CANVAS_WIDTH, CANVAS_HEIGHT);
    }
  };
  clear = function() {
    if (BG_CLEAR) {
      return ctx.clearRect(0, 0, CANVAS_WIDTH, CANVAS_HEIGHT);
    } else {
      ctx.clearRect(0, 0, CANVAS_WIDTH, CANVAS_HEIGHT);
      ctx.fillStyle = BG_COLOR;
      return ctx.fillRect(0, 0, CANVAS_WIDTH, CANVAS_HEIGHT);
    }
  };
  sortfaces = function(poly) {
    var centroids, idx, normals, planesort, ray_origin, zcentroidsort, zsortIndex, _i, _ref, _results;
    centroids = poly.centers();
    normals = poly.normals();
    ray_origin = [0, 0, (persp_z_max * persp_ratio - persp_z_min) / (1 - persp_ratio)];
    planesort = function(a, b) {
      return -dot(sub(ray_origin, a[0]), a[1]) * dot(sub(b[0], a[0]), a[1]);
    };
    zcentroidsort = function(a, b) {
      return a[0][2] - b[0][2];
    };
    zsortIndex = _.zip(centroids, normals, (function() {
      _results = [];
      for (var _i = 0, _ref = poly.face.length - 1; 0 <= _ref ? _i <= _ref : _i >= _ref; 0 <= _ref ? _i++ : _i--){ _results.push(_i); }
      return _results;
    }).apply(this)).sort(zcentroidsort).map(function(x) {
      return x[2];
    });
    poly.face = (function() {
      var _j, _len, _results2;
      _results2 = [];
      for (_j = 0, _len = zsortIndex.length; _j < _len; _j++) {
        idx = zsortIndex[_j];
        _results2.push(poly.face[idx]);
      }
      return _results2;
    })();
    return poly.face_colors = (function() {
      var _j, _len, _results2;
      _results2 = [];
      for (_j = 0, _len = zsortIndex.length; _j < _len; _j++) {
        idx = zsortIndex[_j];
        _results2.push(poly.face_colors[idx]);
      }
      return _results2;
    })();
  };
  drawpoly = function(poly, tvec) {
    var centers, clr, face, face_labels, face_verts, fno, illum, normals, oldxyz, v, v0, vert, vert_labels, vno, x, y, _i, _j, _k, _l, _len, _len2, _len3, _len4, _len5, _m, _ref, _ref10, _ref11, _ref12, _ref13, _ref14, _ref2, _ref3, _ref4, _ref5, _ref6, _ref7, _ref8, _ref9;
    tvec || (tvec = [3, 3, 3]);
    oldxyz = _.map(poly.xyz, function(x) {
      return x;
    });
    poly.xyz = _.map(poly.xyz, function(x) {
      return mv3(globRotM, x);
    });
    sortfaces(poly);
    _ref = enumerate(poly.face);
    for (_i = 0, _len = _ref.length; _i < _len; _i++) {
      _ref2 = _ref[_i], fno = _ref2[0], face = _ref2[1];
      ctx.beginPath();
      v0 = face[face.length - 1];
      _ref3 = perspT(add(tvec, poly.xyz[v0]), persp_z_max, persp_z_min, persp_ratio, perspective_scale), x = _ref3[0], y = _ref3[1];
      ctx.moveTo(x + _2d_x_offset, y + _2d_y_offset);
      for (_j = 0, _len2 = face.length; _j < _len2; _j++) {
        v = face[_j];
        _ref4 = perspT(add(tvec, poly.xyz[v]), persp_z_max, persp_z_min, persp_ratio, perspective_scale), x = _ref4[0], y = _ref4[1];
        ctx.lineTo(x + _2d_x_offset, y + _2d_y_offset);
      }
      clr = poly.face_colors[fno];
      face_verts = (function() {
        var _k, _len3, _results;
        _results = [];
        for (_k = 0, _len3 = face.length; _k < _len3; _k++) {
          v = face[_k];
          _results.push(poly.xyz[v]);
        }
        return _results;
      })();
      illum = dot(normal(face_verts), unit([1, -1, 0]));
      clr = mult((illum / 2.0 + .5) * 0.7 + 0.3, clr);
      if (PaintMode === "fill" || PaintMode === "fillstroke") {
        ctx.fillStyle = "rgba(" + (round(clr[0] * 255)) + ", " + (round(clr[1] * 255)) + ", " + (round(clr[2] * 255)) + ", " + 1.0 + ")";
        ctx.fill();
        ctx.strokeStyle = "rgba(" + (round(clr[0] * 255)) + ", " + (round(clr[1] * 255)) + ", " + (round(clr[2] * 255)) + ", " + 1.0 + ")";
        ctx.stroke();
      }
      if (PaintMode === "fillstroke") {
        ctx.fillStyle = "rgba(" + (round(clr[0] * 255)) + ", " + (round(clr[1] * 255)) + ", " + (round(clr[2] * 255)) + ", " + 1.0 + ")";
        ctx.fill();
        ctx.strokeStyle = "rgba(0,0,0, .3)";
        ctx.stroke();
      }
      if (PaintMode === "stroke") {
        ctx.strokeStyle = "rgba(0,0,0, .8)";
        ctx.stroke();
      }
    }
    centers = _.map(old_centers, function(x) {
      return mv3(globRotM, x);
    });
    normals = _.map(old_normals, function(x) {
      return mv3(globRotM, x);
    });
    face_labels = (function() {
      var _ref5, _results;
      _results = [];
      for (fno = 0, _ref5 = centers.length - 1; 0 <= _ref5 ? fno <= _ref5 : fno >= _ref5; 0 <= _ref5 ? fno++ : fno--) {
        _results.push("" + fno);
      }
      return _results;
    })();
    vert_labels = (function() {
      var _ref5, _results;
      _results = [];
      for (vno = 0, _ref5 = poly.xyz.length - 1; 0 <= _ref5 ? vno <= _ref5 : vno >= _ref5; 0 <= _ref5 ? vno++ : vno--) {
        _results.push("" + vno);
      }
      return _results;
    })();
    _ref5 = enumerate(poly.face);
    for (_k = 0, _len3 = _ref5.length; _k < _len3; _k++) {
      _ref6 = _ref5[_k], fno = _ref6[0], face = _ref6[1];
      ctx.textAlign = "center";
      ctx.fillStyle = "rgba(0,0,0,1)";
      _ref7 = perspT(add(tvec, centers[fno]), persp_z_max, persp_z_min, persp_ratio, perspective_scale), x = _ref7[0], y = _ref7[1];
      ctx.fillText(face_labels[fno], x + _2d_x_offset, y + _2d_y_offset);
    }
    _ref8 = enumerate(poly.face);
    for (_l = 0, _len4 = _ref8.length; _l < _len4; _l++) {
      _ref9 = _ref8[_l], fno = _ref9[0], face = _ref9[1];
      ctx.strokeStyle = "rgba(0,0,0,1)";
      ctx.beginPath();
      _ref10 = perspT(add(tvec, centers[fno]), persp_z_max, persp_z_min, persp_ratio, perspective_scale), x = _ref10[0], y = _ref10[1];
      ctx.moveTo(x + _2d_x_offset, y + _2d_y_offset);
      _ref11 = perspT(add(tvec, add(centers[fno], mult(0.1, normals[fno]))), persp_z_max, persp_z_min, persp_ratio, perspective_scale), x = _ref11[0], y = _ref11[1];
      ctx.lineTo(x + _2d_x_offset, y + _2d_y_offset);
      ctx.stroke();
    }
    _ref12 = enumerate(poly.xyz);
    for (_m = 0, _len5 = _ref12.length; _m < _len5; _m++) {
      _ref13 = _ref12[_m], vno = _ref13[0], vert = _ref13[1];
      _ref14 = perspT(add(tvec, poly.xyz[vno]), persp_z_max, persp_z_min, persp_ratio, perspective_scale), x = _ref14[0], y = _ref14[1];
    }
    return poly.xyz = oldxyz;
  };
  drawShape = function() {
    var i, p, _i, _len, _ref, _ref2, _results;
    clear();
    _ref = enumerate(globPolys);
    _results = [];
    for (_i = 0, _len = _ref.length; _i < _len; _i++) {
      _ref2 = _ref[_i], i = _ref2[0], p = _ref2[1];
      _results.push(drawpoly(p, [0 + 3 * i, 0, 3]));
    }
    return _results;
  };
  animateShape = function() {
    var globtheta, i, p, _i, _len, _ref, _ref2;
    clear();
    globtheta = (2 * Math.PI) / 180.0 * globtime.getSeconds() * 0.1;
    _ref = enumerate(globPolys);
    for (_i = 0, _len = _ref.length; _i < _len; _i++) {
      _ref2 = _ref[_i], i = _ref2[0], p = _ref2[1];
      drawpoly(p, [0 + 3 * i, 0, 3]);
    }
    return setTimeout(animateShape, 100);
  };
  $(function() {
    var c0, c1, c2, c4, specs, urlParams, v, _i, _len, _ref;
    init();
    urlParams = parseurl();
    if ("recipe" in urlParams) {
      specs = [urlParams["recipe"]];
      $("#spec").val(specs);
    } else {
      specs = [randomchoice(DEFAULT_RECIPES)];
      $("#spec").val(specs);
    }
    c0 = cube();
    c1 = cube();
    c2 = cube().rescale(.5);
    _ref = c0.xyz;
    for (_i = 0, _len = _ref.length; _i < _len; _i++) {
      v = _ref[_i];
      c2.translate(mult(1, v));
      c1 = csgUnion(c1, c2);
      c2.translate(mult(-1, v));
    }
    c4 = c1;
    globPolys = [paintPolyhedron(c4)];
    old_centers = _.map(globPolys[0].centers(), function(x) {
      return x;
    });
    old_normals = _.map(globPolys[0].normals(), function(x) {
      return x;
    });
    drawShape();
    try {
      uniteFaces(clone(c4));
    } catch (_e) {}
    $("#spec").change(function(e) {
      specs = $("#spec").val().split(/\s+/g).slice(0, 2);
      globPolys = _.map(specs, function(x) {
        return generatePoly(x);
      });
      return drawShape();
    });
    $("#poly").mousewheel(function(e, delta, deltaX, deltaY) {
      e.preventDefault();
      perspective_scale *= (10 + delta) / 10;
      return drawShape();
    });
    $("#poly").mousedown(function(e) {
      var tmpvec;
      e.preventDefault();
      MOUSEDOWN = true;
      LastMouseX = e.clientX - $(this).offset().left;
      LastMouseY = e.clientY - $(this).offset().top;
      tmpvec = invperspT(LastMouseX, LastMouseY, _2d_x_offset, _2d_y_offset, persp_z_max, persp_z_min, persp_ratio, perspective_scale);
      if (tmpvec[0] * tmpvec[1] * tmpvec[2] * 0 === 0) {
        LastSphVec = tmpvec;
      }
      return globlastRotM = clone(globRotM);
    });
    $("#poly").mouseup(function(e) {
      e.preventDefault();
      return MOUSEDOWN = false;
    });
    $("#poly").mouseleave(function(e) {
      e.preventDefault();
      return MOUSEDOWN = false;
    });
    $("#poly").mousemove(function(e) {
      var MouseX, MouseY, SphVec;
      e.preventDefault();
      if (MOUSEDOWN) {
        MouseX = e.clientX - $(this).offset().left;
        MouseY = e.clientY - $(this).offset().top;
        SphVec = invperspT(MouseX, MouseY, _2d_x_offset, _2d_y_offset, persp_z_max, persp_z_min, persp_ratio, perspective_scale);
        if (SphVec[0] * SphVec[1] * SphVec[2] * 0 === 0 && LastSphVec[0] * LastSphVec[1] * LastSphVec[2] * 0 === 0) {
          globRotM = mm3(getVec2VecRotM(LastSphVec, SphVec), globlastRotM);
        }
        return drawShape();
      }
    });
    $("#strokeonly").click(function(e) {
      PaintMode = "stroke";
      return drawShape();
    });
    $("#fillonly").click(function(e) {
      PaintMode = "fill";
      return drawShape();
    });
    $("#fillandstroke").click(function(e) {
      PaintMode = "fillstroke";
      return drawShape();
    });
    $("#siderot").click(function(e) {
      globRotM = vec_rotm(PI / 2, 0, 1, 0);
      return drawShape();
    });
    $("#toprot").click(function(e) {
      globRotM = vec_rotm(PI / 2, 1, 0, 0);
      return drawShape();
    });
    $("#frontrot").click(function(e) {
      globRotM = rotm(0, 0, 0);
      return drawShape();
    });
    $("#pngsavebutton").click(function(e) {
      var canvas, filename, spec;
      canvas = $("#poly")[0];
      spec = $("#spec").val().split(/\s+/g)[0];
      filename = "polyhedronisme-" + spec + ".png";
      return canvas.toBlob(function(blob) {
        return saveAs(blob, filename);
      });
    });
    $("#objsavebutton").click(function(e) {
      var filename, objtxt, spec;
      objtxt = globPolys[0].toOBJ();
      spec = $("#spec").val().split(/\s+/g)[0];
      filename = "polyhedronisme-" + spec + ".obj";
      return saveText(objtxt, filename);
    });
    return $("#x3dsavebutton").click(function(e) {
      var filename, spec, triangulated, x3dtxt;
      triangulated = triangulate(globPolys[0], true);
      x3dtxt = triangulated.toVRML();
      spec = $("#spec").val().split(/\s+/g)[0];
      filename = "polyhedronisme-" + spec + ".wrl";
      return saveText(x3dtxt, filename);
    });
  });
}).call(this);
