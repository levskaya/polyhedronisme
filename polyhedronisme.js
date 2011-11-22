(function() {
  var BG_CLEAR, BG_COLOR, CANVAS_HEIGHT, CANVAS_WIDTH, COLOR_METHOD, LastMouseX, LastMouseY, MOUSEDOWN, PALETTE, PI, abs, add, adjustXYZ, ambo, animateShape, antiprism, canonicalXYZ, canonicalize, centroid, clear, colorassign, convexarea, copyVecArray, cos, cross, ctx, ctx_linewidth, cube, def_palette, dodecahedron, dot, drawShape, drawpoly, dual, edgeDist, enumerate, faceCenters, faceToEdges, floor, generatePoly, getOps, globPolys, globphi, globtheta, globtime, gyro, hextofloats, icosahedron, init, intersect, kisN, mag, mag2, midName, midpoint, mm3, mult, mv3, normal, octahedron, oneThird, orthogonal, paintPolyhedron, palette, parseurl, perspT, persp_ratio, persp_z_max, persp_z_min, perspective_scale, planarize, polyflag, polyhedron, pow, prism, propellor, pyramid, random, randomchoice, recenter, reciprocal, reciprocalC, reciprocalN, reflect, rotm, round, rwb_palette, rwbg_palette, sin, sortfaces, specreplacements, sqrt, sub, tan, tangentPoint, tangentify, testrig, tetrahedron, toOBJ, topolog, tween, unit, _2d_x_offset, _2d_y_offset, _mult;
  var __bind = function(fn, me){ return function(){ return fn.apply(me, arguments); }; };
  random = Math.random;
  round = Math.round;
  floor = Math.floor;
  sqrt = Math.sqrt;
  sin = Math.sin;
  cos = Math.cos;
  tan = Math.tan;
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
  mag = function(vec) {
    return sqrt(dot(vec, vec));
  };
  mag2 = function(vec) {
    return dot(vec, vec);
  };
  unit = function(vec) {
    return mult(1 / sqrt(mag2(vec)), vec);
  };
  cross = function(d1, d2) {
    return [d1[1] * d2[2] - d1[2] * d2[1], d1[2] * d2[0] - d1[0] * d2[2], d1[0] * d2[1] - d1[1] * d2[0]];
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
  rotm = function(phi, theta, psi) {
    var xy_mat, xz_mat, yz_mat;
    xy_mat = [[cos(phi), -1.0 * sin(phi), 0.0], [sin(phi), cos(phi), 0.0], [0.0, 0.0, 1.0]];
    yz_mat = [[cos(theta), 0, -1.0 * sin(theta)], [0, 1, 0], [sin(theta), 0, cos(theta)]];
    xz_mat = [[1.0, 0, 0], [0, cos(psi), -1.0 * sin(psi)], [0, sin(psi), cos(psi)]];
    return mm3(xz_mat, mm3(yz_mat, xy_mat));
  };
  perspT = function(vec3, max_real_depth, min_real_depth, desired_ratio, desired_length) {
    var scalefactor, z0;
    z0 = (max_real_depth * desired_ratio - min_real_depth) / (1 - desired_ratio);
    scalefactor = desired_length * desired_ratio / (1 - desired_ratio);
    return [scalefactor * vec3[0] / (vec3[2] + z0), scalefactor * vec3[1] / (vec3[2] + z0)];
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
      alledges = _(this.face).map(faceToEdges);
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
    return polyhedron;
  })();
  toOBJ = function(poly) {
    var f, i, norm, objstr, v, _i, _j, _k, _l, _len, _len2, _len3, _len4, _ref, _ref2, _ref3, _ref4;
    objstr = "#Produced by polyHédronisme http://levskaya.github.com/polyhedronisme\n";
    objstr += "group poly\n";
    objstr += "#vertices\n";
    _ref = poly.xyz;
    for (_i = 0, _len = _ref.length; _i < _len; _i++) {
      v = _ref[_i];
      objstr += "v " + v[0] + " " + v[1] + " " + v[2] + "\n";
    }
    objstr += "#normal vector defs \n";
    _ref2 = poly.face;
    for (_j = 0, _len2 = _ref2.length; _j < _len2; _j++) {
      f = _ref2[_j];
      norm = normal((function() {
        var _k, _len3, _results;
        _results = [];
        for (_k = 0, _len3 = f.length; _k < _len3; _k++) {
          v = f[_k];
          _results.push(poly.xyz[v]);
        }
        return _results;
      })());
      objstr += "vn " + norm[0] + " " + norm[1] + " " + norm[2] + "\n";
    }
    objstr += "#face defs \n";
    _ref3 = enumerate(poly.face);
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
      poly.xyz.push([cos(i * theta), sin(i * theta), h]);
    }
    for (i = 0, _ref2 = n - 1; 0 <= _ref2 ? i <= _ref2 : i >= _ref2; 0 <= _ref2 ? i++ : i--) {
      poly.xyz.push([cos(i * theta), sin(i * theta), -h]);
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
    r = r / f;
    h = h / f;
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
    var i, poly, theta, _i, _ref, _ref2, _ref3, _results;
    theta = 2 * PI / n;
    poly = new polyhedron();
    poly.name = "Y" + n;
    for (i = 0, _ref = n - 1; 0 <= _ref ? i <= _ref : i >= _ref; 0 <= _ref ? i++ : i--) {
      poly.xyz.push([cos(i * theta), sin(i * theta), 0.2]);
    }
    poly.xyz.push([0, 0, -2]);
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
      var ctr, f, i, j, poly, v, v0, _ref, _ref2;
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
        while (v !== v0) {
          poly.face[ctr].push(this.verts[v]);
          v = this.flags[i][v];
        }
        ctr++;
      }
      poly.name = "unknown polyhedron";
      return poly;
    };
    return polyflag;
  })();
  kisN = function(poly, n) {
    var centers, f, flag, fname, foundAny, i, newpoly, p, v, v1, v2, _i, _j, _k, _len, _len2, _len3, _ref, _ref2, _ref3, _ref4;
    console.log("Taking kis of " + (n === 0 ? "" : n) + "-sided faces of " + poly.name + "...");
    flag = new polyflag();
    _ref = enumerate(poly.xyz);
    for (_i = 0, _len = _ref.length; _i < _len; _i++) {
      _ref2 = _ref[_i], i = _ref2[0], p = _ref2[1];
      flag.newV("v" + i, p);
    }
    centers = faceCenters(poly);
    foundAny = false;
    _ref3 = enumerate(poly.face);
    for (_j = 0, _len2 = _ref3.length; _j < _len2; _j++) {
      _ref4 = _ref3[_j], i = _ref4[0], f = _ref4[1];
      v1 = "v" + f.slice(-1)[0];
      for (_k = 0, _len3 = f.length; _k < _len3; _k++) {
        v = f[_k];
        v2 = "v" + v;
        if (f.length === n || n === 0) {
          foundAny = true;
          flag.newV("f" + i, centers[i]);
          fname = i + v1;
          flag.newFlag(fname, v1, v2);
          flag.newFlag(fname, v2, "f" + i);
          flag.newFlag(fname, "f" + i, v1);
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
    newpoly.name = "k" + (n === 0 ? "" : n) + poly.name;
    newpoly.xyz = adjustXYZ(newpoly, 3);
    return newpoly;
  };
  midName = function(v1, v2) {
    if (v1 < v2) {
      return v1 + "_" + v2;
    } else {
      return v2 + "_" + v1;
    }
  };
  ambo = function(poly) {
    var f, flag, i, newpoly, v1, v2, v3, _i, _j, _len, _len2, _ref, _ref2, _ref3, _ref4;
    console.log("Taking ambo of " + poly.name + "...");
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
        flag.newFlag("f" + i, midName(v1, v2), midName(v2, v3));
        flag.newFlag("v" + v2, midName(v2, v3), midName(v1, v2));
        _ref4 = [v2, v3], v1 = _ref4[0], v2 = _ref4[1];
      }
    }
    newpoly = flag.topoly();
    newpoly.name = "a" + poly.name;
    newpoly.xyz = adjustXYZ(newpoly, 2);
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
    centers = faceCenters(poly);
    _ref3 = enumerate(poly.face);
    for (_j = 0, _len2 = _ref3.length; _j < _len2; _j++) {
      _ref4 = _ref3[_j], i = _ref4[0], f = _ref4[1];
      flag.newV("f" + i, unit(centers[i]));
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
        flag.newFlag(fname, "f" + i, v1 + "~" + v2);
        flag.newFlag(fname, v1 + "~" + v2, v2 + "~" + v1);
        flag.newFlag(fname, v2 + "~" + v1, "v" + v2);
        flag.newFlag(fname, "v" + v2, v2 + "~" + v3);
        flag.newFlag(fname, v2 + "~" + v3, "f" + i);
        _ref10 = [v2, v3], v1 = _ref10[0], v2 = _ref10[1];
      }
    }
    newpoly = flag.topoly();
    newpoly.name = "g" + poly.name;
    newpoly.xyz = adjustXYZ(newpoly, 3);
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
    newpoly.xyz = adjustXYZ(newpoly, 3);
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
    poly.xyz = adjustXYZ(poly, 1);
    return poly;
  };
  dual = function(poly) {
    var centers, dpoly, f, face, flag, i, j, k, sortF, v1, v2, _i, _len, _ref, _ref2, _ref3, _ref4, _ref5, _ref6, _ref7;
    console.log("Taking dual of " + poly.name + "...");
    flag = new polyflag();
    face = [];
    for (i = 0, _ref = poly.xyz.length - 1; 0 <= _ref ? i <= _ref : i >= _ref; 0 <= _ref ? i++ : i--) {
      face[i] = {};
    }
    for (i = 0, _ref2 = poly.face.length - 1; 0 <= _ref2 ? i <= _ref2 : i >= _ref2; 0 <= _ref2 ? i++ : i--) {
      v1 = poly.face[i][poly.face[i].length - 1];
      for (j = 0, _ref3 = poly.face[i].length - 1; 0 <= _ref3 ? j <= _ref3 : j >= _ref3; 0 <= _ref3 ? j++ : j--) {
        v2 = poly.face[i][j];
        face[v1]["v" + v2] = "" + i;
        v1 = v2;
      }
    }
    centers = faceCenters(poly);
    for (i = 0, _ref4 = poly.face.length - 1; 0 <= _ref4 ? i <= _ref4 : i >= _ref4; 0 <= _ref4 ? i++ : i--) {
      flag.newV("" + i, centers[i]);
    }
    for (i = 0, _ref5 = poly.face.length - 1; 0 <= _ref5 ? i <= _ref5 : i >= _ref5; 0 <= _ref5 ? i++ : i--) {
      v1 = poly.face[i][poly.face[i].length - 1];
      for (j = 0, _ref6 = poly.face[i].length - 1; 0 <= _ref6 ? j <= _ref6 : j >= _ref6; 0 <= _ref6 ? j++ : j--) {
        v2 = poly.face[i][j];
        flag.newFlag(v1, face[v2]["v" + v1], "" + i);
        v1 = v2;
      }
    }
    dpoly = flag.topoly();
    sortF = [];
    _ref7 = dpoly.face;
    for (_i = 0, _len = _ref7.length; _i < _len; _i++) {
      f = _ref7[_i];
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
      c = centroid(coords);
      if (dot(n, centroid) < 0) {
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
    centers = faceCenters(poly);
    for (_i = 0, _len = centers.length; _i < _len; _i++) {
      c = centers[_i];
      c = mult(1 / dot(c, c), c);
    }
    return centers;
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
  reciprocalN = function(poly) {
    var ans, avgEdgeDist, centroid, f, normal, tmp, v1, v2, v3, _i, _j, _len, _len2, _ref, _ref2, _ref3;
    ans = [];
    _ref = poly.face;
    for (_i = 0, _len = _ref.length; _i < _len; _i++) {
      f = _ref[_i];
      centroid = [0, 0, 0];
      normal = [0, 0, 0];
      avgEdgeDist = 0.0;
      _ref2 = f.slice(-2), v1 = _ref2[0], v2 = _ref2[1];
      for (_j = 0, _len2 = f.length; _j < _len2; _j++) {
        v3 = f[_j];
        centroid = add(centroid, poly.xyz[v3]);
        normal = add(normal, orthogonal(poly.xyz[v1], poly.xyz[v2], poly.xyz[v3]));
        avgEdgeDist += edgeDist(poly.xyz[v1], poly.xyz[v2]);
        _ref3 = [v2, v3], v1 = _ref3[0], v2 = _ref3[1];
      }
      centroid = mult(1 / f.length, centroid);
      normal = unit(normal);
      avgEdgeDist = avgEdgeDist / f.length;
      tmp = reciprocal(mult(dot(centroid, normal), normal));
      ans.push(mult((1 + avgEdgeDist) / 2, tmp));
    }
    return ans;
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
  faceCenters = function(poly) {
    var centers, i, j, _ref, _ref2;
    centers = [];
    for (i = 0, _ref = poly.face.length - 1; 0 <= _ref ? i <= _ref : i >= _ref; 0 <= _ref ? i++ : i--) {
      centers[i] = [0, 0, 0];
      for (j = 0, _ref2 = poly.face[i].length - 1; 0 <= _ref2 ? j <= _ref2 : j >= _ref2; 0 <= _ref2 ? j++ : j--) {
        centers[i] = add(centers[i], poly.xyz[poly.face[i][j]]);
      }
      centers[i] = mult(1.0 / poly.face[i].length, centers[i]);
    }
    return centers;
  };
  centroid = function(xyzs) {
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
  specreplacements = [[/P4$/g, "C"], [/e/g, "aa"], [/b/g, "ta"], [/o/g, "jj"], [/m/g, "kj"], [/t(\d*)/g, "dk$1d"], [/j/g, "dad"], [/s/g, "dgd"], [/dd/g, ""], [/aO/g, "aC"], [/aI/g, "aD"], [/gO/g, "gC"], [/gI/g, "gD"]];
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
        case ".":
          poly.xyz = canonicalXYZ(poly, n === 0 ? 5 : n * 5);
          break;
        case "!":
          poly.xyz = canonicalize(poly, n === 0 ? 5 : n * 80);
      }
      ops = ops.slice(0, -1);
    }
    poly.xyz = recenter(poly.xyz, poly.getEdges());
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
    var t, tests, _i, _j, _k, _l, _len, _len2, _len3, _len4, _len5, _len6, _m, _n;
    tests = ["T", "O", "C", "I", "D", "P3", "P4", "A4", "A5", "Y3", "Y4"];
    console.log("===== Test Basic Ops =====");
    console.log("--- primitives ----------------------------------------------------------------------------- ");
    for (_i = 0, _len = tests.length; _i < _len; _i++) {
      t = tests[_i];
      console.log(generatePoly(t));
    }
    console.log("--- kis ----------------------------------------------------------------------------- ");
    for (_j = 0, _len2 = tests.length; _j < _len2; _j++) {
      t = tests[_j];
      console.log(generatePoly("k" + t));
    }
    console.log("--- ambo ----------------------------------------------------------------------------- ");
    for (_k = 0, _len3 = tests.length; _k < _len3; _k++) {
      t = tests[_k];
      console.log(generatePoly("a" + t));
    }
    console.log("--- gyro ----------------------------------------------------------------------------- ");
    for (_l = 0, _len4 = tests.length; _l < _len4; _l++) {
      t = tests[_l];
      console.log(generatePoly("g" + t));
    }
    console.log("--- propellor ----------------------------------------------------------------------------- ");
    for (_m = 0, _len5 = tests.length; _m < _len5; _m++) {
      t = tests[_m];
      console.log(generatePoly("p" + t));
    }
    console.log("--- dual ----------------------------------------------------------------------------- ");
    for (_n = 0, _len6 = tests.length; _n < _len6; _n++) {
      t = tests[_n];
      console.log(generatePoly("d" + t));
    }
    return console.log("===== Done Testing Basic Ops =====");
  };
  ctx = {};
  CANVAS_WIDTH = 600;
  CANVAS_HEIGHT = 300;
  globPolys = {};
  globtheta = 0;
  globphi = 0;
  perspective_scale = 500;
  persp_z_max = 5;
  persp_z_min = 0;
  persp_ratio = 0.8;
  _2d_x_offset = CANVAS_WIDTH / 2;
  _2d_y_offset = CANVAS_HEIGHT / 2;
  globtime = new Date();
  BG_CLEAR = true;
  BG_COLOR = "rgba(255,255,255,1.0)";
  COLOR_METHOD = "area";
  ctx_linewidth = 0.5;
  MOUSEDOWN = false;
  LastMouseX = 0;
  LastMouseY = 0;
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
    var centroids, idx, zsortIndex, _i, _ref, _results;
    centroids = faceCenters(poly);
    zsortIndex = _.zip(centroids, (function() {
      _results = [];
      for (var _i = 0, _ref = poly.face.length - 1; 0 <= _ref ? _i <= _ref : _i >= _ref; 0 <= _ref ? _i++ : _i--){ _results.push(_i); }
      return _results;
    }).apply(this)).sort(function(a, b) {
      return a[0][2] - b[0][2];
    }).map(function(x) {
      return x[1];
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
  drawpoly = function(poly, tvec, rot) {
    var clr, face, face_verts, fno, illum, oldxyz, v, v0, x, y, _i, _j, _len, _len2, _ref, _ref2, _ref3, _ref4;
    tvec || (tvec = [3, 3, 3]);
    rot || (rot = [1, 0, 1]);
    oldxyz = _.map(poly.xyz, function(x) {
      return x;
    });
    poly.xyz = _.map(poly.xyz, function(x) {
      return mv3(rotm(rot[0], rot[1], rot[2]), x);
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
      illum = dot(normal(face_verts), unit([-1, 1, 0]));
      clr = mult((illum / 2.0 + .5) * 0.7 + 0.3, clr);
      ctx.fillStyle = "rgba(" + (round(clr[0] * 255)) + ", " + (round(clr[1] * 255)) + ", " + (round(clr[2] * 255)) + ", " + 1.0 + ")";
      ctx.fill();
      ctx.stroke();
    }
    return poly.xyz = oldxyz;
  };
  $(function() {
    var defspecs, specs, urlParams;
    init();
    defspecs = ["dakD", "oopD", "ajI", ".akY5", ".ooC", "bT"];
    urlParams = parseurl();
    if ("recipe" in urlParams) {
      specs = [urlParams["recipe"]];
      $("#spec").val(specs);
    } else {
      specs = [randomchoice(defspecs)];
      $("#spec").val(specs);
    }
    globPolys = _.map(specs, function(x) {
      return generatePoly(x);
    });
    drawShape();
    $("#spec").change(__bind(function(e) {
      specs = $("#spec").val().split(/\s+/g);
      globPolys = _.map(specs, function(x) {
        return generatePoly(x);
      });
      return drawShape();
    }, this));
    $("#poly").mousewheel(function(e, delta, deltaX, deltaY) {
      event.preventDefault();
      perspective_scale *= (10 + delta) / 10;
      return drawShape();
    });
    $("#poly").mousedown(function(e) {
      event.preventDefault();
      MOUSEDOWN = true;
      LastMouseX = e.clientX - $(this).offset().left;
      return LastMouseY = e.clientY - $(this).offset().top;
    });
    $("#poly").mouseup(function(e) {
      event.preventDefault();
      return MOUSEDOWN = false;
    });
    $("#poly").mouseleave(function(e) {
      event.preventDefault();
      return MOUSEDOWN = false;
    });
    return $("#poly").mousemove(function(e) {
      event.preventDefault();
      if (MOUSEDOWN) {
        globtheta += -(e.clientX - $(this).offset().left - LastMouseX) * (Math.PI / 180);
        globphi += -(e.clientY - $(this).offset().top - LastMouseY) * (Math.PI / 180);
        LastMouseX = e.clientX - $(this).offset().left;
        LastMouseY = e.clientY - $(this).offset().top;
        return drawShape();
      }
    });
  });
  animateShape = function() {
    var i, p, _i, _len, _ref, _ref2;
    clear();
    globtheta = (2 * Math.PI) / 180.0 * globtime.getSeconds() * 0.1;
    _ref = enumerate(globPolys);
    for (_i = 0, _len = _ref.length; _i < _len; _i++) {
      _ref2 = _ref[_i], i = _ref2[0], p = _ref2[1];
      drawpoly(p, [0 + 3 * i, 0, 3], [0, globtheta, globphi]);
    }
    return setTimeout(animateShape, 100);
  };
  drawShape = function() {
    var i, p, _i, _len, _ref, _ref2, _results;
    clear();
    _ref = enumerate(globPolys);
    _results = [];
    for (_i = 0, _len = _ref.length; _i < _len; _i++) {
      _ref2 = _ref[_i], i = _ref2[0], p = _ref2[1];
      _results.push(drawpoly(p, [0 + 3 * i, 0, 3], [0, globtheta, globphi]));
    }
    return _results;
  };
}).call(this);
