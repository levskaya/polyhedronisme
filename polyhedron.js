// Polyhédronisme
//===================================================================================================
//
// A toy for constructing and manipulating polyhedra.
//
// Copyright 2019, Anselm Levskaya
// Released under the MIT License
//

function __range__(left, right, inclusive) {
  let range = [];
  let ascending = left < right;
  let end = !inclusive ? right : ascending ? right + 1 : right - 1;
  for (let i = left; ascending ? i < end : i > end; ascending ? i++ : i--) {
    range.push(i);
  }
  return range;
}

// Polyhedra Functions
//=================================================================================================
//
// Topology stored as set of faces.  Each face is list of n vertex indices
// corresponding to one oriented, n-sided face.  Vertices listed clockwise as seen from outside.

// Generate an array of edges [v1,v2] for the face.
const faceToEdges = function(face) {
  const edges = [];
  let [v1] = face.slice(-1);
  for (let v2 of face) {
    edges.push([v1, v2]);
    v1 = v2;
  }
  return edges;
};

const vertColors = function(poly) {
  const vertcolors=[];
  for (let i = 0; i < poly.faces.length; i++) {
    const face = poly.faces[i];
    for (let v of face) {
      vertcolors[v] = poly.face_classes[i];
    }
  }
  return vertcolors;
};

// Polyhedra Coloring Functions
//=================================================================================================
const rwb_palette = ["#ff7777", "#dddddd", "#889999", "#fff0e5",
                     "#aa3333", "#ff0000", "#ffffff", "#aaaaaa"];
let PALETTE = rwb_palette;  // GLOBAL
const palette = function(n) {
  const k = n % PALETTE.length;
  return hextofloats(PALETTE[k])
};
                     
// converts [h,s,l] float args to [r,g,b] list
function hsl2rgb(h, s, l) {
  let r, g, b;
  if (s == 0) {
    r = g = b = l; // achromatic
  } else {
    const hue2rgb = function(p, q, t) {
      if (t < 0) t += 1;
      if (t > 1) t -= 1;
      if (t < 1/6) return p + (q - p) * 6 * t;
      if (t < 1/2) return q;
      if (t < 2/3) return p + (q - p) * (2/3 - t) * 6;
      return p;
    }
    let q = l < 0.5 ? l * (1 + s) : l + s - l * s;
    let p = 2 * l - q;
    r = hue2rgb(p, q, h + 1/3);
    g = hue2rgb(p, q, h);
    b = hue2rgb(p, q, h - 1/3);
  }
  return [r, g, b];
}

// converts #xxxxxx / #xxx format into list of [r,g,b] floats
const hextofloats = function(hexstr){
  let rgb;
  if (hexstr[0] === "#") {
    hexstr = hexstr.slice(1);
  }
  if (hexstr.length === 3) {
    rgb = hexstr.split('').map(c=> parseInt(c+c, 16)/255);
  } else {
    rgb = hexstr.match(/.{2}/g).map(c=> parseInt(c, 16)/255);
  }
  return rgb;
};

// converts [r,g,b] floats to #xxxxxx form
const floatstohex = function(rgb){
  let r_hex = Number(parseInt(255 * rgb[0], 10)).toString(16);
  let g_hex = Number(parseInt(255 * rgb[1], 10)).toString(16);
  let b_hex = Number(parseInt(255 * rgb[2], 10)).toString(16);
  return "#" + r_hex + g_hex + b_hex;
}

// randomize color palette
const rndcolors = function(){
  let newpalette=[];
  for(let i=0; i<100; i++){
    let h = random();
    let s = 0.5*random() + 0.3;
    let l = 0.5*random() + 0.45;
    let rgb = hsl2rgb(h, s, l);
    newpalette.push(floatstohex(rgb));
  }
  return newpalette;
}


// Color the faces of the polyhedra for display
const paintPolyhedron = function(poly) {
  poly.face_classes = [];
  const colormemory={};

  // memoized color assignment to faces of similar areas
  const colorassign = function(hash, colormemory) {
    //const hash = ar;
    if (hash in colormemory) {
      return colormemory[hash];
    } else {
      const fclr = _.toArray(colormemory).length;
      colormemory[hash] = fclr;
      return fclr;
    }
  };

  for (var f of poly.faces) {
    var clr, face_verts;
    if (COLOR_METHOD === "area") {
      // color by face planar area assuming flatness
      face_verts = f.map(v=>poly.vertices[v])
      clr = colorassign(sigfigs(planararea(face_verts), COLOR_SENSITIVITY), colormemory);
    } else if (COLOR_METHOD === "signature") {
      // color by congruence signature
      face_verts = f.map(v=>poly.vertices[v])
      clr = colorassign(faceSignature(face_verts, COLOR_SENSITIVITY), colormemory);
    } else {
      // color by face-sidedness
      clr = f.length - 3;
    }
    poly.face_classes.push(clr);
  }
  console.log(_.toArray(colormemory).length+" face classes");
  return poly;
};

// z sorts faces of poly
// -------------------------------------------------------------------------
const sortfaces = function(poly) {
  //smallestZ = (x) -> _.sortBy(x,(a,b)->a[2]-b[2])[0]
  //closests = (smallestZ(poly.vertices[v] for v in f) for f in poly.faces)
  let idx;
  const centroids  = poly.centers();
  const normals    = poly.normals();
  const ray_origin = [0,0, ((persp_z_max * persp_ratio) - persp_z_min)/(1-persp_ratio)];

  // sort by binary-space partition: are you on same side as view-origin or not?
  // !!! there is something wrong with this. even triangulated surfaces have artifacts.
  const planesort = (a,b) =>
    //console.log dot(sub(ray_origin,a[0]),a[1]), dot(sub(b[0],a[0]),a[1])
    -dot(sub(ray_origin,a[0]),a[1])*dot(sub(b[0],a[0]),a[1]);

  // sort by centroid z-depth: not correct but more stable heuristic w. weird non-planar "polygons"
  const zcentroidsort = (a, b) => a[0][2]-b[0][2];

  const zsortIndex = _.zip(centroids, normals, __range__(0, poly.faces.length, false))
    //.sort(planesort)
    .sort(zcentroidsort)
    .map(x=> x[2]);

  // sort all face-associated properties
  poly.faces = zsortIndex.map(idx=>poly.faces[idx]);
  poly.face_classes = zsortIndex.map(idx=>poly.face_classes[idx]);
};


class polyhedron {
  // constructor of initially null polyhedron
  constructor(verts, faces, name) {
    // array of faces.  faces.length = # faces
    this.faces = faces || new Array();
    // array of vertex coords.  vertices.length = # of vertices
    this.vertices  = verts || new Array();
    this.name = name  || "null polyhedron";
  }
  
  // return a non-redundant list of the polyhedron's edges
  edges() {
    let e, a, b;
    const uniqEdges = {};
    const faceEdges = this.faces.map(faceToEdges);
    for (let edgeSet of faceEdges) {
      for (e of edgeSet) {
        if (e[0] < e[1]) {
          [a, b] = e;
        } else {
          [b, a] = e;
        }
        uniqEdges[`${a}~${b}`] = e;
      }
    }
    return _.values(uniqEdges);
  }
      
  // get array of face centers
  centers() {
    const centersArray = [];
    for (let face of this.faces) {
      let fcenter = [0, 0, 0];
      // average vertex coords
      for (let vidx of face) {
        fcenter = add(fcenter, this.vertices[vidx]);
      }
      centersArray.push(mult(1.0 / face.length, fcenter));
    }
    // return face-ordered array of centroids
    return centersArray;
  }

  // get array of face normals
  normals() {
    const normalsArray = [];
    for (let face of this.faces) {
      normalsArray.push(normal(face.map(vidx => this.vertices[vidx])));
    }
    return normalsArray;
  }

  // informative string
  data() {
    const nEdges = (this.faces.length + this.vertices.length) - 2; // E = V + F - 2
    return `${this.faces.length} faces, ${nEdges} edges, ${this.vertices.length} vertices`;
  }

  moreData() {
    return `min. edge length ${this.minEdgeLength().toPrecision(2)}; ` +
           `min. face radius ${this.minFaceRadius().toPrecision(2)}`;
  }

  minEdgeLength() {
    let min2 = Number.MAX_VALUE;
    // Compute minimum edge length
    for (let e of this.edges()) {
      // square of edge length
      const d2 = mag2(sub(this.vertices[e[0]], this.vertices[e[1]]));
      if (d2 < min2) {
        min2 = d2;
      }
    }
    // This is normalized if rescaling has happened.
    return sqrt(min2); 
  }
    
  minFaceRadius() {
    let min2 = Number.MAX_VALUE;
    const nFaces = this.faces.length;
    const centers = this.centers();
    for (let f = 0, end = nFaces; f < end; f++) {
      const c = centers[f];
      for (let e of faceToEdges(this.faces[f])) {
        // Check distance from center to each edge.
        const de2 = linePointDist2(this.vertices[e[0]], this.vertices[e[1]], c);
        if (de2 < min2) {
          min2 = de2;
        }
      }
    }
    return sqrt(min2);
  }

  // Export / Formatting Routines --------------------------------------------------

  // produces vanilla OBJ files for import into 3d apps
  toOBJ() {
    let f;
    let v;
    let objstr="#Produced by polyHédronisme http://levskaya.github.com/polyhedronisme\n";
    objstr+=`group ${this.name}\n`;
    objstr+="#vertices\n";
    for (v of this.vertices) {
      objstr += `v ${v[0]} ${v[1]} ${v[2]}\n`;
    }
    objstr += "#normal vector defs \n";
    for (f of this.faces) {
      const norm = normal(f.map(v=>this.vertices[v]))
      objstr += `vn ${norm[0]} ${norm[1]} ${norm[2]}\n`;
    }
    objstr += "#face defs \n";
    for (let i = 0; i < this.faces.length; i++) {
      f = this.faces[i];
      objstr += "f ";
      for (v of f) {
        objstr += `${v+1}//${i+1} `;
      }
      objstr += "\n";
    }
    return objstr;
  }

  toX3D() {
    let v;
    // ShapeWays uses 1unit = 1meter, so reduce to 3cm scale
    const SCALE_FACTOR = .03;
    // opening cruft
    let x3dstr=`\
<?xml version="1.0" encoding ="UTF-8"?>
<X3D profile="Interchange" version="3.0">
<head>
<component name="Rendering" level="3"/>
<meta name="generator" content="Polyhedronisme"/>
<meta name="version" content="0.1.0"/>
</head>
<Scene>
<Shape>
<IndexedFaceSet normalPerVertex="false" coordIndex="\
`;
    // face indices
    for (let f of this.faces) {
      for (v of f) {
        x3dstr+=`${v} `;
      }
      x3dstr+='-1\n';
    }
    x3dstr+='">\n';

    // per-face Color
    x3dstr+='<Color color="';
    for (let cl of vertColors(this)) {//@face_class
      const clr=palette(cl);
      x3dstr+=`${clr[0]} ${clr[1]} ${clr[2]} `;
    }
    x3dstr+='"/>';

    // re-scaled xyz coordinates
    x3dstr+='<Coordinate point="';
    for (v of this.vertices) {
      x3dstr+=`${v[0]*SCALE_FACTOR} ${v[1]*SCALE_FACTOR} ${v[2]*SCALE_FACTOR} `;
    }
    x3dstr+='"/>\n';

      // end cruft
    x3dstr+=`\
</IndexedFaceSet>
</Shape>
</Scene>
</X3D>`;

    return x3dstr;
  }

  toVRML() {
    let v;
    // ShapeWays uses 1unit = 1meter, so reduce to 3cm scale
    const SCALE_FACTOR = .03;
    // opening cruft
    let x3dstr=`\
#VRML V2.0 utf8
#Generated by Polyhedronisme
NavigationInfo {
	type [ "EXAMINE", "ANY" ]
}
Transform {
  scale 1 1 1
  translation 0 0 0
  children
  [
    Shape
    {
      geometry IndexedFaceSet
      {
        creaseAngle .5
        solid FALSE
        coord Coordinate
        {
          point
          [\
`;
    // re-scaled xyz coordinates
    for (v of this.vertices) {
      x3dstr+=`${v[0]*SCALE_FACTOR} ${v[1]*SCALE_FACTOR} ${v[2]*SCALE_FACTOR},`;
    }
    x3dstr=x3dstr.slice(0, +-2 + 1 || undefined);
    x3dstr+=`\
    ]
}
color Color
{
  color
  [\
`;
    // per-face Color
    for (let cl of this.face_classes) {
      const clr=palette(cl);
      x3dstr+=`${clr[0]} ${clr[1]} ${clr[2]} ,`;
    }
    x3dstr=x3dstr.slice(0, +-2 + 1 || undefined);
    x3dstr+=`\
  ]
}
colorPerVertex FALSE
coordIndex
[\
`;
    // face indices
    for (let f of this.faces) {
      for (v of f) {
        x3dstr+=`${v}, `;
      }
      x3dstr+='-1,';
    }
    x3dstr=x3dstr.slice(0, +-2 + 1 || undefined);
    x3dstr+=`\
          ]
      }
      appearance Appearance
      {
        material Material
        {
	       ambientIntensity 0.2
	       diffuseColor 0.9 0.9 0.9
	       specularColor .1 .1 .1
	       shininess .5
        }
      }
    }
  ]
}\
`;
    return x3dstr;
  }
}

//===================================================================================================
// Primitive Polyhedra Seeds
//===================================================================================================

const tetrahedron = function() {
  const poly = new polyhedron();
  poly.name = "T";
  poly.faces = [ [0,1,2], [0,2,3], [0,3,1], [1,3,2] ];
  poly.vertices  = [ [1.0,1.0,1.0], [1.0,-1.0,-1.0], [-1.0,1.0,-1.0], [-1.0,-1.0,1.0] ];
  return poly;
};

const octahedron = function() {
  const poly = new polyhedron();
  poly.name = "O";
  poly.faces = [ [0,1,2], [0,2,3], [0,3,4], [0,4,1], [1,4,5], [1,5,2], [2,5,3], [3,5,4] ];
  poly.vertices  = [ [0,0,1.414], [1.414,0,0], [0,1.414,0], [-1.414,0,0], [0,-1.414,0], [0,0,-1.414] ];
  return poly;
};

const cube = function() {
  const poly = new polyhedron();
  poly.name = "C";
  poly.faces = [ [3,0,1,2], [3,4,5,0], [0,5,6,1], [1,6,7,2], [2,7,4,3], [5,4,7,6] ];
  poly.vertices  = [ [0.707,0.707,0.707], [-0.707,0.707,0.707], [-0.707,-0.707,0.707], [0.707,-0.707,0.707],
                [0.707,-0.707,-0.707], [0.707,0.707,-0.707], [-0.707,0.707,-0.707], [-0.707,-0.707,-0.707] ];
  return poly;
};

const icosahedron = function() {
  const poly = new polyhedron();
  poly.name = "I";
  poly.faces = [ [0,1,2], [0,2,3], [0,3,4], [0,4,5],
    [0,5,1], [1,5,7], [1,7,6], [1,6,2],
    [2,6,8], [2,8,3], [3,8,9], [3,9,4],
    [4,9,10], [4,10,5], [5,10,7], [6,7,11],
    [6,11,8], [7,10,11], [8,11,9], [9,11,10] ];

  poly.vertices = [ [0,0,1.176], [1.051,0,0.526],
    [0.324,1.0,0.525], [-0.851,0.618,0.526],
    [-0.851,-0.618,0.526], [0.325,-1.0,0.526],
    [0.851,0.618,-0.526], [0.851,-0.618,-0.526],
    [-0.325,1.0,-0.526], [-1.051,0,-0.526],
    [-0.325,-1.0,-0.526], [0,0,-1.176] ];
  return poly;
};

const dodecahedron = function() {
   const poly = new polyhedron();
   poly.name = "D";
   poly.faces = [ [0,1,4,7,2], [0,2,6,9,3], [0,3,8,5,1],
      [1,5,11,10,4], [2,7,13,12,6], [3,9,15,14,8],
      [4,10,16,13,7], [5,8,14,17,11], [6,12,18,15,9],
      [10,11,17,19,16], [12,13,16,19,18], [14,15,18,19,17] ];
   poly.vertices = [ [0,0,1.07047], [0.713644,0,0.797878],
      [-0.356822,0.618,0.797878], [-0.356822,-0.618,0.797878],
      [0.797878,0.618034,0.356822], [0.797878,-0.618,0.356822],
      [-0.934172,0.381966,0.356822], [0.136294,1.0,0.356822],
      [0.136294,-1.0,0.356822], [-0.934172,-0.381966,0.356822],
      [0.934172,0.381966,-0.356822], [0.934172,-0.381966,-0.356822],
      [-0.797878,0.618,-0.356822], [-0.136294,1.0,-0.356822],
      [-0.136294,-1.0,-0.356822], [-0.797878,-0.618034,-0.356822],
      [0.356822,0.618,-0.797878], [0.356822,-0.618,-0.797878],
      [-0.713644,0,-0.797878], [0,0,-1.07047] ];
   return poly;
 };

const prism = function(n) {
  let i;
  const theta = (2*PI)/n; // pie angle
  const h = Math.sin(theta/2); // half-edge
  let poly = new polyhedron();
  poly.name = `P${n}`;

  for (i = 0; i < n; i++) { // vertex #'s 0 to n-1 around one face
    poly.vertices.push([-cos(i*theta), -sin(i*theta),  -h]);
  }
  for (i = 0; i < n; i++) { // vertex #'s n to 2n-1 around other
    poly.vertices.push([-cos(i*theta), -sin(i*theta), h]);
  }

  poly.faces.push(__range__(n-1, 0, true));  //top
  poly.faces.push(__range__(n, 2*n, false)); //bottom
  for (i = 0; i < n; i++) { //n square sides
    poly.faces.push([i, (i+1)%n, ((i+1)%n)+n, i+n]);
  }

  poly = adjustXYZ(poly,1);
  return poly;
};

const antiprism = function(n) {
  let i;
  const theta = (2*PI)/n; // pie angle
  let h = sqrt(1-(4/((4+(2*cos(theta/2)))-(2*cos(theta)))));
  let r = sqrt(1-(h*h));
  const f = sqrt((h*h) + pow(r*cos(theta/2),2));
  // correction so edge midpoints (not vertices) on unit sphere
  r = -r/f;
  h = -h/f;
  let poly = new polyhedron();
  poly.name = `A${n}`;

  for (i = 0; i < n; i++) { // vertex #'s 0...n-1 around one face
    poly.vertices.push([r * cos(i*theta), r * sin(i*theta), h]);
  }
  for (i = 0; i < n; i++) { // vertex #'s n...2n-1 around other
    poly.vertices.push([r * cos((i+0.5)*theta), r * sin((i+0.5)*theta), -h]);
  }

  poly.faces.push(__range__(n-1, 0, true));   //top
  poly.faces.push(__range__(n, (2*n)-1, true)); //bottom
  for (i = 0; i <= n-1; i++) { //2n triangular sides
    poly.faces.push([i, (i+1)%n, i+n]);
    poly.faces.push([i, i+n, ((((n+i)-1)%n)+n)]);
  }

  poly = adjustXYZ(poly,1);
  return poly;
};

const pyramid = function(n) {
  let i;
  const theta = (2*PI)/n; // pie angle
  const height = 1;
  let poly = new polyhedron();
  poly.name = `Y${n}`;

  for (i = 0; i < n; i++) { // vertex #'s 0...n-1 around one face
    poly.vertices.push([-cos(i*theta), -sin(i*theta), -0.2]);
  }
  poly.vertices.push([0,0, height]); // apex

  poly.faces.push(__range__(n-1, 0, true)); // base
  for (i = 0; i < n; i++) { // n triangular sides
    poly.faces.push([i, (i+1)%n, n]);
  }

  poly = canonicalXYZ(poly, 3);
  return poly;
};

const cupola = function(n, alpha, height) {
  let i;
  if (n===undefined) { n = 3; }
  if (alpha===undefined) { alpha = 0.0; }

  let poly = new polyhedron();
  poly.name = `U${n}`;

  if (n < 2) {
    return poly;
  }

  let s = 1.0;
  // alternative face/height scaling 
  //let rb = s / 2 / sin(PI / 2 / n - alpha);
  let rb = s / 2 / sin(PI / 2 / n);
  let rt = s / 2 / sin(PI / n);
  if (height===undefined) { 
    height = (rb - rt);
    // set correct height for regularity for n=3,4,5
    if (2 <= n && n <= 5) {
      height = s * sqrt(1 - 1 / 4 / sin(PI/n) / sin(PI/n));
    }
  }
  // init 3N vertices
  for (i = 0; i < 3*n; i++) {
    poly.vertices.push([0,0,0]);
  }
  // fill vertices
  for (i = 0; i < n; i++) {
    poly.vertices[2*i] = [rb * cos(PI*(2*i)/n + PI/2/n+alpha), 
                          rb * sin(PI*(2*i)/n + PI/2/n+alpha),
                          0.0];
    poly.vertices[2*i+1] = [rb * cos(PI*(2*i+1)/n + PI/2/n-alpha), 
                            rb * sin(PI*(2*i+1)/n + PI/2/n-alpha), 
                            0.0];
    poly.vertices[2*n+i] = [rt * cos(2*PI*i/n), 
                            rt * sin(2*PI*i/n), 
                            height];
  }
  
  poly.faces.push(__range__(2*n-1, 0, true)); // base
  poly.faces.push(__range__(2*n, 3*n-1, true)); // top
  for (i = 0; i < n; i++) { // n triangular sides and n square sides
    poly.faces.push([(2*i+1)%(2*n), (2*i+2)%(2*n), 2*n+(i+1)%n]);
    poly.faces.push([2*i, (2*i+1)%(2*n), 2*n+(i+1)%n, 2*n+i]);
  }

  return poly;  
}

const anticupola = function(n, alpha, height) {
  let i;
  if (n===undefined) { n = 3; }
  if (alpha===undefined) { alpha = 0.0; }

  let poly = new polyhedron();
  poly.name = `U${n}`;

  if (n < 3) {
    return poly;
  }

  let s = 1.0;
  // alternative face/height scaling 
  //let rb = s / 2 / sin(PI / 2 / n - alpha);
  let rb = s / 2 / sin(PI / 2 / n);
  let rt = s / 2 / sin(PI / n);
  if (height===undefined) { 
    height = (rb - rt);
  }
  // init 3N vertices
  for (i = 0; i < 3*n; i++) {
    poly.vertices.push([0,0,0]);
  }
  // fill vertices
  for (i = 0; i < n; i++) {
    poly.vertices[2*i] = [rb * cos(PI*(2*i)/n + alpha), 
                          rb * sin(PI*(2*i)/n + alpha),
                          0.0];
    poly.vertices[2*i+1] = [rb * cos(PI*(2*i+1)/n - alpha), 
                            rb * sin(PI*(2*i+1)/n - alpha), 
                            0.0];
    poly.vertices[2*n+i] = [rt * cos(2*PI*i/n), 
                            rt * sin(2*PI*i/n), 
                            height];
  }
  
  poly.faces.push(__range__(2*n-1, 0, true)); // base
  poly.faces.push(__range__(2*n, 3*n-1, true)); // top
  for (i = 0; i < n; i++) { // n triangular sides and n square sides
    poly.faces.push([(2*i)%(2*n), (2*i+1)%(2*n), 2*n+(i)%n]);
    poly.faces.push([2*n+(i+1)%n, (2*i+1)%(2*n), (2*i+2)%(2*n)]);
    poly.faces.push([2*n+(i+1)%n, 2*n+(i)%n, (2*i+1)%(2*n)]);
  }

  return poly;  
}