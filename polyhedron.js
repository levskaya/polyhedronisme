// Polyhédronisme
//===================================================================================================
//
// A toy for constructing and manipulating polyhedra and other meshes
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
//===================================================================================================
//
// Topology stored as set of "faces."  Each face is list of n vertex indices
// corresponding to one n-sided face.  Vertices listed clockwise as seen from outside.

// Generate an array of edges [v1,v2] for the face.
const faceToEdges = function(face) {
  const edges = [];
  let [v1] = face.slice(-1);
  for (let v2 of face) {
    edges.push([v1,v2]);
    v1 = v2;
  }
  return edges;
};

const vertColors = function(poly) {
  const vertcolors=[];
  for (let i = 0; i < poly.face.length; i++) {
    const f = poly.face[i];
    for (let v of f) {
      vertcolors[v] = poly.face_class[i];
    }
  }
  return vertcolors;
};

// Polyhedra Coloring Functions
//===================================================================================================

const rwb_palette  = ["#ff7777","#dddddd","#889999","#fff0e5","#aa3333","#ff0000","#ffffff","#aaaaaa"];

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

let PALETTE = rwb_palette; //GLOBAL
const palette = function(n) {
  const k = n % PALETTE.length;
  return hextofloats(PALETTE[k])
};

const paintPolyhedron = function(poly) {
  // Color the faces of the polyhedra for display
  let v;
  poly.face_class = [];
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

  for (var f of poly.face) {
    var clr, face_verts;
    if (COLOR_METHOD === "area") {
      // color by face planar area assuming flatness
      face_verts = f.map(v=>poly.xyz[v])
      clr = colorassign(sigfigs(planararea(face_verts), COLOR_SENSITIVITY), colormemory);
    } else if (COLOR_METHOD === "signature") {
      // color by congruence signature
      face_verts = f.map(v=>poly.xyz[v])
      clr = colorassign(faceSignature(face_verts, COLOR_SENSITIVITY), colormemory);
    } else {
      // color by face-sidedness
      clr = f.length - 3;
    }
    poly.face_class.push(clr);
  }
  console.log(_.toArray(colormemory).length+" face classes");
  return poly;
};

// z sorts faces of poly
// -------------------------------------------------------------------------
const sortfaces = function(poly) {
  //smallestZ = (x) -> _.sortBy(x,(a,b)->a[2]-b[2])[0]
  //closests = (smallestZ(poly.xyz[v] for v in f) for f in poly.face)
  let idx;
  const centroids  = poly.centers();
  const normals    = poly.normals();
  const ray_origin = [0,0, ((persp_z_max * persp_ratio) - persp_z_min)/(1-persp_ratio)];

  // sort by binary-space partition: are you on same side as view-origin or not?
  // !!! there is something wrong with this. even triangulated surfaces have artifacts.
  const planesort = (a,b)=>
    //console.log dot(sub(ray_origin,a[0]),a[1]), dot(sub(b[0],a[0]),a[1])
    -dot(sub(ray_origin,a[0]),a[1])*dot(sub(b[0],a[0]),a[1]);

  // sort by centroid z-depth: not correct but more stable heuristic w. weird non-planar "polygons"
  const zcentroidsort = (a, b)=> a[0][2]-b[0][2];

  const zsortIndex = _.zip(centroids, normals, __range__(0, poly.face.length, false))
    //.sort(planesort)
    .sort(zcentroidsort)
    .map(x=> x[2]);

  // sort all face-associated properties
  poly.face = zsortIndex.map(idx=>poly.face[idx]);
  poly.face_class = zsortIndex.map(idx=>poly.face_class[idx]);
};


class polyhedron {
  // constructor of initially null polyhedron
  constructor(verts, faces, name) {
    // array of faces.  face.length = # faces
    this.face = faces || new Array();
    // array of vertex coords.  xyz.length = # of vertices
    this.xyz  = verts || new Array();
    this.name = name  || "null polyhedron";
  }

  data() {   // informative string
    const nEdges = (this.face.length + this.xyz.length) - 2; // E = V + F - 2
    return `${this.face.length} faces, ${nEdges} edges, ${this.xyz.length} vertices`;
  }

  moreData() {
    return `min. edge length ${this.minEdgeLength().toPrecision(2)}; min. face radius ${this.minFaceRadius().toPrecision(2)}`;
  }
    
  edges() {
    let e;
    const finalset={};
    const uniqedges=[];
    const alledges = _.map(this.face, faceToEdges);
    for (let edgeset of alledges) {
      for (e of edgeset) {
        var a, b;
        if (e[0] < e[1]) {
          [a, b] = e;
        } else {
          [b, a] = e;
        }
        finalset[a+'~'+b] = e;
      }
    }
    for (let hash in finalset) {
      e = finalset[hash];
      uniqedges.push(e);
    }

    //return edges
    return uniqedges;
  }

  minEdgeLength() {
    let min2 = Number.MAX_VALUE;
    // Compute minimum edge length
    for (let e of this.edges()) {
      // square of edge length
      const d2 = mag2(sub(this.xyz[e[0]], this.xyz[e[1]]));
      if (d2 < min2) {
        min2 = d2;
      }
    }
    // This is normalized if rescaling has happened.
    return sqrt(min2); 
  }
    
  minFaceRadius() {
    let min2 = Number.MAX_VALUE;
    const nFaces = this.face.length;
    const centers = this.centers();
    for (let f = 0, end = nFaces; f < end; f++) {
      const c = centers[f];
      for (let e of faceToEdges(this.face[f])) {
        // Check distance from center to each edge.
        const de2 = linePointDist2(this.xyz[e[0]], this.xyz[e[1]], c);
        if (de2 < min2) {
          min2 = de2;
        }
      }
    }

    return sqrt(min2);
  }
      
  centers() {
    // get array of face centers
    const centers_array = [];
    for (let f of this.face) {
      let fcenter = [0, 0, 0];
      // average vertex coords
      for (let v of f) {
        fcenter = add(fcenter, this.xyz[v]);
      }
      centers_array.push(mult(1.0 / f.length, fcenter));
    }
    // return face-ordered array of centroids
    return centers_array;
  }

  normals() {
    // get array of face normals
    const normals_array = [];
    for (let f of this.face) {
      normals_array.push(normal(f.map((v) => this.xyz[v])));
    }
    return normals_array;
  }

  // Export / Formatting Routines --------------------------------------------------

  // produces vanilla OBJ files for import into 3d apps
  toOBJ() {
    let f;
    let v;
    let objstr="#Produced by polyHédronisme http://levskaya.github.com/polyhedronisme\n";
    objstr+=`group ${this.name}\n`;
    objstr+="#vertices\n";
    for (v of this.xyz) {
      objstr += `v ${v[0]} ${v[1]} ${v[2]}\n`;
    }

    objstr += "#normal vector defs \n";
    for (f of this.face) {
      // const norm = normal((() => {
      //   const result = [];
      //   for (v of f) {           result.push(this.xyz[v]);
      //   }
      //   return result;
      // })());
      const norm = normal(f.map(v=>this.xyz[v]))
      objstr += `vn ${norm[0]} ${norm[1]} ${norm[2]}\n`;
    }

    objstr += "#face defs \n";
    for (let i = 0; i < this.face.length; i++) {
      f = this.face[i];
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
    const SCALE_FACTOR = .01; //ShapeWays uses 1unit = 1meter, so reduce to 1cm scale
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
    for (let f of this.face) {
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
    for (v of this.xyz) {
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
    const SCALE_FACTOR = .01; //ShapeWays uses 1unit = 1meter, so reduce to 1cm scale
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
    for (v of this.xyz) {
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
    for (let cl of this.face_class) {
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
    for (let f of this.face) {
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
  poly.face = [ [0,1,2], [0,2,3], [0,3,1], [1,3,2] ];
  poly.xyz  = [ [1.0,1.0,1.0], [1.0,-1.0,-1.0], [-1.0,1.0,-1.0], [-1.0,-1.0,1.0] ];
  return poly;
};

const octahedron = function() {
  const poly = new polyhedron();
  poly.name = "O";
  poly.face = [ [0,1,2], [0,2,3], [0,3,4], [0,4,1], [1,4,5], [1,5,2], [2,5,3], [3,5,4] ];
  poly.xyz  = [ [0,0,1.414], [1.414,0,0], [0,1.414,0], [-1.414,0,0], [0,-1.414,0], [0,0,-1.414] ];
  return poly;
};

const cube = function() {
  const poly = new polyhedron();
  poly.name = "C";
  poly.face = [ [3,0,1,2], [3,4,5,0], [0,5,6,1], [1,6,7,2], [2,7,4,3], [5,4,7,6] ];
  poly.xyz  = [ [0.707,0.707,0.707], [-0.707,0.707,0.707], [-0.707,-0.707,0.707], [0.707,-0.707,0.707],
                [0.707,-0.707,-0.707], [0.707,0.707,-0.707], [-0.707,0.707,-0.707], [-0.707,-0.707,-0.707] ];
  return poly;
};

const icosahedron = function() {
  const poly = new polyhedron();
  poly.name = "I";
  poly.face = [ [0,1,2], [0,2,3], [0,3,4], [0,4,5],
    [0,5,1], [1,5,7], [1,7,6], [1,6,2],
    [2,6,8], [2,8,3], [3,8,9], [3,9,4],
    [4,9,10], [4,10,5], [5,10,7], [6,7,11],
    [6,11,8], [7,10,11], [8,11,9], [9,11,10] ];

  poly.xyz = [ [0,0,1.176], [1.051,0,0.526],
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
   poly.face = [ [0,1,4,7,2], [0,2,6,9,3], [0,3,8,5,1],
      [1,5,11,10,4], [2,7,13,12,6], [3,9,15,14,8],
      [4,10,16,13,7], [5,8,14,17,11], [6,12,18,15,9],
      [10,11,17,19,16], [12,13,16,19,18], [14,15,18,19,17] ];
   poly.xyz = [ [0,0,1.07047], [0.713644,0,0.797878],
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
    poly.xyz.push([-cos(i*theta), -sin(i*theta),  -h]);
  }
  for (i = 0; i < n; i++) { // vertex #'s n to 2n-1 around other
    poly.xyz.push([-cos(i*theta), -sin(i*theta), h]);
  }

  poly.face.push(__range__(n-1, 0, true));  //top
  poly.face.push(__range__(n, 2*n, false)); //bottom
  for (i = 0; i < n; i++) { //n square sides
    poly.face.push([i, (i+1)%n, ((i+1)%n)+n, i+n]);
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
    poly.xyz.push([r * cos(i*theta), r * sin(i*theta), h]);
  }
  for (i = 0; i < n; i++) { // vertex #'s n...2n-1 around other
    poly.xyz.push([r * cos((i+0.5)*theta), r * sin((i+0.5)*theta), -h]);
  }

  poly.face.push(__range__(n-1, 0, true));   //top
  poly.face.push(__range__(n, (2*n)-1, true)); //bottom
  for (i = 0; i <= n-1; i++) { //2n triangular sides
    poly.face.push([i, (i+1)%n, i+n]);
    poly.face.push([i, i+n, ((((n+i)-1)%n)+n)]);
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
    poly.xyz.push([-cos(i*theta), -sin(i*theta), -0.2]);
  }
  poly.xyz.push([0,0, height]); // apex

  poly.face.push(__range__(n-1, 0, true)); // base
  for (i = 0; i < n; i++) { // n triangular sides
    poly.face.push([i, (i+1)%n, n]);
  }

  poly = canonicalXYZ(poly, 3);
  return poly;
};
