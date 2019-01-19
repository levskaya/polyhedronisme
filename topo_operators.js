// Polyhédronisme
//===================================================================================================
//
// A toy for constructing and manipulating polyhedra and other meshes
//
// Includes implementation of the conway polyhedral operators derived
// from code by mathematician and mathematical sculptor
// George W. Hart http://www.georgehart.com/
//
// Copyright 2019, Anselm Levskaya
// Released under the MIT License

//===================================================================================================
// Polyhedron Flagset Construct
//
// A Flag is an associative triple of a face index and two adjacent vertex indices,
// listed in geometric clockwise order (staring into the normal)
//
// Face_i -> V_i -> V_j
//
// They are a useful abstraction for defining topological transformations of the polyhedral mesh, as
// one can refer to vertices and faces that don't yet exist or haven't been traversed yet in the
// transformation code.
//
// A flag is similar in concept to a directed halfedge in halfedge data structures.
//
class polyflag {
  constructor() {
    this.flags = new Object(); // flags[face][vertex] = next vertex of flag; symbolic triples
    this.verts = new Object(); // XYZ coordinates
    this.xyzs = new Object(); // [symbolic names] holds vertex index
  }

  // Add a new vertex named "name" with coordinates "xyz".
  newV(name, xyz) {
    if (this.verts[name] === undefined) {
      this.verts[name] = 0;
      this.xyzs[name] = xyz;
    }
  }

  newFlag(facename, v1, v2) {
    if (this.flags[facename] === undefined) {
      this.flags[facename] = {};
    }
    this.flags[facename][v1] = v2;
  }

  topoly() {
    let i, v;
    const poly = new polyhedron();

    let ctr=0; // first number the vertices
    for (i in this.verts) {
      v = this.verts[i];
      poly.vertices[ctr]=this.xyzs[i]; // store in array
      this.verts[i] = ctr;
      ctr++;
    }

    ctr=0;
    for (i in this.flags) {
      var v0;
      const f = this.flags[i];
      poly.faces[ctr] = []; // new face
      // grab _any_ vertex as starting point
      for (let j in f) {
        v = f[j];
        v0 = v;
        break;
      }  // need just one
      // build face out of all the edge relations in the flag assoc array
      v = v0; // v moves around face
      poly.faces[ctr].push(this.verts[v]); //record index
      v = this.flags[i][v]; // goto next vertex
      let faceCTR=0;
      while (v !== v0) { // loop until back to start
        poly.faces[ctr].push(this.verts[v]);
        v = this.flags[i][v];
        faceCTR++;
        if (faceCTR>1000) { // necessary to prevent browser hangs on badly formed flagsets!
          console.log("Bad flag spec, have a neverending face:", i, this.flags[i]);
          break;
        }
      }
      ctr++;
    }

    poly.name = "unknown polyhedron";
    return poly;
  }
}


//===================================================================================================
// Polyhedron Operators
//===================================================================================================
// for each vertex of new polyhedron:
//     call newV(Vname, xyz) with a symbolic name and coordinates
// for each flag of new polyhedron:
//     call newFlag(Fname, Vname1, Vname2) with a symbolic name for the new face
//     and the symbolic name for two vertices forming an oriented edge
// ORIENTATION -must- be dealt with properly to make a manifold (correct) mesh.
// Specifically, no edge v1->v2 can ever be crossed in the -same direction- by
// two different faces
// 
// call topoly() to assemble flags into polyhedron structure by following the orbits
// of the vertex mapping stored in the flagset for each new face
// 
// set name as appropriate

// Kis(N)
// ------------------------------------------------------------------------------------------
// Kis (abbreviated from triakis) transforms an N-sided face into an N-pyramid rooted at the
// same base vertices.
// only kis n-sided faces, but n==0 means kis all.
//
const kisN = function(poly, n, apexdist){
  let i;
  if (!n) { n = 0; }
  if (!apexdist) { apexdist = 0.1; }
  console.log(`Taking kis of ${n===0 ? "" : n}-sided faces of ${poly.name}...`);

  const flag = new polyflag();
  for (i = 0; i < poly.vertices.length; i++) {
    // each old vertex is a new vertex
    const p = poly.vertices[i];
    flag.newV(`v${i}`, p);
  }

  const normals = poly.normals();
  const centers = poly.centers();
  let foundAny = false;                 // alert if don't find any
  for (i = 0; i < poly.faces.length; i++) {
    const f = poly.faces[i];
    let v1 = `v${f[f.length-1]}`;
    for (let v of f) {
      const v2 = `v${v}`;
      if ((f.length === n) || (n === 0)) {
        foundAny = true;
        const apex = `apex${i}`;
        const fname = `${i}${v1}`;
        flag.newV(apex, add(centers[i],mult(apexdist,normals[i]))); // new vertices in centers of n-sided face
        flag.newFlag(fname,   v1,   v2); // the old edge of original face
        flag.newFlag(fname,   v2, apex); // up to apex of pyramid
        flag.newFlag(fname, apex,   v1); // and back down again
      } else {
        flag.newFlag(`${i}`, v1, v2);  // same old flag, if non-n
      }
      v1=v2;
    }
  }  // current becomes previous

  if (!foundAny) {
    console.log(`No ${n}-fold components were found.`);
  }

  const newpoly = flag.topoly();
  newpoly.name = `k${n === 0 ? "" : n}${poly.name}`;
  return newpoly;
};


// Ambo
// ------------------------------------------------------------------------------------------
// The best way to think of the ambo operator is as a topological "tween" between a polyhedron
// and its dual polyhedron.  Thus the ambo of a dual polyhedron is the same as the ambo of the
// original. Also called "Rectify".
//
const ambo = function(poly){
  console.log(`Taking ambo of ${poly.name}...`);

  // helper func to insure unique names of midpoints
  const midName = function(v1, v2) { 
    if (v1<v2) { 
      return v1+"_"+v2; 
    } 
    else { 
      return v2+"_"+v1; 
    } 
  };

  const flag = new polyflag();

  // For each face f in the original poly
  for (let i = 0; i < poly.faces.length; i++) {
    const f = poly.faces[i];
    let [v1, v2] = f.slice(-2);
    for (let v3 of f) {
      if (v1 < v2) { // vertices are the midpoints of all edges of original poly
        flag.newV(midName(v1,v2), midpoint(poly.vertices[v1], poly.vertices[v2]));
      }
      // two new flags:
      // One whose face corresponds to the original f:
      flag.newFlag(`orig${i}`,  midName(v1,v2), midName(v2,v3));
      // Another flag whose face  corresponds to (the truncated) v2:
      flag.newFlag(`dual${v2}`, midName(v2,v3), midName(v1,v2));
      // shift over one
      [v1, v2] = [v2, v3];
    }
  }

  const newpoly = flag.topoly();
  newpoly.name = `a${poly.name}`;
  return newpoly;
};


// Gyro
// ----------------------------------------------------------------------------------------------
// This is the dual operator to "snub", i.e dual*Gyro = Snub.  It is a bit easier to implement
// this way.
//
// Snub creates at each vertex a new face, expands and twists it, and adds two new triangles to
// replace each edge.
//
const gyro = function(poly){
  let f, i, v;
  console.log(`Taking gyro of ${poly.name}...`);

  const flag = new polyflag();

  for (i = 0; i < poly.vertices.length; i++) {
    v = poly.vertices[i];
    flag.newV(`v${i}`, unit(v));
  }  // each old vertex is a new vertex

  const centers = poly.centers(); // new vertices in center of each face
  for (i = 0; i < poly.faces.length; i++) {
    f = poly.faces[i];
    flag.newV(`center${i}`, unit(centers[i]));
  }

  for (i = 0; i < poly.faces.length; i++) {
    f = poly.faces[i];
    let [v1, v2] = f.slice(-2);
    for (let j = 0; j < f.length; j++) {
      v = f[j];
      const v3 = v;
      flag.newV(v1+"~"+v2, oneThird(poly.vertices[v1],poly.vertices[v2]));  // new v in face
      const fname = i+"f"+v1;
      flag.newFlag(fname, `center${i}`,      v1+"~"+v2); // five new flags
      flag.newFlag(fname, v1+"~"+v2,  v2+"~"+v1);
      flag.newFlag(fname, v2+"~"+v1,  `v${v2}`);
      flag.newFlag(fname, `v${v2}`,     v2+"~"+v3);
      flag.newFlag(fname, v2+"~"+v3,  `center${i}`);
      [v1, v2] = [v2, v3];
    }
  }                       // shift over one

  const newpoly = flag.topoly();
  newpoly.name = `g${poly.name}`;
  return newpoly;
};


// Propellor
// ------------------------------------------------------------------------------------------
// builds a new 'skew face' by making new points along edges, 1/3rd the distance from v1->v2,
// then connecting these into a new inset face.  This breaks rotational symmetry about the
// faces, whirling them into gyres
//
const propellor = function(poly) {
  let i, v;
  console.log(`Taking propellor of ${poly.name}...`);

  const flag = new polyflag();

  for (i = 0; i < poly.vertices.length; i++) {
    v = poly.vertices[i];
    flag.newV(`v${i}`, unit(v));
  }  // each old vertex is a new vertex

  for (i = 0; i < poly.faces.length; i++) {
    const f = poly.faces[i];
    let [v1, v2] = f.slice(-2);
    for (v of f) {
      const v3 = `${v}`;
      flag.newV(v1+"~"+v2, oneThird(poly.vertices[v1], poly.vertices[v2]));  // new v in face, 1/3rd along edge
      const fname = `${i}f${v2}`;
      flag.newFlag(`v${i}`, v1+"~"+v2,  v2+"~"+v3); // five new flags
      flag.newFlag(fname,   v1+"~"+v2,  v2+"~"+v1);
      flag.newFlag(fname,   v2+"~"+v1,     `v${v2}`);
      flag.newFlag(fname,      `v${v2}`,  v2+"~"+v3);
      flag.newFlag(fname,   v2+"~"+v3,  v1+"~"+v2);
      [v1, v2] = [v2, v3];
    }
  } // shift over one

  const newpoly = flag.topoly();
  newpoly.name = `p${poly.name}`;
  return newpoly;
};


// Reflection
// ------------------------------------------------------------------------------------------
// geometric reflection through origin
const reflect = function(poly) {
  let i;
  console.log(`Taking reflection of ${poly.name}...`);
  // reflect each point through origin
  for (i = 0; i <= poly.vertices.length-1; i++) {
     poly.vertices[i] = mult(-1, poly.vertices[i]);
  }
  // repair clockwise-ness of faces
  for (i = 0; i <= poly.faces.length-1; i++) {
     poly.faces[i] = poly.faces[i].reverse();
  }
  poly.name = `r${poly.name}`;
  return poly;
};


// Dual
// ------------------------------------------------------------------------------------------------
// The dual of a polyhedron is another mesh wherein:
// - every face in the original becomes a vertex in the dual
// - every vertex in the original becomes a face in the dual
//
// So N_faces, N_vertices = N_dualfaces, N_dualvertices
//
// The new vertex coordinates are convenient to set to the original face centroids.
//
const dual = function(poly) {
  let f, i, v1, v2;
  console.log(`Taking dual of ${poly.name}...`);

  const flag = new polyflag();

  const face = []; // make table of face as fn of edge
  for (i = 0; i <= poly.vertices.length-1; i++) {
    face[i] = {};
  } // create empty associative table

  for (i = 0; i < poly.faces.length; i++) {
    f = poly.faces[i];
    v1 = f[f.length-1]; //previous vertex
    for (v2 of f) {
      // THIS ASSUMES that no 2 faces that share an edge share it in the same orientation!
      // which of course never happens for proper manifold meshes, so get your meshes right.
      face[v1][`v${v2}`] = `${i}`;
      v1=v2;
    }
  } // current becomes previous

  const centers = poly.centers();
  for (i = 0; i <= poly.faces.length-1; i++) {
    flag.newV(`${i}`,centers[i]);
  }

  for (i = 0; i < poly.faces.length; i++) {
    f = poly.faces[i];
    v1 = f[f.length-1]; //previous vertex
    for (v2 of f) {
      flag.newFlag(v1, face[v2][`v${v1}`], `${i}`);
      v1=v2;
    }
  } // current becomes previous

  const dpoly = flag.topoly(); // build topological dual from flags

  // match F index ordering to V index ordering on dual
  const sortF = [];
  for (f of dpoly.faces) {
    const k = intersect(poly.faces[f[0]], poly.faces[f[1]], poly.faces[f[2]]);
    sortF[k] = f;
  }
  dpoly.faces = sortF;

  if (poly.name[0] !== "d") {
    dpoly.name = `d${poly.name}`;
  } else {
    dpoly.name = poly.name.slice(1);
  }

  return dpoly;
};


// Chamfer
// ----------------------------------------------------------------------------------------
// A truncation along a polyhedron's edges.
// Chamfering or edge-truncation is similar to expansion, moving faces apart and outward,
// but also maintains the original vertices. Adds a new hexagonal face in place of each
// original edge.
// A polyhedron with e edges will have a chamfered form containing 2e new vertices,
// 3e new edges, and e new hexagonal faces. -- Wikipedia
// See also http://dmccooey.com/polyhedra/Chamfer.html
//
// The dist parameter could control how deeply to chamfer.
// But I'm not sure about implementing that yet.
//
// Q: what is the dual operation of chamfering? I.e.
// if cX = dxdX, and xX = dcdX, what operation is x?

// We could "almost" do this in terms of already-implemented operations:
// cC = t4daC = t4jC, cO = t3daO, cD = t5daD, cI = t3daI
// But it doesn't work for cases like T.

const chamfer = function(poly, dist) {
  console.log(`Taking chamfer of ${poly.name}...`);

  if (!dist) { dist = 0.5; }

  const flag = new polyflag();

  const normals = poly.normals();

  // For each face f in the original poly
  for (let i = 0; i < poly.faces.length; i++) {
    const f = poly.faces[i];
    let v1 = f[f.length-1];
    let v1new = i + "_" + v1;

    for (let v2 of f) {
      // TODO: figure out what distances will give us a planar hex face.
      // Move each old vertex further from the origin.
      flag.newV(v2, mult(1.0 + dist, poly.vertices[v2]));
      // Add a new vertex, moved parallel to normal.
      const v2new = i + "_" + v2;
      flag.newV(v2new, add(poly.vertices[v2], mult(dist*1.5, normals[i])));
      // Four new flags:
      // One whose face corresponds to the original face:
      flag.newFlag(`orig${i}`, v1new, v2new);
      // And three for the edges of the new hexagon:
      const facename = (v1<v2 ? `hex${v1}_${v2}` : `hex${v2}_${v1}`);
      flag.newFlag(facename, v2, v2new);
      flag.newFlag(facename, v2new, v1new);
      flag.newFlag(facename, v1new, v1);
      v1 = v2;
      v1new = v2new;
    }
  }

  const newpoly = flag.topoly();
  newpoly.name = `c${poly.name}`;
  return newpoly;
};


// Whirl
// ----------------------------------------------------------------------------------------------
// Gyro followed by truncation of vertices centered on original faces.
// This create 2 new hexagons for every original edge.
// (https://en.wikipedia.org/wiki/Conway_polyhedron_notation#Operations_on_polyhedra)
//
// Possible extension: take a parameter n that means only whirl n-sided faces.
// If we do that, the flags marked #* below will need to have their other sides
// filled in one way or another, depending on whether the adjacent face is
// whirled or not.

const whirl = function(poly, n) {
  let i, v;
  console.log(`Taking whirl of ${poly.name}...`);
  if (!n) { n = 0; }
  
  const flag = new polyflag();

  // each old vertex is a new vertex
  for (i = 0; i < poly.vertices.length; i++) {
    v = poly.vertices[i];
    flag.newV(`v${i}`, unit(v));
  }

  // new vertices around center of each face
  const centers = poly.centers();
  //for f,i in poly.face
  //  # Whirl: use "center"+i+"~"+v1
  //  flag.newV "center"+i+"~"+v1, unit(centers[i])

  for (i = 0; i < poly.faces.length; i++) {
    const f = poly.faces[i];
    let [v1, v2] = f.slice(-2);
    for (let j = 0; j < f.length; j++) {
      v = f[j];
      const v3 = v;
      // New vertex along edge
      const v1_2 = oneThird(poly.vertices[v1],poly.vertices[v2]);
      flag.newV(v1+"~"+v2, v1_2);
      // New vertices near center of face
      const cv1name = `center${i}~${v1}`;
      const cv2name = `center${i}~${v2}`;
      flag.newV(cv1name, unit(oneThird(centers[i], v1_2))); 
      const fname = i+"f"+v1;
      // New hexagon for each original edge
      flag.newFlag(fname, cv1name,   v1+"~"+v2);
      flag.newFlag(fname, v1+"~"+v2, v2+"~"+v1); //*
      flag.newFlag(fname, v2+"~"+v1, `v${v2}`);  //*
      flag.newFlag(fname, `v${v2}`,  v2+"~"+v3); //*
      flag.newFlag(fname, v2+"~"+v3, cv2name);
      flag.newFlag(fname, cv2name,   cv1name);
      // New face in center of each old face      
      flag.newFlag(`c${i}`, cv1name, cv2name);
      
      [v1, v2] = [v2, v3];
    }
  } // shift over one

  const newpoly = flag.topoly();
  newpoly.name = `w${poly.name}`;
  return newpoly;
};
  
// Quinto

const quinto = function(poly){
  console.log(`Taking quinto of ${poly.name}...`);

  // helper func to insure unique names of midpoints
  const midName = function(v1, v2) { 
    if (v1<v2) { return v1+"_"+v2; } 
    else       { return v2+"_"+v1; } 
  };

  const flag = new polyflag();

  // For each face f in the original poly
  for (let i = 0; i < poly.faces.length; i++) {
    const f = poly.faces[i];
    centroid = calcCentroid(f.map(idx=>poly.vertices[idx]))
    // walk over face vertex-triplets
    let [v1, v2] = f.slice(-2);
    for (let v3 of f) {
      // for each face-corner, we make two new points:
      midpt = midpoint(poly.vertices[v1], poly.vertices[v2])
      innerpt = midpoint(midpt, centroid)
      flag.newV(midName(v1,v2), midpt);
      flag.newV(`inner_${i}_` + midName(v1,v2), innerpt);
      // and add the old corner-vertex
      flag.newV(`${v2}`, poly.vertices[v2]);
    
      // pentagon for each vertex in original face
      flag.newFlag(`f${i}_${v2}`, `inner_${i}_`+midName(v1, v2), midName(v1, v2));
      flag.newFlag(`f${i}_${v2}`, midName(v1, v2), `${v2}`);
      flag.newFlag(`f${i}_${v2}`, `${v2}`, midName(v2, v3));
      flag.newFlag(`f${i}_${v2}`, midName(v2, v3), `inner_${i}_`+midName(v2, v3));
      flag.newFlag(`f${i}_${v2}`, `inner_${i}_`+midName(v2, v3), `inner_${i}_`+midName(v1, v2));

      // inner rotated face of same vertex-number as original
      flag.newFlag(`f_in_${i}`, `inner_${i}_`+midName(v1, v2), `inner_${i}_`+midName(v2, v3));

      // shift over one
      [v1, v2] = [v2, v3];
    }
  }

  const newpoly = flag.topoly();
  newpoly.name = `q${poly.name}`;
  //topolog(newpoly);
  //console.log(newpoly);
  return newpoly;
};


// Triangular Subdivision Operator
// ----------------------------------------------------------------------------------------------
// limited version of the Goldberg-Coxeter u_n operator for triangular meshes
// We subdivide manually here, instead of using the usual flag machinery.
const trisub = function(poly, n) {
  console.log(`Taking trisub of ${poly.name}...`);
  if (!n) { n = 2; }
  
  // No-Op for non-triangular meshes.
  for (let fn = 0; fn < poly.faces.length; fn++) {
    if(poly.faces[fn].length != 3){
      return poly;
    }
  }

  // Calculate redundant set of new vertices for subdivided mesh.
  let newVs=[];
  let vmap={};
  let pos = 0;
  for (let fn = 0; fn < poly.faces.length; fn++) {
    const f = poly.faces[fn];
    let [i1, i2, i3] = f.slice(-3);
    v1 = poly.vertices[i1];
    v2 = poly.vertices[i2];
    v3 = poly.vertices[i3];
    v21 = sub(v2, v1);
    v31 = sub(v3, v1);
    for (let i = 0; i <= n; i++) {
      for (let j = 0; j+i <= n; j++) {
        let v = add(add(v1, mult(i * 1.0 / n, v21)), mult(j * 1.0 / n, v31));
        vmap[`v${fn}-${i}-${j}`] = pos++;
        newVs.push(v);
      }
    }
  }

  // The above vertices are redundant along original edges, 
  // we need to build an index map into a uniqueified list of them.
  // We identify vertices that are closer than a certain epsilon distance.
  const epsilon_close = 1.0e-8;
  let uniqVs = [];
  let newpos = 0;
  let uniqmap = {};
  for (const [i, v] of newVs.entries()) {
    if (i in uniqmap) { continue; } // already mapped
    uniqmap[i] = newpos;
    uniqVs.push(v);
    for(let j = i+1; j < newVs.length; j++) {
      w = newVs[j];
      if (mag(sub(v, w)) < epsilon_close) {
        uniqmap[j] = newpos;
      }
    }
    newpos++;
  }

  let faces = [];
  for (fn = 0; fn < poly.faces.length; fn++) {
    for (let i = 0; i < n; i++) {
      for (let j = 0; j+i < n; j++) {
        faces.push([uniqmap[vmap[`v${fn}-${i}-${j}`]], 
                    uniqmap[vmap[`v${fn}-${i+1}-${j}`]], 
                    uniqmap[vmap[`v${fn}-${i}-${j+1}`]]])
      }
    }
    for (let i = 1; i < n; i++) {
      for (let j = 0; j+i < n; j++) {
        faces.push([uniqmap[vmap[`v${fn}-${i}-${j}`]], 
                    uniqmap[vmap[`v${fn}-${i}-${j+1}`]], 
                    uniqmap[vmap[`v${fn}-${i-1}-${j+1}`]]])
      }
    }
  }

  // Create new polygon out of faces and unique vertices.
  const newpoly = new polyhedron();
  newpoly.name = `u${n}${poly.name}`;
  newpoly.faces = faces;
  newpoly.vertices = uniqVs; 

  return newpoly;
};


// insetN
// ------------------------------------------------------------------------------------------
const insetN = function(poly, n, inset_dist, popout_dist){
  let f, i, v;
  if (!n) { n = 0; }
  if (!inset_dist) { inset_dist = 0.5; }
  if (!popout_dist) { popout_dist = -0.2; }

  console.log(`Taking inset of ${n===0 ? "" : n}-sided faces of ${poly.name}...`);

  const flag = new polyflag();
  for (i = 0; i < poly.vertices.length; i++) {
    // each old vertex is a new vertex
    const p = poly.vertices[i];
    flag.newV(`v${i}`, p);
  }

  const normals = poly.normals();
  const centers = poly.centers();
  for (i = 0; i < poly.faces.length; i++) { //new inset vertex for every vert in face
    f = poly.faces[i];
    if ((f.length === n) || (n === 0)) {
      for (v of f) {
        flag.newV(`f${i}v${v}`, add(tween(poly.vertices[v],centers[i],inset_dist), 
                                    mult(popout_dist,normals[i])));
      }
    }
  }

  let foundAny = false;    // alert if don't find any
  for (i = 0; i < poly.faces.length; i++) {
    f = poly.faces[i];
    let v1 = `v${f[f.length-1]}`;
    for (v of f) {
      const v2 = `v${v}`;
      if ((f.length === n) || (n === 0)) {
        foundAny = true;
        const fname = i + v1;
        flag.newFlag(fname,      v1,       v2);
        flag.newFlag(fname,      v2,       `f${i}${v2}`);
        flag.newFlag(fname, `f${i}${v2}`,  `f${i}${v1}`);
        flag.newFlag(fname, `f${i}${v1}`,  v1);
        //new inset, extruded face
        flag.newFlag(`ex${i}`, `f${i}${v1}`,  `f${i}${v2}`);
      } else {
        flag.newFlag(i, v1, v2);  // same old flag, if non-n
      }
      v1=v2;
    }
  }  // current becomes previous

  if (!foundAny) {
    console.log(`No ${n}-fold components were found.`);
  }

  const newpoly = flag.topoly();
  newpoly.name = `n${n === 0 ? "" : n}${poly.name}`;
  return newpoly;
};


// extrudeN
// ------------------------------------------------------------------------------------------
const extrudeN = function(poly, n){
  let f, i, v;
  if (!n) { n = 0; }
  console.log(`Taking extrusion of ${n===0 ? "" : n}-sided faces of ${poly.name}...`);

  const flag = new polyflag();
  for (i = 0; i < poly.vertices.length; i++) {
    // each old vertex is a new vertex
    const p = poly.vertices[i];
    flag.newV(`v${i}`, p);
  }

  const normals = poly.normals();
  const centers = poly.centers();
  for (i = 0; i < poly.faces.length; i++) { //new inset vertex for every vert in face
    f = poly.faces[i];
    if ((f.length === n) || (n === 0)) {
      for (v of f) {
        //flag.newV "f"+i+"v"+v, add(midpoint(poly.vertices[v],centers[i]),mult(-0.2,normals[i]))
        flag.newV(`f${i}v${v}`, add(poly.vertices[v], mult(0.3,normals[i])));
      }
    }
  }

  let foundAny = false;                 // alert if don't find any
  for (i = 0; i < poly.faces.length; i++) {
    f = poly.faces[i];
    let v1 = `v${f[f.length-1]}`;
    for (v of f) {
      const v2 = `v${v}`;
      if ((f.length === n) || (n === 0)) {
        foundAny = true;
        //fname = i+v1
        flag.newFlag(i+v1,       v1,       v2);
        flag.newFlag(i+v1,       v2, `f${i}${v2}`);
        flag.newFlag(i+v1, `f${i}${v2}`, `f${i}${v1}`);
        flag.newFlag(i+v1, `f${i}${v1}`,       v1);
        //new inset, extruded face
        flag.newFlag(`ex${i}`, `f${i}${v1}`,  `f${i}${v2}`);
      } else {
        flag.newFlag(i, v1, v2);  // same old flag, if non-n
      }
      v1=v2;
    }
  }  // current becomes previous

  if (!foundAny) {
    console.log(`No ${n}-fold components were found.`);
  }

  const newpoly = flag.topoly();
  newpoly.name = `x${n === 0 ? "" : n}${poly.name}`;
  return newpoly;
};


// hollow / skeletonize
// ------------------------------------------------------------------------------------------
const hollow = function(poly, n, inset_dist, thickness){
  let f, i, v;
  if (!n) { n = 0; }
  if (!inset_dist) { inset_dist = 0.5; }
  if (!thickness) { thickness = 0.2; }

  console.log(`Skeletonizing ${n===0 ? "" : n}-sided faces of ${poly.name}...`);

  const dualnormals = dual(poly).normals();
  const normals = poly.normals();
  const centers = poly.centers();

  const flag = new polyflag();
  for (i = 0; i < poly.vertices.length; i++) {
    // each old vertex is a new vertex
    const p = poly.vertices[i];
    flag.newV(`v${i}`, p);
    flag.newV(`downv${i}`,  add(p,mult(-1*thickness,dualnormals[i])));
  }
  // new inset vertex for every vert in face
  for (i = 0; i < poly.faces.length; i++) {
    //if f.length is n or n is 0
    f = poly.faces[i];
    for (v of f) {
      flag.newV(`fin${i}v${v}`, tween(poly.vertices[v],centers[i],inset_dist));
      flag.newV(`findown${i}v${v}`, add(tween(poly.vertices[v],centers[i],inset_dist), 
                                        mult(-1*thickness,normals[i])));
    }
  }

  //foundAny = false    # alert if we don't find any N-ary faces
  for (i = 0; i < poly.faces.length; i++) {
    f = poly.faces[i];
    let v1 = `v${f[f.length-1]}`;
    for (v of f) {
      const v2 = `v${v}`;
      //if f.length is n or n is 0
      const foundAny = true;
      let fname = i + v1;
      flag.newFlag(fname,      v1,       v2);
      flag.newFlag(fname,      v2,       `fin${i}${v2}`);
      flag.newFlag(fname, `fin${i}${v2}`,  `fin${i}${v1}`);
      flag.newFlag(fname, `fin${i}${v1}`,  v1);

      fname = `sides${i}${v1}`;
      flag.newFlag(fname, `fin${i}${v1}`,     `fin${i}${v2}`);
      flag.newFlag(fname, `fin${i}${v2}`,     `findown${i}${v2}`);
      flag.newFlag(fname, `findown${i}${v2}`, `findown${i}${v1}`);
      flag.newFlag(fname, `findown${i}${v1}`, `fin${i}${v1}`);

      fname = `bottom${i}${v1}`;
      flag.newFlag(fname,  `down${v2}`,      `down${v1}`);
      flag.newFlag(fname,  `down${v1}`,      `findown${i}${v1}`);
      flag.newFlag(fname,  `findown${i}${v1}`, `findown${i}${v2}`);
      flag.newFlag(fname,  `findown${i}${v2}`, `down${v2}`);

      //new inset, extruded face
      //flag.newFlag "ex"+i, "fin"+i+v1,  "fin"+i+v2
      //else
      //  flag.newFlag i, v1, v2  # same old flag, if non-n
      v1=v2;
    }
  }  // current becomes previous

  //if not foundAny
  //  console.log "No #{n}-fold components were found."

  const newpoly = flag.topoly();
  //newpoly.name = "h" + (if n is 0 then "" else n) + poly.name
  newpoly.name = `h${poly.name}`;
  return newpoly;
};


// StellaN
// ------------------------------------------------------------------------------------------
const stellaN = function(poly){
  let i;
  console.log(`Taking stella of ${poly.name}...`);

  const centers = poly.centers();  // calculate face centers

  const flag = new polyflag();
  for (i = 0; i < poly.vertices.length; i++) {
    const p = poly.vertices[i];
    flag.newV(`v${i}`, p);
  }      // each old vertex is a new vertex

  // iterate over triplets of faces v1,v2,v3
  for (i = 0; i < poly.faces.length; i++) {
    const f = poly.faces[i];
    let v1 = `v${f[f.length-2]}`;
    let v2 = `v${f[f.length-1]}`;
    let vert1 = poly.vertices[f[f.length-2]];
    let vert2 = poly.vertices[f[f.length-1]];
    for (let v of f) {
      const v3 = `v${v}`;
      const vert3 = poly.vertices[v];
      const v12=v1+"~"+v2; // names for "oriented" midpoints
      const v21=v2+"~"+v1;
      const v23=v2+"~"+v3;

      // on each Nface, N new points inset from edge midpoints towards center = "stellated" points
      flag.newV(v12, midpoint( midpoint(vert1,vert2), centers[i] ));

      // inset Nface made of new, stellated points
      flag.newFlag(`in${i}`,      v12,       v23);

      // new tri face constituting the remainder of the stellated Nface
      flag.newFlag(`f${i}${v2}`,      v23,      v12);
      flag.newFlag(`f${i}${v2}`,       v12,      v2);
      flag.newFlag(`f${i}${v2}`,      v2,      v23);

      // one of the two new triangles replacing old edge between v1->v2
      flag.newFlag(`f${v12}`,     v1,        v21);
      flag.newFlag(`f${v12}`,     v21,       v12);
      flag.newFlag(`f${v12}`,      v12,       v1);

      [v1, v2] = [v2, v3];  // current becomes previous
      [vert1, vert2] = [vert2, vert3];
    }
  }

  const newpoly = flag.topoly();
  newpoly.name = `l${poly.name}`;
  return newpoly;
};
