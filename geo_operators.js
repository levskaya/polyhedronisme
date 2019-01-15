/*
 * decaffeinate suggestions:
 * DS102: Remove unnecessary code created because of implicit returns
 * DS201: Simplify complex destructure assignments
 * DS202: Simplify dynamic range loops
 * DS205: Consider reworking code to avoid use of IIFEs
 * Full docs: https://github.com/decaffeinate/decaffeinate/blob/master/docs/suggestions.md
 */
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
// Canonicalization Algorithms
//===================================================================================================

// Slow Canonicalization Algorithm
// ----------------------------------------------------------------
//
// This algorithm has some convergence problems, what really needs to be done is to
// sum the three forcing factors together as a conherent force and to use a half-decent
// integrator to make sure that it converges well as opposed to the current hack of
// ad-hoc stability multipliers.  Ideally one would implement a conjugate gradient
// descent or similar pretty thing.
//
// Only try to use this on convex polyhedra that have a chance of being canonicalized,
// otherwise it will probably blow up the geometry.  A much trickier / smarter seed-symmetry
// based geometrical regularizer should be used for fancier/weirder polyhedra.

// adjusts vertices on edges such that each edge is tangent to an origin sphere
const tangentify = function(xyzs, edges) {
  // hack to improve convergence
  const STABILITY_FACTOR = 0.1; 
  // copy vertices
  const newVs = copyVecArray(xyzs);
  for (let e of edges) {
    // the point closest to origin
    const t = tangentPoint( newVs[e[0]], newVs[e[1]] );
    // adjustment from sphere
    const c = mult(((STABILITY_FACTOR*1)/2)*(1-sqrt(dot(t,t))), t);
    newVs[e[0]] = add(newVs[e[0]], c);
    newVs[e[1]] = add(newVs[e[1]], c);
  }
  return newVs;
};

// recenters entire polyhedron such that center of mass is at origin
const recenter = function(xyzs, edges) {
  const edgecenters = ((() => {
    const result = [];
    for (let [a,b] of edges) {
      result.push(tangentPoint(xyzs[a], xyzs[b]));
    }
    return result;
  })()); //centers of edges
  let polycenter = [0,0,0];
  // sum centers to find center of gravity
  for (let v of edgecenters) { 
    polycenter = add(polycenter, v);
  }
  polycenter = mult(1/edges.length, polycenter);
  // subtract off any deviation from center
  return _.map(xyzs, x=> sub(x, polycenter));
};

// rescales maximum radius of polyhedron to 1
const rescale = function(xyzs) {
  const polycenter = [0,0,0];
  const maxExtent = _.max(_.map(xyzs, x=>mag(x)));
  const s = 1 / maxExtent;
  return _.map(xyzs, x=>[s*x[0],s*x[1],s*x[2]]);
};

// adjusts vertices in each face to improve its planarity
const planarize = function(xyzs, faces) {
  let v;
  const STABILITY_FACTOR = 0.1; // Hack to improve convergence
  const newVs = copyVecArray(xyzs); // copy vertices
  for (var f of faces) {
    const coords = ((() => {
      const result = [];
      for (v of f) {         
        result.push(xyzs[v]);
      }
      return result;
    })());
    let n = normal(coords); // find avg of normals for each vertex triplet
    const c = calcCentroid(coords); // find planar centroid
    if (dot(n,c) < 0) { // correct sign if needed
      n = mult(-1.0,n);
    }
    for (v of f) {  // project (vertex - centroid) onto normal, subtract off this component
      newVs[v] = add(newVs[v], mult(dot(mult(STABILITY_FACTOR, n), sub(c, xyzs[v])), n) );
    }
  }
  return newVs;
};

// combines above three constraint adjustments in iterative cycle
const canonicalize = function(poly, Niter) {
  if (!Niter) { Niter = 1; }
  console.log(`Canonicalizing ${poly.name}...`);
  const faces = poly.face;
  const edges = poly.edges();
  let newVs = poly.xyz;
  let maxChange=1.0; // convergence tracker
  for (let i = 0, end = Niter, asc = 0 <= end; asc ? i <= end : i >= end; asc ? i++ : i--) {
    const oldVs = copyVecArray(newVs); //copy vertices
    newVs = tangentify(newVs, edges);
    newVs = recenter(newVs, edges);
    newVs = planarize(newVs, faces);
    maxChange = _.max(_.map(_.zip(newVs,oldVs), 
    function(...args){ 
      const [x,y] = args[0]; 
      return mag(sub(x,y));  
    }));
    if (maxChange < 1e-8) {
      break;
    }
  }
  // one should now rescale, but not rescaling here makes for very interesting numerical
  // instabilities that make interesting mutants on multiple applications...
  // more experience will tell what to do
  //newVs = rescale(newVs)
  console.log(`[canonicalization done, last |deltaV|=${maxChange}]`);
  const newpoly = new polyhedron(newVs, poly.face, poly.name);
  console.log("canonicalize" , newpoly);
  return newpoly;
};

// Hacky Canonicalization Algorithm
// --------------------------------------------------------------------
// Using center of gravity of vertices for each face to planarize faces

// get the spherical reciprocals of face centers
const reciprocalC = function(poly) {
  const centers = poly.centers();
  for (let c of centers) {
    c = mult(1.0/dot(c,c), c);
  }
  return centers;
};

// make array of vertices reciprocal to given planes
const reciprocalN = function(poly) {
  const ans = [];
  for (let f of poly.face) { //for each face
    let centroid    = [0,0,0]; // running sum of vertex coords
    let normalV     = [0,0,0]; // running sum of normal vectors
    let avgEdgeDist =    0.0;  // running sum for avg edge distance

    let [v1, v2] = f.slice(-2);
    for (let v3 of f) {
      centroid     = add(centroid, poly.xyz[v3]);
      normalV      = add(normalV, orthogonal(poly.xyz[v1], poly.xyz[v2], poly.xyz[v3]));
      avgEdgeDist += edgeDist(poly.xyz[v1], poly.xyz[v2]);
      [v1, v2] = [v2, v3];
    } // shift over one

    centroid    = mult(1.0/f.length, centroid);
    normalV     = unit(normalV);
    avgEdgeDist = avgEdgeDist / f.length;
    const tmp   = reciprocal(mult(dot(centroid, normalV), normalV)); // based on face
    ans.push(mult((1 + avgEdgeDist) / 2, tmp));
  } // edge correction

  return ans;
};

const canonicalXYZ = function(poly, nIterations) {
  if (!nIterations) { nIterations = 1; }
  const dpoly = dual(poly);
  console.log(`Pseudo-canonicalizing ${poly.name}...`);

  // iteratively reciprocate face normals
  for (let count = 0, end = nIterations; count < end; count++) {
    dpoly.xyz = reciprocalN(poly);
    poly.xyz  = reciprocalN(dpoly);
  }

  return new polyhedron(poly.xyz, poly.face, poly.name);
};


// quick planarization
const adjustXYZ = function(poly, nIterations) {
  if (!nIterations) { nIterations = 1; }
  const dpoly = dual(poly); // v's of dual are in order of arg's f's
  console.log(`Planarizing ${poly.name}...`);

  for (let count = 0, end = nIterations; count < end; count++) {
    // reciprocate face centers
    dpoly.xyz = reciprocalC(poly);
    poly.xyz  = reciprocalC(dpoly);
  }

  return new polyhedron(poly.xyz, poly.face, poly.name);
};


