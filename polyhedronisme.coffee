#===================================================================================================
# Polyhedronisme
#
# A toy for constructing and manipulating polyhedra and other meshes
#
# Includes implementation of the conway polyhedral operators derived
# from code by mathematician and mathematical sculptor
# George W. Hart http://www.georgehart.com/
#
# Copyright 2011, Anselm Levskaya
# Released under the MIT License

#===================================================================================================
# Math / Vector / Matrix Functions
#===================================================================================================

# Math is primal, people.
# import math functions to local namespace
random = Math.random
round = Math.round
floor = Math.floor
sqrt  = Math.sqrt
sin   = Math.sin
cos   = Math.cos
tan   = Math.tan
pow   = Math.pow
abs   = Math.abs
PI    = Math.PI

# for python-style enumerated for-in loops
enumerate = (ar) ->  [i,ar[i]] for i in [0..ar.length-1]

# often useful
randomchoice = (array)->
  n = floor(random()*array.length)
  array[n]

# scalar multiplication
mult = (c,vec) -> [c*vec[0],c*vec[1],c*vec[2]]

# element-wise multiply
_mult = (vec1, vec2) -> [vec1[0]*vec2[0],vec1[1]*vec2[1],vec1[2]*vec2[2]]

# 3d vector addition
add = (vec1, vec2) -> [vec1[0]+vec2[0],vec1[1]+vec2[1],vec1[2]+vec2[2]]

# 3d vector subtraction
sub = (vec1, vec2) -> [vec1[0]-vec2[0],vec1[1]-vec2[1],vec1[2]-vec2[2]]

# 3d dot product
dot = (vec1, vec2) -> vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2]

# vector norm
mag = (vec) -> sqrt(dot(vec,vec))

# vector magnitude squared
mag2 = (vec) -> dot(vec,vec)

# makes vector unit length
unit = (vec) -> mult( 1/sqrt(mag2(vec)), vec)

# cross product d1 x d2
cross = (d1, d2) -> [ d1[1]*d2[2]-d1[2]*d2[1], d1[2]*d2[0]-d1[0]*d2[2],  d1[0]*d2[1]-d1[1]*d2[0] ]

midpoint = (vec1, vec2) -> mult(1/2.0,add(vec1,vec2))

# parametric segment between vec1, vec2 w. parameter t ranging from 0 to 1
tween = (vec1,vec2,t) -> [ (1-t)*vec1[0] + t*vec2[0], (1-t)*vec1[1] + t*vec2[1], (1-t)*vec1[2] + t*vec2[2] ]

# uses above to go one-third of the way along vec1->vec2 line
oneThird = (vec1, vec2) -> tween(vec1, vec2, 1/3.0)

# reflect 3vec in unit sphere, spherical reciprocal
reciprocal = (vec) -> mult( 1.0/mag2(vec), vec)

# point where line v1...v2 tangent to an origin sphere
tangentPoint= (v1,v2) ->
  d = sub v2, v1
  sub v1, mult(dot(d,v1)/mag2(d),d)

# distance of line v1...v2 to origin
edgeDist = (v1,v2) ->
  sqrt mag2(tangentPoint v1, v2)

# find vector orthogonal to plane of 3 pts
# -- do the below algos assume this be normalized or not?
orthogonal = (v1,v2,v3) ->
  # adjacent edge vectors
  d1 = sub v2, v1
  d2 = sub v3, v2
  # cross product
  cross d1, d2

# find element common to 3 sets by brute force search
intersect = (set1, set2, set3) ->
  for s1 in set1
    for s2 in set2
      if s1 is s2
        for s3 in set3
          if s1 is s3
            return s1
  return null

# copies array of arrays by value (deep copy)
copyVecArray = (vecArray)->
  newVecArray = new Array(vecArray.length)
  for i in [0..vecArray.length-1]
    newVecArray[i] = vecArray[i][0..]
  newVecArray

# 3d matrix vector multiply
mv3 = (mat,vec) ->
  #Ghetto custom def of matrix-vector mult
  #example matrix: [[a,b,c],[d,e,f],[g,h,i]]
  [mat[0][0]*vec[0]+mat[0][1]*vec[1]+mat[0][2]*vec[2],
  mat[1][0]*vec[0]+mat[1][1]*vec[1]+mat[1][2]*vec[2],
  mat[2][0]*vec[0]+mat[2][1]*vec[1]+mat[2][2]*vec[2]]

# 3d matrix matrix multiply
mm3 = (A,B) ->
  [[A[0][0]*B[0][0]+A[0][1]*B[1][0]+A[0][2]*B[2][0],
  A[0][0]*B[0][1]+A[0][1]*B[1][1]+A[0][2]*B[2][1],
  A[0][0]*B[0][2]+A[0][1]*B[1][2]+A[0][2]*B[2][2]],
  [A[1][0]*B[0][0]+A[1][1]*B[1][0]+A[1][2]*B[2][0],
  A[1][0]*B[0][1]+A[1][1]*B[1][1]+A[1][2]*B[2][1],
  A[1][0]*B[0][2]+A[1][1]*B[1][2]+A[1][2]*B[2][2]],
  [A[2][0]*B[0][0]+A[2][1]*B[1][0]+A[2][2]*B[2][0],
  A[2][0]*B[0][1]+A[2][1]*B[1][1]+A[2][2]*B[2][1],
  A[2][0]*B[0][2]+A[2][1]*B[1][2]+A[2][2]*B[2][2]]]

# Rotation Matrix
# Totally ghetto, not at all in agreement with euler angles!
rotm = (phi,theta,psi)->
    xy_mat = [
      [cos(phi), -1.0*sin(phi),  0.0],
      [sin(phi),      cos(phi),  0.0],
      [0.0,                0.0,  1.0]]
    yz_mat = [
      [cos(theta), 0, -1.0*sin(theta)],
      [         0, 1,               0],
      [sin(theta), 0,      cos(theta)]]
    xz_mat = [
      [1.0,        0,             0],
      [  0, cos(psi), -1.0*sin(psi)],
      [  0, sin(psi),      cos(psi)]]

    mm3(xz_mat, mm3(yz_mat,xy_mat))

perspT = (vec3, max_real_depth, min_real_depth, desired_ratio, desired_length) ->
  # assumes world's been rotated appropriately such that Z is depth
  # scales perspective such that inside depth regions min_real_depth <--> max_real_depth
  # perspective lengths vary no more than:   desired_ratio
  # with target dimension of roughly length: desired_length
  z0          = (max_real_depth * desired_ratio - min_real_depth)/(1-desired_ratio)
  scalefactor =  desired_length * desired_ratio/(1-desired_ratio)

  # projected [X, Y]
  [scalefactor*vec3[0]/(vec3[2]+z0), scalefactor*vec3[1]/(vec3[2]+z0)]



#===================================================================================================
# Polyhedra Functions
#===================================================================================================
#
# Topology stored as set of "faces."  Each face is list of n vertex indices
# corresponding to one n-sided face.  Vertices listed clockwise as seen from outside.

faceToEdges = (face) ->
  edges = []
  [v1] = face[-1..]
  for v2 in face
    edges.push [v1,v2]
    v1 = v2
  edges

class polyhedron
  constructor: () ->      # constructor of initially null polyhedron
    @face = new Array()   # array of faces.          face.length = # faces
    @xyz  = new Array()   # array of vertex coords.  xyz.length = # of vertices
    @name = "null polyhedron"

  data: () ->   # informative string
    nEdges = @face.length + @xyz.length - 2 # E = V + F - 2
    "(#{@face.length} faces, #{nEdges} edges, #{@xyz.length} vertices)"

  getEdges: () ->
    finalset={}
    uniqedges=[]
    alledges = _(@face).map(faceToEdges)
    for edgeset in alledges
      for e in edgeset
        if e[0] < e[1]
          [a,b] = e
        else
          [b,a] = e
        finalset[a+'~'+b] = e
    for hash,e of finalset
      uniqedges.push e

    #return edges
    uniqedges

# produces vanilla OBJ files for import into 3d apps
toOBJ = (poly) ->
  objstr="#Produced by polyHédronisme http://levskaya.github.com/polyhedronisme\n"
  objstr+="group poly\n"
  objstr+="#vertices\n"
  for v in poly.xyz
    objstr += "v #{v[0]} #{v[1]} #{v[2]}\n"

  objstr += "#normal vector defs \n"
  for f in poly.face
    norm = normal(poly.xyz[v] for v in f)
    objstr += "vn #{norm[0]} #{norm[1]} #{norm[2]}\n"

  objstr += "#face defs \n"
  for [i,f] in enumerate(poly.face)
    objstr += "f "
    for v in f
      objstr += "#{v+1}//#{i+1} "
    objstr += "\n"

  objstr

#===================================================================================================
# Primitive Polyhedra Seeds
#===================================================================================================

tetrahedron = () ->
  poly = new polyhedron()
  poly.name = "T"
  poly.face = [ [0,1,2], [0,2,3], [0,3,1], [1,3,2] ]
  poly.xyz  = [ [1.0,1.0,1.0], [1.0,-1.0,-1.0], [-1.0,1.0,-1.0], [-1.0,-1.0,1.0] ]
  poly

octahedron = () ->
  poly = new polyhedron()
  poly.name = "O"
  poly.face = [ [0,1,2], [0,2,3], [0,3,4], [0,4,1], [1,4,5], [1,5,2], [2,5,3], [3,5,4] ]
  poly.xyz  = [ [0,0,1.414], [1.414,0,0], [0,1.414,0], [-1.414,0,0], [0,-1.414,0], [0,0,-1.414] ]
  poly

cube = () ->
  poly = new polyhedron()
  poly.name = "C"
  poly.face = [ [3,0,1,2], [3,4,5,0], [0,5,6,1], [1,6,7,2], [2,7,4,3], [5,4,7,6] ]
  poly.xyz  = [ [0.707,0.707,0.707], [-0.707,0.707,0.707], [-0.707,-0.707,0.707], [0.707,-0.707,0.707],
                [0.707,-0.707,-0.707], [0.707,0.707,-0.707], [-0.707,0.707,-0.707], [-0.707,-0.707,-0.707] ]
  poly

icosahedron = () ->
  poly = new polyhedron()
  poly.name = "I"
  poly.face = [ [0,1,2], [0,2,3], [0,3,4], [0,4,5],
    [0,5,1], [1,5,7], [1,7,6], [1,6,2],
    [2,6,8], [2,8,3], [3,8,9], [3,9,4],
    [4,9,10], [4,10,5], [5,10,7], [6,7,11],
    [6,11,8], [7,10,11], [8,11,9], [9,11,10] ]

  poly.xyz = [ [0,0,1.176], [1.051,0,0.526],
    [0.324,1.0,0.525], [-0.851,0.618,0.526],
    [-0.851,-0.618,0.526], [0.325,-1.0,0.526],
    [0.851,0.618,-0.526], [0.851,-0.618,-0.526],
    [-0.325,1.0,-0.526], [-1.051,0,-0.526],
    [-0.325,-1.0,-0.526], [0,0,-1.176] ]
  poly

dodecahedron = ->
   poly = new polyhedron()
   poly.name = "D"
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
      [-0.713644,0,-0.797878], [0,0,-1.07047] ]
   poly

prism = (n) ->
  theta = 2*PI/n # pie angle
  h = Math.sin(theta/2) # half-edge
  poly = new polyhedron()
  poly.name = "P#{n}"

  for i in [0..n-1] # vertex #'s 0...n-1 around one face
    poly.xyz.push [cos(i*theta), sin(i*theta),  h]
  for i in [0..n-1] # vertex #'s n...2n-1 around other
    poly.xyz.push [cos(i*theta), sin(i*theta), -h]

  poly.face.push [n-1..0] #top
  poly.face.push [n..2*n-1] #bottom
  for i in [0..n-1] #n square sides
    poly.face.push [i, (i+1)%n, (i+1)%n+n, i+n]

  poly.xyz = adjustXYZ(poly,1)
  poly

antiprism = (n) ->
  theta = 2*PI/n # pie angle
  h = sqrt(1-4/(4+2*cos(theta/2)-2*cos(theta)))
  r = sqrt(1-h*h)
  f = sqrt(h*h + pow(r*cos(theta/2),2) )
  # correction so edge midpoints (not vertices) on unit sphere
  r = r/f
  h = h/f
  poly = new polyhedron()
  poly.name = "A#{n}"

  for i in [0..n-1] # vertex #'s 0...n-1 around one face
    poly.xyz.push [r * cos(i*theta), r * sin(i*theta), h]
  for i in [0..n-1] # vertex #'s n...2n-1 around other
    poly.xyz.push [r * cos((i+0.5)*theta), r * sin((i+0.5)*theta), -h]

  poly.face.push [n-1..0]  #top
  poly.face.push [n..2*n-1] #bottom
  for i in [0..n-1] #2n triangular sides
    poly.face.push [i, (i+1)%n, i+n]
    poly.face.push [i, i+n, ((n+i-1)%n+n)]

  poly.xyz = adjustXYZ(poly,1)
  poly

pyramid = (n) ->
  theta = 2*PI/n # pie angle
  poly = new polyhedron()
  poly.name = "Y#{n}"

  for i in [0..n-1] # vertex #'s 0...n-1 around one face
    poly.xyz.push [cos(i*theta), sin(i*theta), 0.2]
  poly.xyz.push [0,0,-2] # apex

  poly.face.push [n-1..0] # base
  for i in [0..n-1] # n triangular sides
    poly.face.push [i, (i+1)%n, n]

  poly.xyz = canonicalXYZ(poly,3)
  poly



#===================================================================================================
# Flag Construct
#
# A Flag is an associative triple of a face index and two adjacent vertex indices,
# listed in clockwise order
#
# FaceID -> {Vi -> Vj}
#
class polyflag
  constructor: ->
    @flags= new Object() # flags[face][vertex] = next vertex of flag; symbolic triples
    @verts= new Object() # XYZ coordinates
    @xyzs = new Object() # [symbolic names] holds vertex index

  newV: (name, xyz) ->
    if @verts[name] is undefined
      @verts[name] = 0
      @xyzs[name] = xyz

  newFlag: (face, v1, v2) ->
    if @flags[face] is undefined
      @flags[face] = {}
    @flags[face][v1] = v2

  topoly: () ->
    poly = new polyhedron()

    ctr=0 # first number the vertices
    for i,v of @verts
      poly.xyz[ctr]=@xyzs[i] # store in array
      @verts[i] = ctr
      ctr++

    ctr=0
    for i,f of @flags
      poly.face[ctr] = [] # new face
      # grab _any_ vertex as starting point
      for j,v of f
        v0 = v
        break  # need just one
      # build face out of all the edge relations in the flag assoc array
      v = v0 # v moves around face
      poly.face[ctr].push @verts[v] #record index
      v = @flags[i][v] # goto next vertex
      while v isnt v0  # loop until back to start
        poly.face[ctr].push @verts[v]
        v = @flags[i][v]
      ctr++

    poly.name = "unknown polyhedron"
    poly


#===================================================================================================
# Polyhedron Operators
#===================================================================================================
#          for each vertex of new polyhedron:
#              call newV(Vname, xyz) with symbolic name and approx location
#          for each flag of new polyhedron:
#              call newFlag(Fname, Vname1, Vname2)  with symbolic names
#          call topoly() to assemble flags into polyhedron structure
#          canonicalize vertex locations
#          set name as appropriate

# Kis(N) ------------------------------------------------------------------------------------------
# only kis n-sided faces, but n==0 means kiss all.
kisN = (poly, n)->
  console.log "Taking kis of #{if n==0 then "" else n}-sided faces of #{poly.name}..."

  flag = new polyflag()
  for [i,p] in enumerate(poly.xyz)
    # each old vertex is a new vertex
    flag.newV "v#{i}", p

  centers = faceCenters(poly)      # new vertices in centers of n-sided face
  foundAny = false                 # alert if don't find any
  for [i,f] in enumerate(poly.face)
    v1 = "v"+f[-1..][0]
    for v in f
      v2 = "v" + v
      if f.length is n or n is 0
        foundAny = true
        flag.newV "f#{i}", centers[i]
        fname = i + v1
        flag.newFlag fname,      v1,      v2
        flag.newFlag fname,      v2, "f#{i}"
        flag.newFlag fname, "f#{i}",      v1
      else
        flag.newFlag i, v1, v2  # same old flag, if non-n
      v1=v2  # current becomes previous

  if not foundAny
    console.log "No #{n}-fold components were found."

  newpoly = flag.topoly()
  newpoly.name = "k" + (if n is 0 then "" else n) + poly.name
  newpoly.xyz = adjustXYZ(newpoly, 3)
  #newpoly.xyz = canonicalXYZ(newpoly, 3)  # this tends to make results look like shit
  newpoly


# Ambo ------------------------------------------------------------------------------------------
midName = (v1, v2) -> if v1<v2 then v1+"_"+v2 else v2+"_"+v1
ambo = (poly)->
  console.log "Taking ambo of #{poly.name}..."

  flag = new polyflag()

  for [i,f] in enumerate(poly.face)
    [v1, v2] = f[-2..-1]
    for v3 in f
      if v1 < v2
        flag.newV(midName(v1,v2), midpoint(poly.xyz[v1], poly.xyz[v2]))
      # two new flags
      flag.newFlag("f"+i,  midName(v1,v2), midName(v2,v3))
      flag.newFlag("v"+v2, midName(v2,v3), midName(v1,v2))

      # shift over one
      [v1, v2] = [v2, v3]

  newpoly = flag.topoly()
  newpoly.name = "a" + poly.name
  newpoly.xyz = adjustXYZ(newpoly, 2)
  newpoly

# Gyro ----------------------------------------------------------------------------------------------
gyro = (poly)->
  console.log "Taking gyro of #{poly.name}..."

  flag = new polyflag()

  for [i,v] in enumerate(poly.xyz)
    flag.newV "v"+i, unit(v)  # each old vertex is a new vertex

  centers = faceCenters(poly) # new vertices in center of each face
  for [i,f] in enumerate(poly.face)
    flag.newV "f"+i, unit(centers[i])

  for [i,f] in enumerate(poly.face)
    [v1, v2] = f[-2..-1]
    for [j,v] in enumerate(f)
      v3 = v
      flag.newV(v1+"~"+v2, oneThird(poly.xyz[v1],poly.xyz[v2]))  # new v in face
      fname = i + "f" + v1
      flag.newFlag(fname, "f"+i,      v1+"~"+v2) # five new flags
      flag.newFlag(fname, v1+"~"+v2,  v2+"~"+v1)
      flag.newFlag(fname, v2+"~"+v1,  "v"+v2)
      flag.newFlag(fname, "v"+v2,     v2+"~"+v3)
      flag.newFlag(fname, v2+"~"+v3,  "f"+i)
      [v1, v2] = [v2, v3]                       # shift over one

  newpoly = flag.topoly()
  newpoly.name = "g" + poly.name
  newpoly.xyz = adjustXYZ(newpoly, 3)
  newpoly

# Propellor ------------------------------------------------------------------------------------------
propellor = (poly) ->
  console.log "Taking propellor of #{poly.name}..."

  flag = new polyflag()

  for [i,v] in enumerate(poly.xyz)
    flag.newV("v"+i, unit(v))  # each old vertex is a new vertex

  for [i,f] in enumerate(poly.face)
    [v1, v2] = f[-2..-1]
    for v in f
      v3 = "#{v}"
      flag.newV(v1+"~"+v2, oneThird(poly.xyz[v1], poly.xyz[v2]))  # new v in face, 1/3rd along edge
      fname = "#{i}f#{v2}"
      flag.newFlag("v#{i}", v1+"~"+v2,  v2+"~"+v3) # five new flags
      flag.newFlag(fname,   v1+"~"+v2,  v2+"~"+v1)
      flag.newFlag(fname,   v2+"~"+v1,     "v"+v2)
      flag.newFlag(fname,      "v"+v2,  v2+"~"+v3)
      flag.newFlag(fname,   v2+"~"+v3,  v1+"~"+v2)
      [v1, v2] = [v2, v3]                       # shift over one

  newpoly = flag.topoly()
  newpoly.name = "p" + poly.name
  newpoly.xyz  = adjustXYZ(newpoly, 3)
  newpoly

# Reflection ------------------------------------------------------------------------------------------
# compute reflection through origin
reflect = (poly) ->
  console.log "Taking reflection of #{poly.name}..."
  for i in [0..poly.xyz.length-1]
       poly.xyz[i] = mult(-1, poly.xyz[i])         # reflect each point
  for i in [0..poly.face.length-1]
       poly.face[i] = poly.face[i].reverse()       # repair clockwise-ness
  poly.name = "r" + poly.name
  poly.xyz = adjustXYZ(poly, 1)                    # build dual
  poly

# Dual ------------------------------------------------------------------------------------------------
# the dual function computes the dual's topology, matching F and V indices, and invokes simple
# canonicalization for the determination of xyz coordinates
dual = (poly) ->
  console.log "Taking dual of #{poly.name}..."

  flag = new polyflag()

  face = [] # make table of face as fn of edge
  for i in [0..poly.xyz.length-1]
    face[i] = {} # create empty associative table

  for i in [0..poly.face.length-1]
    v1 = poly.face[i][poly.face[i].length-1] #previous vertex
    for j in [0..poly.face[i].length-1]
      v2 = poly.face[i][j] # this vertex
      face[v1]["v#{v2}"] = "#{i}" # fill it. 2nd index is associative
      v1=v2 # current becomes previous

  for i in [0..poly.face.length-1]
    flag.newV("#{i}",[])

  for i in [0..poly.face.length-1]
    v1 = poly.face[i][poly.face[i].length-1]
    for j in [0..poly.face[i].length-1]
      v2 = poly.face[i][j]
      flag.newFlag(v1, face[v2]["v#{v1}"], "#{i}")
      v1=v2

  dpoly = flag.topoly() # build topological dual from flags

  # match F index ordering to V index ordering on dual
  sortF = []
  for f in dpoly.face
    k = intersect(poly.face[f[0]],poly.face[f[1]],poly.face[f[2]])
    sortF[k] = f
  dpoly.face = sortF

  # compute coordinates as dual to those of original poly
  dpoly.xyz = reciprocalN(poly)

  if poly.name[0] isnt "d"
    dpoly.name = "d"+poly.name
  else
    dpoly.name = poly.name[1..]

  dpoly


#===================================================================================================
# Canonicalization Algorithms
#===================================================================================================

# SLOW Canonicalization Algorithm ----------------------------------------------------------------
#
# This algorithm has some convergence problems, what really needs to be done is to
# sum the three forcing factors together as a conherent force and to use a half-decent
# integrator to make sure that it converges well as opposed to the current hack of
# ad-hoc stability multipliers.  Ideally one would implement a conjugate gradient
# descent or similar pretty thing.
#

# adjusts vertices on edges such that each edge is tangent to an origin sphere
tangentify = (xyzs, edges) ->
  STABILITY_FACTOR = 0.1
  newVs = copyVecArray(xyzs) #copy vertices
  for e in edges
    t = tangentPoint( newVs[e[0]], newVs[e[1]] ) #the point closest to origin
    c = mult(STABILITY_FACTOR*1/2*(1-sqrt(dot(t,t))), t) # adjustment from sphere
    newVs[e[0]] = add(newVs[e[0]], c)
    newVs[e[1]] = add(newVs[e[1]], c)
  newVs

# recenters entire polyhedron such that center of mass is at origin
recenter = (xyzs, edges) ->
  edgecenters = (tangentPoint(xyzs[a], xyzs[b]) for [a,b] in edges) #centers of edges
  polycenter = [0,0,0]
  for v in edgecenters # sum centers to find center of gravity
    polycenter = add(polycenter, v)
  polycenter = mult(1/edges.length, polycenter)
  _.map(xyzs, (x)->sub(x, polycenter) ) # subtract off any deviation from center

# adjusts vertices in each face to improve its planarity
planarize = (xyzs, faces) ->
  STABILITY_FACTOR = 0.1
  newVs = copyVecArray(xyzs) # copy vertices
  for f in faces
    coords = (xyzs[v] for v in f)
    n = normal(coords) # find avg of normals for each vertex triplet
    c = centroid(coords) # find planar centroid
    if dot(n,centroid) < 0 # correct sign if needed
      n = mult(-1.0,n)
    for v in f  # project (vertex - centroid) onto normal, subtract off this component
      newVs[v] = add(newVs[v], mult(dot(mult(STABILITY_FACTOR, n), sub(c, xyzs[v])), n) )
  newVs

# combines above three constraint adjustments in iterative cycle
canonicalize = (poly, Niter) ->
  console.log "Canonicalizing #{poly.name}..."
  faces = poly.face
  edges = poly.getEdges()
  newVs = poly.xyz
  maxChange=1.0 # convergence tracker
  for i in [0..Niter]
    oldVs = copyVecArray(newVs) #copy vertices
    newVs = tangentify(newVs, edges)
    newVs = recenter(newVs, edges)
    newVs = planarize(newVs, faces)
    maxChange=_.max(_.map(_.zip(newVs,oldVs), ([x,y])->mag(sub(x,y)) ))
    if maxChange < 1e-8
      break
  console.log "[canonicalization done, last |deltaV|="+maxChange+"]"
  newVs


# Hacky Canonicalization Algorithm--------------------------------------------------------
# Using center of gravity of vertices for each face to planarize faces

# get the spherical reciprocals of face centers
reciprocalC = (poly) ->
  centers = faceCenters(poly)
  for c in centers
    c = mult(1/dot(c,c), c)
  centers

canonicalXYZ = (poly, nIterations) ->
  dpoly = dual(poly)
  console.log "Pseudo-canonicalizing #{poly.name}..."

  # iteratively reciprocate face normals
  for count in [0..nIterations-1]
    dpoly.xyz = reciprocalN(poly)
    poly.xyz  = reciprocalN(dpoly)

  poly.xyz

# make array of vertices reciprocal to given planes
reciprocalN = (poly) ->
  ans = []
  for f in poly.face #for each face
    centroid    = [0,0,0] # running sum of vertex coords
    normal      = [0,0,0] # running sum of normal vectors
    avgEdgeDist =    0.0  # running sum for avg edge distance

    [v1, v2] = f[-2..-1]
    for v3 in f
      centroid = add centroid, poly.xyz[v3]
      normal   = add normal, orthogonal(poly.xyz[v1], poly.xyz[v2], poly.xyz[v3])
      avgEdgeDist += edgeDist(poly.xyz[v1], poly.xyz[v2])
      [v1, v2] = [v2, v3] # shift over one

    centroid    = mult(1/f.length, centroid)
    normal      = unit(normal)
    avgEdgeDist = avgEdgeDist / f.length
    tmp         = reciprocal mult dot(centroid, normal), normal # based on face
    ans.push mult((1 + avgEdgeDist) / 2, tmp) # edge correction

  ans

# quick planarization
adjustXYZ = (poly, nIterations) ->
  dpoly = dual(poly) # v's of dual are in order of arg's f's
  console.log "Planarizing #{poly.name}..."

  for count in [0..nIterations-1]
    # reciprocate face centers
    dpoly.xyz = reciprocalC(poly)
    poly.xyz  = reciprocalC(dpoly)

  poly.xyz

# get array of face centers
faceCenters = (poly) ->
  centers = []
  for i in [0..poly.face.length-1]
    centers[i] = [0,0,0]
    for j in [0..poly.face[i].length-1] #avg vertex coords
      centers[i] = add(centers[i], poly.xyz[poly.face[i][j]]) # add
    centers[i] = mult(1.0/poly.face[i].length, centers[i]) # div by n

  centers

# calculate centroid of array of vertices
centroid = (xyzs) ->
    centroidV = [0,0,0] # running sum of vertex coords
    for v in xyzs
      centroidV = add(centroidV, v)
    mult(1 / xyzs.length, centroidV )

# calculate average normal vector for array of vertices
normal = (xyzs) ->
    normalV = [0,0,0] # running sum of normal vectors
    [v1,v2] = xyzs[-2..-1]
    for v3 in xyzs
      normalV = add(normalV, orthogonal(v1, v2, v3))
      [v1,v2] = [v2,v3] # shift over one
    unit(normalV)

# calculates area planar face by summing over subtriangle areas
#  _Assumes_ Convexity!
convexarea = (xyzs) ->
    area = 0.0
    [v1,v2] = xyzs[0..1]
    for v3 in xyzs[2..]
      #area of sub-triangle
      area += mag( cross(sub(v2,v1), sub(v3,v1)) )
      v2 = v3 # shift over one
    area


#===================================================================================================
# Parser Routines
#===================================================================================================

specreplacements = [[/P4$/g, "C"],  # P4 --> C   (C is prism)
  [/e/g, "aa"],   # e --> aa   (abbr. for explode)
  [/b/g, "ta"],   # b --> ta   (abbr. for bevel)
  [/o/g, "jj"],   # o --> jj   (abbr. for ortho)
  [/m/g, "kj"],   # m --> kj   (abbr. for meta)
  [/t(\d*)/g, "dk$1d"],  # t(n) --> dk(n)d  (dual operations)
  [/j/g, "dad"],  # j --> dad  (dual operations)
  [/s/g, "dgd"],  # s --> dgd  (dual operations)
  [/dd/g, ""],    # dd --> null  (order 2)
# [/ad/g, "a"],   # ad --> a   (a_ = ad_)
# [/gd/g, "g"],   # gd --> g   (g_ = gd_)
  [/aO/g, "aC"],  # aO --> aC  (for uniqueness)
  [/aI/g, "aD"],  # aI --> aD  (for uniqueness)
  [/gO/g, "gC"],  # gO --> gC  (for uniqueness)
  [/gI/g, "gD"]]  # gI --> gD  (for uniqueness)

getOps = (notation) ->
  expanded = notation
  for [orig,equiv] in specreplacements
    expanded = expanded.replace(orig,equiv)
  console.log "#{notation} executed as #{expanded}"
  expanded

# create polyhedron from notation
generatePoly = (notation) ->
  poly = new polyhedron()
  n=0

  ops = getOps(notation)
  if ops.search(/([0-9]+)$/) != -1
    n = 1 * RegExp.lastParen
    ops = ops.slice(0, -RegExp.lastParen.length)

  switch ops.slice(-1)
    when "T" then poly = tetrahedron()
    when "O" then poly = octahedron()
    when "C" then poly = cube()
    when "I" then poly = icosahedron()
    when "D" then poly = dodecahedron()
    when "P" then poly = prism(n)
    when "A" then poly = antiprism(n)
    when "Y" then poly = pyramid(n)
    else return

  while ops != ""
    n=0
    if ops.search(/([0-9]+)$/) != -1
      n = 1 * RegExp.lastParen
      ops = ops.slice(0, -RegExp.lastParen.length)
    switch ops.slice(-1)
      # geometrical operators
      when "d" then poly     = dual(poly)
      when "k" then poly     = kisN(poly, n)
      when "a" then poly     = ambo(poly)
      when "g" then poly     = gyro(poly)
      when "p" then poly     = propellor(poly)
      when "r" then poly     = reflect(poly)
      # canonicalization operators
      when "." then poly.xyz = canonicalXYZ(poly, if n is 0 then 5 else n*5)
      when "!" then poly.xyz = canonicalize(poly, if n is 0 then 5 else n*80)

    ops = ops.slice(0,-1);  # remove last character

  # Recenter polyhedra at origin (rarely needed)
  poly.xyz = recenter(poly.xyz, poly.getEdges())

  # Color the faces of the polyhedra for display
  poly = paintPolyhedron(poly)

  # return the poly object
  poly


# parses URL string for polyhedron recipe, for bookmarking
# should use # href format instead
parseurl = () ->
  urlParams = {}
  a = /\+/g  # Regex for replacing addition symbol with a space
  r = /([^&=]+)=?([^&]*)/g
  d = (s) -> decodeURIComponent(s.replace(a, " "))
  q = window.location.search.substring(1)

  while e=r.exec(q)
    urlParams[d(e[1])] = d(e[2])
  urlParams


#===================================================================================================
# Testing Functions
#===================================================================================================

#report on face topology
topolog = (poly) ->
  str=""
  for f in poly.face
    for v in f
      str+="#{v}->"
    str+="\n"
  console.log str


testrig = ->
  tests=["T","O","C","I","D","P3","P4","A4","A5","Y3","Y4"]
  console.log "===== Test Basic Ops ====="
  console.log "--- primitives ----------------------------------------------------------------------------- "
  for t in tests
    console.log generatePoly t
  console.log "--- kis ----------------------------------------------------------------------------- "
  for t in tests
    console.log generatePoly "k"+t
  console.log "--- ambo ----------------------------------------------------------------------------- "
  for t in tests
    console.log generatePoly "a"+t
  console.log "--- gyro ----------------------------------------------------------------------------- "
  for t in tests
    console.log generatePoly "g"+t
  console.log "--- propellor ----------------------------------------------------------------------------- "
  for t in tests
    console.log generatePoly "p"+t
  console.log "--- dual ----------------------------------------------------------------------------- "
  for t in tests
    console.log generatePoly "d"+t
  console.log "===== Done Testing Basic Ops ====="
#test basic stuff
#testrig()


#===================================================================================================
# GLOBALS
#===================================================================================================
ctx={} # for global access to canvas context
CANVAS_WIDTH  = 600 #canvas dims
CANVAS_HEIGHT = 300 #canvas dims
globPolys={} # constructed polyhedras

globtheta = 0 # rotation and projective mapping parameters
globphi   = 0
perspective_scale = 500
persp_z_max = 5
persp_z_min = 0
persp_ratio = 0.8
_2d_x_offset = CANVAS_WIDTH/2 #300
_2d_y_offset = CANVAS_HEIGHT/2 #140

globtime = new Date() # for animation


BG_CLEAR = true # clear background or colored?
BG_COLOR = "rgba(255,255,255,1.0)" # background color

COLOR_METHOD = "area"

ctx_linewidth = 0.5 # for outline of faces

# Mouse Event Variables
MOUSEDOWN=false
LastMouseX=0
LastMouseY=0


# Polyhedra Coloring Functions
#===================================================================================================

def_palette  = ["#ff3333","#33ff33","#3333ff","#ffff33","#ff33ff","#33ffff","#dddddd","#555555","#dd0000","#00dd00","#0000dd"]
rwb_palette  = ["#ff8888","#dddddd","#777777","#aa3333","#ff0000","#ffffff","#aaaaaa"]
rwbg_palette = ["#ff8888","#ffeeee","#88ff88","#dd7777","#ff2222","#22ff22","#ee4422","#aaaaaa"]

# converts #xxxxxx / #xxx format into list of [r,g,b] floats
hextofloats = (hexstr)->
  if hexstr[0] is "#"
    hexstr = hexstr[1..]
  if hexstr.length is 3
    rgb = hexstr.split('').map(       (c)->parseInt(c+c, 16)/255 )
  else
    rgb = hexstr.match(/.{2}/g).map(  (c)->parseInt(c, 16)/255 )
  rgb

PALETTE = rwb_palette
palette = (n) -> if n < PALETTE.length then hextofloats(PALETTE[n]) else hextofloats(PALETTE[PALETTE.length-1])

#memoized color assignment to faces of similar areas
colorassign = (ar, colormemory) ->
  hash = round(100*ar)
  if hash of colormemory
    return colormemory[hash]
  else
    fclr = palette _.toArray(colormemory).length
    colormemory[hash] = fclr
    return fclr

paintPolyhedron = (poly) ->
  # Color the faces of the polyhedra for display
  poly.face_colors = []
  colormemory={}
  for f in poly.face
    if COLOR_METHOD is "area"
      # color by face area (quick proxy for different kinds of faces) convexarea
      face_verts = (poly.xyz[v] for v in f)
      clr = colorassign(convexarea(face_verts), colormemory)
    else
      # color by face-sidedness
      clr = palette f.length-3

    poly.face_colors.push clr

  poly


#===================================================================================================
# Drawing Functions
#===================================================================================================

#init canvas element -------------------------------------------------------------------------------
init = ->
  canvas = $('#poly')
  canvas.width(CANVAS_WIDTH)
  canvas.height(CANVAS_HEIGHT)

  ctx = canvas[0].getContext("2d")
  ctx.lineWidth = ctx_linewidth

  if BG_CLEAR
    ctx.clearRect 0, 0, CANVAS_WIDTH, CANVAS_HEIGHT
  else
    ctx.clearRect 0, 0, CANVAS_WIDTH, CANVAS_HEIGHT
    ctx.fillStyle = BG_COLOR
    ctx.fillRect 0, 0, CANVAS_WIDTH, CANVAS_HEIGHT

# clear canvas -----------------------------------------------------------------------------------
clear = ->
  if BG_CLEAR
    ctx.clearRect 0, 0, CANVAS_WIDTH, CANVAS_HEIGHT
  else
    ctx.clearRect 0, 0, CANVAS_WIDTH, CANVAS_HEIGHT
    ctx.fillStyle = BG_COLOR
    ctx.fillRect 0, 0, CANVAS_WIDTH, CANVAS_HEIGHT

# z sorts faces of poly -------------------------------------------------------------------------
sortfaces = (poly) ->
  centroids = faceCenters(poly)
  zsortIndex = _.zip(centroids,[0..poly.face.length-1])
    .sort((a,b) -> a[0][2]-b[0][2]) # js sort is lexicographic even for numbers!
    .map((x)->x[1])

  # sort all face-associated properties
  poly.face = (poly.face[idx] for idx in zsortIndex)
  poly.face_colors = (poly.face_colors[idx] for idx in zsortIndex)


#===================================================================================================
# main drawing routine for polyhedra
drawpoly = (poly,tvec,rot) ->
  tvec ||= [3,3,3]
  rot  ||= [1,0,1]

  # rotate poly in 3d
  oldxyz = _.map(poly.xyz, (x)->x)
  poly.xyz = _.map(poly.xyz, (x)->mv3(rotm(rot[0],rot[1],rot[2]),x))
  # z sort faces
  sortfaces(poly)

  for [fno,face] in enumerate(poly.face)
    ctx.beginPath()
    # move to first vertex of face
    v0 = face[face.length-1]
    [x,y] = perspT(add(tvec,poly.xyz[v0]), persp_z_max,persp_z_min,persp_ratio,perspective_scale)
    ctx.moveTo(x+_2d_x_offset, y+_2d_y_offset)
    # loop around face, defining polygon
    for v in face
      [x,y] = perspT(add(tvec,poly.xyz[v]),persp_z_max,persp_z_min,persp_ratio,perspective_scale)
      ctx.lineTo(x+_2d_x_offset, y+_2d_y_offset)

    # use pre-computed colors
    clr = poly.face_colors[fno]

    # shade based on simple cosine illumination factor
    face_verts = (poly.xyz[v] for v in face)
    illum = dot(normal(face_verts), unit([-1,1,0]))
    clr   = mult((illum/2.0+.5)*0.7+0.3,clr)

    ctx.fillStyle = "rgba(#{round(clr[0]*255)}, #{round(clr[1]*255)}, #{round(clr[2]*255)}, #{1.0})"
    ctx.fill()
    ctx.stroke()

  # reset coords, for setting absolute rotation, as poly is passed by ref
  poly.xyz = oldxyz

#===================================================================================================
# Initialization and Basic UI
#===================================================================================================

$( -> #wait for page to load

  init() #init canvas

  defspecs = ["dakD","oopD","ajI",".akY5",".ooC","bT"]   # random grabbag of polyhedra

  urlParams = parseurl() #see if recipe is spec'd in URL
  if "recipe" of urlParams
    specs=[urlParams["recipe"]]
    $("#spec").val(specs)
  else
    specs=[randomchoice(defspecs)]
    $("#spec").val(specs)

  # construct the polyhedra from spec
  globPolys = _.map(specs, (x)->generatePoly(x))

  # draw it
  #animateShape()
  drawShape()

  # when spec changes in input, parse and draw new polyhedra
  $("#spec").change((e) =>
    specs = $("#spec").val().split(/\s+/g)
    globPolys = _.map(specs, (x)->generatePoly(x) )
    #animateShape()
    drawShape()
  )

  # basic manipulation: rotation and scaling of geometry ---------------------------------------------

  # mousewheel changes scale of drawing
  $("#poly").mousewheel( (e,delta, deltaX, deltaY)->
    event.preventDefault()
    perspective_scale*=(10+delta)/10
    drawShape()
  )

  # implement standard trackball routines ------------------------------
  $("#poly").mousedown( (e)->
    event.preventDefault()
    MOUSEDOWN=true
    LastMouseX=e.clientX-$(this).offset().left
    LastMouseY=e.clientY-$(this).offset().top
  )
  $("#poly").mouseup( (e)->
    event.preventDefault()
    MOUSEDOWN=false
  )
  $("#poly").mouseleave( (e)->
    event.preventDefault()
    MOUSEDOWN=false
  )
  $("#poly").mousemove( (e)->
    event.preventDefault()
    if MOUSEDOWN
      #console.log e.clientX-$(this).offset().left-LastMouseX, e.clientY-$(this).offset().top-LastMouseY
      globtheta += -( e.clientX-$(this).offset().left-LastMouseX)*(Math.PI/180)
      globphi   += -( e.clientY-$(this).offset().top-LastMouseY)*(Math.PI/180)

      LastMouseX=e.clientX-$(this).offset().left
      LastMouseY=e.clientY-$(this).offset().top
      drawShape()
  )

)

# loop for animation
animateShape = ->
  clear()
  globtheta=(2*Math.PI)/180.0*globtime.getSeconds()*0.1
  for [i,p] in enumerate(globPolys)
    drawpoly(p,[0+3*i,0,3],[0,globtheta,globphi])
  setTimeout(animateShape, 100)

# just draw polys once
drawShape = ->
  clear()
  for [i,p] in enumerate(globPolys)
    drawpoly(p,[0+3*i,0,3],[0,globtheta,globphi])


#===================================================================================================

