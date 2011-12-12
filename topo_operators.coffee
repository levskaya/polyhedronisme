# Polyhédronisme
#===================================================================================================
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
# Polyhedron Flagset Construct
#
# A Flag is an associative triple of a face index and two adjacent vertex indices,
# listed in geometric clockwise order (staring into the normal)
#
# Face_i -> V_i -> V_j
#
# They are a useful abstraction for defining topological transformations of the polyhedral mesh, as
# one can refer to vertices and faces that don't yet exist or haven't been traversed yet in the
# transformation code.
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
      faceCTR=0
      while v isnt v0 # loop until back to start
        poly.face[ctr].push @verts[v]
        v = @flags[i][v]
        faceCTR++
        if faceCTR>1000 # necessary to prevent browser hangs on badly formed flagsets!
          console.log "Bad flag spec, have a neverending face:", i, @flags[i]
          break
      ctr++

    poly.name = "unknown polyhedron"
    poly


#===================================================================================================
# Polyhedron Operators
#===================================================================================================
#          for each vertex of new polyhedron:
#              call newV(Vname, xyz) with a symbolic name and coordinates
#          for each flag of new polyhedron:
#              call newFlag(Fname, Vname1, Vname2) with a symbolic name for the new face
#              and the symbolic name for two vertices forming an oriented edge
#          ORIENTATION -must- be dealt with properly to make a manifold (correct) mesh.
#          Specifically, no edge v1->v2 can ever be crossed in the -same direction- by
#          two different faces
#
#          call topoly() to assemble flags into polyhedron structure by following the orbits
#          of the vertex mapping stored in the flagset for each new face
#
#          set name as appropriate

# Kis(N)
# ------------------------------------------------------------------------------------------
# Kis (abbreviated from triakis) transforms an N-sided face into an N-pyramid rooted at the
# same base vertices.
# only kis n-sided faces, but n==0 means kiss all.
#
kisN = (poly, n, apexdist)->
  n or= 0
  apexdist or= 0.1
  console.log "Taking kis of #{if n==0 then "" else n}-sided faces of #{poly.name}..."

  flag = new polyflag()
  for p,i in poly.xyz
    # each old vertex is a new vertex
    flag.newV "v#{i}", p

  normals = poly.normals()
  centers = poly.centers()
  foundAny = false                 # alert if don't find any
  for f,i in poly.face
    v1 = "v"+f[f.length-1]
    for v in f
      v2 = "v"+v
      if f.length is n or n is 0
        foundAny = true
        apex = "apex#{i}"
        fname = "#{i}#{v1}"
        flag.newV apex, add(centers[i],mult(apexdist,normals[i])) # new vertices in centers of n-sided face
        flag.newFlag fname,   v1,   v2 # the old edge of original face
        flag.newFlag fname,   v2, apex # up to apex of pyramid
        flag.newFlag fname, apex,   v1 # and back down again
      else
        flag.newFlag "#{i}", v1, v2  # same old flag, if non-n
      v1=v2  # current becomes previous

  if not foundAny
    console.log "No #{n}-fold components were found."

  newpoly = flag.topoly()
  newpoly.name = "k" + (if n is 0 then "" else n) + poly.name
  #newpoly.xyz = adjustXYZ(newpoly, 3)
  #newpoly.xyz = canonicalXYZ(newpoly, 3)  # this tends to make results look like shit
  newpoly


# Ambo
# ------------------------------------------------------------------------------------------
# The best way to think of the ambo operator is as a topological "tween" between a polyhedron
# and its dual polyhedron.  Thus the ambo of a dual polyhedron is the same as the ambo of the
# original. Also called "Rectify".
#
ambo = (poly)->
  console.log "Taking ambo of #{poly.name}..."

  # helper func to insure unique names of midpoints
  midName = (v1, v2) -> if v1<v2 then v1+"_"+v2 else v2+"_"+v1

  flag = new polyflag()

  for f,i in poly.face
    [v1, v2] = f[-2..-1]
    for v3 in f
      if v1 < v2 # vertices are the midpoints of all edges of original poly
        flag.newV(midName(v1,v2), midpoint(poly.xyz[v1], poly.xyz[v2]))
      # two new flags
      flag.newFlag("orig"+i,  midName(v1,v2), midName(v2,v3))
      flag.newFlag("dual"+v2, midName(v2,v3), midName(v1,v2))
      # shift over one
      [v1, v2] = [v2, v3]

  newpoly = flag.topoly()
  newpoly.name = "a" + poly.name
  #newpoly.xyz = adjustXYZ(newpoly, 2)
  newpoly


# Gyro
# ----------------------------------------------------------------------------------------------
# This is the dual operator to "snub", i.e dual*Gyro = Snub.  It is a bit easier to implement
# this way.
#
# Snub creates at each vertex a new face, expands and twists it, and adds two new triangles to
# replace each edge.
#
gyro = (poly)->
  console.log "Taking gyro of #{poly.name}..."

  flag = new polyflag()

  for v,i in poly.xyz
    flag.newV "v"+i, unit(v)  # each old vertex is a new vertex

  centers = poly.centers() # new vertices in center of each face
  for f,i in poly.face
    flag.newV "center"+i, unit(centers[i])

  for f,i in poly.face
    [v1, v2] = f[-2..-1]
    for v,j in f
      v3 = v
      flag.newV(v1+"~"+v2, oneThird(poly.xyz[v1],poly.xyz[v2]))  # new v in face
      fname = i+"f"+v1
      flag.newFlag(fname, "center"+i,      v1+"~"+v2) # five new flags
      flag.newFlag(fname, v1+"~"+v2,  v2+"~"+v1)
      flag.newFlag(fname, v2+"~"+v1,  "v"+v2)
      flag.newFlag(fname, "v"+v2,     v2+"~"+v3)
      flag.newFlag(fname, v2+"~"+v3,  "center"+i)
      [v1, v2] = [v2, v3]                       # shift over one

  newpoly = flag.topoly()
  newpoly.name = "g" + poly.name
  #newpoly.xyz = adjustXYZ(newpoly, 3)
  newpoly


# Propellor
# ------------------------------------------------------------------------------------------
# builds a new 'skew face' by making new points along edges, 1/3rd the distance from v1->v2,
# then connecting these into a new inset face.  This breaks rotational symmetry about the
# faces, whirling them into gyres
#
propellor = (poly) ->
  console.log "Taking propellor of #{poly.name}..."

  flag = new polyflag()

  for v,i in poly.xyz
    flag.newV("v"+i, unit(v))  # each old vertex is a new vertex

  for f,i in poly.face
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
  #newpoly.xyz  = adjustXYZ(newpoly, 3)
  newpoly


# Reflection
# ------------------------------------------------------------------------------------------
# geometric reflection through origin
reflect = (poly) ->
  console.log "Taking reflection of #{poly.name}..."
  for i in [0..poly.xyz.length-1]
       poly.xyz[i] = mult(-1, poly.xyz[i])         # reflect each point through origin
  for i in [0..poly.face.length-1]
       poly.face[i] = poly.face[i].reverse()       # repair clockwise-ness of faces!
  poly.name = "r" + poly.name
  poly


# Dual
# ------------------------------------------------------------------------------------------------
# The dual of a polyhedron is another mesh wherein:
# - every face in the original becomes a vertex in the dual
# - every vertex in the original becomes a face in the dual
#
# So N_faces, N_vertices = N_dualfaces, N_dualvertices
#
# The new vertex coordinates are convenient to set to the original face centroids.
#
dual = (poly) ->
  console.log "Taking dual of #{poly.name}..."

  flag = new polyflag()

  face = [] # make table of face as fn of edge
  for i in [0..poly.xyz.length-1]
    face[i] = {} # create empty associative table

  for f,i in poly.face
    v1 = f[f.length-1] #previous vertex
    for v2 in f
      # THIS ASSUMES that no 2 faces that share an edge share it in the same orientation!
      # which of course never happens for proper manifold meshes, so get your meshes right.
      face[v1]["v#{v2}"] = "#{i}"
      v1=v2 # current becomes previous

  centers = poly.centers()
  for i in [0..poly.face.length-1]
    flag.newV("#{i}",centers[i])

  for f,i in poly.face
    v1 = f[f.length-1] #previous vertex
    for v2 in f
      flag.newFlag(v1, face[v2]["v#{v1}"], "#{i}")
      v1=v2 # current becomes previous

  dpoly = flag.topoly() # build topological dual from flags

  # match F index ordering to V index ordering on dual
  sortF = []
  for f in dpoly.face
    k = intersect(poly.face[f[0]],poly.face[f[1]],poly.face[f[2]])
    sortF[k] = f
  dpoly.face = sortF

  if poly.name[0] isnt "d"
    dpoly.name = "d"+poly.name
  else
    dpoly.name = poly.name[1..]

  dpoly

# insetN
# ------------------------------------------------------------------------------------------
insetN = (poly, n, inset_dist, popout_dist)->
  n or= 0
  inset_dist  or= 0.5
  popout_dist or= -0.2

  console.log "Taking inset of #{if n==0 then "" else n}-sided faces of #{poly.name}..."

  flag = new polyflag()
  for p,i in poly.xyz
    # each old vertex is a new vertex
    flag.newV "v#{i}", p

  normals = poly.normals()
  centers = poly.centers()
  for f,i in poly.face #new inset vertex for every vert in face
    if f.length is n or n is 0
      for v in f
        flag.newV "f"+i+"v"+v, add(tween(poly.xyz[v],centers[i],inset_dist),mult(popout_dist,normals[i]))

  foundAny = false                 # alert if don't find any
  for f,i in poly.face
    v1 = "v"+f[f.length-1]
    for v in f
      v2 = "v"+v
      if f.length is n or n is 0
        foundAny = true
        fname = i + v1
        flag.newFlag fname,      v1,       v2
        flag.newFlag fname,      v2,       "f"+i+v2
        flag.newFlag fname, "f"+i+v2,  "f"+i+v1
        flag.newFlag fname, "f"+i+v1,  v1
        #new inset, extruded face
        flag.newFlag "ex"+i, "f"+i+v1,  "f"+i+v2
      else
        flag.newFlag i, v1, v2  # same old flag, if non-n
      v1=v2  # current becomes previous

  if not foundAny
    console.log "No #{n}-fold components were found."

  newpoly = flag.topoly()
  newpoly.name = "n" + (if n is 0 then "" else n) + poly.name
  #newpoly.xyz = adjustXYZ(newpoly, 3)
  #newpoly.xyz = canonicalXYZ(newpoly, 3)  # this tends to make results look like shit
  newpoly


# ExtrudeN
# ------------------------------------------------------------------------------------------
extrudeN = (poly, n)->
  n or= 0
  console.log "Taking extrusion of #{if n==0 then "" else n}-sided faces of #{poly.name}..."

  flag = new polyflag()
  for p,i in poly.xyz
    # each old vertex is a new vertex
    flag.newV "v#{i}", p

  normals = poly.normals()
  centers = poly.centers()
  for f,i in poly.face #new inset vertex for every vert in face
    if f.length is n or n is 0
      for v in f
        #flag.newV "f"+i+"v"+v, add(midpoint(poly.xyz[v],centers[i]),mult(-0.2,normals[i]))
        flag.newV "f"+i+"v"+v, add(poly.xyz[v], mult(0.3,normals[i]))

  foundAny = false                 # alert if don't find any
  for f,i in poly.face
    v1 = "v"+f[f.length-1]
    for v in f
      v2 = "v"+v
      if f.length is n or n is 0
        foundAny = true
        #fname = i+v1
        flag.newFlag i+v1,       v1,       v2
        flag.newFlag i+v1,       v2, "f"+i+v2
        flag.newFlag i+v1, "f"+i+v2, "f"+i+v1
        flag.newFlag i+v1, "f"+i+v1,       v1
        #new inset, extruded face
        flag.newFlag "ex"+i, "f"+i+v1,  "f"+i+v2
      else
        flag.newFlag i, v1, v2  # same old flag, if non-n
      v1=v2  # current becomes previous

  if not foundAny
    console.log "No #{n}-fold components were found."

  newpoly = flag.topoly()
  newpoly.name = "x" + (if n is 0 then "" else n) + poly.name
  #console.log newpoly
  #newpoly.xyz = adjustXYZ(newpoly, 3)
  #newpoly.xyz = canonicalXYZ(newpoly, 3)  # this tends to make results look like shit
  newpoly


# StellaN
# ------------------------------------------------------------------------------------------
stellaN = (poly)->
  console.log "Taking stella of #{poly.name}..."

  centers = poly.centers()  # calculate face centers

  flag = new polyflag()
  for p,i in poly.xyz
    flag.newV "v#{i}", p      # each old vertex is a new vertex

  # iterate over triplets of faces v1,v2,v3
  for f,i in poly.face
    v1 = "v"+f[f.length-2]
    v2 = "v"+f[f.length-1]
    vert1 = poly.xyz[f[f.length-2]]
    vert2 = poly.xyz[f[f.length-1]]
    for v in f
      v3 = "v"+v
      vert3 = poly.xyz[v]
      v12=v1+"~"+v2 # names for "oriented" midpoints
      v21=v2+"~"+v1
      v23=v2+"~"+v3

      # on each Nface, N new points inset from edge midpoints towards center = "stellated" points
      flag.newV v12, midpoint( midpoint(vert1,vert2), centers[i] )

      # inset Nface made of new, stellated points
      flag.newFlag "in#{i}",      v12,       v23

      # new tri face constituting the remainder of the stellated Nface
      flag.newFlag "f#{i}#{v2}",      v23,      v12
      flag.newFlag "f#{i}#{v2}",       v12,      v2
      flag.newFlag "f#{i}#{v2}",      v2,      v23

      # one of the two new triangles replacing old edge between v1->v2
      flag.newFlag "f"+v12,     v1,        v21
      flag.newFlag "f"+v12,     v21,       v12
      flag.newFlag "f"+v12,      v12,       v1

      [v1,v2]=[v2,v3]  # current becomes previous
      [vert1,vert2]=[vert2,vert3]

  newpoly = flag.topoly()
  newpoly.name = "l" + poly.name
  #newpoly.xyz = adjustXYZ(newpoly, 3)
  #newpoly.xyz = canonicalXYZ(newpoly, 3)  # this tends to make results look like shit
  newpoly
