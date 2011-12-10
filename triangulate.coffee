# Polyhédronisme
#===================================================================================================
#
# A toy for constructing and manipulating polyhedra and other meshes
#
# Copyright 2011, Anselm Levskaya
# Released under the MIT License
#

# Polyhedra Functions
#===================================================================================================
#
# Set of routines for transforming N-face meshes into triangular meshes, necessary for exporting
# STL or VRML for 3D Printing.
#

# Ear-based triangulation of 2d faces, takes array of 2d coords in the face ordering
# Returns indices of the new diagonal lines to cut.
#
# assumes planarity of course, so this isn't the ideal algo for making aesthetically pleasing
# "flattening" choices in distorted polyhedral planes.
#
getDiagonals = (verts)->
  limiter = 999
  diagonals = []
  ear = []
  facelen = verts.length

  XOR = (x, y) -> (x or y) and not (x and y)
  Area2     = (Va,Vb,Vc)   -> (Vb[0]-Va[0])*(Vc[1]-Va[1]) - (Vc[0]-Va[0])*(Vb[1]-Va[1])
  Left      = (Va, Vb, Vc) ->  Area2(Va, Vb, Vc) > 0
  LeftOn    = (Va, Vb, Vc) ->  Area2(Va, Vb, Vc) >= 0
  Collinear = (Va, Vb, Vc) ->  Area2(Va, Vb, Vc) is 0

  Between   = (Va, Vb, Vc) ->
    return false if Collinear(Va, Vb, Vc)
    unless Va[0] is Vb[0]
      (Va[0] <= Vc[0]) and (Vc[0] <= Vb[0]) or (Va[0] >= Vc[0]) and (Vc[0] >= Vb[0])
    else
      (Va[1] <= Vc[1]) and (Vc[1] <= Vb[1]) or (Va[1] >= Vc[1]) and (Vc[1] >= Vb[1])

  IntersectProp = (Va, Vb, Vc, Vd) ->
    return false if Collinear(Va, Vb, Vc) or Collinear(Va, Vb, Vd) or Collinear(Vc, Vd, Va) or Collinear(Vc, Vd, Vb)
    XOR(Left(Va, Vb, Vc), Left(Va, Vb, Vd)) and XOR(Left(Vc, Vd, Va), Left(Vc, Vd, Vb))

  Intersect = (Va, Vb, Vc, Vd) ->
    if IntersectProp(Va, Vb, Vc, Vd)
      true
    else
      if Between(Va, Vb, Vc) or Between(Va, Vb, Vd) or Between(Vc, Vd, Va) or Between(Vc, Vd, Vb)
        true
      else
        false

  InCone = (a, b) ->
    a1 = (a+1+facelen)%facelen
    a0 = (a-1+facelen)%facelen
    if LeftOn(verts[a], verts[a1], verts[a0])
      return (Left(verts[a], verts[b], verts[a0]) and Left(verts[b], verts[a], verts[a1]))
    not (LeftOn(verts[a], verts[b], verts[a1]) and LeftOn(verts[b], verts[a], verts[a0]))

  Diagonalie = (a, b) ->
    c = 0
    loop
      c1 = (c+1+facelen)%facelen
      if (c isnt a) and (c1 isnt a) and (c isnt b) and (c1 isnt b) and IntersectProp(verts[a], verts[b], verts[c], verts[c1])
        return false
      c  = (c+1+facelen)%facelen
      break unless c isnt 0
    true

  Diagonal = (a, b) -> InCone(a, b) and InCone(b, a) and Diagonalie(a, b)

  v1 = 0
  loop
    v2 = (v1+1+facelen)%facelen#v1.next
    v0 = (v1-1+facelen)%facelen#v1.prev
    ear[v1] = Diagonal(v0, v2)
    v1 = (v1+1+facelen)%facelen
    break if v1 is 0

  origIdx = [0..facelen-1]
  n = facelen#verts.length
  z = limiter
  head = 0 #??
  while z > 0 and n > 3
    z -= 1
    v2 = head
    y = limiter
    loop
      y -= 1
      broke = false
      if ear[v2]
        v3 = (v2+1+facelen)%facelen#v2.next
        v4 = (v3+1+facelen)%facelen#v3.next
        v1 = (v2-1+facelen)%facelen#v2.prev
        v0 = (v1-1+facelen)%facelen#v1.prev
        diagonals.push [ origIdx[v1], origIdx[v3] ]
        ear[v1] = Diagonal(v0, v3)
        ear[v3] = Diagonal(v1, v4)
        #v1.next = v3
        #v3.prev = v1
        verts   = verts[0..v2].concat(verts[v3..])
        origIdx = origIdx[0..v2].concat(origIdx[v3..])
        if v0>v2 then v0 -= 1
        if v1>v2 then v1 -= 1
        if v3>v2 then v3 -= 1
        if v4>v2 then v4 -= 1
        facelen--
        head = v3
        n--
        broke = true
      v2 = (v2+1+facelen)%facelen#v2.next
      break unless y > 0 and not broke and v2 isnt head

  #return diagonals
  diagonals

# equates triplets of numbers if they can be rotated into identity
triEq = (tri1,tri2)->
    if ((tri1[0] is tri2[0]) and (tri1[1] is tri2[1]) and (tri1[2] is tri2[2]))\
    or  (tri1[0] is tri2[1]) and (tri1[1] is tri2[2]) and (tri1[2] is tri2[0])\
    or  (tri1[0] is tri2[2]) and (tri1[1] is tri2[0]) and (tri1[2] is tri2[1])
      true
    else
      false

# god-awful but working hack to turn diagonals into triangles
# switch to an edge-matching algo, it would be 10x simpler
diagsToTris = (f,diags)->
  edges = []
  redges = []
  # get edges from faces as assoc arrays
  for [v1,v2] in ([i,(i+1)%f.length] for i in [0..f.length-1])
    edges[v1]  = [v2]
    redges[v2] = [v1]
  for d in diags # push the diagonals into the assoc arrays in both directions!
    edges[d[0]].push d[1]
    edges[d[1]].push d[0]
    redges[d[0]].push d[1]
    redges[d[1]].push d[0]
  tris=[]
  for d in diags  #orig N-face, N-2 triangles from the N-3 diagonals
    for e1 in edges[d[1]] # edge after diag
      for e2 in redges[d[0]] # edge before diag
        if e1 is e2 # if they meet we have a triangle!
          tris.push [d[0],d[1],e1]
    for e1 in edges[d[0]] # same as above for other dir along diagonal
      for e2 in redges[d[1]]
        if e1 is e2
          tris.push [d[1],d[0],e1]
  # unfortunately the above duplicates triangles, so filter out repeats
  uniques = [tris.pop()]
  for tri in tris
    already_present = false
    for extant_tri in uniques
      if triEq tri, extant_tri
        already_present=true
        break
    if not already_present then uniques.push tri

  uniques

# driver routine, projects 3d face to 2d, get diagonals then triangles,
# then builds new polyhedron out of them, preserving original face colors
triangulate = (poly, colors)->
  colors = colors || false
  console.log "Triangulating faces of #{poly.name}..."

  newpoly = new polyhedron()
  newpoly.xyz = clone poly.xyz
  newpoly.face_colors = [ ]
  # iterate over triplets of faces v1,v2,v3
  for f,i in poly.face
    if f.length > 3
      TwoDface = project2dface(poly.xyz[v] for v in f)
      diags = getDiagonals(TwoDface)
      tris  = diagsToTris(f,diags)
      for tri,j in tris
        newpoly.face.push [ f[tri[0]], f[tri[1]], f[tri[2]] ]
        if colors then newpoly.face_colors.push poly.face_colors[i]
    else
      newpoly.face.push [ f[0], f[1], f[2] ]
      if colors then newpoly.face_colors.push poly.face_colors[i]
  newpoly.name = poly.name # don't change the name for export
  newpoly
