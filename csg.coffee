# Polyhédronisme
#===================================================================================================
#
# Constructive Solid Geometry Operations
#
# Copyright 2011, Anselm Levskaya
# Released under the MIT License

#converts [x,y,z] to object that CSG.js uses
arrayToCSGVertex = (pos,norm)->
  pos  or= [0,0,0]
  norm or= unit(pos)
  # return CSG Vertex
  csgvec = new CSG.Vertex(new CSG.Vector(pos[0],pos[1],pos[2]), new CSG.Vector(norm[0],norm[1],norm[2]))
  csgvec

#converts native polyhedral mesh to a BSP tree that CSG.js uses
polyToCSG = (poly)->
  # convert polygons to CSG.js format
  CSGPolys=[]
  for f in poly.face
    CSGPolys.push new CSG.Polygon((arrayToCSGVertex(poly.xyz[v]) for v in f),{})

  #console.log poly, CSGPolys
  # build BSP CSG tree
  csg=CSG.fromPolygons(CSGPolys)
  csg

# get index of vertex vert in list of vertices
vertIndex = (vertlist,vert)->
  for [i,v] in enumerate(vertlist)
    if mag(sub(vert, v))<1E-6 then return i
  return -1

# remove vertex duplicates
uniquifyverts = (vertlist)->
  uniqverts=[]
  uniqverts.push vertlist[0]
  for v in vertlist
    already_present=false
    for u in uniqverts
      if mag(sub(v,u))<1E-6 then already_present = true
    unless already_present
      uniqverts.push v
  uniqverts

# does testvert lie between evert1->evert2 on their line?
# this is a hack and no optimization has gone into it
edge_vert_isect = (evert1,evert2,testvert)->
  dv = sub(evert2, evert1)
  dT = sub(testvert, evert1)
  a  = dot(unit(dv), dT)
  b  = mag(sub(unit(dv),unit(dT)))
  if b<1E-6 and 0<mag(dT)<mag(dv)
    true
  else
    false

# CSG routine creates polygons that intersect internally at edges but does not include these
# edge-intersecting vertices in the larger face, this wreaks havoc on the topological operators
# that assume faces include all vertices with edge-sharing faces
# this fixes that problem, it's probably slow as hell on large meshes.
meshFix = (poly)->
  faces = poly.face
  verts = poly.xyz

  newfaces=[]
  for f in faces
    newf=[]
    edge0=verts[f[f.length-1]]
    for v in f
      edge1=verts[v]
      for [vno,testv] in enumerate(verts)
        if edge_vert_isect(edge0,edge1,testv)
          newf.push vno
      newf.push v
      edge0=edge1
    newfaces.push newf

  newpoly = new polyhedron()
  newpoly.face = newfaces
  newpoly.xyz = clone verts
  newpoly

#convert CSG BSP tree back into a coherent polygon mesh
CSGToPoly = (csgpoly)->
  faces = []
  verts = []
  for csg in csgpoly.toPolygons()
    for v in csg.vertices
      verts.push [v.pos.x,v.pos.y,v.pos.z]
  verts = uniquifyverts(verts)
  for csg in csgpoly.toPolygons()
    faces.push (vertIndex(verts,[v.pos.x,v.pos.y,v.pos.z]) for v in csg.vertices)

  poly = new polyhedron()
  poly.face = faces
  poly.xyz = verts

  # fix edge vertex inclusion problem, return result
  meshFix poly

csgUnion = (polyA,polyB)->
  csgA = polyToCSG(polyA)
  csgB = polyToCSG(polyB)
  #console.log csgA
  csgresult = csgA.union(csgB)
  CSGToPoly csgresult

csgIntersect = (polyA,polyB)->
  csgA = polyToCSG(polyA)
  csgB = polyToCSG(polyB)
  csgresult = csgA.intersect(csgB)
  CSGToPoly csgresult

csgSubtract = (polyA,polyB)->
  csgA = polyToCSG(polyA)
  csgB = polyToCSG(polyB)
  csgresult = csgA.subtract(csgB)
  CSGToPoly csgresult


edgesToFace = (edges)->
  face=[]
  edict={}
  for e in edges
    #console.log "e2f" ,e
    edict[e[0]] =e[1]
  v0=edges[0][0]
  v =edges[0][1]
  #console.log v,v0
  #face.push v0
  face.push v
  itCTR=0
  while v isnt v0
    #console.log v
    v = edict[v]
    face.push v
    itCTR++
    if itCTR>1000 #protect against infinite recursion during mistakes
      console.log "Bad edges to face join, have a neverending face:", edges
      break
  #console.log "e2f",edges,"->",face
  face

uniquifyedges = (edgelist)->
  uniqverts=[]
  uniqverts.push vertlist[0]
  for v in vertlist
    already_present=false
    for u in uniqverts
      if mag(sub(v,u))<1E-6 then already_present = true
    unless already_present
      uniqverts.push v
  uniqverts

faceprint = (face) ->
  str=""
  for v in face
    str+="#{v}->"
  str[0..str.length-3]

joinFaces = (faces)->
  alledges=[]
  alledges=(_.reduce( _.map(faces,faceToEdges) , ((memo,x)->memo.concat(x)), alledges))[0..]
  #console.log faces, alledges.length, alledges[0..]
  newedges=[]
  while alledges.length>0
    e=alledges.pop()
    nuked=false
    if alledges.length > 0
      for [i,e2] in enumerate(alledges)
        if e[0] is e2[1] and e[1] is e2[0]
          nuked=true
          alledges.splice(i,1) #knock out the other annihilating edge
      unless nuked
        newedges.push e
    else
        newedges.push e

  #console.log _.map(faces,faceprint),"joinFaces", _.map(newedges,faceprint)
  edgesToFace newedges

uniteFaces = (poly)->
  face_idx = [0..poly.face.length-1]
  dpoly = dual poly
  #for f in dpoly.face
  #  if f.length<4
  #    console.log faceprint f
  edict={}
  for e in dpoly.getEdges()
    edict[e[0]]=[]
  for e in dpoly.getEdges()
    edict[e[0]].push(e[1])
    edict[e[1]].push(e[0])

  normals = poly.normals()

  facepiles=[]
  facepile=[]
  #totally unoptimized N*(N-1)/2 comps ~N^2!
  Fi = face_idx.pop()
  facepile.push Fi
  unused_faces=[]
  itCTR=0
  while face_idx.length>0 and itCTR<1000000
    itCTR++
    Fi2 = face_idx.pop()
    #if mag(sub(unit(normals[Fi]),unit(normals[Fi2])))<1E-6
      #console.log Fi,edict[Fi],Fi2
    if mag(sub(unit(normals[Fi]),unit(normals[Fi2])))<1E-6 and (edict[Fi].indexOf(Fi2) isnt -1)
      facepile.push Fi2
    else
      unused_faces.unshift Fi2
    if face_idx.length is 0 and unused_faces.length > 0
      face_idx=unused_faces
      unused_faces=[]
      console.log Fi,"pile",facepile.length
      facepiles.push _.map(facepile,(idx)->poly.face[idx])#facepile
      facepile=[]
      Fi = face_idx.pop()
      facepile.push Fi

  #console.log facepiles
  #console.log _.map(facepiles,joinFaces)
  newfaces=_.map(facepiles,joinFaces)
  newpoly = new polyhedron()
  newpoly.xyz = clone poly.xyz
  newpoly.face = newfaces
  newpoly