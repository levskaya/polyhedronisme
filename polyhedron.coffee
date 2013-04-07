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
# Topology stored as set of "faces."  Each face is list of n vertex indices
# corresponding to one n-sided face.  Vertices listed clockwise as seen from outside.

faceToEdges = (face) ->
  edges = []
  [v1] = face[-1..]
  for v2 in face
    edges.push [v1,v2]
    v1 = v2
  edges

vertColors = (poly) ->
  vertcolors=[]
  for f,i in poly.face
    for v in f
      vertcolors[v] = poly.face_class[i]
  vertcolors

# Polyhedra Coloring Functions
#===================================================================================================

rwb_palette  = ["#ff7777","#dddddd","#889999","#fff0e5","#aa3333","#ff0000","#ffffff","#aaaaaa"]

# converts #xxxxxx / #xxx format into list of [r,g,b] floats
hextofloats = (hexstr)->
  if hexstr[0] is "#"
    hexstr = hexstr[1..]
  if hexstr.length is 3
    rgb = hexstr.split('').map(       (c)->parseInt(c+c, 16)/255 )
  else
    rgb = hexstr.match(/.{2}/g).map(  (c)->parseInt(c, 16)/255 )
  rgb

PALETTE = rwb_palette #GLOBAL
palette = (n) ->
  if n < PALETTE.length
    hextofloats(PALETTE[n])
  else
    hextofloats(PALETTE[PALETTE.length-1])

paintPolyhedron = (poly) ->
  # Color the faces of the polyhedra for display
  poly.face_class = []
  colormemory={}

  #memoized color assignment to faces of similar areas
  colorassign = (ar, colormemory) ->
    hash = round(100*ar)
    if hash of colormemory
      return colormemory[hash]
    else
      fclr = _.toArray(colormemory).length #palette _.toArray(colormemory).length
      colormemory[hash] = fclr
      return fclr

  for f in poly.face
    if COLOR_METHOD is "area"
      # color by face area (quick proxy for different kinds of faces) convexarea
      face_verts = (poly.xyz[v] for v in f)
      clr = colorassign(convexarea(face_verts), colormemory)
    else if COLOR_METHOD is "signature"
      face_verts = (poly.xyz[v] for v in f)
      clr = colorassign(faceSignature(face_verts), colormemory)
    else
      # color by face-sidedness
      clr = f.length-3

    poly.face_class.push clr
  console.log _.toArray(colormemory).length+" face classes"
  poly

# z sorts faces of poly
# -------------------------------------------------------------------------
sortfaces = (poly) ->
  #smallestZ = (x) -> _.sortBy(x,(a,b)->a[2]-b[2])[0]
  #closests = (smallestZ(poly.xyz[v] for v in f) for f in poly.face)
  centroids  = poly.centers()
  normals    = poly.normals()
  extents    = poly.extents()
  ray_origin = [0,0, (persp_z_max * persp_ratio - persp_z_min)/(1-persp_ratio)]

  # sort by centroid z-depth: not correct but more stable heuristic w. weird non-planar "polygons"
  zcentroidsort = (a,b) -> a[1][2]-b[1][2]
  zmaxsort = (a,b) -> a[3][2][1]-b[3][2][1]

  zsortIndex = _.zip([0..poly.face.length-1],centroids, normals, extents)
    #.sort(zmaxsort)
    .sort(zcentroidsort)
    .map((x)->x[0])

  # sort all face-associated properties
  poly.face = (poly.face[idx] for idx in zsortIndex)
  poly.face_class = (poly.face_class[idx] for idx in zsortIndex)

# determines if face represented by faceverts2 is on same side or opposite of
# refpoint from faceverts1
facesidedness = (faceverts1, faceverts2, refpoint) ->
  n1 = normal(faceverts1)
  c1 = calcCentroid(faceverts1)
  refside = dot(sub(refpoint,c1),n1)
  samesided=true
  oppsided=true
  for vert in faceverts2
    sidedness = dot(sub(vert,c1),n1)*refside
    # correct for numerical imprecision ~0 is as good as 0
    if abs(sidedness) < 1e-10 then sidedness = 0
    #console.log "sidedness ", sidedness
    if sidedness < 0 then samesided = false
    if sidedness > 0 then oppsided = false
  #console.log "SIDED", samesided, oppsided
  if samesided
    return 1
  else if oppsided
    return -1
  else
    return 0

sortfaces_fancy = (poly) ->
  ray_origin = [0,0, (persp_z_max * persp_ratio - persp_z_min)/(1-persp_ratio)]
  centroids  = poly.centers()
  extents    = poly.extents()
  zminsort = (a,b) -> a[1][2][0]-b[1][2][0]
  #zcentroidsort = (a,b) -> a[2][2]-b[2][2]

  extentsort = (a,b) ->
    [amin,amax] = a[1][2]
    [bmin,bmax] = b[1][2]
    diff = amin-bmin
    if abs(diff)<1e-10 then diff=0
    if diff != 0
      return diff
    diff = amax-bmax
    if abs(diff)<1e-10 then diff=0
    diff

  zsortIndex = _.zip([0..poly.face.length-1], extents, centroids)
    #.sort(zcentroidsort)
    .sort(extentsort)
    .map((x)->x[0])

  facenums = (idx for idx in zsortIndex)
  sortPtr = 0

  console.log facenums
  ii=0
  while sortPtr<facenums.length-1 and ii<500
    ii+=1
    aidx = facenums[sortPtr]
    for bidx,fno2 in facenums
      if fno2<=sortPtr then continue

      a=(poly.xyz[v] for v in poly.face[aidx])
      b=(poly.xyz[v] for v in poly.face[bidx])
      [[AminX,AmaxX],[AminY,AmaxY],[AminZ,AmaxZ]] = calcExtents(a)
      [[BminX,BmaxX],[BminY,BmaxY],[BminZ,BmaxZ]] = calcExtents(b)

      # If z-extents don't overlap, crude ordering OK
      if BmaxZ < AminZ
        console.log "wtf pre-ordering fail", aidx, "=", AminZ, AmaxZ, " ",bidx,"=", BminZ, BmaxZ
      if AmaxZ <= BminZ or BmaxZ <= AminZ
        #console.log "Zsep ", aidx, bidx
        #sortPtr+=1
        #break
        continue
        #console.log fno1, faces.length

      # If X,Y extents don't overlap, OK
      if AmaxX <= BminX or BmaxX <= AminX
        #console.log "Xsep ", aidx, bidx
        continue
      if AmaxY <= BminY or BmaxY <= AminY
        #console.log "Ysep ", aidx, bidx
        continue

      A_B_o = facesidedness(a, b, ray_origin)
      B_A_o = facesidedness(b, a, ray_origin)
      # A further back than B
      # B entirely on same side of A as view origin point
      # A entirely on opposite side of B as view origin point
      if A_B_o == 1 or B_A_o == -1
        #console.log "noswitch",aidx,bidx,A_B_o,B_A_o
        continue

      # B further back than A
      # A entirely on same side of B as view origin point
      # B entirely on opposite side of A as view origin point
      if B_A_o == 1 or A_B_o == -1
        facenums[sortPtr] = bidx
        facenums[fno2] = aidx
        console.log "switch",aidx,bidx,A_B_o,B_A_o
        #console.log "... order test", aidx, "=", AminZ, AmaxZ, " ",bidx,"=", BminZ, BmaxZ

        sortPtr=-1 #HACK to keep sortPtr unchanged below
        break

      #if fno2 == facenums.length-1
      #  console.log "end ", aidx
      #  sortPtr+=1
      #  break

      # All tests failed, emergency bailout
      console.log "Warning: All tests failed on ", aidx, bidx
      #sortPtr+=1
      break
    # no more faces to test against
    #console.log "end ", aidx
    sortPtr+=1

  console.log sortPtr, ii
  #facenums=facenums.reverse()
  poly.face = (poly.face[idx] for idx in facenums)
  poly.face_class = (poly.face_class[idx] for idx in facenums)

# Tests for line intersection in 2D
lineIsect2D = (a,b, sameLineException=false) ->
  [a1,a2] = a
  [b1,b2] = b
  da = sub2D(a2,a1)
  db = sub2D(b2,b1)
  crossp = cross2D(da,db)
  #console.log "crossp = ", crossp
  # Parallel Line Case
  if chop(crossp)==0 # parallel lines
    if cross2D(sub2D(b1,a1),da)==0 #colinear lines
      # Edge Case: collinear intersection
      # project points along a1->a2 line
      # pa1 = 0
      pa2 = dot2D(da,da)
      pb1 = dot2D(sub2D(b1,a1),da)
      pb2 = dot2D(sub2D(b2,a1),da)
      # ensure proper ordering of extent for bounds check
      if pb1>pb2 then [pb1,pb2]=[pb2,pb1]
      # special treatment of identical lines if specified
      if sameLineException and pb1==0 and pb2=pa2
        return false
      if pb2 <= 0 or pb1 >= pa2
        return false #colinear but NO intersection
      else
        return true #colinear intersection
    else
      return false # offset parallels
  # General Case
  du = chop(cross2D(sub2D(b1,a1),da)/crossp)
  dt = chop(cross2D(sub2D(b1,a1),db)/crossp)
  #console.log "du,dt = ", du,dt
  if 0<du<1 and 0<dt<1
    return true # interior intersection
  else
    return false

#complicated by shared colinear lines
polygonIsectTest = (a,b) ->
  # Build Edges from point list
  Aedges = _.zip(a[1..], a[0..a.length-2])
  Aedges.push([a[a.length-1],a[0]])
  Bedges = _.zip(b[1..], b[0..b.length-2])
  Bedges.push([b[b.length-1],b[0]])
  #console.log a,b,Aedges,Bedges
  # Check Edge Intersections
  for Aedge in Aedges
    for Bedge in Bedges
      if lineIsect2D(Aedge, Bedge, true) then return true
  # Check points of A interior to B
  for Apt in a
    cnt=0
    for Bedge in Bedges
      if lineIsect2D(Bedge, [Apt,[10000,0]]) then cnt+=1
    if cnt%2 != 0 then return true
  # Check points of B interior to A
  for Bpt in b
    cnt=0
    for Aedge in Aedges
      if lineIsect2D(Aedge, [Bpt,[10000,0]]) then cnt+=1
    if cnt%2 != 0 then return true
  # No Intersection!
  return false

polygonclip = (polyA, polyB) ->
  # attempt to use clipper lib implementing Vatti clipping algo
  scale_factor = 10000.0
  #subj_polygons = new ClipperLib.Polygons()
  #subj_polygon = new ClipperLib.Polygon()
  #for [x,y] in polyA:
  #  subj_polygon.push(new ClipperLib.IntPoint(x*scale_factor,y*scale_factor))
  #subj_polygons.push(subj_polygon)
  # clip_polygons = new ClipperLib.Polygons()
  # clip_polygon = new ClipperLib.Polygon()
  # for [x,y] in polyA:
  #   clip_polygon.push(new ClipperLib.IntPoint(x*scale_factor,y*scale_factor))
  # clip_polygons.push(clip_polygon)
  subj_polygons = []
  subj_polygon = []
  for [x,y] in polyA
    subj_polygon.push({'X':x*scale_factor,'Y':y*scale_factor})
  subj_polygons.push(subj_polygon)
  clip_polygons = []
  clip_polygon = []
  for [x,y] in polyA
    clip_polygon.push({'X':x*scale_factor,'Y':y*scale_factor})
  clip_polygons.push(clip_polygon)

  cpr = new ClipperLib.Clipper()
  cpr.AddPolygons(subj_polygons, ClipperLib.PolyType.ptSubject)
  cpr.AddPolygons(clip_polygons, ClipperLib.PolyType.ptClip)
  solution_polygons = new ClipperLib.Polygons()
  #cpr.Execute(clipType, solution_polygons, subject_fillType, clip_fillType)
  succeeded = cpr.Execute(0, solution_polygons, 1, 1)

  soln=[]
  for xy in solution_polygons[0]
    soln.push([xy['X']/scale_factor,xy['Y']/scale_factor])
  return soln

sortfaces_fancy2 = (poly,persp_z_max,persp_z_min,persp_ratio,perspective_scale) ->
  ray_origin = [0,0, (persp_z_max * persp_ratio - persp_z_min)/(1-persp_ratio)]
  centroids  = poly.centers()
  extents    = poly.extents()

  facePts = ((poly.xyz[v] for v in f) for f in poly.face)
  perspTfacePts = ((perspT(poly.xyz[v],persp_z_max,persp_z_min,persp_ratio,perspective_scale) for v in f) for f in poly.face)

  chop = (x) -> if abs(x) < 1e-10 then 0 else x
  extentsort = (a,b) ->
    [amin,amax] = a
    [bmin,bmax] = b
    diff = chop(amin-bmin)
    if diff != 0
      return diff
    diff = chop(amax-bmax)
    diff

  zsortIndex = _.zip([0..poly.face.length-1], extents, centroids)
    .sort((a,b) -> extentsort(a[1][2],b[1][2]))
    .map((x)->x[0])

  wtfs=[]

  sortF = (aidx,bidx) ->
    #a=(poly.xyz[v] for v in poly.face[aidx])
    #b=(poly.xyz[v] for v in poly.face[bidx])
    a = facePts[aidx]
    b = facePts[bidx]
    #aC = centroids[aidx]
    #bC = centroids[bidx]
    [[AminX,AmaxX],[AminY,AmaxY],[AminZ,AmaxZ]] =
      calcPerspExtents(a,persp_z_max,persp_z_min,persp_ratio,perspective_scale)
    [[BminX,BmaxX],[BminY,BmaxY],[BminZ,BmaxZ]] =
      calcPerspExtents(b,persp_z_max,persp_z_min,persp_ratio,perspective_scale)

    # If z-extents don't overlap, extent ordering holds
    if chop(AmaxZ-BminZ) <= 0 or chop(BmaxZ-AminZ) <= 0
      #console.log "Zsep ", aidx, bidx
      return chop(AminZ-BminZ)
    # If projected X,Y extents don't overlap
    if chop(AmaxX-BminX) <= 0 or chop(BmaxX-AminX) <= 0
      #console.log "Xsep ", aidx, bidx
      return chop(AminZ-BminZ)
    if chop(AmaxY-BminY) <= 0 or chop(BmaxY-AminY) <= 0
      #console.log "Ysep ", aidx, bidx
      return chop(AminZ-BminZ)

    pA = perspTfacePts[aidx]
    pB = perspTfacePts[bidx]
    if not polygonIsectTest(pA,pB)
      #console.log "no isect", aidx, bidx
      return chop(AminZ-BminZ)
    #else
      #console.log "isect", aidx, bidx

    A_B_o = facesidedness(a, b, ray_origin)
    B_A_o = facesidedness(b, a, ray_origin)
    # B entirely on same side of A as view origin point
    # A entirely on opposite side of B as view origin point
    if A_B_o == 1 or B_A_o == -1
      console.log "noswitch",aidx,bidx,A_B_o,B_A_o
      return -1

    # A entirely on same side of B as view origin point
    # B entirely on opposite side of A as view origin point
    if B_A_o == 1 or A_B_o == -1
      console.log "switch",aidx,bidx,A_B_o,B_A_o
      return 1

    if wtfs.length == 0
      wtfs.push(poly.face[aidx])
      wtfs.push(poly.face[bidx])
      console.log "wtf", aidx, bidx, A_B_o, B_A_o

      pA = perspTfacePts[aidx]
      pB = perspTfacePts[bidx]
      if not polygonIsectTest(pA,pB)
        console.log "no isect", aidx, bidx
      else
        console.log "isect", aidx, bidx

    return chop(AminZ-BminZ)

  zsortIndex.sort(sortF)

  # sort all face-associated properties
  poly.face = (poly.face[idx] for idx in zsortIndex)
  poly.face_class = (poly.face_class[idx] for idx in zsortIndex)
  return wtfs




# Main Polyhedron Class
class polyhedron
  constructor: (verts,faces,name) ->      # constructor of initially null polyhedron
    @face = faces or new Array()   # array of faces.          face.length = # faces
    @xyz  = verts or new Array()   # array of vertex coords.  xyz.length = # of vertices
    @name = name  or "null polyhedron"

  data: () ->   # informative string
    nEdges = @face.length + @xyz.length - 2 # E = V + F - 2
    "(#{@face.length} faces, #{nEdges} edges, #{@xyz.length} vertices)"

  edges: ->
    finalset={}
    uniqedges=[]
    alledges = _.map(@face, faceToEdges)
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

  centers: ->
    # get array of face centers
    centers_array = []
    for f in @face
      fcenter = [0,0,0]
      for v in f #avg vertex coords
        fcenter = add(fcenter, @xyz[v]) # add
      centers_array.push mult(1.0/f.length, fcenter) # div by n
    # return face-ordered array of centroids
    centers_array

  normals: ->
  # get array of face centers
    normals_array = []
    for f in @face
      normals_array.push normal(@xyz[v] for v in f)
    normals_array

  extents: ->
  # get x,y,z bounds for each polygon in polyhedron
    extents_array = []
    for f in @face
      [x,y,z] = @xyz[f[0]]
      [maxX,minX] = [x,x]
      [maxY,minY] = [y,y]
      [maxZ,minZ] = [z,z]
      for v in f
        [x,y,z] = @xyz[v]
        maxX = if maxX < x then x else maxX
        minX = if minX > x then x else minX
        maxY = if maxY < y then y else maxY
        minY = if minY > y then y else minY
        maxZ = if maxZ < z then z else maxZ
        minZ = if minZ > z then z else minZ
      extents_array.push [[minX,maxX],[minY,maxY],[minZ,maxZ]]
    extents_array

  # Export / Formatting Routines --------------------------------------------------

  # produces vanilla OBJ files for import into 3d apps
  toOBJ: () ->
    objstr="#Produced by polyHédronisme http://levskaya.github.com/polyhedronisme\n"
    objstr+="group #{@name}\n"
    objstr+="#vertices\n"
    for v in @xyz
      objstr += "v #{v[0]} #{v[1]} #{v[2]}\n"

    objstr += "#normal vector defs \n"
    for f in @face
      norm = normal(@xyz[v] for v in f)
      objstr += "vn #{norm[0]} #{norm[1]} #{norm[2]}\n"

    objstr += "#face defs \n"
    for f,i in @face
      objstr += "f "
      for v in f
        objstr += "#{v+1}//#{i+1} "
      objstr += "\n"

    objstr

  toX3D: () ->
    SCALE_FACTOR = .01 #ShapeWays uses 1unit = 1meter, so reduce to 1cm scale
    # opening cruft
    x3dstr='''
      <?xml version="1.0" encoding ="UTF-8"?>
      <X3D profile="Interchange" version="3.0">
      <head>
      <component name="Rendering" level="3"/>
      <meta name="generator" content="Polyhedronisme"/>
      <meta name="version" content="0.1.0"/>
      </head>
      <Scene>
      <Shape>
      <IndexedFaceSet normalPerVertex="false" coordIndex="
      '''
    # face indices
    for f in @face
      for v in f
        x3dstr+="#{v} "
      x3dstr+='-1\n'
    x3dstr+='">\n'

    # per-face Color
    x3dstr+='<Color color="'
    for cl in vertColors(this)#@face_class
      clr=palette cl
      x3dstr+="#{clr[0]} #{clr[1]} #{clr[2]} "
    x3dstr+='"/>'

    # re-scaled xyz coordinates
    x3dstr+='<Coordinate point="'
    for v in @xyz
      x3dstr+="#{v[0]*SCALE_FACTOR} #{v[1]*SCALE_FACTOR} #{v[2]*SCALE_FACTOR} "
    x3dstr+='"/>\n'

      # end cruft
    x3dstr+='''
      </IndexedFaceSet>
      </Shape>
      </Scene>
      </X3D>'''

    x3dstr

  toVRML: () ->
    SCALE_FACTOR = .01 #ShapeWays uses 1unit = 1meter, so reduce to 1cm scale
    # opening cruft
    x3dstr='''
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
                      [
            '''
    # re-scaled xyz coordinates
    for v in @xyz
      x3dstr+="#{v[0]*SCALE_FACTOR} #{v[1]*SCALE_FACTOR} #{v[2]*SCALE_FACTOR},"
    x3dstr=x3dstr[0..-2]
    x3dstr+='''
                          ]
                      }
                      color Color
                      {
                        color
                        [
                   '''
    # per-face Color
    for cl in @face_class
      clr=palette cl
      x3dstr+="#{clr[0]} #{clr[1]} #{clr[2]} ,"
    x3dstr=x3dstr[0..-2]
    x3dstr+='''
                      ]
                    }
                    colorPerVertex FALSE
                    coordIndex
                    [
                   '''
    # face indices
    for f in @face
      for v in f
        x3dstr+="#{v}, "
      x3dstr+='-1,'
    x3dstr=x3dstr[0..-2]
    x3dstr+='''
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
                  }
                   '''
    x3dstr

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
    poly.xyz.push [-cos(i*theta), -sin(i*theta),  -h]
  for i in [0..n-1] # vertex #'s n...2n-1 around other
    poly.xyz.push [-cos(i*theta), -sin(i*theta), h]

  poly.face.push [n-1..0]   #top
  poly.face.push [n..2*n-1] #bottom
  for i in [0..n-1] #n square sides
    poly.face.push [i, (i+1)%n, (i+1)%n+n, i+n]

  poly = adjustXYZ(poly,1)
  poly

antiprism = (n) ->
  theta = 2*PI/n # pie angle
  h = sqrt(1-4/(4+2*cos(theta/2)-2*cos(theta)))
  r = sqrt(1-h*h)
  f = sqrt(h*h + pow(r*cos(theta/2),2) )
  # correction so edge midpoints (not vertices) on unit sphere
  r = -r/f
  h = -h/f
  poly = new polyhedron()
  poly.name = "A#{n}"

  for i in [0..n-1] # vertex #'s 0...n-1 around one face
    poly.xyz.push [r * cos(i*theta), r * sin(i*theta), h]
  for i in [0..n-1] # vertex #'s n...2n-1 around other
    poly.xyz.push [r * cos((i+0.5)*theta), r * sin((i+0.5)*theta), -h]

  poly.face.push [n-1..0]   #top
  poly.face.push [n..2*n-1] #bottom
  for i in [0..n-1] #2n triangular sides
    poly.face.push [i, (i+1)%n, i+n]
    poly.face.push [i, i+n, ((n+i-1)%n+n)]

  poly = adjustXYZ(poly,1)
  poly

pyramid = (n) ->
  theta = 2*PI/n # pie angle
  height = 1
  poly = new polyhedron()
  poly.name = "Y#{n}"

  for i in [0..n-1] # vertex #'s 0...n-1 around one face
    poly.xyz.push [-cos(i*theta), -sin(i*theta), -0.2]
  poly.xyz.push [0,0, height] # apex

  poly.face.push [n-1..0] # base
  for i in [0..n-1] # n triangular sides
    poly.face.push [i, (i+1)%n, n]

  poly = canonicalXYZ(poly,3)
  poly

