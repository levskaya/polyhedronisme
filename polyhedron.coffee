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
    for [i,f] in enumerate(@face)
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
    for clr in @face_colors
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


# get array of face centers
# faceCenters = (poly) ->
#   centers = []
#   for i in [0..poly.face.length-1]
#     centers[i] = [0,0,0]
#     for j in [0..poly.face[i].length-1] #avg vertex coords
#       centers[i] = add(centers[i], poly.xyz[poly.face[i][j]]) # add
#     centers[i] = mult(1.0/poly.face[i].length, centers[i]) # div by n

#   centers

# # get array of face centers
# faceNormals = (poly) ->
#   normals = []
#   for f in poly.face
#     normals.push normal(poly.xyz[v] for v in f)

#   normals



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
  height = 1
  poly = new polyhedron()
  poly.name = "Y#{n}"

  for i in [0..n-1] # vertex #'s 0...n-1 around one face
    poly.xyz.push [cos(i*theta), sin(i*theta), 0.2]
  poly.xyz.push [0,0,-1*height] # apex

  poly.face.push [n-1..0] # base
  for i in [0..n-1] # n triangular sides
    poly.face.push [i, (i+1)%n, n]

  poly.xyz = canonicalXYZ(poly,3)
  poly

