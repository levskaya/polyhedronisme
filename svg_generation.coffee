# Polyhédronisme
#===================================================================================================
#
# Routines for making nice SVG exports of polyhedra (2d clipping of partially obscured faces,
# etc. so that the vector representation isn't filled with junk)
#
# Copyright 2011, Anselm Levskaya
# Released under the MIT License


to2Dfaces = (poly,tvec) ->
  tvec ||= [3,3,3]

  poly.xyz = _.map(poly.xyz, (x)->mv3(globRotM,x))

  normals = _.map(poly.normals(), (x)->mv3(globRotM,x))

  # z sort faces
  sortfaces(poly)
  poly.face = poly.face.reverse()

  twoDfaces = []
  facecolors=[]
  for face,fno in poly.face
    twoDface=[]
    # loop around face, defining polygon
    for v in face
      [x,y] = perspT(add(tvec,poly.xyz[v]),persp_z_max,persp_z_min,persp_ratio,perspective_scale)
      twoDface.push [x+_2d_x_offset, y+_2d_y_offset]
    twoDfaces.push twoDface
    console.log "twoDface", _.map(twoDface,(ar)->" ["+ar[0]+","+ar[1]+"] ")

    # use pre-computed colors
    clr = palette poly.face_class[fno]

    # shade based on simple cosine illumination factor
    face_verts = (poly.xyz[v] for v in face)
    illum = dot(normal(face_verts), unit([1,-1,0]))
    clr   = mult((illum/2.0+.5)*0.7+0.3,clr)

    facecolors.push clr

    # if PaintMode is "fill" or PaintMode is "fillstroke"
    #   ctx.fillStyle = "rgba(#{round(clr[0]*255)}, #{round(clr[1]*255)}, #{round(clr[2]*255)}, #{1.0})"
    #   ctx.fill()
    # # make cartoon stroke (=black) / realistic stroke an option (=below)
    #   ctx.strokeStyle = "rgba(#{round(clr[0]*255)}, #{round(clr[1]*255)}, #{round(clr[2]*255)}, #{1.0})"
    #   ctx.stroke()
    # if PaintMode is "fillstroke"
    #   ctx.fillStyle = "rgba(#{round(clr[0]*255)}, #{round(clr[1]*255)}, #{round(clr[2]*255)}, #{1.0})"
    #   ctx.fill()
    #   ctx.strokeStyle = "rgba(0,0,0, .3)"  # light lines, less cartoony, more render-y
    #   ctx.stroke()
    # if PaintMode is "stroke"
    #   ctx.strokeStyle = "rgba(0,0,0, .8)"
    #   ctx.stroke()
  [twoDfaces, facecolors]

# Returns 2 times the signed triangle area. The result is positive if
# abc is ccw, negative if abc is cw, zero if abc is degenerate.
Signed2DTriArea = (a, b, c) -> (a[0] - c[0]) * (b[1] - c[1]) - (a[1] - c[1]) * (b[0] - c[0])

# Test if segments ab and cd overlap. If they do, compute and return
# intersection t value along ab and intersection position p
Test2DSegmentSegment = ( a, b, c, d ) ->
  #float &t, Point &p
  # Sign of areas correspond to which side of ab points c and d are
  a1 = Signed2DTriArea(a, b, d) # Compute winding of abd (+ or -)
  a2 = Signed2DTriArea(a, b, c) # To intersect, must have sign opposite of a1

  # If c and d are on different sides of ab, areas have different signs
  if (a1 isnt 0.0 and a2 isnt 0.0 and a1 * a2 < 0.0)
    # Compute signs for a and b with respect to segment cd
    a3 = Signed2DTriArea(c, d, a) # Compute winding of cda (+ or -)
    # Since area is constant a1 - a2 = a3 - a4, or a4 = a3 + a2 - a1
    a4 = Signed2DTriArea(c, d, b) # Must have opposite sign of a3 float a4 = a3 + a2 - a1;
    # Points a and b on different sides of cd if areas have different signs
    if (a3 * a4 < 0.0)
      # Segments intersect. Find intersection point along L(t) = a + t * (b - a).
      # Given height h1 of an over cd and height h2 of b over cd,
      # t = h1 / (h1 - h2) = (b*h1/2) / (b*h1/2 - b*h2/2) = a3 / (a3 - a4),
      # where b (the base of the triangles cda and cdb, i.e., the length
      # of cd) cancels out.
      t = a3 / (a3 - a4)
      p = [ a[0] + t * (b[0] - a[0]) , a[1] + t * (b[1] - a[1]) ]
      #console.log "ISECT: ",p[0]+" "+p[1]
      return p

  # Segments not intersecting (or collinear)
  false


twoDpolygonarea = (xys) ->
    area = 0.0
    v1 = xys[0]
    for v2 in xys[1..]
      #area of sub-triangle
      area += (v1[0]+v2[0])*(v1[1]-v2[1])
      v1 = v2 # shift over one
    area*0.5

clockwisepoly = (xys) -> twoDpolygonarea(xys)>0

idxof = (vlist, vtest)->
  #console.log "idxof",vlist, vtest
  for v,i in vlist
    if abs(v[0]-vtest[0])+abs(v[1]-vtest[1]) < 1e-6
      return i
  -1

intersectionGraph=(_2dfaceA, _2dfaceB)->
  verts = []
  for v in _2dfaceA
    if idxof(verts,v) is -1
      verts.push v
  for v in _2dfaceB
    if idxof(verts,v) is -1
      verts.push v

  #console.log _2dfaceA
  edgesA = []
  idx0 = idxof(verts, _2dfaceA[_2dfaceA.length-1])
  for v in _2dfaceA
    idx=idxof(verts,v)
    edgesA.push [idx0, idx]
    idx0=idx

  edgesB = []
  idx0 = idxof(verts,_2dfaceB[_2dfaceB.length-1])
  for v in _2dfaceB
    idx=idxof(verts,v)
    edgesB.push [idx0,idx]
    idx0=idx

  newvertsA = ([] for i in [0..edgesA.length-1])
  newvertsB = ([] for i in [0..edgesB.length-1])

  intersectpts=[]

  for eA,i in edgesA
    for eB,j in edgesB
      if p=Test2DSegmentSegment(verts[eA[0]],verts[eA[1]],verts[eB[0]],verts[eB[1]])
        newvertsA[i].push p
        newvertsB[j].push p
        if idxof(verts,p) is -1 then verts.push p
        intersectpts.push idxof(verts,p)

  #console.log "isectpts", intersectpts

  newedgesA = []
  for vA,i in newvertsA
    if newvertsA[i].length > 0
      newvertsA[i].sort((pa,pb)->mag(pa-verts[edgesA[i][0]])-mag(pb-verts[edgesA[i][0]]))#sort new points on edge along edge direction
      #add new subedges
      newedgesA.push [edgesA[i][0], idxof(verts,newvertsA[i][0])]
      for v,k in newvertsA[i][1..]
        newedgesA.push [idxof(verts,newvertsA[i][k]),idxof(verts,newvertsA[i][k+1])]
      newedgesA.push [idxof(verts,newvertsA[i][newvertsA[i].length-1]), edgesA[i][1]]
    else
      #no new subedges, add original
      newedgesA.push edgesA[i]

  newedgesB = []
  for vB,i in newvertsB
    if newvertsB[i].length > 0
      newvertsB[i].sort((pa,pb)->mag(pa-verts[edgesB[i][0]])-mag(pb-verts[edgesB[i][0]]))
      newedgesB.push [edgesB[i][0],idxof(verts,newvertsB[i][0])]
      for v,k in newvertsB[i][1..]
        newedgesB.push [idxof(verts,newvertsB[i][k]),idxof(verts,newvertsB[i][k+1])]
      newedgesB.push [idxof(verts,newvertsB[i][newvertsB[i].length-1]),edgesB[i][1]]
    else
      newedgesB.push edgesB[i]

  #console.log "intersectionGraph",  edgesA, newedgesA, newedgesB

  [verts, newedgesA, newedgesB, intersectpts]


pointInPoly=(_2dfaceA, pt)->
  verts = []
  for v in _2dfaceA
    if idxof(verts,v) is -1
      verts.push v

  edgesA = []
  idx0 = idxof(verts,_2dfaceA[_2dfaceA.length-1])
  for v in _2dfaceA
    idx=idxof(verts,v)
    edgesA.push [idx0,idx]
    idx0=idx

  isects=0
  for eA in edgesA
    if Test2DSegmentSegment(verts[eA[0]],verts[eA[1]],pt, [pt[0],10000])
      isects++

  console.log "pointInpoly called : "+isects, (isects%2!=0)

  if isects%2==0 then false else true


AminusB = (A, B)->
  # get intersection vertices and new edges from A and B overlap
  [verts, newedgesA, newedgesB, intersectpts] = intersectionGraph(A,B.reverse())

  Agraph={}
  for e in newedgesA
    Agraph[e[0]] = e[1]
  Bgraph={}
  for e in newedgesB
    Bgraph[e[0]] = e[1]

  sullied = {}#(false for v in verts)

  if intersectpts.length is 0 #special case, no intersecting lines
    if pointInPoly(B,A[0]) # either B totally occludes A or they are disjoint
      return []
    else
      return [A]

  newpolys=[]
  for pt in intersectpts
    state = "B"
    if sullied[pt+state] then continue
    newpoly=[]
    movingpt = pt
    sullied[movingpt+state]=true
    #newpoly.push movingpt

    itrCTR=0
    while movingpt isnt pt and itrCTR<1000
      if Agraph[pt] and Bgraph[pt] #on an intersection point
        state = if state is "A" then "B" else "A"
      if state is "A"
        movingpt = Agraph[movingpt]
      else
        movingpt = Bgraph[movingpt]
      newpoly.push movingpt
      sullied[movingpt+state]=true

    newpolys.push (verts[idx] for idx in newpoly)

  #console.log "AminusB", newpolys
  #return only correctly oriented polygons
  _.filter(newpolys, clockwisepoly)


getSVG = (poly) ->
  console.log poly.data()

  [twoDfaces, facecolors] = to2Dfaces(poly)

  frontfaces=[]
  frontfaces.push twoDfaces[0]
  for fA in twoDfaces[1..]
    fset=[fA]
    for fB in frontfaces
      fsetnew=[]
      for f in fset
        fsetnew = fsetnew.concat AminusB(f,fB)
      fset=fsetnew

    frontfaces = frontfaces.concat fset

  frontfaces = frontfaces.reverse()

  console.log frontfaces.length, frontfaces

  svg_open="""<?xml version="1.0" encoding="utf-8"?>
  <!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1 Tiny//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11-tiny.dtd">
  <svg version="1.1" baseProfile="tiny" id="Layer_1" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink"
           x="0px" y="0px" width="1024px" height="768px" viewBox="0 0 1024 768" overflow="inherit" xml:space="preserve">

  """
  svg_close="""</svg>"""

  svgtxt=svg_open
  for f in frontfaces
    svgtxt+='<path fill="#25DF00" stroke="#000000" d="'
    svgtxt+='M'+f[0][0]+","+f[0][1]
    for v in f[1..]
      svgtxt+='L'+v[0]+","+v[1]
    svgtxt+='z"/>\n'
  svgtxt+=svg_close

  console.log svgtxt
