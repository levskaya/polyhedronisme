# Polyhédronisme
#===================================================================================================
#
# A toy for constructing and manipulating polyhedra and other meshes
#
# Copyright 2011, Anselm Levskaya
# Released under the MIT License


# GLOBALS
#===================================================================================================
ctx={} # for global access to canvas context
CANVAS_WIDTH  = 500 #canvas dims
CANVAS_HEIGHT = 400 #canvas dims
globPolys={} # constructed polyhedras

globRotM = clone eye3
globlastRotM = clone eye3
#globtheta = 0 # rotation and projective mapping parameters
#globphi   = 0
perspective_scale = 800
persp_z_max = 5
persp_z_min = 0
persp_ratio = 0.8
_2d_x_offset = CANVAS_WIDTH/2 #300
_2d_y_offset = CANVAS_HEIGHT/2 #140

globtime = new Date() # for animation

BG_CLEAR = true # clear background or colored?
BG_COLOR = "rgba(255,255,255,1.0)" # background color
COLOR_METHOD = "area" #"signature"
PaintMode = "fillstroke"
ctx_linewidth = 0.5 # for outline of faces

# Mouse Event Variables
MOUSEDOWN=false
LastMouseX=0
LastMouseY=0
LastSphVec=[1,0,0] #for 3d trackball

# random grabbag of polyhedra
#DEFAULT_RECIPES = [
#  "C2dakD","oC20kkkT","kn4C40A0dA4","opD",
#  "lT","lK5oC","knD","dn6x4K5bT","oox4P7",
#  "n18n18n9n9n9soxY9","khD","lhD"]
DEFAULT_RECIPES = ["T"]

# File-saving objects used to export txt/canvas-png
saveText = (text, filename) ->
  BB = window.BlobBuilder || window.WebKitBlobBuilder || window.MozBlobBuilder
  bb = new BB()
  bb.append(text)
  saveAs(bb.getBlob("text/plain;charset="+document.characterSet), filename)

# parses URL string for polyhedron recipe, for bookmarking
# should use #! href format instead
parseurl = () ->
  urlParams = {}
  a = /\+/g  # Regex for replacing addition symbol with a space
  r = /([^&=]+)=?([^&]*)/g
  d = (s) -> decodeURIComponent(s.replace(a, " "))
  q = window.location.search.substring(1)

  while e=r.exec(q)
    urlParams[d(e[1])] = d(e[2])
  urlParams

# Drawing Functions
#==================================================================================================

# init canvas element
# -------------------------------------------------------------------------------
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

# clear canvas
# -----------------------------------------------------------------------------------
clear = ->
  if BG_CLEAR
    ctx.clearRect 0, 0, CANVAS_WIDTH, CANVAS_HEIGHT
  else
    ctx.clearRect 0, 0, CANVAS_WIDTH, CANVAS_HEIGHT
    ctx.fillStyle = BG_COLOR
    ctx.fillRect 0, 0, CANVAS_WIDTH, CANVAS_HEIGHT


# main drawing routine for polyhedra
#===================================================================================================
drawpoly = (poly,tvec) ->
  tvec ||= [3,3,3]

  #centers = _.map(poly.centers(), (x)->mv3(rotm(rot[0],rot[1],rot[2]),x))
  #oldfaces = ("#{fno}" for fno in [0..centers.length-1])

  # rotate poly in 3d
  oldxyz = _.map(poly.xyz, (x)->x)
  #poly.xyz = _.map(poly.xyz, (x)->mv3(rotm(rot[0],rot[1],rot[2]),x))
  poly.xyz = _.map(poly.xyz, (x)->mv3(globRotM,x))

  # z sort faces
  sortfaces(poly)
  window.polyobj = clone poly #for debugging inspection
  #sortfaces_fancy(poly)
  #wtfs = sortfaces_fancy2(poly,persp_z_max,persp_z_min,persp_ratio,perspective_scale)

  wtfs = [polygonclip((poly.xyz[v] for v in poly.face[0]),(poly.xyz[v] for v in poly.face[2]))]

  #for face culling
  normals = poly.normals()

  for face,fno in poly.face
    if dot(normals[fno],[0,0,1]) < 0 then continue
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
    clr = palette poly.face_class[fno]

    # shade based on simple cosine illumination factor
    face_verts = (poly.xyz[v] for v in face)
    illum = dot(normal(face_verts), unit([1,-1,0]))
    clr   = mult((illum/2.0+.5)*0.7+0.3,clr)

    #console.log dot(normal(face_verts), calcCentroid(face_verts))
    #  clr = "rgba(0,0,0,1)"

    if PaintMode is "fill" or PaintMode is "fillstroke"
      ctx.fillStyle = "rgba(#{round(clr[0]*255)}, #{round(clr[1]*255)}, #{round(clr[2]*255)}, #{1.0})"
      ctx.fill()
    # make cartoon stroke (=black) / realistic stroke an option (=below)
      ctx.strokeStyle = "rgba(#{round(clr[0]*255)}, #{round(clr[1]*255)}, #{round(clr[2]*255)}, #{1.0})"
      ctx.stroke()
    if PaintMode is "fillstroke"
      ctx.fillStyle = "rgba(#{round(clr[0]*255)}, #{round(clr[1]*255)}, #{round(clr[2]*255)}, #{1.0})"
      ctx.fill()
      ctx.strokeStyle = "rgba(0,0,0, .3)"  # light lines, less cartoony, more render-y
      ctx.stroke()
    if PaintMode is "stroke"
      ctx.strokeStyle = "rgba(0,0,0, .8)"
      ctx.stroke()

  for face in wtfs
    ctx.beginPath()
    # move to first vertex of face
    v0 = face[face.length-1]
    [x,y] = v0 #perspT(add(tvec,poly.xyz[v0]), persp_z_max,persp_z_min,persp_ratio,perspective_scale)
    ctx.moveTo(10*x+_2d_x_offset, 10*y+_2d_y_offset)
    # loop around face, defining polygon
    for v in face
      [x,y] = v#perspT(add(tvec,poly.xyz[v]),persp_z_max,persp_z_min,persp_ratio,perspective_scale)
      ctx.lineTo(10*x+_2d_x_offset, 10*y+_2d_y_offset)

    ctx.lineWidth = 1.0
    ctx.strokeStyle = "rgba(0,0,0, 1.0)"
    ctx.stroke()
    ctx.lineWidth = ctx_linewidth

  #for face,fno in poly.face
  #  ctx.textAlign = "center"
  #  ctx.fillStyle = "rgba(0,0,0,1)"
  #  [x,y] = perspT(add(tvec, centers[fno]),persp_z_max,persp_z_min,persp_ratio,perspective_scale)
  #  ctx.fillText(oldfaces[fno],x+_2d_x_offset,y+_2d_y_offset)

  # reset coords, for setting absolute rotation, as poly is passed by ref
  poly.xyz = oldxyz


# draw polyhedra just once
# -----------------------------------------------------------------------------------
drawShape = ->
  clear()
  for p,i in globPolys
    drawpoly(p,[0+3*i,0,3])


# loop for animation
# -----------------------------------------------------------------------------------
animateShape = ->
  clear()
  globtheta=(2*Math.PI)/180.0*globtime.getSeconds()*0.1
  for p,i in globPolys
    drawpoly(p,[0+3*i,0,3])
  setTimeout(animateShape, 100)


# Initialization and Basic UI
#===================================================================================================

$( -> #wait for page to load

  init() #init canvas

  urlParams = parseurl() #see if recipe is spec'd in URL
  if "recipe" of urlParams
    specs=[urlParams["recipe"]]
    $("#spec").val(specs)
  else
    specs=[randomchoice(DEFAULT_RECIPES)]
    $("#spec").val(specs)

  # set initial palette spec
  $("#palette").val( PALETTE.reduce((x,y)->x+" "+y) )

  # construct the polyhedra from spec
  globPolys = _.map(specs, (x)->newgeneratePoly(x))

  # draw it
  drawShape()


  # Event Handlers
  # ----------------------------------------------------

  # when spec changes in input, parse and draw new polyhedra
  $("#spec").change((e) ->
    specs = $("#spec").val().split(/\s+/g)[0..1] #only allow one recipe for now
    globPolys = _.map(specs, (x)->newgeneratePoly(x) )
    #animateShape()
    #window.location.replace("?recipe="+specs[0])
    drawShape()
  )

  # when spec changes in input, parse and draw new polyhedra
  $("#palette").change((e) ->
    PALETTE = $(this).val().split(/\s+/g)
    drawShape()
  )

  # Basic manipulation: rotation and scaling of geometry
  # ----------------------------------------------------

  # mousewheel changes scale of drawing
  $("#poly").mousewheel( (e,delta, deltaX, deltaY)->
    e.preventDefault()
    perspective_scale*=(10+delta)/10
    drawShape()
  )

  # Implement standard trackball routines
  # ---------------------------------------
  $("#poly").mousedown( (e)->
    e.preventDefault()
    MOUSEDOWN=true
    LastMouseX=e.clientX-$(this).offset().left #relative mouse coords
    LastMouseY=e.clientY-($(this).offset().top-$(window).scrollTop())
    #calculate inverse projection of point to sphere
    tmpvec=invperspT(LastMouseX,LastMouseY,_2d_x_offset,_2d_y_offset,persp_z_max,persp_z_min,persp_ratio,perspective_scale)
    if tmpvec[0]*tmpvec[1]*tmpvec[2]*0 is 0  #quick NaN check
      LastSphVec=tmpvec
    globlastRotM = clone globRotM #copy last transform state
    #console.log LastSphVec[0],LastSphVec[1],LastSphVec[2]
  )
  $("#poly").mouseup( (e)->
    e.preventDefault()
    MOUSEDOWN=false
  )
  $("#poly").mouseleave( (e)->
    e.preventDefault()
    MOUSEDOWN=false
  )
  $("#poly").mousemove( (e)->
    e.preventDefault()
    if MOUSEDOWN
      MouseX=e.clientX-$(this).offset().left
      MouseY=e.clientY-($(this).offset().top-$(window).scrollTop())
      SphVec=invperspT(MouseX,MouseY,_2d_x_offset,_2d_y_offset,persp_z_max,persp_z_min,persp_ratio,perspective_scale)

      # quick NaN check
      if SphVec[0]*SphVec[1]*SphVec[2]*0 is 0 and LastSphVec[0]*LastSphVec[1]*LastSphVec[2]*0 is 0
        globRotM = mm3(getVec2VecRotM(LastSphVec,SphVec),globlastRotM)

      drawShape()
  )

  # State control via some buttons
  # ---------------------------------------

  $("#strokeonly").click((e) ->
    PaintMode = "stroke"
    drawShape()
  )
  $("#fillonly").click((e) ->
    PaintMode = "fill"
    drawShape()
  )
  $("#fillandstroke").click((e) ->
    PaintMode = "fillstroke"
    drawShape()
  )

  $("#siderot").click((e) ->
    globRotM = vec_rotm(PI/2,0,1,0)
    drawShape()
  )
  $("#toprot").click((e) ->
    globRotM = vec_rotm(PI/2,1,0,0)
    drawShape()
  )
  $("#frontrot").click((e) ->
    globRotM = rotm(0,0,0)
    drawShape()
  )

  # Export Options
  # ---------------------------------------

  $("#pngsavebutton").click((e)->
    getSVG(clone globPolys[0])
    canvas=$("#poly")[0]
    #this works, but is janky
    #window.location = canvas.toDataURL("image/png")
    spec = $("#spec").val().split(/\s+/g)[0]
    filename = "polyhedronisme-"+spec.replace(/\([^\)]+\)/g, "")+".png"
    canvas.toBlob( (blob)->saveAs(blob, filename) )
  )
  $("#objsavebutton").click((e)->
    objtxt = globPolys[0].toOBJ()
    spec = $("#spec").val().split(/\s+/g)[0]
    filename = "polyhedronisme-"+spec.replace(/\([^\)]+\)/g, "")+".obj"
    saveText(objtxt,filename)
  )
  $("#x3dsavebutton").click((e)->
    triangulated = triangulate(globPolys[0],true) #triangulate to preserve face_colors for 3d printing
    x3dtxt = triangulated.toVRML()
    spec = $("#spec").val().split(/\s+/g)[0]
    #filename = "polyhedronisme-"+spec+".x3d"
    filename = "polyhedronisme-"+spec.replace(/\([^\)]+\)/g, "")+".wrl"
    saveText(x3dtxt,filename)
  )

)
