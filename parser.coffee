# Polyhédronisme
#===================================================================================================
#
# A toy for constructing and manipulating polyhedra and other meshes
#
# Copyright 2011, Anselm Levskaya
# Released under the MIT License

# Parser Routines
#===================================================================================================

# fairly straightforward Parser Expression Grammar spec for simple
# operator-chain-on-base-polyhedra recipes
PEG_parser_spec = '''
/* series of opspecs */
start  = opspec+

/* opspec one of:
 A  - single letter
 A3 - single letter and float
 B(5,4.3,3) - function call format w. float args
*/
opspec =
   let:opcode args:opargs {return {"op":let,"args":args};}
/ let:opcode float:float     {return {"op":let,"args":[float]};}
/ let:opcode                     {return {"op":let,"args":[]};}

/*
parentheses surrounding comma-delimited list of floats i.e.
( 1 , 3.2, 4 ) or (1) or (2,3)
*/
opargs = "("
           num:( float:float ","? {return float} )+
         ")" {return num;}

/* just a letter */
opcode = op:[a-zA-Z] {return op;}

/* standard numerical types */
int   = digits:[0-9-]+   { return parseInt(digits.join(""), 10);  }
float = digits:[0-9.-]+  { return parseFloat(digits.join(""), 10); }
'''
op_parser = PEG.buildParser PEG_parser_spec

#applies func fn to array args i.e. f, [1,2,3] -> f(1,2,3)
dispatch = (fn, args) -> fn.apply(this, args || [])

basemap = {
  "T": tetrahedron
  "O": octahedron
  "C": cube
  "I": icosahedron
  "D": dodecahedron
  "P": prism     #takes integer arg
  "A": antiprism #takes integer arg
  "Y": pyramid   #takes integer arg
  }

opmap = {
  "d": dual
  "k": kisN
  "a": ambo
  "g": gyro
  "p": propellor
  "r": reflect
  "c": chamfer
  "w": whirl
  "n": insetN
  "x": extrudeN
  "l": stellaN
  "z": triangulate
  "K": canonicalXYZ
  "C": canonicalize
  "A": adjustXYZ
  }

#list of basic equivalences, easier to replace before parsing
specreplacements = [
  [/e/g, "aa"],   # e --> aa   (abbr. for explode)
  [/b/g, "ta"],   # b --> ta   (abbr. for bevel)
  [/o/g, "jj"],   # o --> jj   (abbr. for ortho)
  [/m/g, "kj"],   # m --> kj   (abbr. for meta)
  [/t(\d*)/g, "dk$1d"],  # t(n) --> dk(n)d  (dual operations)
  [/j/g, "dad"],  # j --> dad  (dual operations) # Why not j --> da ?
  [/s/g, "dgd"],  # s --> dgd  (dual operations) # Why not s --> dg ?
  [/dd/g, ""],    # dd --> null  (order 2)
  [/ad/g, "a"],   # ad --> a   (a_ = ad_)
  [/gd/g, "g"],   # gd --> g   (g_ = gd_)
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
newgeneratePoly = (notation) ->
  #poly = new polyhedron()

  ops_spec = getOps(notation)
  oplist = op_parser.parse(ops_spec).reverse()

  op = oplist.shift()
  basefunc = basemap[op["op"]]
  baseargs = op["args"]
  poly     = dispatch basefunc, baseargs

  #console.log "base", poly

  for op in oplist
    opfunc = opmap[op["op"]]
    opargs = [poly].concat(op["args"])
    #console.log opargs
    poly   = dispatch opfunc, opargs

  #console.log "final", poly

  # Recenter polyhedra at origin (rarely needed)
  poly.xyz = recenter(poly.xyz, poly.edges())
  poly.xyz = rescale(poly.xyz)

  # Color the faces of the polyhedra for display
  poly = paintPolyhedron(poly)

  # return the poly object
  poly