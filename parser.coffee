# Polyhédronisme
#===================================================================================================
#
# A toy for constructing and manipulating polyhedra and other meshes
#
# Copyright 2011, Anselm Levskaya
# Released under the MIT License

# Parser Routines
#===================================================================================================

specreplacements = [
  [/e/g, "aa"],   # e --> aa   (abbr. for explode)
  [/b/g, "ta"],   # b --> ta   (abbr. for bevel)
  [/o/g, "jj"],   # o --> jj   (abbr. for ortho)
  [/m/g, "kj"],   # m --> kj   (abbr. for meta)
  [/t(\d*)/g, "dk$1d"],  # t(n) --> dk(n)d  (dual operations)
  [/j/g, "dad"],  # j --> dad  (dual operations)
  [/s/g, "dgd"],  # s --> dgd  (dual operations)
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
      when "K" then poly.xyz = canonicalXYZ(poly, if n is 0 then 1 else n)
      when "C" then poly.xyz = canonicalize(poly, if n is 0 then 1 else n)
      when "A" then poly.xyz =    adjustXYZ(poly, if n is 0 then 1 else n)
      # experimental
      when "n" then poly     = insetN(poly, n)
      when "x" then poly     = extrudeN(poly, n)
      when "l" then poly     = stellaN(poly, n)
      when "z" then poly     = triangulate(poly,false)

    ops = ops.slice(0,-1);  # remove last character

  # Recenter polyhedra at origin (rarely needed)
  poly.xyz = recenter(poly.xyz, poly.getEdges())
  poly.xyz = rescale(poly.xyz)

  # Color the faces of the polyhedra for display
  poly = paintPolyhedron(poly)

  # return the poly object
  poly
