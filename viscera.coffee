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
# Random Offcuts, maybe interesting operators, etc.
#===================================================================================================

# Really Starry StellaN
# ------------------------------------------------------------------------------------------
#
# This is the bad dog version that makes non-convex star-shaped faces
#
stellaN = (poly)->
  console.log "Taking stella of #{poly.name}..."

  flag = new polyflag()
  for [i,p] in enumerate(poly.xyz)
    # each old vertex is a new vertex
    flag.newV "v#{i}", p

  centers = faceCenters(poly)

  for [i,f] in enumerate(poly.face)
    v1 = "v"+f[f.length-2]
    v2 = "v"+f[f.length-1]
    vert1 = poly.xyz[f[f.length-2]]
    vert2 = poly.xyz[f[f.length-1]]
    for v in f
      v3 = "v"+v
      vert3 = poly.xyz[v]
      #if f.length is n or n is 0
      foundAny = true
      v12=v1+"~"+v2
      v21=v2+"~"+v1
      flag.newV v12, midpoint( midpoint(vert1,vert2), centers[i] )
      #old, stellated face
      flag.newFlag "f#{i}",      v1,       v12
      flag.newFlag "f#{i}",      v12,       v2
      #new 1/2 of the two new triangles
      flag.newFlag "f"+v12,      v1,       v12
      flag.newFlag "f"+v12,     v12,       v21
      flag.newFlag "f"+v12,     v21,        v1
      v1=v2  # current becomes previous
      vert1=vert2

  newpoly = flag.topoly()
  newpoly.name = "*" + poly.name
  console.log newpoly
  #newpoly.xyz = adjustXYZ(newpoly, 3)
  #newpoly.xyz = canonicalXYZ(newpoly, 3)  # this tends to make results look like shit
  newpoly

