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


# Math / Vector / Matrix Functions
#===================================================================================================

# Math is primal, people.
# import math functions to local namespace
random = Math.random
round = Math.round
floor = Math.floor
sqrt  = Math.sqrt
sin   = Math.sin
cos   = Math.cos
tan   = Math.tan
atan  = Math.atan
pow   = Math.pow
abs   = Math.abs
PI    = Math.PI

# for python-style enumerated for-in loops
# --should use "for [i,x] in AR then do (i,x)->" idiom instead
enumerate = (ar) ->  [i,ar[i]] for i in [0..ar.length-1]

# often useful
randomchoice = (array)->
  n = floor(random()*array.length)
  array[n]

# scalar multiplication
mult = (c,vec) -> [c*vec[0],c*vec[1],c*vec[2]]

# element-wise multiply
_mult = (vec1, vec2) -> [vec1[0]*vec2[0],vec1[1]*vec2[1],vec1[2]*vec2[2]]

# 3d vector addition
add = (vec1, vec2) -> [vec1[0]+vec2[0],vec1[1]+vec2[1],vec1[2]+vec2[2]]

# 3d vector subtraction
sub = (vec1, vec2) -> [vec1[0]-vec2[0],vec1[1]-vec2[1],vec1[2]-vec2[2]]

# 3d dot product
dot = (vec1, vec2) -> vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2]

# vector norm
mag = (vec) -> sqrt(dot(vec,vec))

# vector magnitude squared
mag2 = (vec) -> dot(vec,vec)

# makes vector unit length
unit = (vec) -> mult( 1/sqrt(mag2(vec)), vec)

# cross product d1 x d2
cross = (d1, d2) -> [ d1[1]*d2[2]-d1[2]*d2[1], d1[2]*d2[0]-d1[0]*d2[2],  d1[0]*d2[1]-d1[1]*d2[0] ]

midpoint = (vec1, vec2) -> mult(1/2.0,add(vec1,vec2))

# parametric segment between vec1, vec2 w. parameter t ranging from 0 to 1
tween = (vec1,vec2,t) -> [ (1-t)*vec1[0] + t*vec2[0], (1-t)*vec1[1] + t*vec2[1], (1-t)*vec1[2] + t*vec2[2] ]

# uses above to go one-third of the way along vec1->vec2 line
oneThird = (vec1, vec2) -> tween(vec1, vec2, 1/3.0)

# reflect 3vec in unit sphere, spherical reciprocal
reciprocal = (vec) -> mult( 1.0/mag2(vec), vec)

# point where line v1...v2 tangent to an origin sphere
tangentPoint= (v1,v2) ->
  d = sub v2, v1
  sub v1, mult(dot(d,v1)/mag2(d),d)

# distance of line v1...v2 to origin
edgeDist = (v1,v2) ->
  sqrt mag2(tangentPoint v1, v2)

# find vector orthogonal to plane of 3 pts
# -- do the below algos assume this be normalized or not?
orthogonal = (v1,v2,v3) ->
  # adjacent edge vectors
  d1 = sub v2, v1
  d2 = sub v3, v2
  # cross product
  cross d1, d2

# find first element common to 3 sets by brute force search
intersect = (set1, set2, set3) ->
  for s1 in set1
    for s2 in set2
      if s1 is s2
        for s3 in set3
          if s1 is s3
            return s1
  return null # oh noes!

# copies array of arrays by value (deep copy)
copyVecArray = (vecArray)->
  newVecArray = new Array(vecArray.length)
  for i in [0..vecArray.length-1]
    newVecArray[i] = vecArray[i][0..]
  newVecArray

# 3d matrix vector multiply
mv3 = (mat,vec) ->
  #Ghetto custom def of matrix-vector mult
  #example matrix: [[a,b,c],[d,e,f],[g,h,i]]
  [mat[0][0]*vec[0]+mat[0][1]*vec[1]+mat[0][2]*vec[2],
  mat[1][0]*vec[0]+mat[1][1]*vec[1]+mat[1][2]*vec[2],
  mat[2][0]*vec[0]+mat[2][1]*vec[1]+mat[2][2]*vec[2]]

# 3d matrix matrix multiply
mm3 = (A,B) ->
  [[A[0][0]*B[0][0]+A[0][1]*B[1][0]+A[0][2]*B[2][0],
  A[0][0]*B[0][1]+A[0][1]*B[1][1]+A[0][2]*B[2][1],
  A[0][0]*B[0][2]+A[0][1]*B[1][2]+A[0][2]*B[2][2]],
  [A[1][0]*B[0][0]+A[1][1]*B[1][0]+A[1][2]*B[2][0],
  A[1][0]*B[0][1]+A[1][1]*B[1][1]+A[1][2]*B[2][1],
  A[1][0]*B[0][2]+A[1][1]*B[1][2]+A[1][2]*B[2][2]],
  [A[2][0]*B[0][0]+A[2][1]*B[1][0]+A[2][2]*B[2][0],
  A[2][0]*B[0][1]+A[2][1]*B[1][1]+A[2][2]*B[2][1],
  A[2][0]*B[0][2]+A[2][1]*B[1][2]+A[2][2]*B[2][2]]]

# Rotation Matrix
# Totally ghetto, not at all in agreement with euler angles!
# use quaternions instead
rotm = (phi,theta,psi)->
    xy_mat = [
      [cos(phi), -1.0*sin(phi),  0.0],
      [sin(phi),      cos(phi),  0.0],
      [0.0,                0.0,  1.0]]
    yz_mat = [
      [cos(theta), 0, -1.0*sin(theta)],
      [         0, 1,               0],
      [sin(theta), 0,      cos(theta)]]
    xz_mat = [
      [1.0,        0,             0],
      [  0, cos(psi), -1.0*sin(psi)],
      [  0, sin(psi),      cos(psi)]]

    mm3(xz_mat, mm3(yz_mat,xy_mat))

# Perspective Transform
# assumes world's been rotated appropriately such that Z is depth
# scales perspective such that inside depth regions min_real_depth <--> max_real_depth
# perspective lengths vary no more than:   desired_ratio
# with target dimension of roughly length: desired_length
perspT = (vec3, max_real_depth, min_real_depth, desired_ratio, desired_length) ->
  z0          = (max_real_depth * desired_ratio - min_real_depth)/(1-desired_ratio)
  scalefactor =  desired_length * desired_ratio/(1-desired_ratio)

  # projected [X, Y]
  [scalefactor*vec3[0]/(vec3[2]+z0), scalefactor*vec3[1]/(vec3[2]+z0)]

