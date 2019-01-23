A polyhedral recipe looks like: 

[__op__][__op__] ... [__op__][__base__] _no spaces, just a string of characters_

where [__base__] is one of

&nbsp;| &nbsp;
-------|--------------------------
__T__ | Tetrahedron
__C__ | Cube
__O__ | Octahedron
__I__ | Icosahedron
__D__ | Dodecahedron
__Y__<i>N</i> | N-sided Pyramid
__P__<i>N</i> | N-sided Prism
__A__<i>N</i> | N-sided Anti-prism
__U__<i>N</i> | 2N-sided Cupola (N&ge;3, regular for 3,4,5)
__V__<i>N</i> | 2N-sided Anticupola (N&ge;3)
__J__<i>N</i> | Nth Johnson Solid (1 to 92)
  
and __op__ is one of these [polyhedron-building operators][1]

&nbsp;| &nbsp;
-------|--------------------------
  __d__ | [dual][11]
  __a__ | ambo
  __k__<i>N</i> | [kis][2] on N-sided faces (if no N, then general kis)
  __g__ | gyro
  __r__ | reflect
  __e__ | [expand][3] (=_aa_)
  __b__ | bevel (=_ta_) 
  __o__ | ortho (=_jj_) 
  __m__ | meta (=_k3j_) 
  __t__<i>N</i> | truncate vertices of degree N (=_dkNd_) if no N, then truncate all vertices) 
  __j__ | join (=_dad_) 
  __s__ | snub (=_dgd_) 
  __p__ | propellor
  __c__ | chamfer
  __w__ | whirl
  __q__ | quinto
  __l__ | loft
  <!--
  __h__ | half (caution: requires even-sided faces, and can produce digons)
  __n__ | needle
  __z__ | zip
  -->

Also, some newer, experimental operators:

&nbsp;| &nbsp;
:-----|:-------------------------
__P__ | from _Perspectiva Corporum Regularium_
__n__<i>N</i> | insetN 
__x__<i>N</i> | extrudeN 
__Z__ | triangulate
__H__(_inset_, _depth_) | hollow - useful for 3D printing, makes a hollow-faced shell version of the polyhedron,this applies the hollowing operator on all faces, insetting by _inset_ (scaled from 0 to 1), and  with a shell thickness of _depth_
__u__<i>N</i> | limited version of the Goldberg-Coxeter u_n operator (for triangular meshes only)


 There are more complicated, parameterized forms for __k__ and __n__: 

&nbsp;| &nbsp;
:-----|:-------------------------
__n__(_n_,_inset_,_depth_) | this applies the inset operator on _n_-sided faces, insetting by _inset_ scaled from 0 to 1, and extruding in or out along the normal by _depth_ (can be negative)
__k__(_n_,_depth_) | this applies the kis operator on _n_-sided faces, setting the pyramidal height out or in along the normal by _depth_ (can be negative)
  
Note that for most of the above operations, while the _topology_ of the result is uniquely specified, a great variety of _geometry_ is possible. For maximum flexibility, the above operators do not enforce convexity of the polyhedron, or planarity of the faces, at each step. If these properties are desired in the final result, the following geometric "refinement" operators can be used. These operators are for canonicalizing the polyhedral shape, and are mainly intended for making the more traditional, convex polyhedra more symmetric:

&nbsp;| &nbsp;
:-----|:-------------------------
__A__<i>N</i> | convex spherical <b>A</b>djustment. Iterates N times. May give more pleasing symmetry, but can be nstable for certain shapes. Usually an N of 20-40 is enough.
__C__<i>N</i> | proper <b>C</b>anonicalization, iteratively refines shape N times.  Flattens faces. A typical N is 200 or 300.

_Remember that these can blow up the geometry of nonconvex polyhedra._

### 3D Printing

You can export these shapes in forms appropriate for 3D printing by
shapeways. Export in VRML2 format to preserve face colors if you want
to use their colored fused-sand process.  You'll may want to rescale 
the exported geometry to a different size.

### More Information

 - [Wikipedia on Conway Polyhedral Notation][1]
 - [George W. Hart's Polyhedral Site][4]
 - [Source code at Github][5]

### Related Sites

- [Interactive Polyhedra Viewer][6] by [@tesseralis][7]
- [David McCooey's Visual Polyhedra][8]
- [Antiprism Polyhedral Modelling Software][9] - can do much more than conway operators!
- [Stella Polyhedral Modelling Software][10]

### Thanks

  George Hart - for his original pages, artworks and software characterizing higher polyhedra.
  Lars Huttar - for adding several new operators, and helping improve this site.
  Lei Willems - for inventing quinto.
  Everyone else - for all of your kind words and suggestions!

[1]:http://en.wikipedia.org/wiki/Conway_polyhedron_notation
[2]:https://en.wikipedia.org/wiki/Kleetope
[3]:https://en.wikipedia.org/wiki/Expansion_%28geometry%29
[4]:http://www.georgehart.com/
[5]:http://github.com/levskaya/polyhedronisme
[6]:https://polyhedra.tessera.li/
[7]:https://www.tessera.li/
[8]:http://dmccooey.com/polyhedra/index.html
[9]:http://www.antiprism.com/index.html
[10]:https://www.software3d.com/Stella.php
[11]:https://en.wikipedia.org/wiki/Dual_polyhedron