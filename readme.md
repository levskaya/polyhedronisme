polyHédronisme
===============

This is a project built to demonstrate the
[Conway Polyhedral Operators][1] and their extensions by
[George W. Hart][2]. It was heavily inspired by some [old code][3] of
George Hart's that did much the same in VRML.  Additional operators are being added to increase the range of procedurally generated forms possible.

One derives any of a vast family of polyhedral shapes by specifying a core geometric 'seed' and then modifying it with topological and geometrical operators that derive more complex shapes from this base set.

It should runs on any modern canvas-enabled browser.

It is written in coffeescript.

Background
----------

This project was inspired by the renaissance tome [Perspectiva Corporum Regularium][4] featuring engravings by Jost Amman after designs and drawings by Wenzel Jamnitzer.

### Future
- 2d vector export to SVG
- proper touch event handling for mobile
- geometric distortion operators, i.e.:
```
G x**2,y**2,z**2
poly.xyz = _.map(poly.xyz, ([x,y,z]) -> [x*x,y*y,z*z])
```
- more advanced (optional) parser for more complicated recipes
- Simple boolean operations on like-numbered faces: "cut and glue"
- Perhaps more full boolean operations on triangulated meshes +
  detriangulation to make complicated joins of, i.e. rotated cubes
  into a coherent mesh
- WebGL frontend to replace the hand-rolled, slow canvas renderer, for
  those browsers that can use it.

* * *
Text CC-BY, Code MIT License
2011 Anselm Levskaya

[1]: http://en.wikipedia.org/wiki/Conway_polyhedron_notation "Conway Operators"
[2]: http://www.georgehart.com/ "George W. Hart"
[3]: http://www.georgehart.com/virtual-polyhedra/conway_notation.html
[4]: http://bibliodyssey.blogspot.com/2009/08/jamnitzer-perspectiva.html

