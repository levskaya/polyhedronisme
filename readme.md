polyHédronisme
===============

This is a project built to demonstrate the
[Conway Polyhedral Operators][1] and their extensions by
[George W. Hart][2]. It was heavily inspired by some [old code][3] of
George Hart's that did much the same in VRML.  Additional operators
are being added to increase the range of procedurally generated forms
possible.

One derives any of a vast family of polyhedral shapes by specifying a
core geometric 'seed' and then modifying it with topological and
geometrical operators that derive more complex shapes from this base
set.

It should runs on any modern canvas-enabled browser.

This project was also inspired by the renaissance tome
[Perspectiva Corporum Regularium][4] featuring engravings by Jost
Amman after designs and drawings by Wenzel Jamnitzer.

It is written in coffeescript. Uses jQuery, underscore, [PEG.js][5] parser for parsing recipes, and
[Eli Grey's][eli] BlobBuilder, FileSaver, canvas-toBlob for saving files.

### Future
- Switch rendering to Three.JS
- WebGL frontend to replace the hand-rolled, slow canvas renderer, for
  those browsers that can use it.
- 2d vector export to SVG
- integration of [csg.js][6] lib for complex mesh joins
- geometric distortion operators, i.e.:
```
G x**2,y**2,z**2 => poly.xyz = _.map(poly.xyz, ([x,y,z]) -> [x*x,y*y,z*z])
```
- Simple boolean operations on like-numbered faces: "cut and glue"
- proper touch event handling for mobile

* * *
Text CC-BY, Code MIT License
2011 Anselm Levskaya

[1]: http://en.wikipedia.org/wiki/Conway_polyhedron_notation "Conway Operators"
[2]: http://www.georgehart.com/ "George W. Hart"
[3]: http://www.georgehart.com/virtual-polyhedra/conway_notation.html
[4]: http://bibliodyssey.blogspot.com/2009/08/jamnitzer-perspectiva.html
[5]: http://pegjs.majda.cz/
[6]: https://github.com/evanw/csg.js
[eli]: https://github.com/eligrey
