# csg.js

![](http://evanw.github.com/csg.js/image.png)

Constructive Solid Geometry (CSG) is a modeling technique that uses Boolean operations like union and intersection to combine 3D solids. This library implements CSG operations on meshes elegantly and concisely using BSP trees, and is meant to serve as an easily understandable implementation of the algorithm. All edge cases involving overlapping coplanar polygons in both solids are correctly handled.

Example usage:

    var cube = CSG.cube();
    var sphere = CSG.sphere({ radius: 1.3 });
    var polygons = cube.subtract(sphere).toPolygons();

# Documentation

[Detailed documentation](http://evanw.github.com/csg.js/docs/) can be automatically generated using [Docco](http://jashkenas.github.com/docco/).

# Demos

* [All CSG operations](http://evanw.github.com/csg.js/)
* [Coplanar test cases](http://evanw.github.com/csg.js/coplanar.html)
* [More test cases](http://evanw.github.com/csg.js/more.html)

# Implementation Details

All CSG operations are implemented in terms of two functions, `clipTo()` and `invert()`, which remove parts of a BSP tree inside another BSP tree and swap solid and empty space, respectively. To find the union of `a` and `b`, we want to remove everything in `a` inside `b` and everything in `b` inside `a`, then combine polygons from `a` and `b` into one solid:

    a.clipTo(b);
    b.clipTo(a);
    a.build(b.allPolygons());

The only tricky part is handling overlapping coplanar polygons in both trees. The code above keeps both copies, but we need to keep them in one tree and remove them in the other tree. To remove them from `b` we can clip the inverse of `b` against `a`. The code for union now looks like this:

    a.clipTo(b);
    b.clipTo(a);
    b.invert();
    b.clipTo(a);
    b.invert();
    a.build(b.allPolygons());

Subtraction and intersection naturally follow from set operations. If union is `A | B`, subtraction is `A - B = ~(~A | B)` and intersection is `A & B = ~(~A | ~B)` where `~` is the complement operator.

# License

Copyright (c) 2011 Evan Wallace (http://madebyevan.com/), under the [MIT license](http://www.opensource.org/licenses/mit-license.php).
