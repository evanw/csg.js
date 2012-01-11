var _CSGDEBUG=true;

/*
Todo: lots!

Polygon reduction: after CSG operations we still have lots of adjacent coplanar
polygon fragments, even though we keep track of hierarchical splits. 
These fragments should be joined into one or more convex polygons.
Probably this can be done efficiently using a sweep based algorithm.
To prevent problems due to rounding errors we should keep track for each
side of a polygon by which plane it was split. If we know that two sides
of two polygons were split by the same plane then we can determine exactly
if a point is a T junction. Otherwise we can only find this out by interpolating
line segments, yielding small rounding errors.

CSG.Polygon.prototype.expand: we should rotate the generated cylinders and spheres
such that they join nicely with the extruded plane even if a low precision is used.

Add a invertedClipTo() function so we can get rid of the inverts() here:
    b.invert();
    b.clipTo(a);
    b.invert();

splitByPlane() is built around the original splitPolygon() function, there's 
some juggling of arrays. Can be integrated more tightly and efficiently.

The flip() and invert() methods change the original object, while all other
methods return a newly created object. If we change this we can get rid
of some clone()s, and multiple polygons can for example reference the same
Plane object instead of keeping cloned copies.

-joostn
*/



// Constructive Solid Geometry (CSG) is a modeling technique that uses Boolean
// operations like union and intersection to combine 3D solids. This library
// implements CSG operations on meshes elegantly and concisely using BSP trees,
// and is meant to serve as an easily understandable implementation of the
// algorithm. All edge cases involving overlapping coplanar polygons in both
// solids are correctly handled.
// 
// Example usage:
// 
//     var cube = CSG.cube();
//     var sphere = CSG.sphere({ radius: 1.3 });
//     var polygons = cube.subtract(sphere).toPolygons();
// 
// ## Implementation Details
// 
// All CSG operations are implemented in terms of two functions, `clipTo()` and
// `invert()`, which remove parts of a BSP tree inside another BSP tree and swap
// solid and empty space, respectively. To find the union of `a` and `b`, we
// want to remove everything in `a` inside `b` and everything in `b` inside `a`,
// then combine polygons from `a` and `b` into one solid:
// 
//     a.clipTo(b);
//     b.clipTo(a);
//     a.build(b.allPolygons());
// 
// The only tricky part is handling overlapping coplanar polygons in both trees.
// The code above keeps both copies, but we need to keep them in one tree and
// remove them in the other tree. To remove them from `b` we can clip the
// inverse of `b` against `a`. The code for union now looks like this:
// 
//     a.clipTo(b);
//     b.clipTo(a);
//     b.invert();
//     b.clipTo(a);
//     b.invert();
//     a.build(b.allPolygons());
// 
// Subtraction and intersection naturally follow from set operations. If
// union is `A | B`, subtraction is `A - B = ~(~A | B)` and intersection is
// `A & B = ~(~A | ~B)` where `~` is the complement operator.
// 
// ## License
// 
// Copyright (c) 2011 Evan Wallace (http://madebyevan.com/), under the MIT license.
// Parts Copyright (c) 2012 Joost Nieuwenhuijse (joost@newhouse.nl) under the MIT license.

// # class CSG

// Holds a binary space partition tree representing a 3D solid. Two solids can
// be combined using the `union()`, `subtract()`, and `intersect()` methods.

CSG = function() {
  this.polygons = [];
};

// Construct a CSG solid from a list of `CSG.Polygon` instances.
CSG.fromPolygons = function(polygons) {
  var csg = new CSG();
  csg.polygons = polygons;
  return csg;
};

CSG.prototype = {
  clone: function() {
    var csg = new CSG();
    csg.polygons = this.polygons.map(function(p) { return p.clone(); });
    return csg;
  },

  toPolygons: function() {
    return this.polygons;
  },

  // Return a new CSG solid representing space in either this solid or in the
  // solid `csg`. Neither this solid nor the solid `csg` are modified.
  // 
  //     A.union(B)
  // 
  //     +-------+            +-------+
  //     |       |            |       |
  //     |   A   |            |       |
  //     |    +--+----+   =   |       +----+
  //     +----+--+    |       +----+       |
  //          |   B   |            |       |
  //          |       |            |       |
  //          +-------+            +-------+
  // 
  union: function(csg) {
    var a = new CSG.Tree(this.polygons);
    var b = new CSG.Tree(csg.polygons);
    a.clipTo(b);
    b.clipTo(a);
    b.invert();
    b.clipTo(a);
    b.invert();
    a.addPolygons(b.allPolygons());
    return CSG.fromPolygons(a.allPolygons());
  },

  // Return a new CSG solid representing space in this solid but not in the
  // solid `csg`. Neither this solid nor the solid `csg` are modified.
  // 
  //     A.subtract(B)
  // 
  //     +-------+            +-------+
  //     |       |            |       |
  //     |   A   |            |       |
  //     |    +--+----+   =   |    +--+
  //     +----+--+    |       +----+
  //          |   B   |
  //          |       |
  //          +-------+
  // 
  subtract: function(csg) {
    var a = new CSG.Tree(this.polygons);
    var b = new CSG.Tree(csg.polygons);
    a.invert();
    a.clipTo(b);
    b.clipTo(a);
    b.invert();
    b.clipTo(a);
    b.invert();
    a.addPolygons(b.allPolygons());
    a.invert();
    return CSG.fromPolygons(a.allPolygons());
  },

  // Return a new CSG solid representing space both this solid and in the
  // solid `csg`. Neither this solid nor the solid `csg` are modified.
  // 
  //     A.intersect(B)
  // 
  //     +-------+
  //     |       |
  //     |   A   |
  //     |    +--+----+   =   +--+
  //     +----+--+    |       +--+
  //          |   B   |
  //          |       |
  //          +-------+
  // 
  intersect: function(csg) {
    var a = new CSG.Tree(this.polygons);
    var b = new CSG.Tree(csg.polygons);
    a.invert();
    b.clipTo(a);
    b.invert();
    a.clipTo(b);
    b.clipTo(a);
    a.addPolygons(b.allPolygons());
    a.invert();
    return CSG.fromPolygons(a.allPolygons());
  },

  // Return a new CSG solid with solid and empty space switched. This solid is
  // not modified.
  inverse: function() {
    var csg = this.clone();
    csg.polygons.map(function(p) { p.flip(); });
    return csg;
  },
  
  // Affine transformation of CSG object. Returns a new CSG object
  transform: function(matrix4x4) {
    var newpolygons = this.polygons.map(function(p) { return p.transform(matrix4x4); } );
    return CSG.fromPolygons(newpolygons);  
  },
  
  translate: function(v) {
    return this.transform(CSG.Matrix4x4.translation(v));
  },
  
  scale: function(f) {
    return this.transform(CSG.Matrix4x4.scaling(f));
  },
  
  rotateX: function(deg) {
    return this.transform(CSG.Matrix4x4.rotationX(deg));
  },
  
  rotateY: function(deg) {
    return this.transform(CSG.Matrix4x4.rotationY(deg));
  },
  
  rotateZ: function(deg) {
    return this.transform(CSG.Matrix4x4.rotationZ(deg));
  },
  
  toStlString: function() {
    var result="solid csg.js\n";
    this.polygons.map(function(p){ result += p.toStlString(); });
    result += "endsolid csg.js\n";
    return result;
  },
  
  toString: function() {
    var result = "";
    this.polygons.map(function(p){ result += p.toString(); });
    return result;
  },

};

// Construct an axis-aligned solid cuboid. Optional parameters are `center` and
// `radius`, which default to `[0, 0, 0]` and `[1, 1, 1]`. The radius can be
// specified using a single number or a list of three numbers, one for each axis.
// 
// Example code:
// 
//     var cube = CSG.cube({
//       center: [0, 0, 0],
//       radius: 1
//     });
CSG.cube = function(options) {
  options = options || {};
  var c = new CSG.Vector(options.center || [0, 0, 0]);
  var r = !options.radius ? [1, 1, 1] : options.radius.length ?
           options.radius : [options.radius, options.radius, options.radius];
  return CSG.fromPolygons([
    [[0, 4, 6, 2], [-1, 0, 0]],
    [[1, 3, 7, 5], [+1, 0, 0]],
    [[0, 1, 5, 4], [0, -1, 0]],
    [[2, 6, 7, 3], [0, +1, 0]],
    [[0, 2, 3, 1], [0, 0, -1]],
    [[4, 5, 7, 6], [0, 0, +1]]
  ].map(function(info) {
    return new CSG.Polygon(info[0].map(function(i) {
      var pos = new CSG.Vector(
        c.x + r[0] * (2 * !!(i & 1) - 1),
        c.y + r[1] * (2 * !!(i & 2) - 1),
        c.z + r[2] * (2 * !!(i & 4) - 1)
      );
      return new CSG.Vertex(pos, new CSG.Vector(info[1]));
    }));
  }));
};

// Construct a solid sphere. Optional parameters are `center`, `radius`,
// `slices`, and `stacks`, which default to `[0, 0, 0]`, `1`, `16`, and `8`.
// The `slices` and `stacks` parameters control the tessellation along the
// longitude and latitude directions.
// 
// Example usage:
// 
//     var sphere = CSG.sphere({
//       center: [0, 0, 0],
//       radius: 1,
//       slices: 16,
//       stacks: 8
//     });
CSG.sphere = function(options) {
  options = options || {};
  var c = new CSG.Vector(options.center || [0, 0, 0]);
  var r = options.radius || 1;
  var slices = options.slices || 16;
  var stacks = options.stacks || 8;
  var polygons = [], vertices;
  function vertex(theta, phi) {
    theta *= Math.PI * 2;
    phi *= Math.PI;
    var dir = new CSG.Vector(
      Math.cos(theta) * Math.sin(phi),
      Math.cos(phi),
      Math.sin(theta) * Math.sin(phi)
    );
    vertices.push(new CSG.Vertex(c.plus(dir.times(r)), dir));
  }
  for (var i = 0; i < slices; i++) {
    for (var j = 0; j < stacks; j++) {
      vertices = [];
      vertex(i / slices, j / stacks);
      if (j > 0) vertex((i + 1) / slices, j / stacks);
      if (j < stacks - 1) vertex((i + 1) / slices, (j + 1) / stacks);
      vertex(i / slices, (j + 1) / stacks);
      polygons.push(new CSG.Polygon(vertices));
    }
  }
  return CSG.fromPolygons(polygons);
};

// Construct a solid cylinder. Optional parameters are `start`, `end`,
// `radius`, and `slices`, which default to `[0, -1, 0]`, `[0, 1, 0]`, `1`, and
// `16`. The `slices` parameter controls the tessellation.
// 
// Example usage:
// 
//     var cylinder = CSG.cylinder({
//       start: [0, -1, 0],
//       end: [0, 1, 0],
//       radius: 1,
//       slices: 16
//     });
CSG.cylinder = function(options) {
  options = options || {};
  var s = new CSG.Vector(options.start || [0, -1, 0]);
  var e = new CSG.Vector(options.end || [0, 1, 0]);
  var ray = e.minus(s);
  var r = options.radius || 1;
  var slices = options.slices || 16;
  var axisZ = ray.unit(), isY = (Math.abs(axisZ.y) > 0.5);
  var axisX = new CSG.Vector(isY, !isY, 0).cross(axisZ).unit();
  var axisY = axisX.cross(axisZ).unit();
  var start = new CSG.Vertex(s, axisZ.negated());
  var end = new CSG.Vertex(e, axisZ.unit());
  var polygons = [];
  function point(stack, slice, normalBlend) {
    var angle = slice * Math.PI * 2;
    var out = axisX.times(Math.cos(angle)).plus(axisY.times(Math.sin(angle)));
    var pos = s.plus(ray.times(stack)).plus(out.times(r));
    var normal = out.times(1 - Math.abs(normalBlend)).plus(axisZ.times(normalBlend));
    return new CSG.Vertex(pos, normal);
  }
  for (var i = 0; i < slices; i++) {
    var t0 = i / slices, t1 = (i + 1) / slices;
    polygons.push(new CSG.Polygon([start, point(0, t0, -1), point(0, t1, -1)]));
    polygons.push(new CSG.Polygon([point(0, t1, 0), point(0, t0, 0), point(1, t0, 0), point(1, t1, 0)]));
    polygons.push(new CSG.Polygon([end, point(1, t1, 1), point(1, t0, 1)]));
  }
  return CSG.fromPolygons(polygons);
};

// # class Vector

// Represents a 3D vector.
// 
// Example usage:
// 
//     new CSG.Vector(1, 2, 3);
//     new CSG.Vector([1, 2, 3]);
//     new CSG.Vector({ x: 1, y: 2, z: 3 });

CSG.Vector = function(x, y, z) {
  if (arguments.length == 3) {
    this.x = x;
    this.y = y;
    this.z = z;
  } else if ('x' in x) {
    this.x = x.x;
    this.y = x.y;
    this.z = x.z;
  } else {
    this.x = x[0];
    this.y = x[1];
    this.z = x[2];
  }
};

CSG.Vector.prototype = {
  clone: function() {
    return new CSG.Vector(this.x, this.y, this.z);
  },

  negated: function() {
    return new CSG.Vector(-this.x, -this.y, -this.z);
  },

  plus: function(a) {
    return new CSG.Vector(this.x + a.x, this.y + a.y, this.z + a.z);
  },

  minus: function(a) {
    return new CSG.Vector(this.x - a.x, this.y - a.y, this.z - a.z);
  },

  times: function(a) {
    return new CSG.Vector(this.x * a, this.y * a, this.z * a);
  },

  dividedBy: function(a) {
    return new CSG.Vector(this.x / a, this.y / a, this.z / a);
  },

  dot: function(a) {
    return this.x * a.x + this.y * a.y + this.z * a.z;
  },

  lerp: function(a, t) {
    return this.plus(a.minus(this).times(t));
  },

  length: function() {
    return Math.sqrt(this.dot(this));
  },

  unit: function() {
    return this.dividedBy(this.length());
  },

  cross: function(a) {
    return new CSG.Vector(
      this.y * a.z - this.z * a.y,
      this.z * a.x - this.x * a.z,
      this.x * a.y - this.y * a.x
    );
  },

  equals: function(a) {
    return (this.x == a.x) && (this.y == a.y) && (this.z == a.z);
  },
  
  // Right multiply by a 4x4 matrix (the vector is interpreted as a row vector)
  // Returns a new CSG.Vector
  multiply4x4: function(matrix4x4) {
    return matrix4x4.rightMultiply1x3Vector(this);
  },
  
  toStlString: function() {
    return this.x+" "+this.y+" "+this.z;
  },
  
  toString: function() {
    return "("+this.x+", "+this.y+", "+this.z+")";
  },
  
};

// # class Vertex

// Represents a vertex of a polygon. Use your own vertex class instead of this
// one to provide additional features like texture coordinates and vertex
// colors. Custom vertex classes need to provide a `pos` property and `clone()`,
// `flip()`, and `interpolate()` methods that behave analogous to the ones
// defined by `CSG.Vertex`. This class provides `normal` so convenience
// functions like `CSG.sphere()` can return a smooth vertex normal, but `normal`
// is not used anywhere else.

CSG.Vertex = function(pos, normal) {
  this.pos = new CSG.Vector(pos);
  this.normal = new CSG.Vector(normal);
};

CSG.Vertex.prototype = {
  clone: function() {
    return new CSG.Vertex(this.pos.clone(), this.normal.clone());
  },

  // Invert all orientation-specific data (e.g. vertex normal). Called when the
  // orientation of a polygon is flipped.
  flip: function() {
    this.normal = this.normal.negated();
  },

  // Create a new vertex between this vertex and `other` by linearly
  // interpolating all properties using a parameter of `t`. Subclasses should
  // override this to interpolate additional properties.
  interpolate: function(other, t) {
    return new CSG.Vertex(
      this.pos.lerp(other.pos, t),
      this.normal.lerp(other.normal, t)
    );
  },
  
  // Affine transformation of vertex. Returns a new CSG.Vertex
  transform: function(matrix4x4) {
    var newpos = this.pos.multiply4x4(matrix4x4);
    var posPlusNormal = this.pos.plus(this.normal);
    var newPosPlusNormal = posPlusNormal.multiply4x4(matrix4x4);
    var newnormal = newPosPlusNormal.minus(newpos).unit();
    return new CSG.Vertex(newpos, newnormal);  
  },
  
  toStlString: function() {
    return "vertex "+this.pos.toStlString()+"\n";
  },
  
  toString: function() {
    return this.pos.toString();
  },
};

// # class Plane

// Represents a plane in 3D space.

CSG.Plane = function(normal, w) {
  this.normal = normal;
  this.w = w;
};

// `CSG.Plane.EPSILON` is the tolerance used by `splitPolygon()` to decide if a
// point is on the plane.
CSG.Plane.EPSILON = 1e-5;

CSG.Plane.fromPoints = function(a, b, c) {
  var n = b.minus(a).cross(c.minus(a)).unit();
  return new CSG.Plane(n, n.dot(a));
};

CSG.Plane.prototype = {
  clone: function() {
    return new CSG.Plane(this.normal.clone(), this.w);
  },

  flip: function() {
    this.normal = this.normal.negated();
    this.w = -this.w;
  },

  equals: function(n) {
    return this.normal.equals(n.normal) && this.w == n.w;
  },

  // Split `polygon` by this plane if needed, then put the polygon or polygon
  // fragments in the appropriate lists. Coplanar polygons go into either
  // `coplanarFront` or `coplanarBack` depending on their orientation with
  // respect to this plane. Polygons in front or in back of this plane go into
  // either `front` or `back`.
  splitPolygon: function(polygon, coplanarFront, coplanarBack, front, back) {
    var COPLANAR = 0;
    var FRONT = 1;
    var BACK = 2;
    var SPANNING = 3;

    // Classify each point as well as the entire polygon into one of the above
    // four classes.
    var polygonType = 0;
    var types = [];
    for (var i = 0; i < polygon.vertices.length; i++) {
      var t = this.normal.dot(polygon.vertices[i].pos) - this.w;
      var type = (t < -CSG.Plane.EPSILON) ? BACK : (t > CSG.Plane.EPSILON) ? FRONT : COPLANAR;
      polygonType |= type;
      types.push(type);
    }

    // Put the polygon in the correct list, splitting it when necessary.
    switch (polygonType) {
      case COPLANAR:
        (this.normal.dot(polygon.plane.normal) > 0 ? coplanarFront : coplanarBack).push(polygon);
        break;
      case FRONT:
        front.push(polygon);
        break;
      case BACK:
        back.push(polygon);
        break;
      case SPANNING:
        var f = [], b = [];
        for (var i = 0; i < polygon.vertices.length; i++) {
          var j = (i + 1) % polygon.vertices.length;
          var ti = types[i], tj = types[j];
          var vi = polygon.vertices[i], vj = polygon.vertices[j];
          if (ti != BACK) f.push(vi);
          if (ti != FRONT) b.push(ti != BACK ? vi.clone() : vi);
          if ((ti | tj) == SPANNING) {
            var t = (this.w - this.normal.dot(vi.pos)) / this.normal.dot(vj.pos.minus(vi.pos));
            var v = vi.interpolate(vj, t);
            f.push(v);
            b.push(v.clone());
          }
        }
        if (f.length >= 3) front.push(new CSG.Polygon(f, polygon.shared, polygon.plane.clone()));
        if (b.length >= 3) back.push(new CSG.Polygon(b, polygon.shared, polygon.plane.clone()));
        break;
    }
  },
  
  toString: function() {
    return "[normal: "+this.normal.toString()+", w: "+this.w+"]";
  },
};


// # class Polygon

// Represents a convex polygon. The vertices used to initialize a polygon must
// be coplanar and form a convex loop. They do not have to be `CSG.Vertex`
// instances but they must behave similarly (duck typing can be used for
// customization).
// 
// Each convex polygon has a `shared` property, which is shared between all
// polygons that are clones of each other or were split from the same polygon.
// This can be used to define per-polygon properties (such as surface color).
// 
// The plane of the polygon is calculated from the vertex coordinates
// To avoid unnecessary recalculation, the plane can alternatively be
// passed as the third argument 
CSG.Polygon = function(vertices, shared, plane) {
  this.vertices = vertices;
  var numvertices = vertices.length;

  this.shared = shared;
  if(arguments.length >= 3)
  {
    this.plane = plane;
  }
  else
  {
    this.plane = CSG.Plane.fromPoints(vertices[0].pos, vertices[1].pos, vertices[2].pos);
  }
  
  if(_CSGDEBUG)
  {
    this.checkIfConvex();
  }
};

CSG.Polygon.prototype = {
  clone: function() {
    var vertices = this.vertices.map(function(v) { return v.clone(); });
    return new CSG.Polygon(vertices, this.shared, this.plane.clone() );
  },
  
  // check whether the polygon is convex (it should be, otherwise we will get unexpected results)
  checkIfConvex: function() {
    if(! CSG.Polygon.verticesConvex(this.vertices, this.plane.normal))
    {
      throw new Error("Not convex!");
    }
  },

  flip: function() {
    this.vertices.reverse().map(function(v) { v.flip(); });
    this.plane.flip();
  },
  
  // Affine transformation of polygon. Returns a new CSG.Polygon
  transform: function(matrix4x4) {
    var newvertices = this.vertices.map(function(v) { return v.transform(matrix4x4); } );
    return new CSG.Polygon(newvertices, this.shared);
  },
  
  toStlString: function() {
    var result="";
    if(this.vertices.length >= 3) // should be!
    {
      // STL requires triangular polygons. If our polygon has more vertices, create
      // multiple triangles:
      var firstVertexStl = this.vertices[0].toStlString();
      for(var i=0; i < this.vertices.length-2; i++)
      {
        result += "facet normal "+this.plane.normal.toStlString()+"\nouter loop\n";
        result += firstVertexStl;
        result += this.vertices[i+1].toStlString();
        result += this.vertices[i+2].toStlString();
        result += "endloop\nendfacet\n";    
      } 
    }
    return result;
  },
  
  toString: function() {
    var result = "Polygon plane: "+this.plane.toString()+"\n";
    this.vertices.map(function(vertex) {
      result += "  "+vertex.toString()+"\n";
    });
    return result;
  },  
};

CSG.Polygon.verticesConvex = function(vertices, planenormal) {
  var numvertices = vertices.length;
  if(numvertices > 2)
  {
    var prevprevpos=vertices[numvertices-2].pos;
    var prevpos=vertices[numvertices-1].pos;
    for(var i=0; i < numvertices; i++)
    {
      var pos=vertices[i].pos;
      if(!CSG.Polygon.isConvexPoint(prevprevpos, prevpos, pos, planenormal))
      {
        return false;
      }
      prevprevpos=prevpos;
      prevpos=pos;
    }
  }
  return true;
};

CSG.Polygon.verticesStrictlyConvex = function(vertices, planenormal) {
  var numvertices = vertices.length;
  if(numvertices > 2)
  {
    var prevprevpos=vertices[numvertices-2].pos;
    var prevpos=vertices[numvertices-1].pos;
    for(var i=0; i < numvertices; i++)
    {
      var pos=vertices[i].pos;
      if(!CSG.Polygon.isStrictlyConvexPoint(prevprevpos, prevpos, pos, planenormal))
      {
        return false;
      }
      prevprevpos=prevpos;
      prevpos=pos;
    }
  }
  return true;
};

// Create a polygon from the given points
CSG.Polygon.createFromPoints = function(points, shared) {
  // initially set a dummy vertex normal:
  var dummynormal = new CSG.Vector(0, 0, 0);
  var vertices = [];
  points.map( function(p) {
    var vec = new CSG.Vector(p);
    var vertex = new CSG.Vertex(vec, dummynormal);
    vertices.push(vertex); 
  });            
  var polygon = new CSG.Polygon(vertices, shared);
  // now set the vertex normals to the polygon normal:
  var normal = polygon.plane.normal;
  polygon.vertices.map( function(v) {
    v.normal = normal;
  });
  return polygon;
};

// calculate whether three points form a convex corner 
//  prevpoint, point, nextpoint: the 3 coordinates (CSG.Vector instances)
//  normal: the normal vector of the plane
CSG.Polygon.isConvexPoint = function(prevpoint, point, nextpoint, normal) {
  var crossproduct=point.minus(prevpoint).cross(nextpoint.minus(point));
  var crossdotnormal=crossproduct.dot(normal);
  return (crossdotnormal >= 0);
};

CSG.Polygon.isStrictlyConvexPoint = function(prevpoint, point, nextpoint, normal) {
  var crossproduct=point.minus(prevpoint).cross(nextpoint.minus(point));
  var crossdotnormal=crossproduct.dot(normal);
  return (crossdotnormal >= 1e-5);
};

// # class PolygonTreeNode

// This class manages hierarchical splits of polygons
// At the top is a root node which doesn hold a polygon, only child PolygonTreeNodes
// Below that are zero or more 'top' nodes; each holds a polygon. The polygons can be in different planes 
// splitByPlane() splits a node by a plane. If the plane intersects the polygon, two new child nodes
// are created holding the splitted polygon.
// getPolygons() retrieves the polygon from the tree. If for PolygonTreeNode the polygon is split but 
// the two split parts (child nodes) are still intact, then the unsplit polygon is returned.
// This ensures that we can safely split a polygon into many fragments. If the fragments are untouched,
//  getPolygons() will return the original unsplit polygon instead of the fragments.
// remove() removes a polygon from the tree. Once a polygon is removed, the parent polygons are invalidated 
// since they are no longer intact. 

// constructor creates the root node:
CSG.PolygonTreeNode = function() {
  this.parent = null;
  this.children = [];
  this.polygon = null;
};

CSG.PolygonTreeNode.prototype = {
  // fill the tree with polygons. Should be called on the root node only; child nodes must
  // always be a derivate (split) of the parent node.
  addPolygons: function(polygons) {
    if(!this.isRootNode()) throw new Error("Assertion failed");  // new polygons can only be added to root node; children can only be splitted polygons
    var _this = this;
    polygons.map(function(polygon) {
      _this.addChild(polygon);
    });
  },
  
  // remove a node
  // - the siblings become toplevel nodes
  // - the parent is removed recursively
  remove: function() {
    if(this.isRootNode()) throw new Error("Assertion failed");  // can't remove root node
    if(this.children.length) throw new Error("Assertion failed"); // we shouldn't remove nodes with children
    
    // remove ourselves from the parent's children list:
    var parentschildren = this.parent.children;
    var i = parentschildren.indexOf(this);
    if(i < 0) throw new Error("Assertion failed");
    parentschildren.splice(i,1);
    
    // invalidate the parent's polygon, and of all parents above it:
    this.parent.recursivelyInvalidatePolygon();
  },

  isRootNode: function() {
    return !this.parent;
  },  

  // invert all polygons in the tree. Call on the root node
  invert: function() {
    if(!this.isRootNode()) throw new Error("Assertion failed");  // can only call this on the root node
    this.invertSub();
  },

  getPolygon: function () {
    if(!this.polygon) throw new Error("Assertion failed");  // doesn't have a polygon, which means that it has been broken down
    return this.polygon;
  },

  getPolygons: function (result) {
    if(this.polygon)
    {
      // the polygon hasn't been broken yet. We can ignore the children and return our polygon:
      result.push(this.polygon);
    }
    else
    {
      // our polygon has been split up and broken, so gather all subpolygons from the children:
      var childpolygons = [];
      this.children.map(function(child) {
        child.getPolygons(childpolygons);
      });
//      if( (this.parent) && (childpolygons.length == 2) )
//      {
//        var joinedpolygon = CSG.PolygonTreeNode.tryJoiningPolygons(childpolygons[0], childpolygons[1]);
//        if(joinedpolygon)
//        {
//          childpolygons=[joinedpolygon];
//        }
//      }

//      if(this.parent && (!this.parent.parent) )
//      {
//        // it's a top node:
//        CSG.PolygonTreeNode.tryJoiningManyPolygons(childpolygons);
//      }
      childpolygons.map(function(p) {
        result.push(p);
      });
    }
  },

  // split the node by a plane; add the resulting nodes to the frontnodes and backnodes array  
  // If the plane doesn't intersect the polygon, the 'this' object is added to one of the arrays
  // If the plane does intersect the polygon, two new child nodes are created for the front and back fragments,
  //  and added to both arrays. 
  splitByPlane: function(plane, coplanarfrontnodes, coplanarbacknodes, frontnodes, backnodes) {
    if(this.children.length > 0)
    {
      // if we have children, split the children
      this.children.map(function(child) {
        child.splitByPlane(plane, coplanarfrontnodes, coplanarbacknodes, frontnodes, backnodes);
      });
    }
    else
    {
      // no children. Split the polygon:
      if(this.polygon)
      {
        var coplanarfrontpolygons = [], coplanarbackpolygons = [];
        var frontpolygons = [], backpolygons = [];
        plane.splitPolygon(this.polygon, coplanarfrontpolygons, coplanarbackpolygons, frontpolygons, backpolygons);
        var numpolygons=coplanarfrontpolygons.length + coplanarbackpolygons.length + frontpolygons.length + backpolygons.length; 
        if(numpolygons > 2) throw new Error("Assertion failed");
        if(numpolygons == 1)
        {
          // the polygon was not split:
          if(frontpolygons.length == 1)
          {
            frontnodes.push(this);
          }
          else if(backpolygons.length == 1)
          {
            backnodes.push(this);
          }
          else if(coplanarfrontpolygons.length == 1)
          {
            coplanarfrontnodes.push(this);
          }
          else // if(coplanarbackpolygons.length == 1)
          {
            coplanarbacknodes.push(this);
          }
        }
        else if(numpolygons == 2)
        {
          // the polygon was split. Split our node by creating two children:
          if(frontpolygons.length != 1) throw new Error("Assertion failed");
          if(backpolygons.length != 1) throw new Error("Assertion failed");
          var frontnode = this.addChild(frontpolygons[0].clone());
          var backnode = this.addChild(backpolygons[0].clone());
          frontnodes.push(frontnode);
          backnodes.push(backnode);          
        }
      }
    }
  },
  
 
  // PRIVATE methods from here:

  // add child to a node
  // this should be called whenever the polygon is split
  // a child should be created for every fragment of the split polygon 
  // returns the newly created child
  addChild: function(polygon) {
    var newchild = new CSG.PolygonTreeNode();
    newchild.parent = this;
    newchild.polygon = polygon;
    this.children.push(newchild);
    return newchild;
  },

  invertSub: function() {
    if(this.polygon)
    {
      this.polygon.flip();
    }
    this.children.map(function(child) {
      child.invertSub();
    });
  },
  
  recursivelyInvalidatePolygon: function() {
    this.polygon = null;
    if(this.parent)
    {
      this.parent.recursivelyInvalidatePolygon();
    }
  },
  
};


// Experimental and terribly inefficient:
// Try joining polygons into convex polygons
// polygons: array of CSG.Polygon; will be modified if polygons are merged 
CSG.PolygonTreeNode.tryJoiningManyPolygons = function(polygons) {
  var numpolygons=polygons.length;
  if(numpolygons >= 2)
  {
    var i1=0;
    while(true)
    {
      if(i1 >= numpolygons-1) break;
      var i2=i1+1;
      while(true)
      {
        if(i2 >= numpolygons) break;
        var joined = CSG.PolygonTreeNode.tryJoiningPolygons(polygons[i1], polygons[i2]); 
        if(joined)
        {
          polygons.splice(i1,1,joined);
          polygons.splice(i2,1);
          --numpolygons;
          i2 = i1+1;
          continue;
        }
        ++i2;
      }
      ++i1;
    }
  }
};

// Experimental and terribly inefficient:
// Try joining two polygons into a convex polygon.
// If succesful returns the joined polygon (a  CSG.Polygon), null otherwise
// polygon1, polygon2: CSG.Polygon instances 
CSG.PolygonTreeNode.tryJoiningPolygons = function(polygon1, polygon2) {
  var newvertices = null;
  if(!polygon1.plane.equals(polygon2.plane)) throw new Error("Assertion failed");
  var numpoints1 = polygon1.vertices.length;
  var numpoints2 = polygon2.vertices.length;
  
  if( (numpoints1 >= 3) && (numpoints2 >= 3) )
  {
    var prevpoint1=polygon1.vertices[numpoints1-1].pos;    
    for(var p1=0; p1 < numpoints1; p1++)
    {
      var point1=polygon1.vertices[p1].pos;
      var prevpoint2=polygon2.vertices[numpoints2-1].pos;    
      for(var p2=0; p2 < numpoints2; p2++)
      {
        var point2=polygon2.vertices[p2].pos;
        if(point2.equals(prevpoint1) && point1.equals(prevpoint2))
        {
          // we have matching vertices!
          newvertices = [];
          var i=p1;
          var count=numpoints1 - 1;
          while(count-- > 0)
          {
            newvertices.push(polygon1.vertices[i].clone());
            ++i;
            if(i >= numpoints1) i=0;
          }
          i = p2;
          count = numpoints2 - 1; 
          while(count-- > 0)
          {
            newvertices.push(polygon2.vertices[i].clone());
            ++i;
            if(i >= numpoints2) i=0;
          }
        } 
      }
      prevpoint1 == point1;    
    }
  }
  var result = null;
  if(newvertices)
  {
    // we have succesfully joined the polygons, but we need to check if the resulting polygon is convex:
    var isconvex = CSG.Polygon.verticesStrictlyConvex(newvertices, polygon1.plane.normal);
    if(isconvex)
    {
      result = new CSG.Polygon(newvertices, polygon1.shared, polygon1.plane.clone()); 
    }
  }
  return result;
}


// # class Tree
// This is the root of a BSP tree
// We are using this separate class for the root of the tree, to hold the PolygonTreeNode root
// The actual tree is kept in this.rootnode
CSG.Tree = function(polygons) {
  this.polygonTree = new CSG.PolygonTreeNode();
  this.rootnode = new CSG.Node();
  if (polygons) this.addPolygons(polygons);
};

CSG.Tree.prototype = {
  invert: function() {
    this.polygonTree.invert();
    this.rootnode.invert();
  },
  
  // Remove all polygons in this BSP tree that are inside the other BSP tree
  // `tree`.
  clipTo: function(tree) {
    this.rootnode.clipTo(tree);
  },

  allPolygons: function() {
    var result = [];
    this.polygonTree.getPolygons(result);
    return result;
  },

  addPolygons: function(polygons) {
    var _this = this;
    polygons.map(function(p) {
      _this.addPolygon(p.clone());
    });
  },

  addPolygon: function(polygon) {
    var polygontreenode=this.polygonTree.addChild(polygon);
    this.rootnode.addPolygonTreeNode(polygontreenode);
  },  
};

// # class Node

// Holds a node in a BSP tree. A BSP tree is built from a collection of polygons
// by picking a polygon to split along.
// Polygons are not stored directly in the tree, but in PolygonTreeNodes, stored in
// this.polygontreenodes. Those PolygonTreeNodes are children of the owning
// CSG.Tree.polygonTree
// This is not a leafy BSP tree since there is
// no distinction between internal and leaf nodes.

CSG.Node = function() {
  this.plane = null;
  this.front = null;
  this.back = null;
  this.polygontreenodes = [];
};

CSG.Node.prototype = {
  // Convert solid space to empty space and empty space to solid space.
  invert: function() {
    this.plane.flip();
    if (this.front) this.front.invert();
    if (this.back) this.back.invert();
    var temp = this.front;
    this.front = this.back;
    this.back = temp;
  },

  // clip polygontreenodes to our plane
  // returns a new array of PolygonTreeNodes with the clipped polygons
  // also calls remove() for all clipped PolygonTreeNodes, so the polygon tree is modified!
  clipPolygons: function(polygontreenodes) {
    var clippednodes;
    if(this.plane)
    {
      var backnodes = [];
      var frontnodes = [];
      var plane = this.plane;
      polygontreenodes.map(function(node) {
        node.splitByPlane(plane, frontnodes, backnodes, frontnodes, backnodes);
      });
      var clippedfrontnodes = [];
      var clippedbacknodes = [];
      if(this.front)
      {
        clippedfrontnodes = this.front.clipPolygons(frontnodes);
      }
      else
      {
        clippedfrontnodes = frontnodes;
      }
      if(this.back)
      {
        clippedbacknodes = this.back.clipPolygons(backnodes);
      }
      else
      {
        // there's nothing behind this plane. Delete the nodes behind this plane:
        backnodes.map( function(node) {
          node.remove();
        });
      }
      clippednodes = clippedfrontnodes.concat(clippedbacknodes);
    }
    else
    {
      clippednodes=polygontreenodes.slice();
    }
    return clippednodes;
  },

  // Remove all polygons in this BSP tree that are inside the other BSP tree
  // `tree`.
  clipTo: function(tree) {
    var origpolygontreenodes = this.polygontreenodes;
    if(origpolygontreenodes.length > 0)
    {
      var clippedtreenodes = tree.rootnode.clipPolygons(origpolygontreenodes);
      this.polygontreenodes = clippedtreenodes;
    }
    if (this.front) this.front.clipTo(tree);
    if (this.back) this.back.clipTo(tree);
  },
  
  addPolygonTreeNode: function(polygontreenode) {
    if(!this.plane)
    {
      this.plane = polygontreenode.getPolygon().plane.clone();
    }
    var frontnodes = [];
    var backnodes = [];
    polygontreenode.splitByPlane(this.plane, this.polygontreenodes, this.polygontreenodes, frontnodes, backnodes);
    if(frontnodes.length > 0)
    {
      if (!this.front) this.front = new CSG.Node();
      this.front.addPolygonTreeNode(frontnodes[0]);
    }
    if(backnodes.length > 0)
    {
      if (!this.back) this.back = new CSG.Node();
      this.back.addPolygonTreeNode(backnodes[0]);
    }
  },
};

//////////

// # class Matrix4x4:
// Represents a 4x4 matrix. Elements are specified in row order
CSG.Matrix4x4 = function(elements) {
  if (arguments.length >= 1) {
    this.elements=elements;
  }
  else
  {
    // if no arguments passed: create unity matrix  
    this.elements=[1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1];
  }
}

CSG.Matrix4x4.prototype = {
  plus: function(m) {
    var r=[];
    for(var i=0; i < 16; i++)
    {
      r[i]=this.elements[i]+m.elements[i];
    }
    return new CSG.Matrix4x4(r);
  },
  
  minus: function(m) {
    var r=[];
    for(var i=0; i < 16; i++)
    {
      r[i]=this.elements[i]-m.elements[i];
    }
    return new CSG.Matrix4x4(r);
  },

  // right multiply by another 4x4 matrix:
  multiply: function(m) {
    // cache elements in local variables, for speedup:
    var this0=this.elements[0];
    var this1=this.elements[1];
    var this2=this.elements[2];
    var this3=this.elements[3];
    var this4=this.elements[4];
    var this5=this.elements[5];
    var this6=this.elements[6];
    var this7=this.elements[7];
    var this8=this.elements[8];
    var this9=this.elements[9];
    var this10=this.elements[10];
    var this11=this.elements[11];
    var this12=this.elements[12];
    var this13=this.elements[13];
    var this14=this.elements[14];
    var this15=this.elements[15];
    var m0=m.elements[0];
    var m1=m.elements[1];
    var m2=m.elements[2];
    var m3=m.elements[3];
    var m4=m.elements[4];
    var m5=m.elements[5];
    var m6=m.elements[6];
    var m7=m.elements[7];
    var m8=m.elements[8];
    var m9=m.elements[9];
    var m10=m.elements[10];
    var m11=m.elements[11];
    var m12=m.elements[12];
    var m13=m.elements[13];
    var m14=m.elements[14];
    var m15=m.elements[15];
    
    var result=[];
    result[0] = this0*m0 + this1*m4 + this2*m8 + this3*m12;
    result[1] = this0*m1 + this1*m5 + this2*m9 + this3*m13;
    result[2] = this0*m2 + this1*m6 + this2*m10 + this3*m14;
    result[3] = this0*m3 + this1*m7 + this2*m11 + this3*m15;
    result[4] = this4*m0 + this5*m4 + this6*m8 + this7*m12;
    result[5] = this4*m1 + this5*m5 + this6*m9 + this7*m13;
    result[6] = this4*m2 + this5*m6 + this6*m10 + this7*m14;
    result[7] = this4*m3 + this5*m7 + this6*m11 + this7*m15;
    result[8] = this8*m0 + this9*m4 + this10*m8 + this11*m12;
    result[9] = this8*m1 + this9*m5 + this10*m9 + this11*m13;
    result[10] = this8*m2 + this9*m6 + this10*m10 + this11*m14;
    result[11] = this8*m3 + this9*m7 + this10*m11 + this11*m15;
    result[12] = this12*m0 + this13*m4 + this14*m8 + this15*m12;
    result[13] = this12*m1 + this13*m5 + this14*m9 + this15*m13;
    result[14] = this12*m2 + this13*m6 + this14*m10 + this15*m14;
    result[15] = this12*m3 + this13*m7 + this14*m11 + this15*m15;
    return new CSG.Matrix4x4(result);
  },
  
  clone: function() {
    var elements = this.elements.map(function(p) { return p; }); 
    return new CSG.Matrix4x4(elements);
  },
  
  // Multiply a CSG.Vector (interpreted as 1 row, 3 column) by this matrix 
  // Fourth element is taken as 1
  rightMultiply1x3Vector: function(v) {
    var v0 = v.x;
    var v1 = v.y;
    var v2 = v.z;
    var v3 = 1;    
    var x = v0*this.elements[0] + v1*this.elements[1] + v2*this.elements[2] + v3*this.elements[3];    
    var y = v0*this.elements[4] + v1*this.elements[5] + v2*this.elements[6] + v3*this.elements[7];    
    var z = v0*this.elements[8] + v1*this.elements[9] + v2*this.elements[10] + v3*this.elements[11];    
    var w = v0*this.elements[12] + v1*this.elements[13] + v2*this.elements[14] + v3*this.elements[15];
    // scale such that fourth element becomes 1:
    if(w != 1)
    {
      var invw=1.0/w;
      x *= invw;
      y *= invw;
      z *= invw;
    }
    return new CSG.Vector(x,y,z);       
  },
  
  // Multiply a CSG.Vector2D (interpreted as 1 row, 2 column) by this matrix 
  // Fourth element is taken as 1
  rightMultiply1x2Vector: function(v) {
    var v0 = v.x;
    var v1 = v.y;
    var v2 = 0;
    var v3 = 1;    
    var x = v0*this.elements[0] + v1*this.elements[1] + v2*this.elements[2] + v3*this.elements[3];    
    var y = v0*this.elements[4] + v1*this.elements[5] + v2*this.elements[6] + v3*this.elements[7];    
    var z = v0*this.elements[8] + v1*this.elements[9] + v2*this.elements[10] + v3*this.elements[11];    
    var w = v0*this.elements[12] + v1*this.elements[13] + v2*this.elements[14] + v3*this.elements[15];
    // scale such that fourth element becomes 1:
    if(w != 1)
    {
      var invw=1.0/w;
      x *= invw;
      y *= invw;
      z *= invw;
    }
    return new CSG.Vector2D(x,y);       
  },
};

// return the unity matrix
CSG.Matrix4x4.unity = function() {
  return new CSG.Matrix4x4(); 
};

// Create a rotation matrix for rotating around the x axis
CSG.Matrix4x4.rotationX = function(degrees) {
  var radians = degrees * Math.PI * (1.0/180.0);
  var cos = Math.cos(radians);
  var sin = Math.sin(radians);
  var els = [
    1, 0, 0, 0,
    0, cos, -sin, 0,
    0, sin, cos, 0,
    0, 0, 0, 1
  ];
  return new CSG.Matrix4x4(els);
};

// Create a rotation matrix for rotating around the y axis
CSG.Matrix4x4.rotationY = function(degrees) {
  var radians = degrees * Math.PI * (1.0/180.0);
  var cos = Math.cos(radians);
  var sin = Math.sin(radians);
  var els = [
    cos, 0, sin, 0,
    0, 1, 0, 0,
    -sin, 0, cos, 0,
    0, 0, 0, 1
  ];
  return new CSG.Matrix4x4(els);
};

// Create a rotation matrix for rotating around the z axis
CSG.Matrix4x4.rotationZ = function(degrees) {
  var radians = degrees * Math.PI * (1.0/180.0);
  var cos = Math.cos(radians);
  var sin = Math.sin(radians);
  var els = [
    cos, -sin, 0, 0,
    sin, cos, 0, 0,
    0, 0, 1, 0,
    0, 0, 0, 1
  ];
  return new CSG.Matrix4x4(els);
};

// Create an affine matrix for translation:
CSG.Matrix4x4.translation = function(v) {
  // parse as CSG.Vector, so we can pass an array or a CSG.Vector
  var vec = new CSG.Vector(v);
  var els = [
    1, 0, 0, vec.x,
    0, 1, 0, vec.y,
    0, 0, 1, vec.z,
    0, 0, 0, 1
  ];
  return new CSG.Matrix4x4(els);
};

// Create an affine matrix for scaling:
CSG.Matrix4x4.scaling = function(v) {
  // parse as CSG.Vector, so we can pass an array or a CSG.Vector
  var vec = new CSG.Vector(v);
  var els = [
    vec.x, 0, 0, 0,
    0, vec.y, 0, 0,
    0, 0, vec.z, 0,
    0, 0, 0, 1
  ];
  return new CSG.Matrix4x4(els);
};

///////////////////////////////////////////////////

// # class Vector2D:
// Represents a 2 element vector
CSG.Vector2D = function(x, y) {
  if (arguments.length == 2) {
    this.x = x;
    this.y = y;
  } else if ('x' in x) {
    this.x = x.x;
    this.y = x.y;
  } else {
    this.x = x[0];
    this.y = x[1];
  }
};

CSG.Vector2D.prototype = {
  // extend to a 3D vector by adding a z coordinate:
  toVector3D: function(z) {
    return new CSG.Vector(this.x, this.y, z);
  },
  
  clone: function() {
    return new CSG.Vector(this.x, this.y);
  },

  negated: function() {
    return new CSG.Vector(-this.x, -this.y);
  },

  plus: function(a) {
    return new CSG.Vector(this.x + a.x, this.y + a.y);
  },

  minus: function(a) {
    return new CSG.Vector(this.x - a.x, this.y - a.y);
  },

  times: function(a) {
    return new CSG.Vector(this.x * a, this.y * a);
  },

  dividedBy: function(a) {
    return new CSG.Vector(this.x / a, this.y / a);
  },

  dot: function(a) {
    return this.x * a.x + this.y * a.y;
  },

  lerp: function(a, t) {
    return this.plus(a.minus(this).times(t));
  },

  length: function() {
    return Math.sqrt(this.dot(this));
  },

  unit: function() {
    return this.dividedBy(this.length());
  },

  // Right multiply by a 4x4 matrix (the vector is interpreted as a row vector)
  // Returns a new CSG.Vector2D
  multiply4x4: function(matrix4x4) {
    return matrix4x4.rightMultiply1x2Vector(this);
  },
};

// A polygon in 2D space:
CSG.Polygon2D = function(points, shared) {
  var vectors = [];
  if(arguments.length >= 1) {
    points.map( function(p) {
      vectors.push(new CSG.Vector2D(p) );
    });    
  }
  this.points = vectors;
  this.shared = shared;
};

CSG.Polygon2D.prototype = {
  // Matrix transformation of polygon. Returns a new CSG.Polygon2D
  transform: function(matrix4x4) {
    var newpoints = this.points.map(function(p) { return p.multiply4x4(matrix4x4); } );
    return new CSG.Polygon2D(newpoints, this.shared);
  },
  
  translate: function(v) {
    v=new CSG.Vector2D(v);
    return this.transform(CSG.Matrix4x4.translation(v.toVector3D(0)));
  },
  
  scale: function(f) {
    f=new CSG.Vector2D(f);
    return this.transform(CSG.Matrix4x4.scaling(f.toVector3D(1)));
  },
  
  rotate: function(deg) {
    return this.transform(CSG.Matrix4x4.rotationZ(deg));
  },    
  
  // convert into a CSG.Polygon; set z coordinate to the given value
  toPolygon3D: function(z) {
    var points3d=[];
    this.points.map( function(p) {
      var vec3d = p.toVector3D(z);      
      points3d.push(vec3d);
    });
    var polygon = CSG.Polygon.createFromPoints(points3d, this.shared);
    polygon.checkIfConvex();
    return polygon;
  },
  
  // extruded=shape2d.extrude({offset: [0,0,10], twistangle: 360, twiststeps: 100});
  // linear extrusion of 2D polygon, with optional twist
  // The 2d polygon is placed in in z=0 plane and extruded into direction <offset> (a CSG.Vector)
  // The final face is rotated <twistangle> degrees. Rotation is done around the origin of the 2d shape (i.e. x=0, y=0)
  // twiststeps determines the resolution of the twist (should be >= 1)  
  // returns a CSG object
  extrude: function(params) {
    // parse parameters:
    if(!params) params={};
    var offsetvector;
    if("offset" in params)
    {
      offsetvector = new CSG.Vector(params.offset); // reparse as a CSG.Vector
    }
    else
    {
      offsetvector = new CSG.Vector(0,0,1);
    }
    
    var twistangle=0;
    if("twistangle" in params)
    {
      twistangle = params.twistangle;
    }
    
    var twiststeps=10;
    if("twiststeps" in params)
    {
      twiststeps = params.twiststeps;
    }
    if(twistangle == 0) twiststeps=1;
    if(twiststeps < 1) twiststeps=1;

    // create the polygons:        
    var newpolygons = [];
    
    // bottom face polygon:
    var bottomfacepolygon = this.toPolygon3D(0);
    var direction = bottomfacepolygon.plane.normal.dot(offsetvector);
    if(direction > 0)
    {
      bottomfacepolygon.flip();
    }
    newpolygons.push(bottomfacepolygon);
    
    var getTwistedPolygon = function(twiststep) {
      var fraction = (twiststep + 1) / twiststeps;
      var rotation = twistangle * fraction;
      var offset = offsetvector.times(fraction);
      var transformmatrix = CSG.Matrix4x4.rotationZ(rotation).multiply( CSG.Matrix4x4.translation(offset) );
      var polygon = bottomfacepolygon.transform(transformmatrix);      
      return polygon;
    };

    // create the side face polygons:
    var numvertices = bottomfacepolygon.vertices.length;
    var prevlevelpolygon = bottomfacepolygon;
    for(var twiststep=0; twiststep < twiststeps; ++twiststep)
    {
      var levelpolygon = getTwistedPolygon(twiststep);
      for(var i=0; i < numvertices; i++)
      {
        var sidefacepoints = [];
        var nexti = (i < (numvertices-1))? i+1:0;
        sidefacepoints.push(prevlevelpolygon.vertices[i].pos);
        sidefacepoints.push(levelpolygon.vertices[i].pos);
        sidefacepoints.push(levelpolygon.vertices[nexti].pos);
        sidefacepoints.push(prevlevelpolygon.vertices[nexti].pos);
        var sidefacepolygon=CSG.Polygon.createFromPoints(sidefacepoints, this.shared);
        newpolygons.push(sidefacepolygon);
      }
      if(twiststep == (twiststeps -1) )
      {
        // last level; add the top face polygon:
        levelpolygon.flip(); // flip so that the normal points outwards
        newpolygons.push(levelpolygon);
      }
      prevlevelpolygon = levelpolygon;
    }

    return CSG.fromPolygons(newpolygons);
  }
};

CSG.Polygon.prototype.extrude = function(offsetvector) {
  var newpolygons = [];

  var polygon1=this.clone();
  var direction = polygon1.plane.normal.dot(offsetvector);
  if(direction > 0)
  {
    polygon1.flip();
  }
  newpolygons.push(polygon1);
  var polygon2=polygon1.translate(offsetvector);
  var numvertices=this.vertices.length;
  for(var i=0; i < numvertices; i++)
  {
    var sidefacepoints = [];
    var nexti = (i < (numvertices-1))? i+1:0;
    sidefacepoints.push(polygon1.vertices[i].pos);
    sidefacepoints.push(polygon2.vertices[i].pos);
    sidefacepoints.push(polygon2.vertices[nexti].pos);
    sidefacepoints.push(polygon1.vertices[nexti].pos);
    var sidefacepolygon=CSG.Polygon.createFromPoints(sidefacepoints);
    newpolygons.push(sidefacepolygon);
  }
  polygon2.flip();
  newpolygons.push(polygon2);
  return CSG.fromPolygons(newpolygons);
};

CSG.Polygon.prototype.translate = function(offset) {
  return this.transform(CSG.Matrix4x4.translation(offset));
};

// Expand the polygon with a certain radius
// This extrudes the face of the polygon and adds rounded corners 
// Returns a CSG object (not a polygon anymore!)
// resolution: number of points per 360 degree for the rounded corners
CSG.Polygon.prototype.expand = function(radius, resolution) {
  resolution=resolution || 8;
  var result=new CSG();
  
  // expand each side of the polygon. The expansion of a line is a cylinder with
  // two spheres at the end:
  var numvertices=this.vertices.length;
  for(var i=0; i < numvertices; i++)
  {
    var previ = (i == 0) ? (numvertices-1):i-1;
    var p1 = this.vertices[previ].pos;
    var p2 = this.vertices[i].pos;
    var cylinder=CSG.cylinder({start: p1, end: p2, radius: radius, slices: resolution});
    result = result.union(cylinder);
    var sphere = CSG.sphere({center: p1, radius: radius, slices: resolution, stacks: resolution});
    result = result.union(sphere);
  }
  var extrudevector=this.plane.normal.unit().times(2*radius);
  var translatedpolygon = this.translate(extrudevector.times(-0.5));
  var extrudedface = translatedpolygon.extrude(extrudevector);  
  result=result.union(extrudedface);
  return result;
};

// Expand the solid
// resolution: number of points per 360 degree for the rounded corners
CSG.prototype.expand = function(radius, resolution) {
  var result=this.clone();
  this.polygons.map(function(p) {
    var expanded=p.expand(radius, resolution);
    result=result.union(expanded);
  });
  return result;
};

// Contract the solid
// resolution: number of points per 360 degree for the rounded corners
CSG.prototype.contract = function(radius, resolution) {
  var result=this.clone();
  this.polygons.map(function(p) {
    var expanded=p.expand(radius, resolution);
    result=result.subtract(expanded);
  });
  return result;
};

CSG.roundedCube = function(cuberadius, roundradius, resolution) {
  resolution = resolution || 8;
  cuberadius=new CSG.Vector(cuberadius);
  var innercuberadius=cuberadius.clone();
  innercuberadius.x -= roundradius;
  innercuberadius.y -= roundradius;
  innercuberadius.z -= roundradius;
  var result = CSG.cube({radius: [cuberadius.x, innercuberadius.y, innercuberadius.z]});
  result = result.union( CSG.cube({radius: [innercuberadius.x, cuberadius.y, innercuberadius.z]}));
  result = result.union( CSG.cube({radius: [innercuberadius.x, innercuberadius.y, cuberadius.z]}));
  for(var level=0; level < 2; level++)
  {
    var z = innercuberadius.z;
    if(level == 1) z = -z;
    var p1 = new CSG.Vector(innercuberadius.x, innercuberadius.y, z);
    var p2 = new CSG.Vector(innercuberadius.x, -innercuberadius.y, z);
    var p3 = new CSG.Vector(-innercuberadius.x, -innercuberadius.y, z);
    var p4 = new CSG.Vector(-innercuberadius.x, innercuberadius.y, z);
    var sphere = CSG.sphere({center: p1, radius: roundradius, slices: resolution, stacks: resolution});
    result = result.union(sphere);
    sphere = CSG.sphere({center: p2, radius: roundradius, slices: resolution, stacks: resolution});
    result = result.union(sphere);
    sphere = CSG.sphere({center: p3, radius: roundradius, slices: resolution, stacks: resolution});
    result = result.union(sphere);
    sphere = CSG.sphere({center: p4, radius: roundradius, slices: resolution, stacks: resolution});
    result = result.union(sphere);
    var cylinder = CSG.cylinder({start:p1, end: p2, radius: roundradius, slices: resolution});
    result = result.union(cylinder);
    cylinder = CSG.cylinder({start:p2, end: p3, radius: roundradius, slices: resolution});
    result = result.union(cylinder);
    cylinder = CSG.cylinder({start:p3, end: p4, radius: roundradius, slices: resolution});
    result = result.union(cylinder);
    cylinder = CSG.cylinder({start:p4, end: p1, radius: roundradius, slices: resolution});
    result = result.union(cylinder);
    if(level == 0) {
      var d = new CSG.Vector(0, 0, -2*z);
      cylinder = CSG.cylinder({start:p1, end: p1.plus(d), radius: roundradius, slices: resolution});
      result = result.union(cylinder);
      cylinder = CSG.cylinder({start:p2, end: p2.plus(d), radius: roundradius, slices: resolution});
      result = result.union(cylinder);
      cylinder = CSG.cylinder({start:p3, end: p3.plus(d), radius: roundradius, slices: resolution});
      result = result.union(cylinder);
      cylinder = CSG.cylinder({start:p4, end: p4.plus(d), radius: roundradius, slices: resolution});
      result = result.union(cylinder);
    }
  }
  return result;
}

