var angleX = 20;
var angleY = 20;
var shader;
var meshes = [];

function setup() {
  function setPolygonColor(solid, color) {
    solid.toPolygons().map(function(polygon) {
      polygon.shared = color;
    });
  }

  // Test coplanar cases
  var a = CSG.cube();
  var b = CSG.cylinder({ slices: 8, start: new CSG.Vector(-1, 0, 0), end: new CSG.Vector(1, 0, 0) });
  b = b.union(CSG.cylinder({ slices: 8, start: new CSG.Vector(0, -1, 0), end: new CSG.Vector(0, 1, 0) }));
  b = b.union(CSG.cylinder({ slices: 8, start: new CSG.Vector(0, 0, -1), end: new CSG.Vector(0, 0, 1) }));
  var c = CSG.cube({ radius: 0.5 });
  var d = CSG.cylinder({ radius: 0.5, slices: 8, start: new CSG.Vector(-1, 0, 0), end: new CSG.Vector(-0.5, 0, 0) });
  d = d.union(CSG.cylinder({ radius: 0.5, slices: 8, start: new CSG.Vector(0.5, 0, 0), end: new CSG.Vector(1, 0, 0) }));
  d = d.union(CSG.cylinder({ radius: 0.5, slices: 8, start: new CSG.Vector(0, -1, 0), end: new CSG.Vector(0, -0.5, 0) }));
  d = d.union(CSG.cylinder({ radius: 0.5, slices: 8, start: new CSG.Vector(0, 0.5, 0), end: new CSG.Vector(0, 1, 0) }));
  d = d.union(CSG.cylinder({ radius: 0.5, slices: 8, start: new CSG.Vector(0, 0, -1), end: new CSG.Vector(0, 0, -0.5) }));
  d = d.union(CSG.cylinder({ radius: 0.5, slices: 8, start: new CSG.Vector(0, 0, 0.5), end: new CSG.Vector(0, 0, 1) }));
  setPolygonColor(a, new CSG.Vector(1, 0, 0));
  setPolygonColor(b, new CSG.Vector(0, 0, 1));
  setPolygonColor(c, new CSG.Vector(0, 1, 0));
  setPolygonColor(d, new CSG.Vector(1, 1, 0));
  var operations = [
    a,
    b,
    a.inverse(),
    b.inverse(),
    a.union(b),
    b.union(a),
    a.subtract(b),
    b.subtract(a),
    a.intersect(b),
    b.intersect(a),
    c,
    d,
    c.inverse(),
    d.inverse(),
    c.union(d),
    d.union(c),
    c.subtract(d),
    d.subtract(c),
    c.intersect(d),
    d.intersect(c)
  ];

  // Build meshes from operations
  for (var i = 0; i < operations.length; i++) {
    var polygons = operations[i].toPolygons();
    var indexer = new Indexer();
    mesh = new Mesh({ normals: true, colors: true });
    polygons.map(function(polygon) {
      var indices = polygon.vertices.map(function(vertex) {
        vertex.color = polygon.shared;
        return indexer.add(vertex);
      });
      for (var i = 2; i < indices.length; i++) {
        mesh.triangles.push([indices[0], indices[i - 1], indices[i]]);
      }
    });
    mesh.vertices = indexer.unique.map(function(v) { return [v.pos.x, v.pos.y, v.pos.z]; });
    mesh.normals = indexer.unique.map(function(v) { return [v.normal.x, v.normal.y, v.normal.z]; });
    mesh.colors = indexer.unique.map(function(v) { return [v.color.x, v.color.y, v.color.z]; });
    mesh.compile();
    meshes.push(mesh);
  }

  // Use diffuse and specular lighting
  shader = new Shader('\
    varying vec3 color;\
    varying vec3 normal;\
    varying vec3 vertex;\
    varying vec3 light;\
    void main() {\
      const vec3 lightDir = vec3(1.0, 2.0, 3.0) / 3.741657386773941;\
      light = (gl_ModelViewMatrix * vec4(lightDir, 0.0)).xyz;\
      color = gl_Color;\
      normal = (gl_ModelViewMatrix * vec4(gl_Normal, 0.0)).xyz;\
      vertex = (gl_ModelViewMatrix * vec4(gl_Vertex, 1.0)).xyz;\
      gl_Position = gl_ModelViewProjectionMatrix * vec4(gl_Vertex, 1.0);\
    }\
  ', '\
    varying vec3 color;\
    varying vec3 normal;\
    varying vec3 vertex;\
    varying vec3 light;\
    void main() {\
      vec3 n = normalize(normal);\
      float diffuse = max(0.0, dot(light, n));\
      float specular = pow(max(0.0, dot(reflect(light, n), normalize(vertex))), 32.0) * sqrt(diffuse);\
      gl_FragColor = vec4(mix(color * (0.2 + 0.8 * diffuse), vec3(1.0), specular), 1.0);\
    }\
  ');

  gl.clearColor(0.75, 0.75, 0.75, 0);
  gl.fullscreen();
  gl.enable(gl.CULL_FACE);
  gl.autoDraw = false;
  draw();
}

function mouseDragged() {
  angleX += deltaMouseY;
  angleY += deltaMouseX;
  angleX = Math.max(-90, Math.min(90, angleX));
  draw();
}

function draw() {
  gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

  var scale = 5;
  var countX = Math.ceil(Math.sqrt(meshes.length));
  var countY = Math.ceil(meshes.length / countX);
  for (var i = 0; i < meshes.length; i++) {
    var mesh = meshes[i];
    var x = Math.floor(i % countX) / (countX - 1) * +2 - 1;
    var y = Math.floor(i / countX) / (countY - 1) * -2 + 1;
    gl.loadIdentity();
    gl.translate(x * scale * 1.5, y * scale, -3.5 * scale);
    gl.rotate(angleX, 1, 0, 0);
    gl.rotate(angleY, 0, 1, 0);

    gl.enable(gl.POLYGON_OFFSET_FILL);
    gl.polygonOffset(1, 1);
    shader.draw(mesh);
    gl.disable(gl.POLYGON_OFFSET_FILL);

    gl.enable(gl.BLEND);
    gl.disable(gl.DEPTH_TEST);
    gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);
    gl.color(0, 0, 0, 0.2);
    gl.begin(gl.LINES);
    for (var j = 0; j < mesh.triangles.length; j++) {
      var tri = mesh.triangles[j];
      var a = Vector.fromArray(mesh.vertices[tri[0]]);
      var b = Vector.fromArray(mesh.vertices[tri[1]]);
      var c = Vector.fromArray(mesh.vertices[tri[2]]);
      gl.vertex(a); gl.vertex(b);
      gl.vertex(b); gl.vertex(c);
      gl.vertex(c); gl.vertex(a);
    }
    gl.end();
    gl.disable(gl.BLEND);
    gl.enable(gl.DEPTH_TEST);
  }
}
