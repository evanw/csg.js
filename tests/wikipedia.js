var angleX = 20;
var angleY = 20;
var shader;
var mesh;

function setup() {
  function setPolygonColor(solid, color) {
    solid.toPolygons().map(function(polygon) {
      polygon.shared = color;
    });
  }

  // Build CSG example from Wikipedia
  var cube = CSG.cube();
  var sphere = CSG.sphere({ radius: 1.35, stacks: 12 });
  var cyl1 = CSG.cylinder({ radius: 0.7, start: new CSG.Vector(-2, 0, 0), end: new CSG.Vector(2, 0, 0) });
  var cyl2 = CSG.cylinder({ radius: 0.7, start: new CSG.Vector(0, -2, 0), end: new CSG.Vector(0, 2, 0) });
  var cyl3 = CSG.cylinder({ radius: 0.7, start: new CSG.Vector(0, 0, -2), end: new CSG.Vector(0, 0, 2) });
  setPolygonColor(cube, new CSG.Vector(1, 0, 0));
  setPolygonColor(sphere, new CSG.Vector(0, 0, 1));
  setPolygonColor(cyl1, new CSG.Vector(0, 1, 0));
  setPolygonColor(cyl2, new CSG.Vector(0, 1, 0));
  setPolygonColor(cyl3, new CSG.Vector(0, 1, 0));
  var polygons = cube.intersect(sphere).subtract(cyl1.union(cyl2).union(cyl3)).toPolygons();

  // Build mesh from result
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
  gl.loadIdentity();
  gl.translate(0, 0, -5);
  gl.rotate(angleX, 1, 0, 0);
  gl.rotate(angleY, 0, 1, 0);

  gl.enable(gl.POLYGON_OFFSET_FILL);
  gl.polygonOffset(1, 1);
  shader.draw(mesh);
  gl.disable(gl.POLYGON_OFFSET_FILL);

  gl.enable(gl.BLEND);
  gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);
  gl.color(0, 0, 0, 0.4);
  gl.begin(gl.LINES);
  for (var i = 0; i < mesh.triangles.length; i++) {
    var tri = mesh.triangles[i];
    var a = Vector.fromArray(mesh.vertices[tri[0]]);
    var b = Vector.fromArray(mesh.vertices[tri[1]]);
    var c = Vector.fromArray(mesh.vertices[tri[2]]);
    gl.vertex(a); gl.vertex(b);
    gl.vertex(b); gl.vertex(c);
    gl.vertex(c); gl.vertex(a);
  }
  gl.end();
  gl.disable(gl.BLEND);
}
