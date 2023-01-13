int a = 200, b = 250;
size(560, 315, P3D);
camera(-30, -150, 275, -40, -75, 100, 0, 1, 0);
// camera(300, -150, 100, 0, 0, 0, 0, 1, 0);
perspective();
noStroke();
pointLight(255, 255, 255, -100, -250, 300);
background(69, 150, 243);

// floor
/*fill(255, 228, 108);
beginShape(QUADS);
vertex(-a, 0, -b);
vertex(-a, 0, b);
vertex(a, 0, b);
vertex(a, 0, -b);
endShape();*/

// sphere 1
/*translate(-55, -100, 100);
fill(150);
sphere(60);
translate(55, 100, -100);*/

// sphere 2
/*translate(15, -60, 5);
fill(250);
sphere(45);
translate(-15, 60, -5);*/


// bunny
ArrayList<float[]> points = new ArrayList();;
ArrayList<int[]> triangles = new ArrayList();;

String[] lines = loadStrings("bun_zipper_res4.ply");
//lines = loadStrings("bun_zipper.ply");
String[] line;
for (String s : lines) {
  line = s.split(" ");
  if (line.length == 5) {
    points.add(new float[3]);
    points.get(points.size() - 1)[0] = Float.parseFloat(line[0]);
    points.get(points.size() - 1)[1] = Float.parseFloat(line[1]);
    points.get(points.size() - 1)[2] = Float.parseFloat(line[2]);
  } else {
    triangles.add(new int[3]);
    triangles.get(triangles.size() - 1)[0] = Integer.parseInt(line[1]);
    triangles.get(triangles.size() - 1)[1] = Integer.parseInt(line[2]);
    triangles.get(triangles.size() - 1)[2] = Integer.parseInt(line[3]);
  }
}
fill(255);
float[] point;
int scale = 1000;
float y_trans = -60;
float x_trans = -12;
beginShape(TRIANGLES);
for (int i = 0; i < triangles.size(); i++) {
  point = points.get(triangles.get(i)[0]);
  vertex(scale * point[0] + x_trans, -(scale * point[1] + y_trans), scale * point[2]);
  point = points.get(triangles.get(i)[1]);
  vertex(scale * point[0] + x_trans, -(scale * point[1] + y_trans), scale * point[2]);
  point = points.get(triangles.get(i)[2]);
  vertex(scale * point[0] + x_trans, -(scale * point[1] + y_trans), scale * point[2]);
}
endShape();
//save("kd_II.png");
