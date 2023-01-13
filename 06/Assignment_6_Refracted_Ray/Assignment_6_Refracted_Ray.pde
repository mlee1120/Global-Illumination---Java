void setup() {
  float theta = (89.0 * PI / 180.0);
  Quaternion q = new Quaternion(cos(theta / 2), 0.0, 0.0, sin(theta / 2));
  Quaternion qp = new Quaternion(0.0, 0.0, -1.0, 0.0);
  Quaternion newp = q.product(qp).product(q.inverse());
  println(newp.x + " " + newp.y + " " + newp.z);
  float x = 1.0 * cos(PI / 2 - theta);
  float y = -1.0 * cos(theta);
  float[] v = new float[]{x, y, 0};
  println(v[0] + " " + v[1] + " " + v[2]);
}
float[] normalize(float[] v) {
  float temp = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  return new float[]{v[0] / temp, v[1] / temp, v[2] / temp};
}
