boolean aaa = true;

// light position
ArrayList<float[]> lights;

// camera position
float[] eye;

// camera look at
float[] lookAt;

// right and up vector of the near clipping plane
float[] vector_right;
float[] vector_up;

// top left pixel position on the near clipping plane
float[] upper_left;

// intersecion with which object (multiple when performing supersampling)
ArrayList<Integer> choices;

// the intersecion coordinate (multiple when performing supersampling)
ArrayList<float[]> inter;

// if the pixel is shadowed or not (multiple when performing supersampling)
ArrayList<ArrayList<Boolean>> shadows;

// the color of a pixel (multiple when performing supersampling)
ArrayList<float[]> colors;

// number of supersampling
int sample = 1;

// size of the tiles
int procedural = 25;

// bunny
ArrayList<float[]> vertices;
ArrayList<int[]> triangles;

// to update the canvas or not
boolean toDraw;

float[] intersection;
int index;

float renderstart, renderend;

// this function is executed once for setting up some parameters
void setup() {
  // canvas size
  size(560, 315);

  // lights' positions
  lights = new ArrayList();
  lights.add(new float[]{-100.0, 250.0, 300.0});
  //lights.add(new float[]{300.0, 1000.0, 0.0});

  // camera position
  eye = new float[]{-30.0, 150.0, 275.0};

  // camera look at
  lookAt = new float[]{-40.0, 75.0, 100.0};

  // finding the near clipping plane (upper left pixel position)
  findNear();
  choices = new ArrayList();
  inter = new ArrayList();
  shadows = new ArrayList();
  colors = new ArrayList();
  vertices = new ArrayList();
  triangles = new ArrayList();
  toDraw = true;
  readBunny();
}

void draw() {
  if (toDraw) {
    renderstart = millis();
    for (int i = 0; i < width; i++) {
      for (int j = 0; j < height; j++) {
        choices.clear();
        inter.clear();
        shadows.clear();
        colors.clear();

        for (int k = 0; k < sample; k++) {
          if (k == 0) raytrace(i, j, 0, 0);
          else raytrace(i, j, random(-0.5, 0.5), random(-0.5, 0.5));
          calculate_color();
        }
        float r = 0.0, g = 0.0, b = 0.0;
        for (int k = 0; k < sample; k++) {
          r += colors.get(k)[0];
          g += colors.get(k)[1];
          b += colors.get(k)[2];
        }
        set(i, j, color(r / sample, g / sample, b / sample));
      }
    }
    renderend = millis();
    println(renderend - renderstart);
    //save("no_kd.png");
    toDraw = false;
  }
}

/**
 * This function calculates the position of the top left pixel on the near clipping
 * plane and also the up and right vectors of the near clipping plane.
 */
void findNear() {
  // field of view y (in degree)
  float fovy = 60.0;

  // distance from the camera to the near clipping plane (default value of the perspective() function in processing)
  float near = (height/2.0) / tan(fovy / 2.0 * PI / 180.0) / 10.0;

  // distance from camera to the far clipping plane (default value of the perspective() function in processing)
  //float far = (height/2.0) / tan(fovy / 2.0 * PI / 180.0) * 10.0;

  // the vecter from eye to lookAt
  float[] vector_camera_eye = new float[]{lookAt[0] - eye[0], lookAt[1] - eye[1], lookAt[2] - eye[2]};
  vector_camera_eye = normalize(vector_camera_eye);

  // the center of the near clipping plan
  float[]near_center = new float[]{eye[0] + near * vector_camera_eye[0], eye[1] + near * vector_camera_eye[1], eye[2] + near * vector_camera_eye[2]};

  // the right vector of the near clipping plane (cross product of the vecter from eye to lookAt and the camera up vector)
  vector_right = crossProduct(vector_camera_eye, new float[]{0.0, 1.0, 0.0});

  // normalize the right vector
  vector_right = normalize(vector_right);

  // the up vector of the near clipping plane (cross product of the right vector and the vecter from eye to lookAt)
  vector_up = crossProduct(vector_right, vector_camera_eye);

  // normalize the up vector
  vector_up = normalize(vector_up);

  // auxiliary variable for finding the top center of the near clipping plan
  float n1 = tan(fovy / 2.0 * PI / 180.0) * near;

  // top center of the near clipping plan
  float[] near_topcenter = new float[]{near_center[0] + n1 * vector_up[0], near_center[1] + n1 * vector_up[1], near_center[2] + n1 * vector_up[2]};

  // calculate the pixel size accordingly (height == width)
  float pixel_size;
  if (height % 2 == 1) pixel_size = sqrt((near_topcenter[0] - near_center[0]) * (near_topcenter[0] - near_center[0]) + (near_topcenter[1] - near_center[1]) * (near_topcenter[1] - near_center[1]) + (near_topcenter[2] - near_center[2]) * (near_topcenter[2] - near_center[2])) / (height / 2);
  else pixel_size = sqrt((near_topcenter[0] - near_center[0]) * (near_topcenter[0] - near_center[0]) + (near_topcenter[1] - near_center[1]) * (near_topcenter[1] - near_center[1]) + (near_topcenter[2] - near_center[2]) * (near_topcenter[2] - near_center[2])) / ((height / 2) - 0.5);

  // normalize vector_right wrt the pixel size
  vector_right[0] = vector_right[0] * pixel_size;
  vector_right[1] = vector_right[1] * pixel_size;
  vector_right[2] = vector_right[2] * pixel_size;

  // normalize vector_up wrt the pixel size
  vector_up[0] = vector_up[0] * pixel_size;
  vector_up[1] = vector_up[1] * pixel_size;
  vector_up[2] = vector_up[2] * pixel_size;

  // calculate upper left pixel position accordingly
  if (width % 2 == 1) upper_left = new float[]{near_topcenter[0] - (width / 2) * vector_right[0], near_topcenter[1] - (width / 2) * vector_right[1], near_topcenter[2] - (width / 2) * vector_right[2]};
  else upper_left = new float[]{near_topcenter[0] - (width / 2 - 0.5) * vector_right[0], near_topcenter[1] - (width / 2 - 0.5) * vector_right[1], near_topcenter[2] - (width / 2 - 0.5) * vector_right[2]};
}


/**
 * This function calculates and returns the cross product of two given vectors.
 *
 * @param v1 - the first vector
 * @param v2 - the other vector
 * @return   - the cross product
 */
float[] crossProduct(float[] v1, float[] v2) {
  float[] result = new float[3];
  result[0] = v1[1] * v2[2] - v1[2] * v2[1];
  result[1] = v1[2] * v2[0] - v1[0] * v2[2];
  result[2] = v1[0] * v2[1] - v1[1] * v2[0];
  return result;
}

float dotProduct(float[] v1, float[] v2) {
  return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

/*
 * This function calculates the normalized vector of the given one and returns it.
 *
 * @param v - the given vector
 * @return  - the normalized vector
 */
float[] normalize(float[] v) {
  float temp = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  v[0] = v[0] / temp;
  v[1] = v[1] / temp;
  v[2] = v[2] / temp;
  return v;
}

/**
 * This function perfroms raytracing to find the intersection.
 * sphere1: (x+55)^2 + (y-100)^2 + (z-100)^2 = 60^2
 * sphere2: (x-15)^2 + (y-60)^2 + (z-5)^2 = 45^2
 * floor: y = 0 with -200 <= x <= 200 && -250 <= z <= 250
 *
 * @param i  - varaible for locating the pixel position
 * @param j  - varaible for locating the pixel position
 * @param dx - varaible for translating the pixel position
 * @param dy - varaible for translating the pixel position
 */
void raytrace(int i, int j, float dx, float dy) {
  choices.add(0);
  inter.add(new float[3]);
  shadows.add(new ArrayList());
  for (int k = 0; k < lights.size(); k++) {
    shadows.get(shadows.size() - 1).add(false);
  }
  // the pixel position
  float[] pixel_position = new float[3];
  // the vector of the ray tracing path
  float[] vector = new float[3];


  // vector for checking shadow
  float[] vTemp = new float[3];
  // coeffecients for finding intersections
  float A, B, C, D;
  // auxiliary parameter for finding intersections
  float coefficient;

  pixel_position[0] = upper_left[0] + i * vector_right[0] - j * vector_up[0] + dx * vector_right[0] - dy * vector_up[0];
  pixel_position[1] = upper_left[1] + i * vector_right[1] - j * vector_up[1] + dx * vector_right[1] - dy * vector_up[1];
  pixel_position[2] = upper_left[2] + i * vector_right[2] - j * vector_up[2] + dx * vector_right[2] - dy * vector_up[2];
  vector[0] = pixel_position[0] - eye[0];
  vector[1] = pixel_position[1] - eye[1];
  vector[2] = pixel_position[2] - eye[2];


  // three vertices of a triangle
  float[] vertex0, vertex1, vertex2;
  // vector from vertex1 to vertex0 and vector from vertex1 to vertex2
  float[] vector0 = new float[3];
  float[] vector1 = new float[3];
  // normal of the triangle
  float[] norm = new float[3];
  float[] temp = new float[3];
  // auxiliary parameters for checking if the intersection is inside the triangle
  float uu, vv, a0, b0, c0, a1, b1, c1, a2, b2, c2;
  boolean first = true;
  float d1, d2;
  intersection = new float[3];
  index = -1;
  for (int k = 0; k < triangles.size(); k++) {
    vertex0 = vertices.get(triangles.get(k)[0]);
    vertex1 = vertices.get(triangles.get(k)[1]);
    vertex2 = vertices.get(triangles.get(k)[2]);
    vector0[0] = vertex0[0] - vertex1[0];
    vector0[1] = vertex0[1] - vertex1[1];
    vector0[2] = vertex0[2] - vertex1[2];
    vector1[0] = vertex2[0] - vertex1[0];
    vector1[1] = vertex2[1] - vertex1[1];
    vector1[2] = vertex2[2] - vertex1[2];
    vector0 = normalize(vector0);
    vector1 = normalize(vector1);
    norm = normalize(crossProduct(vector0, vector1));
    // formula of the triangle => Ax + By + Cz = D
    A = norm[0];
    B = norm[1];
    C = norm[2];
    D = A * vertex0[0] + B * vertex0[1] + C * vertex0[2];
    coefficient = (D - A * pixel_position[0] - B * pixel_position[1] - C * pixel_position[2]) / (A * vector[0] + B * vector[1] + C * vector[2]);
    temp[0] = pixel_position[0] + vector[0] * coefficient;
    temp[1] = pixel_position[1] + vector[1] * coefficient;
    temp[2] = pixel_position[2] + vector[2] * coefficient;
    a0 = vertex0[0] - vertex2[0];
    b0 = vertex1[0] - vertex2[0];
    c0 = temp[0] - vertex2[0];
    a1 = vertex0[1] - vertex2[1];
    b1 = vertex1[1] - vertex2[1];
    c1 = temp[1] - vertex2[1];
    a2 = vertex0[2] - vertex2[2];
    b2 = vertex1[2] - vertex2[2];
    c2 = temp[2] - vertex2[2];
    vv = (c0/a0 - c1/a1) / (b0/a0 - b1/a1);
    uu = (c0 - b0 * vv) / a0;
    if (vv >= 0.0 && vv <= 1.0 && uu >= 0.0 && uu <= 1.0 && uu + vv <= 1.0) {
      if (first) {
        intersection[0] = temp[0];
        intersection[1] = temp[1];
        intersection[2] = temp[2];
        index = k;
        first = false;
      } else {
        d1 = sqrt(pow(temp[0] - pixel_position[0], 2) + pow(temp[1] - pixel_position[1], 2) + pow(temp[2] - pixel_position[2], 2));
        d2 = sqrt(pow(intersection[0] - pixel_position[0], 2) + pow(intersection[1] - pixel_position[1], 2) + pow(intersection[2] - pixel_position[2], 2));
        if (d1 < d2) {
          intersection[0] = temp[0];
          intersection[1] = temp[1];
          intersection[2] = temp[2];
          index = k;
        }
      }
    }
  }
}


void calculate_color() {
  float ka= 0.45, kd = 0.45, ks = 0.2, specExp = 10.0;
  float[] ambient, diffuse, specular, c = new float[]{0.0, 0.0, 0.0};
  float[] S, N, V, R;
  float dotTemp, specDot;
  float[] vertex0, vertex1, vertex2, vector0 = new float[3], vector1 = new float[3];
  if (index == -1) colors.add(new float[]{69.0, 150.0, 243.0});
  else {
    vertex0 = vertices.get(triangles.get(index)[0]);
    vertex1 = vertices.get(triangles.get(index)[1]);
    vertex2 = vertices.get(triangles.get(index)[2]);
    vector0[0] = vertex0[0] - vertex1[0];
    vector0[1] = vertex0[1] - vertex1[1];
    vector0[2] = vertex0[2] - vertex1[2];
    vector1[0] = vertex2[0] - vertex1[0];
    vector1[1] = vertex2[1] - vertex1[1];
    vector1[2] = vertex2[2] - vertex1[2];
    vector0 = normalize(vector0);
    vector1 = normalize(vector1);
    N = normalize(crossProduct(vector1, vector0));
    V = new float[]{eye[0] - intersection[0], eye[1] - intersection[1], eye[2] - intersection[2]};
    V = normalize(V);
    if (dotProduct(N, V) < 0) N = normalize(crossProduct(vector0, vector1));

    ambient = new float[]{255.0 * ka, 255.0 * ka, 255.0 * ka};
    c[0] += ambient[0];
    c[1] += ambient[1];
    c[2] += ambient[2];
    S = new float[]{lights.get(0)[0] - intersection[0], lights.get(0)[1] - intersection[1], lights.get(0)[2] - intersection[2]};
    S = normalize(S);
    dotTemp = S[0] * N[0] + S[1] * N[1] + S[2] * N[2];
    R = new float[]{-S[0] + 2 * N[0] * dotTemp, -S[1] + 2 * N[1] * dotTemp, -S[2] + 2 * N[2] * dotTemp};
    R = normalize(R);
    if (dotTemp >= 0) diffuse = new float[]{255.0 * dotTemp * kd, 255.0 * dotTemp * kd, 255.0 * dotTemp * kd};
    else diffuse = new float[]{0.0, 0.0, 0.0};
    dotTemp = R[0] * V[0] + R[1] * V[1] + R[2] * V[2];
    if (dotTemp >= 0) specDot = pow(dotTemp, specExp);
    else specDot = 0.0;
    specular = new float[]{255.0 * specDot * ks, 255.0 * specDot * ks, 255.0 * specDot * ks};
    c[0] += diffuse[0];
    c[1] += diffuse[1];
    c[2] += diffuse[2];
    c[0] += specular[0];
    c[1] += specular[1];
    c[2] += specular[2];
    colors.add(c);
  }
}

void readBunny() {
  String[] lines = loadStrings("bun_zipper_res4.ply");
  lines = loadStrings("bun_zipper.ply");
  String[] line;
  float scale = 1000.0;
  float y_trans = -60.0;
  float x_trans = -12.0;
  for (String s : lines) {
    line = s.split(" ");
    if (line.length == 5) {
      vertices.add(new float[3]);
      vertices.get(vertices.size() - 1)[0] = Float.parseFloat(line[0]) * scale + x_trans;
      vertices.get(vertices.size() - 1)[1] = Float.parseFloat(line[1]) * scale + y_trans;
      vertices.get(vertices.size() - 1)[2] = Float.parseFloat(line[2]) * scale;
    } else {
      triangles.add(new int[3]);
      triangles.get(triangles.size() - 1)[0] = Integer.parseInt(line[1]);
      triangles.get(triangles.size() - 1)[1] = Integer.parseInt(line[2]);
      triangles.get(triangles.size() - 1)[2] = Integer.parseInt(line[3]);
    }
  }
}
