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

// all objects in the scene (for kd tree)
ArrayList<Primitive> objects;

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

// to update the canvas or not
boolean toDraw;

// max reflection time
int max_depth = 5;

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
  eye = new float[]{30.0, 150.0, 275.0};
  //eye = new float[]{300.0, 150.0, 100.0};

  // camera look at
  lookAt = new float[]{-40.0, 75.0, 100.0};
  //lookAt = new float[]{0.0, 0.0, 0.0};

  // create all primitives
  objects = new ArrayList();
  objects.add(new Primitive("background", null, null, 0.0, 0.0, 0.0, 0.0));
  objects.add(new Primitive("sphere", new float[]{-55.0, 100.0, 100.0}, new float[]{60.0}, 0.5, 0.6, 0.3, 30.0));
  objects.add(new Primitive("sphere", new float[]{15.0, 60.0, 5.0}, new float[]{45.0}, 0.5, 0.5, 0.3, 30.0));
  objects.add(new Primitive("plane", new float[]{0.0, 0.0, 0.0}, new float[]{400.0, 500.0}, 0.5, 0.5, 0.3, 30.0));

  // finding the near clipping plane (upper left pixel position)
  findNear();
  choices = new ArrayList();
  inter = new ArrayList();
  shadows = new ArrayList();
  colors = new ArrayList();
  toDraw = true;
}

void draw() {
  if (toDraw) {
    // move camera
    //eye[0] += .5;
    //eye[2] += 1;
    //findNear();

    for (int i = 0; i < width; i++) {
      for (int j = 0; j < height; j++) {
        choices.clear();
        inter.clear();
        shadows.clear();
        colors.clear();

        // multisampling
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
    //save("Assignment_5_Reflection.png");
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

float mod(float f1, float f2) {
  return f1 - f2 * floor(f1 / f2);
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

  // check sphere1
  A = vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2];
  B = 2 * (pixel_position[0] + 55) * vector[0] + 2 * (pixel_position[1] - 100) * vector[1] + 2 * (pixel_position[2] - 100) * vector[2];
  C = pow(pixel_position[0] + 55, 2) + pow(pixel_position[1] - 100, 2) + pow(pixel_position[2] - 100, 2) - 60 * 60;
  D = B * B - 4 * A * C;
  if (D >= 0) {
    coefficient = (-B - sqrt(D)) / (2 * A);
    if (coefficient <= 0) coefficient = (-B + sqrt(D)) / (2 * A);
    if (coefficient > 0) {
      choices.set(choices.size() - 1, 1);
      inter.get(inter.size() - 1)[0] = pixel_position[0] + vector[0] * coefficient;
      inter.get(inter.size() - 1)[1] = pixel_position[1] + vector[1] * coefficient;
      inter.get(inter.size() - 1)[2] = pixel_position[2] + vector[2] * coefficient;

      // check shadow
      for (int k = 0; k < lights.size(); k++) {
        vTemp[0] = lights.get(k)[0] - inter.get(inter.size() - 1)[0];
        vTemp[1] = lights.get(k)[1] - inter.get(inter.size() - 1)[1];
        vTemp[2] = lights.get(k)[2] - inter.get(inter.size() - 1)[2];
        A = vTemp[0] * vTemp[0] + vTemp[1] * vTemp[1] + vTemp[2] * vTemp[2];
        B = 2 * (inter.get(inter.size() - 1)[0] + 55) * vTemp[0] + 2 * (inter.get(inter.size() - 1)[1] - 100) * vTemp[1] + 2 * (inter.get(inter.size() - 1)[2] - 100) * vTemp[2];
        C = pow(inter.get(inter.size() - 1)[0] + 55, 2) + pow(inter.get(inter.size() - 1)[1] - 100, 2) + pow(inter.get(inter.size() - 1)[2] - 100, 2) - 60 * 60;
        D = B * B - 4 * A * C;
        if (D >= 0) {
          if ((-B - sqrt(D)) / (2 * A) > 0 && (-B + sqrt(D)) / (2 * A) > 0) shadows.get(shadows.size() - 1).set(k, true);
        }
      }
    }
  }

  // check sphere2
  if (choices.get(choices.size() - 1) == 0) {
    A = vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2];
    B = 2 * (pixel_position[0] - 15) * vector[0] + 2 * (pixel_position[1] - 60) * vector[1] + 2 * (pixel_position[2] - 5) * vector[2];
    C = (pixel_position[0] - 15) * (pixel_position[0] - 15) + (pixel_position[1] - 60) * (pixel_position[1] - 60) + (pixel_position[2] - 5) * (pixel_position[2] - 5) - 45 * 45;
    D = B * B - 4 * A * C;
    if (D >= 0) {
      coefficient = (-B - sqrt(D)) / (2 * A);
      if (coefficient <= 0) coefficient = (-B + sqrt(D)) / (2 * A);
      if (coefficient > 0) {
        choices.set(choices.size() - 1, 2);
        inter.get(inter.size() - 1)[0] = pixel_position[0] + vector[0] * coefficient;
        inter.get(inter.size() - 1)[1] = pixel_position[1] + vector[1] * coefficient;
        inter.get(inter.size() - 1)[2] = pixel_position[2] + vector[2] * coefficient;

        // check shadow
        for (int k = 0; k < lights.size(); k++) {
          vTemp[0] = lights.get(k)[0] - inter.get(inter.size() - 1)[0];
          vTemp[1] = lights.get(k)[1] - inter.get(inter.size() - 1)[1];
          vTemp[2] = lights.get(k)[2] - inter.get(inter.size() - 1)[2];
          A = vTemp[0] * vTemp[0] + vTemp[1] * vTemp[1] + vTemp[2] * vTemp[2];
          B = 2 * (inter.get(inter.size() - 1)[0] - 15) * vTemp[0] + 2 * (inter.get(inter.size() - 1)[1] - 60) * vTemp[1] + 2 * (inter.get(inter.size() - 1)[2] - 5) * vTemp[2];
          C = pow(inter.get(inter.size() - 1)[0] - 15, 2) + pow(inter.get(inter.size() - 1)[1] - 60, 2) + pow(inter.get(inter.size() - 1)[2] - 5, 2) - 45 * 45;
          D = B * B - 4 * A * C;
          if (D >=0) {
            if ((-B - sqrt(D)) / (2 * A) > 0 && (-B + sqrt(D)) / (2 * A) > 0) shadows.get(shadows.size() - 1).set(k, true);
          }

          A = vTemp[0] * vTemp[0] + vTemp[1] * vTemp[1] + vTemp[2] * vTemp[2];
          B = 2 * (inter.get(inter.size() - 1)[0] + 55) * vTemp[0] + 2 * (inter.get(inter.size() - 1)[1] - 100) * vTemp[1] + 2 * (inter.get(inter.size() - 1)[2] - 100) * vTemp[2];
          C = pow(inter.get(inter.size() - 1)[0] + 55, 2) + pow(inter.get(inter.size() - 1)[1] - 100, 2) + pow(inter.get(inter.size() - 1)[2] - 100, 2) - 60 * 60;
          D = B * B - 4 * A * C;
          if (D >= 0) {
            if ((-B - sqrt(D)) / (2 * A) > 0 && (-B + sqrt(D)) / (2 * A) > 0) shadows.get(shadows.size() - 1).set(k, true);
          }
        }
      }
    }
  }

  // check floor
  if (choices.get(choices.size() - 1) == 0) {
    coefficient = -pixel_position[1] / vector[1];
    if (coefficient > 0) {
      inter.get(inter.size() - 1)[0] = pixel_position[0] + vector[0] * coefficient;
      inter.get(inter.size() - 1)[1] = pixel_position[1] + vector[1] * coefficient;
      inter.get(inter.size() - 1)[2] = pixel_position[2] + vector[2] * coefficient;
      if (inter.get(inter.size() - 1)[0] <= 200 && inter.get(inter.size() - 1)[0] >= -200 && inter.get(inter.size() - 1)[2] <= 250 && inter.get(inter.size() - 1)[2] >= -250) {
        choices.set(choices.size() - 1, 3);
      }

      // check shadow
      for (int k = 0; k < lights.size(); k++) {
        vTemp[0] = lights.get(k)[0] - inter.get(inter.size() - 1)[0];
        vTemp[1] = lights.get(k)[1] - inter.get(inter.size() - 1)[1];
        vTemp[2] = lights.get(k)[2] - inter.get(inter.size() - 1)[2];
        A = vTemp[0] * vTemp[0] + vTemp[1] * vTemp[1] + vTemp[2] * vTemp[2];
        B = 2 * (inter.get(inter.size() - 1)[0] - 15) * vTemp[0] + 2 * (inter.get(inter.size() - 1)[1] - 60) * vTemp[1] + 2 * (inter.get(inter.size() - 1)[2] - 5) * vTemp[2];
        C = pow(inter.get(inter.size() - 1)[0] - 15, 2)  + pow(inter.get(inter.size() - 1)[1] - 60, 2) + pow(inter.get(inter.size() - 1)[2] - 5, 2) - 45 * 45;
        D = B * B - 4 * A * C;
        if (D >=0) {
          if ((-B - sqrt(D)) / (2 * A) > 0 && (-B + sqrt(D)) / (2 * A) > 0) shadows.get(shadows.size() - 1).set(k, true);
        }

        A = vTemp[0] * vTemp[0] + vTemp[1] * vTemp[1] + vTemp[2] * vTemp[2];
        B = 2 * (inter.get(inter.size() - 1)[0] + 55) * vTemp[0] + 2 * (inter.get(inter.size() - 1)[1] - 100) * vTemp[1] + 2 * (inter.get(inter.size() - 1)[2] - 100) * vTemp[2];
        C = pow(inter.get(inter.size() - 1)[0] + 55, 2) + pow(inter.get(inter.size() - 1)[1] - 100, 2) + pow(inter.get(inter.size() - 1)[2] - 100, 2) - 60 * 60;
        D = B * B - 4 * A * C;
        if (D >= 0) {
          if ((-B - sqrt(D)) / (2 * A) > 0 && (-B + sqrt(D)) / (2 * A) > 0) shadows.get(shadows.size() - 1).set(k, true);
        }
      }
    }
  }
}


void calculate_color() {
  float ka, kd, ks, specExp; // reflection coefficients
  float[] ambient, diffuse, specular, c = new float[]{0.0, 0.0, 0.0};
  float[] S, N, V, R;
  float dotTemp, specDot;
  float u, v, dx, dz;

  // draw scene
  switch(choices.get(choices.size() - 1)) {
  case 0:
    colors.add(new float[]{69.0, 150.0, 243.0});
    break;
  case 1:
    ka = 0.5;
    kd = 0.6;
    ks = 0.3;
    specExp = 30.0;

    // normal vector
    N = new float[]{inter.get(inter.size() - 1)[0] + 55.0, inter.get(inter.size() - 1)[1] - 100.0, inter.get(inter.size() - 1)[2] - 100.0};
    // normalize N
    N = normalize(N);

    // destination (to the camera)
    V = new float[]{eye[0] - inter.get(inter.size() - 1)[0], eye[1] - inter.get(inter.size() - 1)[1], eye[2] - inter.get(inter.size() - 1)[2]};
    // normalize V
    V = normalize(V);

    ambient = new float[]{150.0 * ka, 150.0 * ka, 150.0 * ka};
    c[0] += ambient[0];
    c[1] += ambient[1];
    c[2] += ambient[2];
    for (int i = 0; i < lights.size(); i++) {
      // direction of incoming light
      S = new float[]{lights.get(i)[0] - inter.get(inter.size() - 1)[0], lights.get(i)[1] - inter.get(inter.size() - 1)[1], lights.get(i)[2] - inter.get(inter.size() - 1)[2]};
      // normalize S
      S = normalize(S);

      // reflection vector
      dotTemp = S[0] * N[0] + S[1] * N[1] + S[2] * N[2];
      R = new float[]{-S[0] + 2 * N[0] * dotTemp, -S[1] + 2 * N[1] * dotTemp, -S[2] + 2 * N[2] * dotTemp};

      dotTemp = N[0] * S[0] + N[1] * S[1] + N[2] * S[2];
      if (dotTemp >= 0) diffuse = new float[]{150.0 * dotTemp * kd, 150.0 * dotTemp * kd, 150.0 * dotTemp * kd};
      else diffuse = new float[]{0.0, 0.0, 0.0};
      dotTemp = R[0] * V[0] + R[1] * V[1] + R[2] * V[2];
      if (dotTemp >= 0) specDot = pow(dotTemp, specExp);
      else specDot = 0.0;
      specular = new float[]{255.0 * specDot * ks, 255.0 * specDot * ks, 255.0 * specDot * ks};
      if (!shadows.get(shadows.size() - 1).get(i)) {
        c[0] += diffuse[0];
        c[1] += diffuse[1];
        c[2] += diffuse[2];
        c[0] += specular[0];
        c[1] += specular[1];
        c[2] += specular[2];
      }
      colors.add(c);
    }
    break;
  case 2:
    ka = 0.5;
    kd = 0.5;
    ks = 0.3;
    specExp = 30.0;

    // normal vector
    N = new float[]{inter.get(inter.size() - 1)[0] - 15.0, inter.get(inter.size() - 1)[1] - 60.0, inter.get(inter.size() - 1)[2] - 5.0};
    // normalize N
    N = normalize(N);

    // destination (to the camera)
    V = new float[]{eye[0] - inter.get(inter.size() - 1)[0], eye[1] - inter.get(inter.size() - 1)[1], eye[2] - inter.get(inter.size() - 1)[2]};
    // normalize V
    V = normalize(V);

    ambient = new float[]{250.0 * ka, 250.0 * ka, 250.0 * ka};
    c[0] += ambient[0];
    c[1] += ambient[1];
    c[2] += ambient[2];
    for (int i = 0; i < lights.size(); i++) {
      // direction of incoming light
      S = new float[]{lights.get(i)[0] - inter.get(inter.size() - 1)[0], lights.get(i)[1] - inter.get(inter.size() - 1)[1], lights.get(i)[2] - inter.get(inter.size() - 1)[2]};
      // normalize S
      S = normalize(S);

      // reflection vector
      dotTemp = S[0] * N[0] + S[1] * N[1] + S[2] * N[2];
      R = new float[]{-S[0] + 2 * N[0] * dotTemp, -S[1] + 2 * N[1] * dotTemp, -S[2] + 2 * N[2] * dotTemp};

      dotTemp = N[0] * S[0] + N[1] * S[1] + N[2] * S[2];
      if (dotTemp >= 0) diffuse = new float[]{250.0 * dotTemp * kd, 250.0 * dotTemp * kd, 250.0 * dotTemp * kd};
      else diffuse = new float[]{0.0, 0.0, 0.0};
      dotTemp = R[0] * V[0] + R[1] * V[1] + R[2] * V[2];
      if (dotTemp >= 0) specDot = pow(dotTemp, specExp);
      else specDot = 0.0;
      specular = new float[]{255.0 * specDot * ks, 255.0 * specDot * ks, 255.0 * specDot * ks};
      if (!shadows.get(shadows.size() - 1).get(i)) {
        c[0] += diffuse[0];
        c[1] += diffuse[1];
        c[2] += diffuse[2];
        c[0] += specular[0];
        c[1] += specular[1];
        c[2] += specular[2];
      }
      // calculating reflection
      dotTemp = V[0] * N[0] + V[1] * N[1] + V[2] * N[2];
      float[] RR = new float[]{-V[0] + 2 * N[0] * dotTemp, -V[1] + 2 * N[1] * dotTemp, -V[2] + 2 * N[2] * dotTemp};
      float[] reflection = reflect(max_depth, 0.3, inter.get(inter.size() - 1), RR);
      c[0] += reflection[0];
      c[1] += reflection[1];
      c[2] += reflection[2];
      colors.add(c);
    }
    break;
  case 3:
    ka = 0.5;
    kd = 0.5;
    ks = 0.3;
    specExp = 30.0;

    // normal vector
    N = new float[]{0.0, 1.0, 0.0};
    // normalize N
    N = normalize(N);

    // destination (to the camera)
    V = new float[]{eye[0] - inter.get(inter.size() - 1)[0], eye[1] - inter.get(inter.size() - 1)[1], eye[2] - inter.get(inter.size() - 1)[2]};
    // normalize V
    V = normalize(V);

    u = (inter.get(inter.size() - 1)[0] - (-200.0)) / (200.0 - (-200.0));
    v = (inter.get(inter.size() - 1)[2] - (-250.0)) / (250.0 - (-250.0));

    // checkerboard size in uv
    dx = procedural / 400.0;
    dz = procedural / 500.0;
    if ((ceil(u / dx) + ceil(v / dz)) % 2 == 0) ambient = new float[]{255.0 * ka, 255.0 * ka, 0.0 * ka};
    else ambient = new float[]{255.0 * ka, 0.0 * ka, 0.0 * ka};

    //ambient = new float[]{255.0 * ka, 228.0 * ka, 108.0 * ka};
    c[0] += ambient[0];
    c[1] += ambient[1];
    c[2] += ambient[2];

    for (int i = 0; i < lights.size(); i++) {
      // direction of incoming light
      S = new float[]{lights.get(i)[0] - inter.get(inter.size() - 1)[0], lights.get(i)[1] - inter.get(inter.size() - 1)[1], lights.get(i)[2] - inter.get(inter.size() - 1)[2]};
      // normalize S
      S = normalize(S);

      // reflection vector
      dotTemp = S[0] * N[0] + S[1] * N[1] + S[2] * N[2];
      R = new float[]{-S[0] + 2 * N[0] * dotTemp, -S[1] + 2 * N[1] * dotTemp, -S[2] + 2 * N[2] * dotTemp};

      dotTemp = N[0] * S[0] + N[1] * S[1] + N[2] * S[2];
      if (dotTemp >= 0) {
        if ((ceil(u / dx) + ceil(v / dz)) % 2 == 0) diffuse = new float[]{255.0 * dotTemp * kd, 255.0 * dotTemp * kd, 0.0 * dotTemp * kd};
        else diffuse = new float[]{255.0 * dotTemp * kd, 0.0 * dotTemp * kd, 0.0 * dotTemp * kd};
      } else diffuse = new float[]{0.0, 0.0, 0.0};
      dotTemp = R[0] * V[0] + R[1] * V[1] + R[2] * V[2];
      if (dotTemp >= 0) specDot = pow(dotTemp, specExp);
      else specDot = 0.0;
      specular = new float[]{255.0 * specDot * ks, 255.0 * specDot * ks, 255.0 * specDot * ks};
      if (!shadows.get(shadows.size() - 1).get(i)) {
        c[0] += diffuse[0];
        c[1] += diffuse[1];
        c[2] += diffuse[2];
        c[0] += specular[0];
        c[1] += specular[1];
        c[2] += specular[2];
      }
      colors.add(c);
    }
    break;
  }
}

float[] reflect(int depth, float kr, float[] point, float[] direction) {
  float[] result = new float[]{0.0, 0.0, 0.0};
  if (depth != 0) {
    depth -= 1;
    int choiceR = 0;
    // coeffecients for finding reflecting intersections
    float A, B, C, D;
    // auxiliary parameter for finding reflecting intersections
    float coefficient;
    // intersection
    float[] intersectionR = new float[3];
    // if the intersection is under shadows
    ArrayList<Boolean> shadowsR = new ArrayList();
    // vector for checking shadow
    float[] vTemp = new float[3];


    // check sphere1
    A = direction[0] * direction[0] + direction[1] * direction[1] + direction[2] * direction[2];
    B = 2 * (point[0] + 55) * direction[0] + 2 * (point[1] - 100) * direction[1] + 2 * (point[2] - 100) * direction[2];
    C = pow(point[0] + 55, 2) + pow(point[1] - 100, 2) + pow(point[2] - 100, 2) - 60 * 60;
    D = B * B - 4 * A * C;
    if (D >= 0) {
      coefficient = (-B - sqrt(D)) / (2 * A);
      if (coefficient <= 0) coefficient = (-B + sqrt(D)) / (2 * A);
      if (coefficient > 0.001) {
        choiceR = 1;
        intersectionR[0] = point[0] + direction[0] * coefficient;
        intersectionR[1] = point[1] + direction[1] * coefficient;
        intersectionR[2] = point[2] + direction[2] * coefficient;

        // check shadow
        for (int i = 0; i < lights.size(); i++) {
          vTemp[0] = lights.get(i)[0] - intersectionR[0];
          vTemp[1] = lights.get(i)[1] - intersectionR[1];
          vTemp[2] = lights.get(i)[2] - intersectionR[2];
          A = vTemp[0] * vTemp[0] + vTemp[1] * vTemp[1] + vTemp[2] * vTemp[2];
          B = 2 * (intersectionR[0] + 55) * vTemp[0] + 2 * (intersectionR[1] - 100) * vTemp[1] + 2 * (intersectionR[2] - 100) * vTemp[2];
          C = pow(intersectionR[0] + 55, 2) + pow(intersectionR[1] - 100, 2) + pow(intersectionR[2] - 100, 2) - 60 * 60;
          D = B * B - 4 * A * C;
          if (D >= 0) {
            if ((-B - sqrt(D)) / (2 * A) > 0 && (-B + sqrt(D)) / (2 * A) > 0) shadowsR.add(true);
            else shadowsR.add(false);
          } else shadowsR.add(false);
        }
      }
    }

    // check sphere2
    if (choiceR == 0) {
      A = direction[0] * direction[0] + direction[1] * direction[1] + direction[2] * direction[2];
      B = 2 * (point[0] - 15) * direction[0] + 2 * (point[1] - 60) * direction[1] + 2 * (point[2] - 5) * direction[2];
      C = (point[0] - 15) * (point[0] - 15) + (point[1] - 60) * (point[1] - 60) + (point[2] - 5) * (point[2] - 5) - 45 * 45;
      D = B * B - 4 * A * C;
      if (D >= 0) {
        coefficient = (-B - sqrt(D)) / (2 * A);
        if (coefficient <= 0) coefficient = (-B + sqrt(D)) / (2 * A);
        // should be 0, but there might be error
        if (coefficient > 0.001) {
          choiceR = 2;
          intersectionR[0] = point[0] + direction[0] * coefficient;
          intersectionR[1] = point[1] + direction[1] * coefficient;
          intersectionR[2] = point[2] + direction[2] * coefficient;

          // check shadow
          for (int i = 0; i < lights.size(); i++) {
            vTemp[0] = lights.get(i)[0] - intersectionR[0];
            vTemp[1] = lights.get(i)[1] - intersectionR[1];
            vTemp[2] = lights.get(i)[2] - intersectionR[2];
            A = vTemp[0] * vTemp[0] + vTemp[1] * vTemp[1] + vTemp[2] * vTemp[2];
            B = 2 * (intersectionR[0] - 15) * vTemp[0] + 2 * (intersectionR[1] - 60) * vTemp[1] + 2 * (intersectionR[2] - 5) * vTemp[2];
            C = pow(intersectionR[0] - 15, 2) + pow(intersectionR[1] - 60, 2) + pow(intersectionR[2] - 5, 2) - 45 * 45;
            D = B * B - 4 * A * C;
            if (D >=0) {
              if ((-B - sqrt(D)) / (2 * A) > 0 && (-B + sqrt(D)) / (2 * A) > 0) shadowsR.add(true);
              else shadowsR.add(false);
            } else shadowsR.add(false);

            if (!shadowsR.get(shadowsR.size() - 1)) {
              A = vTemp[0] * vTemp[0] + vTemp[1] * vTemp[1] + vTemp[2] * vTemp[2];
              B = 2 * (intersectionR[0] + 55) * vTemp[0] + 2 * (intersectionR[1] - 100) * vTemp[1] + 2 * (intersectionR[2] - 100) * vTemp[2];
              C = pow(intersectionR[0] + 55, 2) + pow(intersectionR[1] - 100, 2) + pow(intersectionR[2] - 100, 2) - 60 * 60;
              D = B * B - 4 * A * C;
              if (D >= 0) {
                if ((-B - sqrt(D)) / (2 * A) > 0 && (-B + sqrt(D)) / (2 * A) > 0) shadowsR.set(shadowsR.size() - 1, true);
              }
            }
          }
        }
      }
    }

    // check floor
    if (choiceR == 0) {
      coefficient = -point[1] / direction[1];
      if (coefficient >= 0) {
        intersectionR[0] = point[0] + direction[0] * coefficient;
        intersectionR[1] = point[1] + direction[1] * coefficient;
        intersectionR[2] = point[2] + direction[2] * coefficient;
        if (intersectionR[0] <= 200 && intersectionR[0] >= -200 && intersectionR[2] <= 250 && intersectionR[2] >= -250) {
          choiceR = 3;
        }


        // check shadow
        for (int i = 0; i < lights.size(); i++) {
          vTemp[0] = lights.get(i)[0] - intersectionR[0];
          vTemp[1] = lights.get(i)[1] - intersectionR[1];
          vTemp[2] = lights.get(i)[2] - intersectionR[2];
          A = vTemp[0] * vTemp[0] + vTemp[1] * vTemp[1] + vTemp[2] * vTemp[2];
          B = 2 * (intersectionR[0] - 15) * vTemp[0] + 2 * (intersectionR[1] - 60) * vTemp[1] + 2 * (intersectionR[2] - 5) * vTemp[2];
          C = pow(intersectionR[0] - 15, 2)  + pow(intersectionR[1] - 60, 2) + pow(intersectionR[2] - 5, 2) - 45 * 45;
          D = B * B - 4 * A * C;
          if (D >=0) {
            if ((-B - sqrt(D)) / (2 * A) > 0 && (-B + sqrt(D)) / (2 * A) > 0) shadowsR.add(true);
            else shadowsR.add(false);
          } else shadowsR.add(false);
          if (!shadowsR.get(shadowsR.size() - 1)) {
            A = vTemp[0] * vTemp[0] + vTemp[1] * vTemp[1] + vTemp[2] * vTemp[2];
            B = 2 * (intersectionR[0] + 55) * vTemp[0] + 2 * (intersectionR[1] - 100) * vTemp[1] + 2 * (intersectionR[2] - 100) * vTemp[2];
            C = pow(intersectionR[0] + 55, 2) + pow(intersectionR[1] - 100, 2) + pow(intersectionR[2] - 100, 2) - 60 * 60;
            D = B * B - 4 * A * C;
            if (D >= 0) {
              if ((-B - sqrt(D)) / (2 * A) > 0 && (-B + sqrt(D)) / (2 * A) > 0) shadowsR.set(shadowsR.size() - 1, true);
            }
          }
        }
      }
    }



    float ka, kd, ks, specExp; // reflection coefficients
    float[] ambient, diffuse, specular, c = new float[]{0.0, 0.0, 0.0};
    float[] S, N, V, R;
    float dotTemp, specDot;
    float u, v, dx, dz;
    float[] cr = new float[]{0.0, 0.0, 0.0};
    float[] RR;

    // draw scene
    switch(choiceR) {
    case 0:
      c = new float[]{69.0, 150.0, 243.0};
      break;
    case 1:
      ka = 0.5;
      kd = 0.6;
      ks = 0.3;
      specExp = 30.0;

      // normal vector
      N = new float[]{intersectionR[0] + 55.0, intersectionR[1] - 100.0, intersectionR[2] - 100.0};
      // normalize N
      N = normalize(N);

      // destination (to the camera)
      V = new float[]{point[0] - intersectionR[0], point[1] - intersectionR[1], point[2] - intersectionR[2]};
      // normalize V
      V = normalize(V);

      ambient = new float[]{150.0 * ka, 150.0 * ka, 150.0 * ka};
      c[0] += ambient[0];
      c[1] += ambient[1];
      c[2] += ambient[2];
      for (int i = 0; i < lights.size(); i++) {
        // direction of incoming light
        S = new float[]{lights.get(i)[0] - intersectionR[0], lights.get(i)[1] - intersectionR[1], lights.get(i)[2] - intersectionR[2]};
        // normalize S
        S = normalize(S);

        // reflection vector
        dotTemp = S[0] * N[0] + S[1] * N[1] + S[2] * N[2];
        R = new float[]{-S[0] + 2 * N[0] * dotTemp, -S[1] + 2 * N[1] * dotTemp, -S[2] + 2 * N[2] * dotTemp};

        dotTemp = N[0] * S[0] + N[1] * S[1] + N[2] * S[2];
        if (dotTemp >= 0) diffuse = new float[]{150.0 * dotTemp * kd, 150.0 * dotTemp * kd, 150.0 * dotTemp * kd};
        else diffuse = new float[]{0.0, 0.0, 0.0};
        dotTemp = R[0] * V[0] + R[1] * V[1] + R[2] * V[2];
        if (dotTemp >= 0) specDot = pow(dotTemp, specExp);
        else specDot = 0.0;
        specular = new float[]{255.0 * specDot * ks, 255.0 * specDot * ks, 255.0 * specDot * ks};
        if (!shadowsR.get(i)) {
          c[0] += diffuse[0];
          c[1] += diffuse[1];
          c[2] += diffuse[2];
          c[0] += specular[0];
          c[1] += specular[1];
          c[2] += specular[2];
        }
      }
      // calculating reflection
      dotTemp = V[0] * N[0] + V[1] * N[1] + V[2] * N[2];
      RR = new float[]{-V[0] + 2 * N[0] * dotTemp, -V[1] + 2 * N[1] * dotTemp, -V[2] + 2 * N[2] * dotTemp};
      cr = reflect(depth, 0.0, intersectionR, RR);
      c[0] += cr[0];
      c[1] += cr[1];
      c[2] += cr[2];
      break;
    case 2:
      ka = 0.5;
      kd = 0.5;
      ks = 0.3;
      specExp = 30.0;

      // normal vector
      N = new float[]{intersectionR[0] - 15.0, intersectionR[1] - 60.0, intersectionR[2] - 5.0};
      // normalize N
      N = normalize(N);

      // destination (to the camera)
      V = new float[]{point[0] - intersectionR[0], point[1] - intersectionR[1], point[2] - intersectionR[2]};
      // normalize V
      V = normalize(V);

      ambient = new float[]{250.0 * ka, 250.0 * ka, 250.0 * ka};
      c[0] += ambient[0];
      c[1] += ambient[1];
      c[2] += ambient[2];
      for (int i = 0; i < lights.size(); i++) {
        // direction of incoming light
        S = new float[]{lights.get(i)[0] - intersectionR[0], lights.get(i)[1] - intersectionR[1], lights.get(i)[2] - intersectionR[2]};
        // normalize S
        S = normalize(S);

        // reflection vector
        dotTemp = S[0] * N[0] + S[1] * N[1] + S[2] * N[2];
        R = new float[]{-S[0] + 2 * N[0] * dotTemp, -S[1] + 2 * N[1] * dotTemp, -S[2] + 2 * N[2] * dotTemp};

        dotTemp = N[0] * S[0] + N[1] * S[1] + N[2] * S[2];
        if (dotTemp >= 0) diffuse = new float[]{250.0 * dotTemp * kd, 250.0 * dotTemp * kd, 250.0 * dotTemp * kd};
        else diffuse = new float[]{0.0, 0.0, 0.0};
        dotTemp = R[0] * V[0] + R[1] * V[1] + R[2] * V[2];
        if (dotTemp >= 0) specDot = pow(dotTemp, specExp);
        else specDot = 0.0;
        specular = new float[]{255.0 * specDot * ks, 255.0 * specDot * ks, 255.0 * specDot * ks};
        if (!shadowsR.get(i)) {
          c[0] += diffuse[0];
          c[1] += diffuse[1];
          c[2] += diffuse[2];
          c[0] += specular[0];
          c[1] += specular[1];
          c[2] += specular[2];
        }
      }
      // calculating reflection
      dotTemp = V[0] * N[0] + V[1] * N[1] + V[2] * N[2];
      RR = new float[]{-V[0] + 2 * N[0] * dotTemp, -V[1] + 2 * N[1] * dotTemp, -V[2] + 2 * N[2] * dotTemp};
      cr = reflect(depth, 0.3, intersectionR, RR);
      c[0] += cr[0];
      c[1] += cr[1];
      c[2] += cr[2];
      break;
    case 3:
      ka = 0.5;
      kd = 0.5;
      ks = 0.3;
      specExp = 30.0;

      // normal vector
      N = new float[]{0.0, 1.0, 0.0};
      // normalize N
      N = normalize(N);

      // destination (to the camera)
      V = new float[]{point[0] - intersectionR[0], point[1] - intersectionR[1], point[2] - intersectionR[2]};
      // normalize V
      V = normalize(V);

      u = (intersectionR[0] - (-200.0)) / (200.0 - (-200.0));
      v = (intersectionR[2] - (-250.0)) / (250.0 - (-250.0));

      // checkerboard size in uv
      dx = procedural / 400.0;
      dz = procedural / 500.0;
      if ((ceil(u / dx) + ceil(v / dz)) % 2 == 0) ambient = new float[]{255.0 * ka, 255.0 * ka, 0.0 * ka};
      else ambient = new float[]{255.0 * ka, 0.0 * ka, 0.0 * ka};

      //ambient = new float[]{255.0 * ka, 228.0 * ka, 108.0 * ka};
      c[0] += ambient[0];
      c[1] += ambient[1];
      c[2] += ambient[2];

      for (int i = 0; i < lights.size(); i++) {
        // direction of incoming light
        S = new float[]{lights.get(i)[0] - intersectionR[0], lights.get(i)[1] - intersectionR[1], lights.get(i)[2] - intersectionR[2]};
        // normalize S
        S = normalize(S);

        // reflection vector
        dotTemp = S[0] * N[0] + S[1] * N[1] + S[2] * N[2];
        R = new float[]{-S[0] + 2 * N[0] * dotTemp, -S[1] + 2 * N[1] * dotTemp, -S[2] + 2 * N[2] * dotTemp};

        dotTemp = N[0] * S[0] + N[1] * S[1] + N[2] * S[2];
        if (dotTemp >= 0) {
          if ((ceil(u / dx) + ceil(v / dz)) % 2 == 0) diffuse = new float[]{255.0 * dotTemp * kd, 255.0 * dotTemp * kd, 0.0 * dotTemp * kd};
          else diffuse = new float[]{255.0 * dotTemp * kd, 0.0 * dotTemp * kd, 0.0 * dotTemp * kd};
        } else diffuse = new float[]{0.0, 0.0, 0.0};
        dotTemp = R[0] * V[0] + R[1] * V[1] + R[2] * V[2];
        if (dotTemp >= 0) specDot = pow(dotTemp, specExp);
        else specDot = 0.0;
        specular = new float[]{255.0 * specDot * ks, 255.0 * specDot * ks, 255.0 * specDot * ks};
        if (!shadowsR.get(i)) {
          c[0] += diffuse[0];
          c[1] += diffuse[1];
          c[2] += diffuse[2];
          c[0] += specular[0];
          c[1] += specular[1];
          c[2] += specular[2];
        }
      }
      // calculating reflection
      dotTemp = V[0] * N[0] + V[1] * N[1] + V[2] * N[2];
      RR = new float[]{-V[0] + 2 * N[0] * dotTemp, -V[1] + 2 * N[1] * dotTemp, -V[2] + 2 * N[2] * dotTemp};
      cr = reflect(depth, 0.0, intersectionR, RR);
      c[0] += cr[0];
      c[1] += cr[1];
      c[2] += cr[2];
      break;
    }
    result[0] = kr * c[0];
    result[1] = kr * c[1];
    result[2] = kr * c[2];
  }
  return result;
}
