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

// buffer
ArrayList<ArrayList<float[]>> buffer;

// number of supersampling
int supersample = 100;

// size of the tiles
int procedural = 25;

// to update the canvas or not
boolean toDraw;

// max reflection time
int max_depth = 5;

float softshadow = 0.72;

boolean inside = false;

boolean blinn = true;

float photon_strength = 0.06;
int photon = 100000;
int bounce = 72;
// KD tree
KdTree tree;
HashMap<String, Integer> photons;
ArrayList<float[]> photons2;
float start, end;


// this function is executed once for setting up some parameters
void setup() {
  // canvas size
  size(600, 600);
  //size(900, 900);

  // lights' positions
  lights = new ArrayList();
  lights.add(new float[]{275.0, 549.9, 280.0});
  //lights.add(new float[]{210.0, 549.9, 190.0});
  //lights.add(new float[]{340.0, 549.9, 190.0});
  //lights.add(new float[]{340.0, 549.9, 370.0});
  //lights.add(new float[]{210.0, 549.9, 370.0});


  // camera position
  eye = new float[]{275.0, 275.0, -477.0};

  // camera look at
  lookAt = new float[]{275.0, 275.0, 0.0};

  // create all primitives
  objects = new ArrayList();
  objects.add(new Primitive("background", null, null, new float[]{0.0, 0.0, 0.0}, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0));
  // floor
  objects.add(new Primitive("plane", new float[]{275.0, 0.0, 280.0}, new float[]{550.0, 0.0, 560.0, 1.0}, new float[]{255.0, 255.0, 255.0}, 0.5, 0.5, 0.3, 30.0, 0.0, 0.0, 1.0));
  // ceiling
  objects.add(new Primitive("plane", new float[]{275.0, 550.0, 280.0}, new float[]{550.0, 0.0, 560.0, -1.0}, new float[]{255.0, 255.0, 255.0}, 0.5, 0.5, 0.3, 30.0, 0.0, 0.0, 1.0));
  // back wall
  objects.add(new Primitive("plane", new float[]{275.0, 275.0, 560.0}, new float[]{550.0, 550.0, 0.0, -1.0}, new float[]{255.0, 255.0, 255.0}, 0.5, 0.5, 0.3, 30.0, 0.0, 0.0, 1.0));
  // right wall
  objects.add(new Primitive("plane", new float[]{0.0, 275.0, 280.0}, new float[]{0.0, 550.0, 560.0, 1.0}, new float[]{0.0, 255.0, 0.0}, 0.5, 0.5, 0.3, 30.0, 0.0, 0.0, 1.0));
  // left wall
  objects.add(new Primitive("plane", new float[]{550.0, 275.0, 280.0}, new float[]{0.0, 550.0, 560.0, -1.0}, new float[]{255.0, 0.0, 0.0}, 0.5, 0.5, 0.3, 30.0, 0.0, 0.0, 1.0));
  // refracted ball
  objects.add(new Primitive("sphere", new float[]{162.0, 90.0, 200.0}, new float[]{90.0}, new float[]{15.0, 15.0, 15.0}, 0.5, 0.6, 0.3, 30.0, 0.0, 1.0, 1.33));
  // reflected ball
  objects.add(new Primitive("sphere", new float[]{400.0, 90.0, 360.0}, new float[]{90.0}, new float[]{15.0, 15.0, 15.0}, 0.5, 0.6, 0.3, 30.0, 0.6, 0.0, 1.08));
  // rectangular light source
  objects.add(new Primitive("light", new float[]{275.0, 549.95, 280.0}, new float[]{130.0, 0.0, 90.0, -1.0}, new float[]{255.0, 0.0, 0.0}, 0.5, 0.5, 0.3, 30.0, 0.0, 0.0, 1.0));

  // finding the near clipping plane (upper left pixel position)
  findNear();
  buffer = new ArrayList();
  for (int i = 0; i < width; i++) {
    buffer.add(new ArrayList());
    for (int j = 0; j < height; j++) {
      buffer.get(i).add(new float[3]);
    }
  }
  toDraw = true;

  // photon mapping
  photons = new HashMap();
  photons2 = new ArrayList();
  start = millis();
  scatter();
  end = millis();
  println("scatter: " + ((end - start) / 1000) + " s");
  start = millis();
  tree = new KdTree();
  tree.root = new Node();
  buildTree(tree.root, photons2, 0, photons2.size() - 1, "x");
  end = millis();
  println("tree: " + ((end - start) / 1000) + " s");
}

void draw() {
  if (toDraw) {
    float[] pixel_position = new float[3];
    float[] ray_direction = new float[3];
    float dx, dy;
    float[] temp, intersection;
    float[] rgb;
    Primitive p;
    float max_color = 0;
    Object[] shadows;
    Object[] temp2;
    // spawn rays
    start = millis();
    for (int i = 0; i < width; i++) {
      for (int j = 0; j < height; j++) {
        rgb = new float[]{0.0, 0.0, 0.0};
        // multisampling
        for (int k = 0; k < supersample; k++) {
          if (k == 0) {
            dx = 0.0;
            dy = 0.0;
          } else {
            dx = random(-0.5, 0.5);
            dy = random(-0.5, 0.5);
          }
          pixel_position[0] = upper_left[0] + i * vector_right[0] - j * vector_up[0] + dx * vector_right[0] - dy * vector_up[0];
          pixel_position[1] = upper_left[1] + i * vector_right[1] - j * vector_up[1] + dx * vector_right[1] - dy * vector_up[1];
          pixel_position[2] = upper_left[2] + i * vector_right[2] - j * vector_up[2] + dx * vector_right[2] - dy * vector_up[2];
          ray_direction[0] = pixel_position[0] - eye[0];
          ray_direction[1] = pixel_position[1] - eye[1];
          ray_direction[2] = pixel_position[2] - eye[2];

          temp2 = check_intersection(pixel_position, ray_direction);
          intersection = (float[]) temp2[0];
          p = (Primitive) temp2[1];

          // check shadow
          shadows = check_shadows(intersection);
          inside = false;
          temp = calculate_color(eye, p, intersection, shadows);
          rgb[0] += temp[0];
          rgb[1] += temp[1];
          rgb[2] += temp[2];
        }
        buffer.get(i).set(j, rgb);
        if (rgb[0] / supersample > max_color) max_color = rgb[0] / supersample;
        if (rgb[1] / supersample > max_color) max_color = rgb[1] / supersample;
        if (rgb[2] / supersample > max_color) max_color = rgb[2] / supersample;
      }
    }
    //max_color = 255;
    println("max color: " + max_color);
    for (int i = 0; i < width; i++) {
      for (int j = 0; j < height; j++) {
        set(i, j, color(buffer.get(i).get(j)[0] / supersample / max_color * 255, buffer.get(i).get(j)[1] / supersample / max_color * 255, buffer.get(i).get(j)[2] / supersample / max_color * 255));
      }
    }
    end = millis();
    println("render: " + ((end - start) / 1000) + " s");
    //save("Assignment_Advanced_Photon_Mapping_100000.png");
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

/**
 * This function calculates and returns the dot product of two given vectors.
 *
 * @param v1 - the first vector
 * @param v2 - the other vector
 * @return   - the dot product
 */
float dotProduct(float[] v1, float[] v2) {
  return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

/**
 * This function calculates and returns the reflected vector of a given vector according to a normal vector.
 *
 * @param v - the given vector
 * @param n - the normal vector
 * @return  - the reflected vector
 */
float[] reflection(float[] v, float[] n) {
  float dotTemp = dotProduct(v, n);
  return new float[]{-v[0] + 2 * n[0] * dotTemp, -v[1] + 2 * n[1] * dotTemp, -v[2] + 2 * n[2] * dotTemp};
}

float[] refraction(float[] v, float[] n, float n_in, float n_re) {
  float theta_in = acos(dotProduct(v, n));
  if (theta_in == PI / 2.0) return v;
  else {
    // refraction
    if (n_in < n_re || theta_in <= asin(n_re / n_in)) {
      float theta_re = asin(n_in * sin(theta_in) / n_re);
      float[] axis = normalize(crossProduct(n, v));
      Quaternion q = new Quaternion(cos((theta_re + PI) / 2.0), axis[0] * sin((theta_re + PI) / 2.0), axis[1] * sin((theta_re + PI) / 2.0), axis[2] * sin((theta_re + PI) / 2.0));
      Quaternion qp = new Quaternion(0.0, n[0], n[1], n[2]);
      Quaternion new_qp = q.product(qp).product(q.inverse());
      inside = !inside;
      return new float[]{new_qp.x, new_qp.y, new_qp.z};
    }
    // total reflection
    else return null;
  }
}

/**
 * This function calculates and returns the distance between two given points.
 *
 * @param v1 - the first point
 * @param v2 - the other point
 * @return   - the distance
 */
float calculate_distance(float[] v1, float[] v2) {
  return sqrt(pow(v1[0] - v2[0], 2) + pow(v1[1] - v2[1], 2) + pow(v1[2] - v2[2], 2));
}


/**
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
 * This function calculates and returns the remainder of one number divided by another.
 *
 * @param f1 the dividend
 * @param f2 the divisor
 * @return the remainder
 */
float mod(float f1, float f2) {
  return f1 - f2 * floor(f1 / f2);
}

/**
 * This function perfroms raytracing to find the intersection.
 * sphere1: (x+55)^2 + (y-100)^2 + (z-100)^2 = 60^2
 * sphere2: (x-15)^2 + (y-60)^2 + (z-5)^2 = 45^2
 * floor: y = 0 with -200 <= x <= 200 && -250 <= z <= 250
 *
 * @param dx - varaible for translating the pixel position
 * @param dy - varaible for translating the pixel position
 */
float[] raytrace(float[] origin, float[] ray, Primitive p) {
  float[] intersection = new float[]{1000, 1000, 1000};
  // coeffecients for finding intersections
  float A, B, C, D;
  // auxiliary parameter for finding intersections
  float coefficient;

  if (p.type.equals("sphere")) {
    A = ray[0] * ray[0] + ray[1] * ray[1] + ray[2] * ray[2];
    B = 2 * (origin[0] - p.center[0]) * ray[0] + 2 * (origin[1] - p.center[1]) * ray[1] + 2 * (origin[2] - p.center[2]) * ray[2];
    C = pow(origin[0] - p.center[0], 2) + pow(origin[1] - p.center[1], 2) + pow(origin[2] - p.center[2], 2) - pow(p.geometry[0], 2);
    D = B * B - 4 * A * C;
    if (D >= 0) {
      coefficient = (-B - sqrt(D)) / (2 * A);
      if (coefficient <= 0) coefficient = (-B + sqrt(D)) / (2 * A);
      if (coefficient > 0) {
        intersection[0] = origin[0] + ray[0] * coefficient;
        intersection[1] = origin[1] + ray[1] * coefficient;
        intersection[2] = origin[2] + ray[2] * coefficient;
      }
    }
  } else if (p.type.equals("plane") || p.type.equals("light")) {
    if (p.geometry[0] == 0.0) {
      coefficient = (p.center[0] - origin[0]) / ray[0];
      if (coefficient > 0) {
        float[] temp = new float[3];
        temp[0] = origin[0] + ray[0] * coefficient;
        temp[1] = origin[1] + ray[1] * coefficient;
        temp[2] = origin[2] + ray[2] * coefficient;
        if (temp[1] <= p.geometry[1] / 2 + p.center[1] && temp[1] >= -p.geometry[1] / 2 + p.center[1] && temp[2] <= p.geometry[2] / 2 + p.center[2] && temp[2] >= -p.geometry[2] / 2 + p.center[2]) {
          intersection[0] = origin[0] + ray[0] * coefficient;
          intersection[1] = origin[1] + ray[1] * coefficient;
          intersection[2] = origin[2] + ray[2] * coefficient;
        }
      }
    } else if (p.geometry[1] == 0.0) {
      coefficient = (p.center[1] - origin[1]) / ray[1];
      if (coefficient > 0) {
        float[] temp = new float[3];
        temp[0] = origin[0] + ray[0] * coefficient;
        temp[1] = origin[1] + ray[1] * coefficient;
        temp[2] = origin[2] + ray[2] * coefficient;
        if (temp[0] <= p.geometry[0] / 2 + p.center[0] && temp[0] >= -p.geometry[0] / 2 + p.center[0] && temp[2] <= p.geometry[2] / 2  + p.center[2] && temp[2] >= -p.geometry[2] / 2 + p.center[2]) {
          intersection[0] = origin[0] + ray[0] * coefficient;
          intersection[1] = origin[1] + ray[1] * coefficient;
          intersection[2] = origin[2] + ray[2] * coefficient;
        }
      }
    } else {
      coefficient = (p.center[2] - origin[2]) / ray[2];
      if (coefficient > 0) {
        float[] temp = new float[3];
        temp[0] = origin[0] + ray[0] * coefficient;
        temp[1] = origin[1] + ray[1] * coefficient;
        temp[2] = origin[2] + ray[2] * coefficient;
        if (temp[0] <= p.geometry[0] / 2 + p.center[0] && temp[0] >= -p.geometry[0] / 2 + p.center[0] && temp[1] <= p.geometry[1] / 2  + p.center[1] && temp[1] >= -p.geometry[1] / 2 + p.center[1]) {
          intersection[0] = origin[0] + ray[0] * coefficient;
          intersection[1] = origin[1] + ray[1] * coefficient;
          intersection[2] = origin[2] + ray[2] * coefficient;
        }
      }
    }
  }
  return intersection;
}

float[] calculate_color(float[] origin, Primitive p, float[] intersection, Object[] shadow) {
  ArrayList<Boolean> shadows = (ArrayList<Boolean>) shadow[0];
  ArrayList<Boolean> softshadows = (ArrayList<Boolean>) shadow[1];
  float[] final_color = new float[]{0.0, 0.0, 0.0};
  float[] ambient, diffuse, specular;
  float[] S, N, V, R, H;
  float dotTemp, specDot;
  int ph = searchTree(tree.root, intersection);
  //int ph = 0;
  if (p.type.equals("background")) {
    final_color = p.rgb;
  } else if (p.type.equals("sphere")) {
    ambient = new float[]{p.rgb[0] * p.ka, p.rgb[1] * p.ka, p.rgb[2] * p.ka};
    final_color[0] += ambient[0];
    final_color[1] += ambient[1];
    final_color[2] += ambient[2];
    N = normalize(new float[]{intersection[0] - p.center[0], intersection[1] - p.center[1], intersection[2] - p.center[2]});
    V = normalize(new float[]{origin[0] - intersection[0], origin[1] - intersection[1], origin[2] - intersection[2]});
    for (int i = 0; i < lights.size(); i++) {
      if (!shadows.get(i) || softshadows.get(i)) {
        S = normalize(new float[]{lights.get(i)[0] - intersection[0], lights.get(i)[1] - intersection[1], lights.get(i)[2] - intersection[2]});
        R = reflection(S, N);
        H = normalize(new float[]{S[0] + N[0], S[1] + N[1], S[2] + N[2]});
        dotTemp = dotProduct(S, N);
        if (dotTemp >= 0) diffuse = new float[]{p.rgb[0] * dotTemp * p.kd, p.rgb[1] * dotTemp * p.kd, p.rgb[2] * dotTemp * p.kd};
        else diffuse = new float[]{0.0, 0.0, 0.0};
        if (blinn) dotTemp = dotProduct(H, N);
        else dotTemp = dotProduct(R, V);
        if (dotTemp >= 0) specDot = pow(dotTemp, p.specExp);
        else specDot = 0.0;
        specular = new float[]{255.0 * specDot * p.ks, 255.0 * specDot * p.ks, 255.0 * specDot * p.ks};
        float transparent = 1.0;
        if (softshadows.get(i)) transparent = softshadow;
        final_color[0] += diffuse[0] * transparent;
        final_color[1] += diffuse[1] * transparent;
        final_color[2] += diffuse[2] * transparent;
        final_color[0] += specular[0] * transparent;
        final_color[1] += specular[1] * transparent;
        final_color[2] += specular[2] * transparent;
      }
    }
    if (p.kr != 0) {
      float[] reflected_color = reflect(intersection, reflection(V, N));
      final_color[0] += p.kr * reflected_color[0];
      final_color[1] += p.kr * reflected_color[1];
      final_color[2] += p.kr * reflected_color[2];
    }
    if (p.kt != 0) {
      if (inside) {
        N[0] = -N[0];
        N[1] = -N[1];
        N[2] = -N[2];
      }

      float[] refracted_color;
      // for checking if total reflection happens
      boolean temp = inside;
      float[] refracted_vector;
      if (inside) refracted_vector = refraction(V, N, p.n, 1.0);
      else refracted_vector = refraction(V, N, 1.0, p.n);
      if (temp != inside) {
        V = normalize(new float[]{origin[0] - intersection[0], origin[1] - intersection[1], origin[2] - intersection[2]});
        refracted_color = refract(intersection, refracted_vector);
        /*final_color[0] = (1.0 - p.kt) * final_color[0] + p.kt * refracted_color[0];
         final_color[1] = (1.0 - p.kt) * final_color[1] + p.kt * refracted_color[1];
         final_color[2] = (1.0 - p.kt) * final_color[2] + p.kt * refracted_color[2];*/
        /*final_color[0] = p.kt * refracted_color[0];
         final_color[1] = p.kt * refracted_color[1];
         final_color[2] = p.kt * refracted_color[2];*/
        final_color[0] += p.kt * refracted_color[0];
        final_color[1] += p.kt * refracted_color[1];
        final_color[2] += p.kt * refracted_color[2];
      }
    }
    final_color[0] *= 1.0 + ph * photon_strength;
    final_color[1] *= 1.0 + ph * photon_strength;
    final_color[2] *= 1.0 + ph * photon_strength;
  } else if (p.type.equals("plane")) {
    ambient = new float[]{p.rgb[0] * p.ka, p.rgb[1] * p.ka, p.rgb[2] * p.ka};
    final_color[0] += ambient[0];
    final_color[1] += ambient[1];
    final_color[2] += ambient[2];
    V = normalize(new float[]{origin[0] - intersection[0], origin[1] - intersection[1], origin[2] - intersection[2]});
    if (p.geometry[0] == 0.0) N = normalize(new float[]{p.geometry[3], 0.0, 0.0});
    else if (p.geometry[1] == 0.0) N = normalize(new float[]{0.0, p.geometry[3], 0.0});
    else N = normalize(new float[]{0.0, 0.0, p.geometry[3]});
    for (int i = 0; i < lights.size(); i++) {
      if (!shadows.get(i) || softshadows.get(i)) {
        S = normalize(new float[]{lights.get(i)[0] - intersection[0], lights.get(i)[1] - intersection[1], lights.get(i)[2] - intersection[2]});
        R = reflection(S, N);
        H = normalize(new float[]{S[0] + N[0], S[1] + N[1], S[2] + N[2]});
        dotTemp = dotProduct(S, N);
        if (dotTemp >= 0) diffuse = new float[]{p.rgb[0] * dotTemp * p.kd, p.rgb[1] * dotTemp * p.kd, p.rgb[2] * dotTemp * p.kd};
        else diffuse = new float[]{0.0, 0.0, 0.0};
        if (blinn) dotTemp = dotProduct(H, N);
        else dotTemp = dotProduct(R, V);
        if (dotTemp >= 0) specDot = pow(dotTemp, p.specExp);
        else specDot = 0.0;
        specular = new float[]{255.0 * specDot * p.ks, 255.0 * specDot * p.ks, 255.0 * specDot * p.ks};
        float transparent = 1.0;
        if (softshadows.get(i)) transparent = softshadow;
        final_color[0] += diffuse[0] * transparent;
        final_color[1] += diffuse[1] * transparent;
        final_color[2] += diffuse[2] * transparent;
        final_color[0] += specular[0] * transparent;
        final_color[1] += specular[1] * transparent;
        final_color[2] += specular[2] * transparent;
      }
    }
    if (p.kr != 0) {
      float[] reflected_color = reflect(intersection, reflection(V, N));
      final_color[0] += p.kr * reflected_color[0];
      final_color[1] += p.kr * reflected_color[1];
      final_color[2] += p.kr * reflected_color[2];
    }
    final_color[0] *= 1.0 + ph * photon_strength;
    final_color[1] *= 1.0 + ph * photon_strength;
    final_color[2] *= 1.0 + ph * photon_strength;
  } else if (p.type.equals("light")) {
    final_color[0] = 450;
    final_color[1] = 450;
    final_color[2] = 450;
  }
  //if (final_color[0] > 255.0) final_color[0] = 255.0;
  //if (final_color[1] > 255.0) final_color[1] = 255.0;
  //if (final_color[2] > 255.0) final_color[2] = 255.0;
  //final_color[0] = final_color[0] * ph / 255;
  //final_color[1] = final_color[1] * ph / 255;
  //final_color[2] = final_color[2] * ph / 255;
  return final_color;
}

float[] reflect(float[] origin, float[] ray) {
  origin[0] += ray[0];
  origin[1] += ray[1];
  origin[2] += ray[2];

  Object[] temp1 = check_intersection(origin, ray);
  float[] intersection = (float[]) temp1[0];
  Primitive p = (Primitive) temp1[1];

  Object[] shadows = check_shadows(intersection);

  float[] temp2 = calculate_color(origin, p, intersection, shadows);
  float[] reflected_color = new float[]{0.0, 0.0, 0.0};
  reflected_color[0] += temp2[0];
  reflected_color[1] += temp2[1];
  reflected_color[2] += temp2[2];
  return reflected_color;
}

float[] refract(float[] origin, float[] ray) {
  origin[0] += ray[0];
  origin[1] += ray[1];
  origin[2] += ray[2];

  Object[] temp1 = check_intersection(origin, ray);
  float[] intersection = (float[]) temp1[0];
  Primitive p = (Primitive) temp1[1];

  Object[] shadows = check_shadows(intersection);

  float[] temp2 = calculate_color(origin, p, intersection, shadows);
  float[] refracted_color = new float[]{0.0, 0.0, 0.0};
  refracted_color[0] += temp2[0];
  refracted_color[1] += temp2[1];
  refracted_color[2] += temp2[2];
  return refracted_color;
}

Object[] check_shadows(float[] intersection) {
  ArrayList<Boolean> shadows = new ArrayList();
  ArrayList<Boolean> transparent = new ArrayList();
  float[] light_direction = new float[3];
  float[] intersection_shift = new float[3];
  float[] temp1;
  for (int l = 0; l < lights.size(); l++) {
    if (intersection[0] != 1000) {
      light_direction[0] = lights.get(l)[0] - intersection[0];
      light_direction[1] = lights.get(l)[1] - intersection[1];
      light_direction[2] = lights.get(l)[2] - intersection[2];
      light_direction = normalize(light_direction);
      intersection_shift = new float[3];
      intersection_shift[0] = intersection[0] + light_direction[0];
      intersection_shift[1] = intersection[1] + light_direction[1];
      intersection_shift[2] = intersection[2] + light_direction[2];
      shadows.add(false);
      transparent.add(false);
      for (int m = 1; m < objects.size(); m++) {
        temp1 = raytrace(intersection_shift, light_direction, objects.get(m));
        if (temp1[0] != 1000 && calculate_distance(intersection, temp1) < calculate_distance(intersection, lights.get(l))) {
          if (objects.get(m).kt != 0.0) transparent.set(transparent.size() - 1, true);
          else transparent.set(transparent.size() - 1, false);
          shadows.set(shadows.size() - 1, true);
        }
      }
    } else shadows.add(false);
  }
  return new Object[]{shadows, transparent};
}

Object[] check_intersection(float[] origin, float[] ray) {
  float[] intersection = new float[]{1000, 1000, 1000};
  float[] temp;
  Primitive p = objects.get(0);
  float distance1, distance2;

  for (int l = 1; l < objects.size(); l++) {
    temp = raytrace(origin, ray, objects.get(l));
    distance1 = calculate_distance(origin, intersection);
    distance2 = calculate_distance(origin, temp);
    if (distance1 > distance2) {
      intersection = temp;
      p = objects.get(l);
    }
  }
  return new Object[]{intersection, p};
}

void scatter() {
  float[] origin = new float[3];
  float[] direction = new float[3];
  float[] axis = new float[3];
  float theta, phi;
  Object[] temp1;
  float[] intersection;
  Primitive p;
  float[] N;
  Quaternion q, qp1, qp2;
  float[] refracted_vector;
  boolean temp;
  float[] temp2 = new float[3];
  for (int i = 0; i < photon; i++) {
    origin[0] = random(210.0, 250.0);
    origin[1] = 549.9;
    origin[2] = random(235.0, 235.0);
    theta = random(0.1, 1.0);
    phi = random(0.0, 3.0);
    do {
      axis[0] = random(0.0, 1.0);
      axis[1] = random(0.0, 1.0);
      axis[2] = random(0.0, 1.0);
    } while (axis[0] == 0.0 && axis[1] == 0.0 && axis[2] == 0.0);
    axis = normalize(axis);
    q = new Quaternion(cos(theta * PI / 4.0), sin(theta * PI / 4.0) * axis[0], sin(theta * PI / 4.0) * axis[1], sin(theta * PI / 4.0) * axis[2]);
    qp1 = new Quaternion(0.0, 0.0, -1.0, 0.0);
    qp2 = q.product(qp1).product(q.inverse());
    q = new Quaternion(cos(phi * PI), sin(phi * PI) * 0.0, sin(phi * PI) * 1.0, sin(phi * PI) * 0.0);
    qp1 = new Quaternion(0.0, qp2.x, qp2.y, qp2.z);
    qp2 = q.product(qp1).product(q.inverse());
    direction[0] = qp2.x;
    direction[1] = qp2.y;
    direction[2] = qp2.z;

    for (int j = 0; j < bounce; j++) {
      temp1 = check_intersection(origin, direction);
      intersection = (float[]) temp1[0];
      p = (Primitive) temp1[1];
      if (p.type.equals("background")) break;
      else if (p.type.equals("plane")) {
        temp2[0] = intersection[0] - mod(intersection[0], 0.001);
        temp2[1] = intersection[1] - mod(intersection[1], 0.001);
        temp2[2] = intersection[2] - mod(intersection[2], 0.001);

        if (photons.containsKey(temp2[0] + " " + temp2[1] + " " + temp2[2])) photons.put(temp2[0] + " " + temp2[1] + " " + temp2[2], photons.get(temp2[0] + " " + temp2[1] + " " + temp2[2]) + 1);
        else {
          photons.put(temp2[0] + " " + temp2[1] + " " + temp2[2], 1);
          photons2.add(new float[]{temp2[0], temp2[1], temp2[2]});
        }
        if (p.geometry[0] == 0.0) N = new float[]{p.geometry[3], 0.0, 0.0};
        else if (p.geometry[1] == 0.0) N = new float[]{0.0, p.geometry[3], 0.0};
        else N = new float[]{0.0, 0.0, p.geometry[3]};
        direction = reflection(new float[]{-direction[0], -direction[1], -direction[2]}, N);
      } else if (p.type.equals("sphere")) {
        if (p.kt == 0.0) {
          temp2[0] = intersection[0] - mod(intersection[0], 0.001);
          temp2[1] = intersection[1] - mod(intersection[1], 0.001);
          temp2[2] = intersection[2] - mod(intersection[2], 0.001);
          if (photons.containsKey(temp2[0] + " " + temp2[1] + " " + temp2[2])) photons.put(temp2[0] + " " + temp2[1] + " " + temp2[2], photons.get(temp2[0] + " " + temp2[1] + " " + temp2[2]) + 1);
          else {
            photons.put(temp2[0] + " " + temp2[1] + " " + temp2[2], 1);
            photons2.add(new float[]{temp2[0], temp2[1], temp2[2]});
          }
        }
        if (p.kt != 0.0) {
          N = normalize(new float[]{intersection[0] - p.center[0], intersection[1] - p.center[1], intersection[2] - p.center[2]});
          if (inside) {
            N[0] = -N[0];
            N[1] = -N[1];
            N[2] = -N[2];
          }
          // for checking if total reflection happens
          temp = inside;
          if (inside) refracted_vector = refraction(new float[]{-direction[0], -direction[1], -direction[2]}, N, p.n, 1.0);
          else refracted_vector = refraction(new float[]{-direction[0], -direction[1], -direction[2]}, N, 1.0, p.n);
          if (temp != inside) direction = refracted_vector;
          else direction = reflection(new float[]{-direction[0], -direction[1], -direction[2]}, N);
        } else {
          N = normalize(new float[]{intersection[0] - p.center[0], intersection[1] - p.center[1], intersection[2] - p.center[2]});
          direction = reflection(new float[]{-direction[0], -direction[1], -direction[2]}, N);
        }
      } else if (p.type.equals("light")) {
        N = new float[]{0.0, -1.0, 0.0};
        direction = reflection(new float[]{-direction[0], -direction[1], -direction[2]}, N);
      }
      origin = intersection;
      origin[0] += direction[0];
      origin[1] += direction[1];
      origin[2] += direction[2];
    }
  }
}

ArrayList<float[]> mergeSort(ArrayList<float[]> aList, int start, int end, String xyz) {
  ArrayList<float[]>  sortedList = new ArrayList();
  if (start == end) sortedList.add(aList.get(start));
  else {
    int middle = (start + end) / 2;
    ArrayList<float[]> leftList = mergeSort(aList, start, middle, xyz);
    ArrayList<float[]> rightList = mergeSort(aList, middle + 1, end, xyz);
    int leftIndex = 0, rightIndex = 0;
    int sort;
    if (xyz.equals("x")) sort = 0;
    else if (xyz.equals("y")) sort = 1;
    else sort = 2;
    while (leftIndex < leftList.size() && rightIndex < rightList.size()) {
      if (leftList.get(leftIndex)[sort] < rightList.get(rightIndex)[sort]) {
        sortedList.add(leftList.get(leftIndex));
        leftIndex++;
      } else {
        sortedList.add(rightList.get(rightIndex));
        rightIndex++;
      }
    }
    while (leftIndex < leftList.size()) {
      sortedList.add(leftList.get(leftIndex));
      leftIndex++;
    }
    while (rightIndex < rightList.size()) {
      sortedList.add(rightList.get(rightIndex));
      rightIndex++;
    }
  }
  return sortedList;
}

void buildTree(Node aNode, ArrayList<float[]> aList, int start, int end, String xyz) {
  if (start == end) {
    aNode.name = "yes";
    aNode.position = aList.get(start);
    if (xyz.equals("x")) aNode.sort = 0;
    else if (xyz.equals("y")) aNode.sort = 1;
    else aNode.sort = 2;
    aNode.leftChild = new Node();
    aNode.rightChild = new Node();
    tree.size++;
  } else {
    // sort accordingly
    ArrayList<float[]> tempList = mergeSort(aList, start, end, xyz);
    int middle = (tempList.size() - 1) / 2;
    int sort;
    if (xyz.equals("x")) sort = 0;
    else if (xyz.equals("y")) sort = 1;
    else sort = 2;
    while (middle < tempList.size() - 1) {
      if (tempList.get(middle)[sort] == tempList.get(middle + 1)[sort]) {
        middle += 1;
      } else break;
    }

    aNode.leftChild = new Node();
    aNode.rightChild = new Node();

    if (middle == tempList.size() - 1) {
      aNode.name = "no";
      if (xyz.equals("x")) buildTree(aNode.leftChild, tempList, 0, middle, "y");
      else if (xyz.equals("y")) buildTree(aNode.leftChild, tempList, 0, middle, "z");
      else if (xyz.equals("z")) buildTree(aNode.leftChild, tempList, 0, middle, "x");
    } else if (middle == 0) {
      tree.size++;
      aNode.name = "yes";
      aNode.position = tempList.get(middle);
      if (xyz.equals("x")) aNode.sort = 0;
      else if (xyz.equals("y")) aNode.sort = 1;
      else aNode.sort = 2;
      if (xyz.equals("x")) buildTree(aNode.rightChild, tempList, middle + 1, tempList.size() - 1, "y");
      else if (xyz.equals("y")) buildTree(aNode.rightChild, tempList, middle + 1, tempList.size() - 1, "z");
      else if (xyz.equals("z")) buildTree(aNode.rightChild, tempList, middle + 1, tempList.size() - 1, "x");
    } else {
      tree.size++;
      aNode.name = "yes";
      aNode.position = tempList.get(middle);
      if (xyz.equals("x")) aNode.sort = 0;
      else if (xyz.equals("y")) aNode.sort = 1;
      else aNode.sort = 2;
      if (xyz.equals("x")) {
        buildTree(aNode.leftChild, tempList, 0, middle - 1, "y");
        buildTree(aNode.rightChild, tempList, middle + 1, tempList.size() - 1, "y");
      } else if (xyz.equals("y")) {
        buildTree(aNode.leftChild, tempList, 0, middle - 1, "z");
        buildTree(aNode.rightChild, tempList, middle + 1, tempList.size() - 1, "z");
      } else if (xyz.equals("z")) {
        buildTree(aNode.leftChild, tempList, 0, middle - 1, "x");
        buildTree(aNode.rightChild, tempList, middle + 1, tempList.size() - 1, "x");
      }
    }
  }
}

int searchTree(Node current, float[] intersection) {
  int result = 0;
  float distance = 1.0;
  if (current.name.equals("yes")) {
    if (calculate_distance(current.position, intersection) < distance) result += photons.get(current.position[0] + " " + current.position[1] + " " + current.position[2]);
    if (intersection[current.sort] > current.position[current.sort]) {
      if (calculate_distance(current.leftChild.position, intersection) < distance) return result + searchTree(current.leftChild, intersection) + searchTree(current.rightChild, intersection);
      else return result + searchTree(current.rightChild, intersection);
    } else {
      if (calculate_distance(current.rightChild.position, intersection) < distance) return result + searchTree(current.leftChild, intersection) + searchTree(current.rightChild, intersection);
      else return result + searchTree(current.leftChild, intersection);
    }
  } else if (current.name.equals("no")) {
    return searchTree(current.leftChild, intersection) + searchTree(current.rightChild, intersection);
  } else return result;
}
