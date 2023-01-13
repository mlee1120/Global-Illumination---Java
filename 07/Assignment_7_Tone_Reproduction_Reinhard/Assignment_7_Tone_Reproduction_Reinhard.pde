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
int supersample = 1;

// size of the tiles
int procedural = 25;

// to update the canvas or not
boolean toDraw;

// max reflection time
int max_depth = 5;

float softshadow = 0.72;

boolean inside = false;

String range = "high";

// this function is executed once for setting up some parameters
void setup() {
  // canvas size
  size(560, 315);
  //size(1280, 720);
  // lights' positions
  lights = new ArrayList();
  lights.add(new float[]{-100.0, 250.0, 300.0});
  //lights.add(new float[]{300.0, 1000.0, 0.0});

  // camera position
  eye = new float[]{-30.0, 150.0, 275.0};
  //eye = new float[]{30.0, 150.0, 275.0};

  // camera look at
  lookAt = new float[]{-40.0, 75.0, 100.0};
  //lookAt = new float[]{-55.0, 0.0, 100.0};
  //lookAt = new float[]{0.0, 0.0, 0.0};

  // create all primitives
  objects = new ArrayList();
  objects.add(new Primitive("background", null, null, new float[]{69.0, 150.0, 243.0}, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0));
  objects.add(new Primitive("sphere", new float[]{-55.0, 100.0, 100.0}, new float[]{60.0}, new float[]{150.0, 150.0, 150.0}, 0.5, 0.6, 0.3, 30.0, 0.0, 0.6, 1.08));
  objects.add(new Primitive("sphere", new float[]{15.0, 60.0, 5.0}, new float[]{45.0}, new float[]{250.0, 250.0, 250.0}, 0.5, 0.5, 0.3, 30.0, 0.3, 0.0, 1.0));
  objects.add(new Primitive("plane", new float[]{0.0, 0.0, 0.0}, new float[]{400.0, 0.0, 500.0, 1.0}, new float[]{255.0, 255.0, 0.0, 255.0, 0.0, 0.0}, 0.5, 0.5, 0.3, 30.0, 0.0, 0.0, 1.0));

  // finding the near clipping plane (upper left pixel position)
  findNear();
  buffer = new ArrayList();
  for (int i = 0; i < width; i++) {
    buffer.add(new ArrayList());
    for (int j = 0; j < height; j++) {
      buffer.get(i).add(new float[4]);
    }
  }
  toDraw = true;
}

void draw() {
  if (toDraw) {
    float[] pixel_position = new float[3];
    float[] ray_direction = new float[3];
    float dx, dy;
    float[] temp, intersection;
    float[] rgbi;
    Primitive p;
    float max_rgb = 0.0;
    float max_intensity = 0.0;
    float logmean_intensity = 0.0;
    Object[] shadows;
    Object[] temp2;
    // spawn rays
    for (int i = 0; i < width; i++) {
      for (int j = 0; j < height; j++) {
        rgbi = new float[]{0.0, 0.0, 0.0, 0.0};
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
          rgbi[0] += temp[0] / supersample;
          rgbi[1] += temp[1] / supersample;
          rgbi[2] += temp[2] / supersample;
        }
        rgbi[3] = 0.27 * rgbi[0] + 0.67 * rgbi[1] + 0.06 * rgbi[2];
        buffer.get(i).set(j, rgbi);
        if (rgbi[0] > max_rgb) max_rgb = rgbi[0];
        if (rgbi[1] > max_rgb) max_rgb = rgbi[1];
        if (rgbi[2] > max_rgb) max_rgb = rgbi[2];
        if (rgbi[3] > max_intensity) max_intensity = rgbi[3];
        logmean_intensity += log(1.0 + rgbi[3]);
      }
    }
    logmean_intensity = exp(logmean_intensity / width / height);
    println(logmean_intensity);
    println(max_intensity);
    
    float Ldmax = 255.0;
    float upperbound, lowerbound;
    if (range.equals("low")) {
      upperbound = max_rgb * 0.7;
      lowerbound = max_rgb * 0.5;
    } else if (range.equals("mid")) {
      upperbound = max_rgb * 0.8;
      lowerbound = max_rgb * 0.2;
    } else {
      upperbound = max_rgb * 1.0;
      lowerbound = max_rgb * 0.0;
    }
    for (int i = 0; i < width; i++) {
      for (int j = 0; j < height; j++) {
        buffer.get(i).get(j)[0] = (buffer.get(i).get(j)[0] - lowerbound) / (upperbound - lowerbound) * Ldmax;
        buffer.get(i).get(j)[1] = (buffer.get(i).get(j)[1] - lowerbound) / (upperbound - lowerbound) * Ldmax;
        buffer.get(i).get(j)[2] = (buffer.get(i).get(j)[2] - lowerbound) / (upperbound - lowerbound) * Ldmax;
      }
    }
    
    float a = 0.18;
    ArrayList<ArrayList<float[]>> scaled_luminance = new ArrayList();
    for (int i = 0; i < width; i++) {
      scaled_luminance.add(new ArrayList());
      for (int j = 0; j < height; j++) {
        scaled_luminance.get(i).add(new float[]{buffer.get(i).get(j)[0] * a / logmean_intensity, buffer.get(i).get(j)[1] * a / logmean_intensity, buffer.get(i).get(j)[2] * a / logmean_intensity});
      }
    }

    Ldmax *= 3.6;
    
    for (int i = 0; i < width; i++) {
      for (int j = 0; j < height; j++) {
        set(i, j, color(Ldmax * scaled_luminance.get(i).get(j)[0] / (1 + scaled_luminance.get(i).get(j)[0]), Ldmax * scaled_luminance.get(i).get(j)[1] / (1 + scaled_luminance.get(i).get(j)[1]), Ldmax * scaled_luminance.get(i).get(j)[2] / (1 + scaled_luminance.get(i).get(j)[2])));
      }
    }
    //save("Assignment_7_Tone_Reproduction_Reinhard_High_II.png");
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
  } else if (p.type.equals("plane")) {
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
  float[] S, N, V, R;
  float dotTemp, specDot;
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
        dotTemp = dotProduct(S, N);
        if (dotTemp >= 0) diffuse = new float[]{p.rgb[0] * dotTemp * p.kd, p.rgb[1] * dotTemp * p.kd, p.rgb[2] * dotTemp * p.kd};
        else diffuse = new float[]{0.0, 0.0, 0.0};
        dotTemp = dotProduct(R, V);
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
  } else if (p.type.equals("plane")) {
    float u, v, dxyz1, dxyz2;

    if (p.geometry[0] == 0) {
      u = (intersection[2] - (-p.geometry[2] / 2)) / (p.geometry[2] / 2 - (-p.geometry[2] / 2));
      v = (intersection[1] - (-p.geometry[1] / 2)) / (p.geometry[1] / 2 - (-p.geometry[1] / 2));
      // checkerboard size in uv
      dxyz1 = procedural / p.geometry[2];
      dxyz2 = procedural / p.geometry[1];
    } else if (p.geometry[1] == 0) {
      u = (intersection[0] - (-p.geometry[0] / 2)) / (p.geometry[0] / 2 - (-p.geometry[0] / 2));
      v = (intersection[2] - (-p.geometry[2] / 2)) / (p.geometry[2] / 2 - (-p.geometry[2] / 2));
      // checkerboard size in uv
      dxyz1 = procedural / p.geometry[0];
      dxyz2 = procedural / p.geometry[2];
    } else {
      u = (intersection[0] - (-p.geometry[0] / 2)) / (p.geometry[0] / 2 - (-p.geometry[0] / 2));
      v = (intersection[1] - (-p.geometry[1] / 2)) / (p.geometry[1] / 2 - (-p.geometry[1] / 2));
      // checkerboard size in uv
      dxyz1 = procedural / p.geometry[0];
      dxyz2 = procedural / p.geometry[1];
    }


    if ((ceil(u / dxyz1) + ceil(v / dxyz2)) % 2 == 0) ambient = new float[]{p.rgb[0] * p.ka, p.rgb[1] * p.ka, p.rgb[2] * p.ka};
    else ambient = new float[]{p.rgb[3] * p.ka, p.rgb[4] * p.ka, p.rgb[5] * p.ka};
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
        dotTemp = dotProduct(S, N);
        if (dotTemp >= 0) {
          if ((ceil(u / dxyz1) + ceil(v / dxyz2)) % 2 == 0) diffuse = new float[]{p.rgb[0] * dotTemp * p.kd, p.rgb[1] * dotTemp * p.kd, p.rgb[2] * dotTemp * p.kd};
          else diffuse = new float[]{p.rgb[3] * dotTemp * p.kd, p.rgb[4] * dotTemp * p.kd, p.rgb[5] * dotTemp * p.kd};
        } else diffuse = new float[]{0.0, 0.0, 0.0};
        dotTemp = dotProduct(R, V);
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
  }
  //if (final_color[0] > 255.0) final_color[0] = 255.0;
  //if (final_color[1] > 255.0) final_color[1] = 255.0;
  //if (final_color[2] > 255.0) final_color[2] = 255.0;
  return final_color;
}

float[] reflect(float[] origin, float[] ray) {
  origin[0] += 2*ray[0];
  origin[1] += 2*ray[1];
  origin[2] += 2*ray[2];

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
  origin[0] += 2*ray[0];
  origin[1] += 2*ray[1];
  origin[2] += 2*ray[2];

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
