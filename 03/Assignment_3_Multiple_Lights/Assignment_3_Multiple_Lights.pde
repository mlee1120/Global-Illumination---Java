size(560, 315);

// lights' positions
float[] light1 = new float[]{-100.0, 250.0, 300.0};
float[] light2 = new float[]{300.0, 1000.0, 0.0};

// finding the near clipping plane (upper left pixel position)
float[] eye = new float[]{-30.0, 150.0, 275.0}; // camera position
float[] lookAt = new float[]{-40.0, 75.0, 100.0}; // camera look at


float fovy = 60.0; // degree
float near = (height/2.0) / tan(fovy / 2.0 * PI / 180.0) / 10.0; // distance from camera to the near clipping plane
float far = (height/2.0) / tan(fovy / 2.0 * PI / 180.0) * 10.0; // distance from camera to the far clipping plane
float vx = lookAt[0] - eye[0], vy = lookAt[1] - eye[1], vz = lookAt[2] - eye[2]; // the vecter from eye to lookAt
float n1 = near / sqrt(vx * vx + vy * vy + vz * vz); // auxiliary variable for finding near center
float[] near_center = new float[]{eye[0] + n1 * vx, eye[1] + n1 * vy, eye[2] + n1 * vz}; // the center of the near clipping plan
float[] vector_right = new float[]{-vz, 0.0, vx}; // the right vector of the near clipping plane (cross product of the vecter from eye to lookAt and the camera up vector)
float[] vector_up = new float[]{vector_right[1] * vz - vector_right[2] * vy, vector_right[2] * vx - vector_right[0] * vz, vector_right[0] * vy - vector_right[1] * vx}; // the up vector of the near clipping plane (cross product of the right vector and the vecter from eye to lookAt)
float n2 = tan(fovy / 2.0 * PI / 180.0) * near / sqrt(vector_up[0] * vector_up[0] + vector_up[1] * vector_up[1] + vector_up[2] * vector_up[2]); // auxiliary variable for finding the top center of the near clipping plan
float[] near_topcenter = new float[]{eye[0] + n1 * vx + n2 * vector_up[0], eye[1] + n1 * vy + n2 * vector_up[1], eye[2] + n1 * vz + n2 * vector_up[2]}; // top center of the near clipping plan
float pixel_size = sqrt((near_topcenter[0] - near_center[0]) * (near_topcenter[0] - near_center[0]) + (near_topcenter[1] - near_center[1]) * (near_topcenter[1] - near_center[1]) + (near_topcenter[2] - near_center[2]) * (near_topcenter[2] - near_center[2])) / (height / 2); // assume height value is odd (if even, we have to modify the code)
float n3 = pixel_size / sqrt(vector_right[0] * vector_right[0] + vector_right[1] * vector_right[1] + vector_right[2] * vector_right[2]); // auxiliary variable for normalizing right vector
vector_right[0] = vector_right[0] * n3; // normalize vector_right (wrt the pixel size)
vector_right[1] = vector_right[1] * n3;
vector_right[2] = vector_right[2] * n3;
float n4 = pixel_size / sqrt(vector_up[0] * vector_up[0] + vector_up[1] * vector_up[1] + vector_up[2] * vector_up[2]); // auxiliary variable for normalizing up vector
vector_up[0] = vector_up[0] * n4; // normalize vector_up (wrt the pixel size)
vector_up[1] = vector_up[1] * n4;
vector_up[2] = vector_up[2] * n4;
float n5 = width / 2 - 0.5; // auxiliary variable for finding upper left pixel (assume width value is even => if odd, we have to modify the code)
float[] upper_left = new float[]{near_topcenter[0] - n5 * vector_right[0], near_topcenter[1] - n5 * vector_right[1], near_topcenter[2] - n5 * vector_right[2]}; // upper left pixel position


float x, y, z; // pixel position
float A, B, C, D; // coeffecients for finding intersections
float[] v = new float[3]; // vector of the ray tracing path
float[] vTemp = new float[3]; // vector for checking shadow
float xx, yy, zz; // floor intersection
int choice; // intersection with which object
float[] intersection = new float[3]; // intersections
float coefficient; // auxiliary parameter for finding intersections
float temp; // for normalization
float ka, kd, ks, specExp; // reflection coefficients
float[] ambient, diffuse1, specular1, c, diffuse2, specular2;
float[] S, N, V, R;
float dotTemp, specDot;
boolean shadow1, shadow2;


/* ray tracing
 sphere1: (x+55)^2 + (y-100)^2 + (z-100)^2 = 60^2
 sphere2: (x-15)^2 + (y-60)^2 + (z-5)^2 = 45^2
 floor: y = 0 with -200 <= x <= 200 && -250 <= z <= 250
 */
for (int i = 0; i < width; i++) {
  for (int j = 0; j < height; j++) {
    choice = 0;
    shadow1 = false;
    shadow2 = false;
    x = upper_left[0] + i * vector_right[0] - j * vector_up[0];
    y = upper_left[1] + i * vector_right[1] - j * vector_up[1];
    z = upper_left[2] + i * vector_right[2] - j * vector_up[2];
    v[0] = x - eye[0];
    v[1] = y - eye[1];
    v[2] = z - eye[2];

    // check sphere1
    A = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
    B = 2 * (x + 55) * v[0] + 2 * (y - 100) * v[1] + 2 * (z - 100) * v[2];
    C = (x + 55) * (x + 55) + (y - 100) * (y - 100) + (z - 100) * (z - 100) - 60 * 60;
    D = B * B - 4 * A * C;
    if (D >= 0) {
      choice = 1;
      coefficient = (-B - sqrt(D)) / (2 * A);
      if (coefficient <= 0) coefficient = (B - sqrt(D)) / (2 * A);
      intersection[0] = x + v[0] * coefficient;
      intersection[1] = y + v[1] * coefficient;
      intersection[2] = z + v[2] * coefficient;

      // check shadow1
      vTemp[0] = light1[0] - intersection[0];
      vTemp[1] = light1[1] - intersection[1];
      vTemp[2] = light1[2] - intersection[2];
      A = vTemp[0] * vTemp[0] + vTemp[1] * vTemp[1] + vTemp[2] * vTemp[2];
      B = 2 * (intersection[0] + 55) * vTemp[0] + 2 * (intersection[1] - 100) * vTemp[1] + 2 * (intersection[2] - 100) * vTemp[2];
      C = (intersection[0] + 55) * (intersection[0] + 55) + (intersection[1] - 100) * (intersection[1] - 100) + (intersection[2] - 100) * (intersection[2] - 100) - 60 * 60;
      D = B * B - 4 * A * C;
      if (D >= 0) {
        if ((-B - sqrt(D)) / (2 * A) > 0 && (-B + sqrt(D)) / (2 * A) > 0) shadow1 = true;
      }

      // check shadow2
      vTemp[0] = light2[0] - intersection[0];
      vTemp[1] = light2[1] - intersection[1];
      vTemp[2] = light2[2] - intersection[2];
      A = vTemp[0] * vTemp[0] + vTemp[1] * vTemp[1] + vTemp[2] * vTemp[2];
      B = 2 * (intersection[0] + 55) * vTemp[0] + 2 * (intersection[1] - 100) * vTemp[1] + 2 * (intersection[2] - 100) * vTemp[2];
      C = (intersection[0] + 55) * (intersection[0] + 55) + (intersection[1] - 100) * (intersection[1] - 100) + (intersection[2] - 100) * (intersection[2] - 100) - 60 * 60;
      D = B * B - 4 * A * C;
      if (D >= 0) {
        if ((-B - sqrt(D)) / (2 * A) > 0 && (-B + sqrt(D)) / (2 * A) > 0) shadow2 = true;
      }
    }

    // check sphere2
    if (choice == 0) {
      A = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
      B = 2 * (x - 15) * v[0] + 2 * (y - 60) * v[1] + 2 * (z - 5) * v[2];
      C = (x - 15) * (x - 15) + (y - 60) * (y - 60) + (z - 5) * (z - 5) - 45 * 45;
      D = B * B - 4 * A * C;
      if (D >= 0) choice = 2;
      coefficient = (-B - sqrt(D)) / (2 * A);
      if (coefficient <= 0) coefficient = (B - sqrt(D)) / (2 * A);
      intersection[0] = x + v[0] * coefficient;
      intersection[1] = y + v[1] * coefficient;
      intersection[2] = z + v[2] * coefficient;

      // check shadow1
      vTemp[0] = light1[0] - intersection[0];
      vTemp[1] = light1[1] - intersection[1];
      vTemp[2] = light1[2] - intersection[2];
      A = vTemp[0] * vTemp[0] + vTemp[1] * vTemp[1] + vTemp[2] * vTemp[2];
      B = 2 * (intersection[0] - 15) * vTemp[0] + 2 * (intersection[1] - 60) * vTemp[1] + 2 * (intersection[2] - 5) * vTemp[2];
      C = (intersection[0] - 15) * (intersection[0] - 15) + (intersection[1] - 60) * (intersection[1] - 60) + (intersection[2] - 5) * (intersection[2] - 5) - 45 * 45;
      D = B * B - 4 * A * C;
      if (D >=0) {
        if ((-B - sqrt(D)) / (2 * A) > 0 && (-B + sqrt(D)) / (2 * A) > 0) shadow1 = true;
      }

      A = vTemp[0] * vTemp[0] + vTemp[1] * vTemp[1] + vTemp[2] * vTemp[2];
      B = 2 * (intersection[0] + 55) * vTemp[0] + 2 * (intersection[1] - 100) * vTemp[1] + 2 * (intersection[2] - 100) * vTemp[2];
      C = (intersection[0] + 55) * (intersection[0] + 55) + (intersection[1] - 100) * (intersection[1] - 100) + (intersection[2] - 100) * (intersection[2] - 100) - 60 * 60;
      D = B * B - 4 * A * C;
      if (D >= 0) {
        if ((-B - sqrt(D)) / (2 * A) > 0 && (-B + sqrt(D)) / (2 * A) > 0) shadow1 = true;
      }

      // check shadow2
      vTemp[0] = light2[0] - intersection[0];
      vTemp[1] = light2[1] - intersection[1];
      vTemp[2] = light2[2] - intersection[2];
      A = vTemp[0] * vTemp[0] + vTemp[1] * vTemp[1] + vTemp[2] * vTemp[2];
      B = 2 * (intersection[0] - 15) * vTemp[0] + 2 * (intersection[1] - 60) * vTemp[1] + 2 * (intersection[2] - 5) * vTemp[2];
      C = (intersection[0] - 15) * (intersection[0] - 15) + (intersection[1] - 60) * (intersection[1] - 60) + (intersection[2] - 5) * (intersection[2] - 5) - 45 * 45;
      D = B * B - 4 * A * C;
      if (D >=0) {
        if ((-B - sqrt(D)) / (2 * A) > 0 && (-B + sqrt(D)) / (2 * A) > 0) shadow2 = true;
      }

      A = vTemp[0] * vTemp[0] + vTemp[1] * vTemp[1] + vTemp[2] * vTemp[2];
      B = 2 * (intersection[0] + 55) * vTemp[0] + 2 * (intersection[1] - 100) * vTemp[1] + 2 * (intersection[2] - 100) * vTemp[2];
      C = (intersection[0] + 55) * (intersection[0] + 55) + (intersection[1] - 100) * (intersection[1] - 100) + (intersection[2] - 100) * (intersection[2] - 100) - 60 * 60;
      D = B * B - 4 * A * C;
      if (D >= 0) {
        if ((-B - sqrt(D)) / (2 * A) > 0 && (-B + sqrt(D)) / (2 * A) > 0) shadow2 = true;
      }
    }


    // check floor
    if (choice == 0) {
      D = -y / v[1];
      if (D >= 0) {
        xx = x + v[0] * D;
        yy = y + v[1] * D;
        zz = z + v[2] * D;
        if (xx <= 200 && xx >= -200 && zz <= 250 && zz >= -250) {
          choice = 3;
          intersection[0] = xx;
          intersection[1] = yy;
          intersection[2] = zz;
        }
      }

      // check shadow1
      vTemp[0] = light1[0] - intersection[0];
      vTemp[1] = light1[1] - intersection[1];
      vTemp[2] = light1[2] - intersection[2];
      A = vTemp[0] * vTemp[0] + vTemp[1] * vTemp[1] + vTemp[2] * vTemp[2];
      B = 2 * (intersection[0] - 15) * vTemp[0] + 2 * (intersection[1] - 60) * vTemp[1] + 2 * (intersection[2] - 5) * vTemp[2];
      C = (intersection[0] - 15) * (intersection[0] - 15) + (intersection[1] - 60) * (intersection[1] - 60) + (intersection[2] - 5) * (intersection[2] - 5) - 45 * 45;
      D = B * B - 4 * A * C;
      if (D >=0) {
        if ((-B - sqrt(D)) / (2 * A) > 0 && (-B + sqrt(D)) / (2 * A) > 0) shadow1 = true;
      }

      A = vTemp[0] * vTemp[0] + vTemp[1] * vTemp[1] + vTemp[2] * vTemp[2];
      B = 2 * (intersection[0] + 55) * vTemp[0] + 2 * (intersection[1] - 100) * vTemp[1] + 2 * (intersection[2] - 100) * vTemp[2];
      C = (intersection[0] + 55) * (intersection[0] + 55) + (intersection[1] - 100) * (intersection[1] - 100) + (intersection[2] - 100) * (intersection[2] - 100) - 60 * 60;
      D = B * B - 4 * A * C;
      if (D >= 0) {
        if ((-B - sqrt(D)) / (2 * A) > 0 && (-B + sqrt(D)) / (2 * A) > 0) shadow1 = true;
      }

      // check shadow2
      vTemp[0] = light2[0] - intersection[0];
      vTemp[1] = light2[1] - intersection[1];
      vTemp[2] = light2[2] - intersection[2];
      A = vTemp[0] * vTemp[0] + vTemp[1] * vTemp[1] + vTemp[2] * vTemp[2];
      B = 2 * (intersection[0] - 15) * vTemp[0] + 2 * (intersection[1] - 60) * vTemp[1] + 2 * (intersection[2] - 5) * vTemp[2];
      C = (intersection[0] - 15) * (intersection[0] - 15) + (intersection[1] - 60) * (intersection[1] - 60) + (intersection[2] - 5) * (intersection[2] - 5) - 45 * 45;
      D = B * B - 4 * A * C;
      if (D >=0) {
        if ((-B - sqrt(D)) / (2 * A) > 0 && (-B + sqrt(D)) / (2 * A) > 0) shadow2 = true;
      }

      A = vTemp[0] * vTemp[0] + vTemp[1] * vTemp[1] + vTemp[2] * vTemp[2];
      B = 2 * (intersection[0] + 55) * vTemp[0] + 2 * (intersection[1] - 100) * vTemp[1] + 2 * (intersection[2] - 100) * vTemp[2];
      C = (intersection[0] + 55) * (intersection[0] + 55) + (intersection[1] - 100) * (intersection[1] - 100) + (intersection[2] - 100) * (intersection[2] - 100) - 60 * 60;
      D = B * B - 4 * A * C;
      if (D >= 0) {
        if ((-B - sqrt(D)) / (2 * A) > 0 && (-B + sqrt(D)) / (2 * A) > 0) shadow2 = true;
      }
    }

    // draw scene
    switch(choice) {
    case 0:
      set(i, j, color(69, 150, 243));
      break;
    case 1:
      // direction of incoming light
      S = new float[]{light1[0] - intersection[0], light1[1] - intersection[1], light1[2] - intersection[2]};
      // normalize S
      temp = sqrt(S[0] * S[0] + S[1] * S[1] + S[2] * S[2]);
      S[0] = S[0] / temp;
      S[1] = S[1] / temp;
      S[2] = S[2] / temp;

      // normal vector
      N = new float[]{intersection[0] + 55.0, intersection[1] - 100.0, intersection[2] - 100.0};
      // normalize N
      temp = sqrt(N[0] * N[0] + N[1] * N[1] + N[2] * N[2]);
      N[0] = N[0] / temp;
      N[1] = N[1] / temp;
      N[2] = N[2] / temp;

      // destination (to the camera)
      V = new float[]{eye[0] - intersection[0], eye[1] - intersection[1], eye[2] - intersection[2]};
      // normalize V
      temp = sqrt(V[0] * V[0] + V[1] * V[1] + V[2] * V[2]);
      V[0] = V[0] / temp;
      V[1] = V[1] / temp;
      V[2] = V[2] / temp;

      // reflection vector
      dotTemp = S[0] * N[0] + S[1] * N[1] + S[2] * N[2];
      R = new float[]{-S[0] + 2 * N[0] * dotTemp, -S[1] + 2 * N[1] * dotTemp, -S[2] + 2 * N[2] * dotTemp};
      // normalize R
      temp = sqrt(R[0] * R[0] + R[1] * R[1] + R[2] * R[2]);
      R[0] = R[0] / temp;
      R[1] = R[1] / temp;
      R[2] = R[2] / temp;

      ka = 0.5;
      kd = 0.6;
      ks = 0.3;
      specExp = 30.0;

      ambient = new float[]{150.0 * ka, 150.0 * ka, 150.0 * ka};
      dotTemp = N[0] * S[0] + N[1] * S[1] + N[2] * S[2];
      if (dotTemp >= 0) diffuse1 = new float[]{150.0 * dotTemp * kd, 150.0 * dotTemp * kd, 150.0 * dotTemp * kd};
      else diffuse1 = new float[]{0.0, 0.0, 0.0};
      dotTemp = R[0] * V[0] + R[1] * V[1] + R[2] * V[2];
      if (dotTemp >= 0) specDot = pow(dotTemp, specExp);
      else specDot = 0.0;
      specular1 = new float[]{255.0 * specDot * ks, 255.0 * specDot * ks, 255.0 * specDot * ks};

      // light source 2
      S = new float[]{light2[0] - intersection[0], light2[1] - intersection[1], light2[2] - intersection[2]};
      // normalize S
      temp = sqrt(S[0] * S[0] + S[1] * S[1] + S[2] * S[2]);
      S[0] = S[0] / temp;
      S[1] = S[1] / temp;
      S[2] = S[2] / temp;

      dotTemp = S[0] * N[0] + S[1] * N[1] + S[2] * N[2];
      R = new float[]{-S[0] + 2 * N[0] * dotTemp, -S[1] + 2 * N[1] * dotTemp, -S[2] + 2 * N[2] * dotTemp};
      // normalize R
      temp = sqrt(R[0] * R[0] + R[1] * R[1] + R[2] * R[2]);
      R[0] = R[0] / temp;
      R[1] = R[1] / temp;
      R[2] = R[2] / temp;

      dotTemp = N[0] * S[0] + N[1] * S[1] + N[2] * S[2];
      if (dotTemp >= 0) diffuse2 = new float[]{150.0 * dotTemp * kd, 150.0 * dotTemp * kd, 150.0 * dotTemp * kd};
      else diffuse2 = new float[]{0.0, 0.0, 0.0};
      dotTemp = R[0] * V[0] + R[1] * V[1] + R[2] * V[2];
      if (dotTemp >= 0) specDot = pow(dotTemp, specExp);
      else specDot = 0.0;
      specular2 = new float[]{255.0 * specDot * ks, 255.0 * specDot * ks, 255.0 * specDot * ks};

      if (!shadow1 && !shadow2) c = new float[]{ambient[0] + diffuse1[0] + specular1[0] + diffuse2[0] + specular2[0], ambient[1] + diffuse1[1] + specular1[1] + diffuse2[1] + specular2[1], ambient[2] + diffuse1[2] + specular1[2] + diffuse2[2] + specular2[2]};
      else if (!shadow1 && shadow2) c = new float[]{ambient[0] + diffuse1[0] + specular1[0], ambient[1] + diffuse1[1] + specular1[1], ambient[2] + diffuse1[2] + specular1[2]};
      else if (shadow1 && !shadow2) c = new float[]{ambient[0] + diffuse2[0] + specular2[0], ambient[1] + diffuse2[1] + specular2[1], ambient[2] + diffuse2[2] + specular2[2]};
      else c = new float[]{ambient[0], ambient[1], ambient[2]};
      set(i, j, color(c[0], c[1], c[2]));
      break;
    case 2:
      // direction of incoming light
      S = new float[]{light1[0] - intersection[0], light1[1] - intersection[1], light1[2] - intersection[2]};
      // normalize S
      temp = sqrt(S[0] * S[0] + S[1] * S[1] + S[2] * S[2]);
      S[0] = S[0] / temp;
      S[1] = S[1] / temp;
      S[2] = S[2] / temp;

      // normal vector
      N = new float[]{intersection[0] - 15.0, intersection[1] - 60.0, intersection[2] - 5.0};
      // normalize N
      temp = sqrt(N[0] * N[0] + N[1] * N[1] + N[2] * N[2]);
      N[0] = N[0] / temp;
      N[1] = N[1] / temp;
      N[2] = N[2] / temp;

      // destination (to the camera)
      V = new float[]{eye[0] - intersection[0], eye[1] - intersection[1], eye[2] - intersection[2]};
      // normalize V
      temp = sqrt(V[0] * V[0] + V[1] * V[1] + V[2] * V[2]);
      V[0] = V[0] / temp;
      V[1] = V[1] / temp;
      V[2] = V[2] / temp;

      // reflection vector
      dotTemp = S[0] * N[0] + S[1] * N[1] + S[2] * N[2];
      R = new float[]{-S[0] + 2 * N[0] * dotTemp, -S[1] + 2 * N[1] * dotTemp, -S[2] + 2 * N[2] * dotTemp};
      // normalize R
      temp = sqrt(R[0] * R[0] + R[1] * R[1] + R[2] * R[2]);
      R[0] = R[0] / temp;
      R[1] = R[1] / temp;
      R[2] = R[2] / temp;

      ka = 0.5;
      kd = 0.5;
      ks = 0.3;
      specExp = 30.0;

      ambient = new float[]{250.0 * ka, 250.0 * ka, 250.0 * ka};
      dotTemp = N[0] * S[0] + N[1] * S[1] + N[2] * S[2];
      if (dotTemp >= 0) diffuse1 = new float[]{250.0 * dotTemp * kd, 250.0 * dotTemp * kd, 250.0 * dotTemp * kd};
      else diffuse1 = new float[]{0.0, 0.0, 0.0};
      dotTemp = R[0] * V[0] + R[1] * V[1] + R[2] * V[2];
      if (dotTemp >= 0) specDot = pow(dotTemp, specExp);
      else specDot = 0.0;
      specular1 = new float[]{255.0 * specDot * ks, 255.0 * specDot * ks, 255.0 * specDot * ks};

      // light source 2
      S = new float[]{light2[0] - intersection[0], light2[1] - intersection[1], light2[2] - intersection[2]};
      // normalize S
      temp = sqrt(S[0] * S[0] + S[1] * S[1] + S[2] * S[2]);
      S[0] = S[0] / temp;
      S[1] = S[1] / temp;
      S[2] = S[2] / temp;

      dotTemp = S[0] * N[0] + S[1] * N[1] + S[2] * N[2];
      R = new float[]{-S[0] + 2 * N[0] * dotTemp, -S[1] + 2 * N[1] * dotTemp, -S[2] + 2 * N[2] * dotTemp};
      // normalize R
      temp = sqrt(R[0] * R[0] + R[1] * R[1] + R[2] * R[2]);
      R[0] = R[0] / temp;
      R[1] = R[1] / temp;
      R[2] = R[2] / temp;

      dotTemp = N[0] * S[0] + N[1] * S[1] + N[2] * S[2];
      if (dotTemp >= 0) diffuse2 = new float[]{250.0 * dotTemp * kd, 250.0 * dotTemp * kd, 250.0 * dotTemp * kd};
      else diffuse2 = new float[]{0.0, 0.0, 0.0};
      dotTemp = R[0] * V[0] + R[1] * V[1] + R[2] * V[2];
      if (dotTemp >= 0) specDot = pow(dotTemp, specExp);
      else specDot = 0.0;
      specular2 = new float[]{255.0 * specDot * ks, 255.0 * specDot * ks, 255.0 * specDot * ks};
      
      if (!shadow1 && !shadow2) c = new float[]{ambient[0] + diffuse1[0] + specular1[0] + diffuse2[0] + specular2[0], ambient[1] + diffuse1[1] + specular1[1] + diffuse2[1] + specular2[1], ambient[2] + diffuse1[2] + specular1[2] + diffuse2[2] + specular2[2]};
      else if (!shadow1 && shadow2) c = new float[]{ambient[0] + diffuse1[0] + specular1[0], ambient[1] + diffuse1[1] + specular1[1], ambient[2] + diffuse1[2] + specular1[2]};
      else if (shadow1 && !shadow2) c = new float[]{ambient[0] + diffuse2[0] + specular2[0], ambient[1] + diffuse2[1] + specular2[1], ambient[2] + diffuse2[2] + specular2[2]};
      else c = new float[]{ambient[0], ambient[1], ambient[2]};


      set(i, j, color(c[0], c[1], c[2]));
      break;
    case 3:
      // direction of incoming light
      S = new float[]{light1[0] - intersection[0], light1[1] - intersection[1], light1[2] - intersection[2]};
      // normalize S
      temp = sqrt(S[0] * S[0] + S[1] * S[1] + S[2] * S[2]);
      S[0] = S[0] / temp;
      S[1] = S[1] / temp;
      S[2] = S[2] / temp;

      // normal vector
      N = new float[]{0.0, 1.0, 0.0};
      // normalize N
      temp = sqrt(N[0] * N[0] + N[1] * N[1] + N[2] * N[2]);
      N[0] = N[0] / temp;
      N[1] = N[1] / temp;
      N[2] = N[2] / temp;

      // destination (to the camera)
      V = new float[]{eye[0] - intersection[0], eye[1] - intersection[1], eye[2] - intersection[2]};
      // normalize V
      temp = sqrt(V[0] * V[0] + V[1] * V[1] + V[2] * V[2]);
      V[0] = V[0] / temp;
      V[1] = V[1] / temp;
      V[2] = V[2] / temp;

      // reflection vector
      dotTemp = S[0] * N[0] + S[1] * N[1] + S[2] * N[2];
      R = new float[]{-S[0] + 2 * N[0] * dotTemp, -S[1] + 2 * N[1] * dotTemp, -S[2] + 2 * N[2] * dotTemp};
      // normalize R
      temp = sqrt(R[0] * R[0] + R[1] * R[1] + R[2] * R[2]);
      R[0] = R[0] / temp;
      R[1] = R[1] / temp;
      R[2] = R[2] / temp;

      ka = 0.5;
      kd = 0.5;
      ks = 0.3;
      specExp = 30.0;

      ambient = new float[]{255.0 * ka, 228.0 * ka, 108.0 * ka};
      dotTemp = N[0] * S[0] + N[1] * S[1] + N[2] * S[2];
      if (dotTemp >= 0) diffuse1 = new float[]{255.0 * dotTemp * kd, 228.0 * dotTemp * kd, 108.0 * dotTemp * kd};
      else diffuse1 = new float[]{0.0, 0.0, 0.0};
      dotTemp = R[0] * V[0] + R[1] * V[1] + R[2] * V[2];
      if (dotTemp >= 0) specDot = pow(dotTemp, specExp);
      else specDot = 0.0;
      specular1 = new float[]{255.0 * specDot * ks, 255.0 * specDot * ks, 255.0 * specDot * ks};

      // light source 2
      S = new float[]{light2[0] - intersection[0], light2[1] - intersection[1], light2[2] - intersection[2]};
      // normalize S
      temp = sqrt(S[0] * S[0] + S[1] * S[1] + S[2] * S[2]);
      S[0] = S[0] / temp;
      S[1] = S[1] / temp;
      S[2] = S[2] / temp;

      dotTemp = S[0] * N[0] + S[1] * N[1] + S[2] * N[2];
      R = new float[]{-S[0] + 2 * N[0] * dotTemp, -S[1] + 2 * N[1] * dotTemp, -S[2] + 2 * N[2] * dotTemp};
      // normalize R
      temp = sqrt(R[0] * R[0] + R[1] * R[1] + R[2] * R[2]);
      R[0] = R[0] / temp;
      R[1] = R[1] / temp;
      R[2] = R[2] / temp;

      dotTemp = N[0] * S[0] + N[1] * S[1] + N[2] * S[2];
      if (dotTemp >= 0) diffuse2 = new float[]{255.0 * dotTemp * kd, 228.0 * dotTemp * kd, 108.0 * dotTemp * kd};
      else diffuse2 = new float[]{0.0, 0.0, 0.0};
      dotTemp = R[0] * V[0] + R[1] * V[1] + R[2] * V[2];
      if (dotTemp >= 0) specDot = pow(dotTemp, specExp);
      else specDot = 0.0;
      specular2 = new float[]{255.0 * specDot * ks, 255.0 * specDot * ks, 255.0 * specDot * ks};

      if (!shadow1 && !shadow2) c = new float[]{ambient[0] + diffuse1[0] + specular1[0] + diffuse2[0] + specular2[0], ambient[1] + diffuse1[1] + specular1[1] + diffuse2[1] + specular2[1], ambient[2] + diffuse1[2] + specular1[2] + diffuse2[2] + specular2[2]};
      else if (!shadow1 && shadow2) c = new float[]{ambient[0] + diffuse1[0] + specular1[0], ambient[1] + diffuse1[1] + specular1[1], ambient[2] + diffuse1[2] + specular1[2]};
      else if (shadow1 && !shadow2) c = new float[]{ambient[0] + diffuse2[0] + specular2[0], ambient[1] + diffuse2[1] + specular2[1], ambient[2] + diffuse2[2] + specular2[2]};
      else c = new float[]{ambient[0], ambient[1], ambient[2]};



      set(i, j, color(c[0], c[1], c[2]));
      break;
    }
  }
}

//save("Assignment_3_Multiple_Lights.png");
