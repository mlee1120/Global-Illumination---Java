size(560, 315);

// light position
float[] light = new float[]{-100.0, 250.0, 300.0};

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
int choice1, choice2, choice3, choice4; // intersection with which object
float[] intersection1 = new float[3]; // intersections
float[] intersection2 = new float[3];
float[] intersection3 = new float[3];
float[] intersection4 = new float[3];
float coefficient; // auxiliary parameter for finding intersections
float temp; // for normalization
float ka, kd, ks, specExp; // reflection coefficients
float[] ambient, diffuse, specular, c;
float[] S, N, V, R;
float dotTemp, specDot;
boolean shadow1, shadow2, shadow3, shadow4;
float[] color1 = new float[3];
float[] color2 = new float[3];
float[] color3 = new float[3];
float[] color4 = new float[3];


/* ray tracing
 sphere1: (x+55)^2 + (y-100)^2 + (z-100)^2 = 60^2
 sphere2: (x-15)^2 + (y-60)^2 + (z-5)^2 = 45^2
 floor: y = 0 with -200 <= x <= 200 && -250 <= z <= 250
 */
for (int i = 0; i < width; i++) {
  for (int j = 0; j < height; j++) {
    choice1 = 0;
    choice2 = 0;
    choice3 = 0;
    choice4 = 0;
    shadow1 = false;
    shadow2 = false;
    shadow3 = false;
    shadow4 = false;

    // upper left
    x = upper_left[0] + i * vector_right[0] - j * vector_up[0] - (vector_right[0] / 4.0) + (vector_up[0] / 4.0);
    y = upper_left[1] + i * vector_right[1] - j * vector_up[1] - (vector_right[1] / 4.0) + (vector_up[1] / 4.0);
    z = upper_left[2] + i * vector_right[2] - j * vector_up[2] - (vector_right[2] / 4.0) + (vector_up[2] / 4.0);
    v[0] = x - eye[0];
    v[1] = y - eye[1];
    v[2] = z - eye[2];

    // check sphere1
    A = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
    B = 2 * (x + 55) * v[0] + 2 * (y - 100) * v[1] + 2 * (z - 100) * v[2];
    C = (x + 55) * (x + 55) + (y - 100) * (y - 100) + (z - 100) * (z - 100) - 60 * 60;
    D = B * B - 4 * A * C;
    if (D >= 0) {
      choice1 = 1;
      coefficient = (-B - sqrt(D)) / (2 * A);
      if (coefficient <= 0) coefficient = (B - sqrt(D)) / (2 * A);
      intersection1[0] = x + v[0] * coefficient;
      intersection1[1] = y + v[1] * coefficient;
      intersection1[2] = z + v[2] * coefficient;

      // check shadow
      vTemp[0] = light[0] - intersection1[0];
      vTemp[1] = light[1] - intersection1[1];
      vTemp[2] = light[2] - intersection1[2];
      A = vTemp[0] * vTemp[0] + vTemp[1] * vTemp[1] + vTemp[2] * vTemp[2];
      B = 2 * (intersection1[0] + 55) * vTemp[0] + 2 * (intersection1[1] - 100) * vTemp[1] + 2 * (intersection1[2] - 100) * vTemp[2];
      C = (intersection1[0] + 55) * (intersection1[0] + 55) + (intersection1[1] - 100) * (intersection1[1] - 100) + (intersection1[2] - 100) * (intersection1[2] - 100) - 60 * 60;
      D = B * B - 4 * A * C;
      if (D >= 0) {
        if ((-B - sqrt(D)) / (2 * A) > 0 && (-B + sqrt(D)) / (2 * A) > 0) shadow1 = true;
      }
    }

    // check sphere2
    if (choice1 == 0) {
      A = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
      B = 2 * (x - 15) * v[0] + 2 * (y - 60) * v[1] + 2 * (z - 5) * v[2];
      C = (x - 15) * (x - 15) + (y - 60) * (y - 60) + (z - 5) * (z - 5) - 45 * 45;
      D = B * B - 4 * A * C;
      if (D >= 0) choice1 = 2;
      coefficient = (-B - sqrt(D)) / (2 * A);
      if (coefficient <= 0) coefficient = (B - sqrt(D)) / (2 * A);
      intersection1[0] = x + v[0] * coefficient;
      intersection1[1] = y + v[1] * coefficient;
      intersection1[2] = z + v[2] * coefficient;

      // check shadow
      vTemp[0] = light[0] - intersection1[0];
      vTemp[1] = light[1] - intersection1[1];
      vTemp[2] = light[2] - intersection1[2];
      A = vTemp[0] * vTemp[0] + vTemp[1] * vTemp[1] + vTemp[2] * vTemp[2];
      B = 2 * (intersection1[0] - 15) * vTemp[0] + 2 * (intersection1[1] - 60) * vTemp[1] + 2 * (intersection1[2] - 5) * vTemp[2];
      C = (intersection1[0] - 15) * (intersection1[0] - 15) + (intersection1[1] - 60) * (intersection1[1] - 60) + (intersection1[2] - 5) * (intersection1[2] - 5) - 45 * 45;
      D = B * B - 4 * A * C;
      if (D >=0) {
        if ((-B - sqrt(D)) / (2 * A) > 0 && (-B + sqrt(D)) / (2 * A) > 0) shadow1 = true;
      }

      A = vTemp[0] * vTemp[0] + vTemp[1] * vTemp[1] + vTemp[2] * vTemp[2];
      B = 2 * (intersection1[0] + 55) * vTemp[0] + 2 * (intersection1[1] - 100) * vTemp[1] + 2 * (intersection1[2] - 100) * vTemp[2];
      C = (intersection1[0] + 55) * (intersection1[0] + 55) + (intersection1[1] - 100) * (intersection1[1] - 100) + (intersection1[2] - 100) * (intersection1[2] - 100) - 60 * 60;
      D = B * B - 4 * A * C;
      if (D >= 0) {
        if ((-B - sqrt(D)) / (2 * A) > 0 && (-B + sqrt(D)) / (2 * A) > 0) shadow1 = true;
      }
    }


    // check floor
    if (choice1 == 0) {
      D = -y / v[1];
      if (D >= 0) {
        xx = x + v[0] * D;
        yy = y + v[1] * D;
        zz = z + v[2] * D;
        if (xx <= 200 && xx >= -200 && zz <= 250 && zz >= -250) {
          choice1 = 3;
          intersection1[0] = xx;
          intersection1[1] = yy;
          intersection1[2] = zz;
        }
      }

      // check shadow
      vTemp[0] = light[0] - intersection1[0];
      vTemp[1] = light[1] - intersection1[1];
      vTemp[2] = light[2] - intersection1[2];
      A = vTemp[0] * vTemp[0] + vTemp[1] * vTemp[1] + vTemp[2] * vTemp[2];
      B = 2 * (intersection1[0] - 15) * vTemp[0] + 2 * (intersection1[1] - 60) * vTemp[1] + 2 * (intersection1[2] - 5) * vTemp[2];
      C = (intersection1[0] - 15) * (intersection1[0] - 15) + (intersection1[1] - 60) * (intersection1[1] - 60) + (intersection1[2] - 5) * (intersection1[2] - 5) - 45 * 45;
      D = B * B - 4 * A * C;
      if (D >=0) {
        if ((-B - sqrt(D)) / (2 * A) > 0 && (-B + sqrt(D)) / (2 * A) > 0) shadow1 = true;
      }

      A = vTemp[0] * vTemp[0] + vTemp[1] * vTemp[1] + vTemp[2] * vTemp[2];
      B = 2 * (intersection1[0] + 55) * vTemp[0] + 2 * (intersection1[1] - 100) * vTemp[1] + 2 * (intersection1[2] - 100) * vTemp[2];
      C = (intersection1[0] + 55) * (intersection1[0] + 55) + (intersection1[1] - 100) * (intersection1[1] - 100) + (intersection1[2] - 100) * (intersection1[2] - 100) - 60 * 60;
      D = B * B - 4 * A * C;
      if (D >= 0) {
        if ((-B - sqrt(D)) / (2 * A) > 0 && (-B + sqrt(D)) / (2 * A) > 0) shadow1 = true;
      }
    }

    // upper right
    x = upper_left[0] + i * vector_right[0] - j * vector_up[0] + (vector_right[0] / 4.0) + (vector_up[0] / 4.0);
    y = upper_left[1] + i * vector_right[1] - j * vector_up[1] + (vector_right[1] / 4.0) + (vector_up[1] / 4.0);
    z = upper_left[2] + i * vector_right[2] - j * vector_up[2] + (vector_right[2] / 4.0) + (vector_up[2] / 4.0);
    v[0] = x - eye[0];
    v[1] = y - eye[1];
    v[2] = z - eye[2];

    // check sphere1
    A = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
    B = 2 * (x + 55) * v[0] + 2 * (y - 100) * v[1] + 2 * (z - 100) * v[2];
    C = (x + 55) * (x + 55) + (y - 100) * (y - 100) + (z - 100) * (z - 100) - 60 * 60;
    D = B * B - 4 * A * C;
    if (D >= 0) {
      choice2 = 1;
      coefficient = (-B - sqrt(D)) / (2 * A);
      if (coefficient <= 0) coefficient = (B - sqrt(D)) / (2 * A);
      intersection2[0] = x + v[0] * coefficient;
      intersection2[1] = y + v[1] * coefficient;
      intersection2[2] = z + v[2] * coefficient;

      // check shadow
      vTemp[0] = light[0] - intersection2[0];
      vTemp[1] = light[1] - intersection2[1];
      vTemp[2] = light[2] - intersection2[2];
      A = vTemp[0] * vTemp[0] + vTemp[1] * vTemp[1] + vTemp[2] * vTemp[2];
      B = 2 * (intersection2[0] + 55) * vTemp[0] + 2 * (intersection2[1] - 100) * vTemp[1] + 2 * (intersection2[2] - 100) * vTemp[2];
      C = (intersection2[0] + 55) * (intersection2[0] + 55) + (intersection2[1] - 100) * (intersection2[1] - 100) + (intersection2[2] - 100) * (intersection2[2] - 100) - 60 * 60;
      D = B * B - 4 * A * C;
      if (D >= 0) {
        if ((-B - sqrt(D)) / (2 * A) > 0 && (-B + sqrt(D)) / (2 * A) > 0) shadow2 = true;
      }
    }

    // check sphere2
    if (choice2 == 0) {
      A = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
      B = 2 * (x - 15) * v[0] + 2 * (y - 60) * v[1] + 2 * (z - 5) * v[2];
      C = (x - 15) * (x - 15) + (y - 60) * (y - 60) + (z - 5) * (z - 5) - 45 * 45;
      D = B * B - 4 * A * C;
      if (D >= 0) choice2 = 2;
      coefficient = (-B - sqrt(D)) / (2 * A);
      if (coefficient <= 0) coefficient = (B - sqrt(D)) / (2 * A);
      intersection2[0] = x + v[0] * coefficient;
      intersection2[1] = y + v[1] * coefficient;
      intersection2[2] = z + v[2] * coefficient;

      // check shadow
      vTemp[0] = light[0] - intersection2[0];
      vTemp[1] = light[1] - intersection2[1];
      vTemp[2] = light[2] - intersection2[2];
      A = vTemp[0] * vTemp[0] + vTemp[1] * vTemp[1] + vTemp[2] * vTemp[2];
      B = 2 * (intersection2[0] - 15) * vTemp[0] + 2 * (intersection2[1] - 60) * vTemp[1] + 2 * (intersection2[2] - 5) * vTemp[2];
      C = (intersection2[0] - 15) * (intersection2[0] - 15) + (intersection2[1] - 60) * (intersection2[1] - 60) + (intersection2[2] - 5) * (intersection2[2] - 5) - 45 * 45;
      D = B * B - 4 * A * C;
      if (D >=0) {
        if ((-B - sqrt(D)) / (2 * A) > 0 && (-B + sqrt(D)) / (2 * A) > 0) shadow2 = true;
      }

      A = vTemp[0] * vTemp[0] + vTemp[1] * vTemp[1] + vTemp[2] * vTemp[2];
      B = 2 * (intersection2[0] + 55) * vTemp[0] + 2 * (intersection2[1] - 100) * vTemp[1] + 2 * (intersection2[2] - 100) * vTemp[2];
      C = (intersection2[0] + 55) * (intersection2[0] + 55) + (intersection2[1] - 100) * (intersection2[1] - 100) + (intersection2[2] - 100) * (intersection2[2] - 100) - 60 * 60;
      D = B * B - 4 * A * C;
      if (D >= 0) {
        if ((-B - sqrt(D)) / (2 * A) > 0 && (-B + sqrt(D)) / (2 * A) > 0) shadow2 = true;
      }
    }


    // check floor
    if (choice2 == 0) {
      D = -y / v[1];
      if (D >= 0) {
        xx = x + v[0] * D;
        yy = y + v[1] * D;
        zz = z + v[2] * D;
        if (xx <= 200 && xx >= -200 && zz <= 250 && zz >= -250) {
          choice2 = 3;
          intersection2[0] = xx;
          intersection2[1] = yy;
          intersection2[2] = zz;
        }
      }

      // check shadow
      vTemp[0] = light[0] - intersection2[0];
      vTemp[1] = light[1] - intersection2[1];
      vTemp[2] = light[2] - intersection2[2];
      A = vTemp[0] * vTemp[0] + vTemp[1] * vTemp[1] + vTemp[2] * vTemp[2];
      B = 2 * (intersection2[0] - 15) * vTemp[0] + 2 * (intersection2[1] - 60) * vTemp[1] + 2 * (intersection2[2] - 5) * vTemp[2];
      C = (intersection2[0] - 15) * (intersection2[0] - 15) + (intersection2[1] - 60) * (intersection2[1] - 60) + (intersection2[2] - 5) * (intersection2[2] - 5) - 45 * 45;
      D = B * B - 4 * A * C;
      if (D >=0) {
        if ((-B - sqrt(D)) / (2 * A) > 0 && (-B + sqrt(D)) / (2 * A) > 0) shadow2 = true;
      }

      A = vTemp[0] * vTemp[0] + vTemp[1] * vTemp[1] + vTemp[2] * vTemp[2];
      B = 2 * (intersection2[0] + 55) * vTemp[0] + 2 * (intersection2[1] - 100) * vTemp[1] + 2 * (intersection2[2] - 100) * vTemp[2];
      C = (intersection2[0] + 55) * (intersection2[0] + 55) + (intersection2[1] - 100) * (intersection2[1] - 100) + (intersection2[2] - 100) * (intersection2[2] - 100) - 60 * 60;
      D = B * B - 4 * A * C;
      if (D >= 0) {
        if ((-B - sqrt(D)) / (2 * A) > 0 && (-B + sqrt(D)) / (2 * A) > 0) shadow2 = true;
      }
    }

    // lower left
    x = upper_left[0] + i * vector_right[0] - j * vector_up[0] - (vector_right[0] / 4.0) - (vector_up[0] / 4.0);
    y = upper_left[1] + i * vector_right[1] - j * vector_up[1] - (vector_right[1] / 4.0) - (vector_up[1] / 4.0);
    z = upper_left[2] + i * vector_right[2] - j * vector_up[2] - (vector_right[2] / 4.0) - (vector_up[2] / 4.0);
    v[0] = x - eye[0];
    v[1] = y - eye[1];
    v[2] = z - eye[2];

    // check sphere1
    A = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
    B = 2 * (x + 55) * v[0] + 2 * (y - 100) * v[1] + 2 * (z - 100) * v[2];
    C = (x + 55) * (x + 55) + (y - 100) * (y - 100) + (z - 100) * (z - 100) - 60 * 60;
    D = B * B - 4 * A * C;
    if (D >= 0) {
      choice3 = 1;
      coefficient = (-B - sqrt(D)) / (2 * A);
      if (coefficient <= 0) coefficient = (B - sqrt(D)) / (2 * A);
      intersection3[0] = x + v[0] * coefficient;
      intersection3[1] = y + v[1] * coefficient;
      intersection3[2] = z + v[2] * coefficient;

      // check shadow
      vTemp[0] = light[0] - intersection3[0];
      vTemp[1] = light[1] - intersection3[1];
      vTemp[2] = light[2] - intersection3[2];
      A = vTemp[0] * vTemp[0] + vTemp[1] * vTemp[1] + vTemp[2] * vTemp[2];
      B = 2 * (intersection3[0] + 55) * vTemp[0] + 2 * (intersection3[1] - 100) * vTemp[1] + 2 * (intersection3[2] - 100) * vTemp[2];
      C = (intersection3[0] + 55) * (intersection3[0] + 55) + (intersection3[1] - 100) * (intersection3[1] - 100) + (intersection3[2] - 100) * (intersection3[2] - 100) - 60 * 60;
      D = B * B - 4 * A * C;
      if (D >= 0) {
        if ((-B - sqrt(D)) / (2 * A) > 0 && (-B + sqrt(D)) / (2 * A) > 0) shadow3 = true;
      }
    }

    // check sphere2
    if (choice3 == 0) {
      A = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
      B = 2 * (x - 15) * v[0] + 2 * (y - 60) * v[1] + 2 * (z - 5) * v[2];
      C = (x - 15) * (x - 15) + (y - 60) * (y - 60) + (z - 5) * (z - 5) - 45 * 45;
      D = B * B - 4 * A * C;
      if (D >= 0) choice3 = 2;
      coefficient = (-B - sqrt(D)) / (2 * A);
      if (coefficient <= 0) coefficient = (B - sqrt(D)) / (2 * A);
      intersection3[0] = x + v[0] * coefficient;
      intersection3[1] = y + v[1] * coefficient;
      intersection3[2] = z + v[2] * coefficient;

      // check shadow
      vTemp[0] = light[0] - intersection3[0];
      vTemp[1] = light[1] - intersection3[1];
      vTemp[2] = light[2] - intersection3[2];
      A = vTemp[0] * vTemp[0] + vTemp[1] * vTemp[1] + vTemp[2] * vTemp[2];
      B = 2 * (intersection3[0] - 15) * vTemp[0] + 2 * (intersection3[1] - 60) * vTemp[1] + 2 * (intersection3[2] - 5) * vTemp[2];
      C = (intersection3[0] - 15) * (intersection3[0] - 15) + (intersection3[1] - 60) * (intersection3[1] - 60) + (intersection3[2] - 5) * (intersection3[2] - 5) - 45 * 45;
      D = B * B - 4 * A * C;
      if (D >=0) {
        if ((-B - sqrt(D)) / (2 * A) > 0 && (-B + sqrt(D)) / (2 * A) > 0) shadow3 = true;
      }

      A = vTemp[0] * vTemp[0] + vTemp[1] * vTemp[1] + vTemp[2] * vTemp[2];
      B = 2 * (intersection3[0] + 55) * vTemp[0] + 2 * (intersection3[1] - 100) * vTemp[1] + 2 * (intersection3[2] - 100) * vTemp[2];
      C = (intersection3[0] + 55) * (intersection3[0] + 55) + (intersection3[1] - 100) * (intersection3[1] - 100) + (intersection3[2] - 100) * (intersection3[2] - 100) - 60 * 60;
      D = B * B - 4 * A * C;
      if (D >= 0) {
        if ((-B - sqrt(D)) / (2 * A) > 0 && (-B + sqrt(D)) / (2 * A) > 0) shadow3 = true;
      }
    }

    // check floor
    if (choice3 == 0) {
      D = -y / v[1];
      if (D >= 0) {
        xx = x + v[0] * D;
        yy = y + v[1] * D;
        zz = z + v[2] * D;
        if (xx <= 200 && xx >= -200 && zz <= 250 && zz >= -250) {
          choice3 = 3;
          intersection3[0] = xx;
          intersection3[1] = yy;
          intersection3[2] = zz;
        }
      }

      // check shadow
      vTemp[0] = light[0] - intersection3[0];
      vTemp[1] = light[1] - intersection3[1];
      vTemp[2] = light[2] - intersection3[2];
      A = vTemp[0] * vTemp[0] + vTemp[1] * vTemp[1] + vTemp[2] * vTemp[2];
      B = 2 * (intersection3[0] - 15) * vTemp[0] + 2 * (intersection3[1] - 60) * vTemp[1] + 2 * (intersection3[2] - 5) * vTemp[2];
      C = (intersection3[0] - 15) * (intersection3[0] - 15) + (intersection3[1] - 60) * (intersection3[1] - 60) + (intersection3[2] - 5) * (intersection3[2] - 5) - 45 * 45;
      D = B * B - 4 * A * C;
      if (D >=0) {
        if ((-B - sqrt(D)) / (2 * A) > 0 && (-B + sqrt(D)) / (2 * A) > 0) shadow3 = true;
      }

      A = vTemp[0] * vTemp[0] + vTemp[1] * vTemp[1] + vTemp[2] * vTemp[2];
      B = 2 * (intersection3[0] + 55) * vTemp[0] + 2 * (intersection3[1] - 100) * vTemp[1] + 2 * (intersection3[2] - 100) * vTemp[2];
      C = (intersection3[0] + 55) * (intersection3[0] + 55) + (intersection3[1] - 100) * (intersection3[1] - 100) + (intersection3[2] - 100) * (intersection3[2] - 100) - 60 * 60;
      D = B * B - 4 * A * C;
      if (D >= 0) {
        if ((-B - sqrt(D)) / (2 * A) > 0 && (-B + sqrt(D)) / (2 * A) > 0) shadow3 = true;
      }
    }

    // lower right
    x = upper_left[0] + i * vector_right[0] - j * vector_up[0] + (vector_right[0] / 4.0) - (vector_up[0] / 4.0);
    y = upper_left[1] + i * vector_right[1] - j * vector_up[1] + (vector_right[1] / 4.0) - (vector_up[1] / 4.0);
    z = upper_left[2] + i * vector_right[2] - j * vector_up[2] + (vector_right[2] / 4.0) - (vector_up[2] / 4.0);
    v[0] = x - eye[0];
    v[1] = y - eye[1];
    v[2] = z - eye[2];

    // check sphere1
    A = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
    B = 2 * (x + 55) * v[0] + 2 * (y - 100) * v[1] + 2 * (z - 100) * v[2];
    C = (x + 55) * (x + 55) + (y - 100) * (y - 100) + (z - 100) * (z - 100) - 60 * 60;
    D = B * B - 4 * A * C;
    if (D >= 0) {
      choice4 = 1;
      coefficient = (-B - sqrt(D)) / (2 * A);
      if (coefficient <= 0) coefficient = (B - sqrt(D)) / (2 * A);
      intersection4[0] = x + v[0] * coefficient;
      intersection4[1] = y + v[1] * coefficient;
      intersection4[2] = z + v[2] * coefficient;

      // check shadow
      vTemp[0] = light[0] - intersection4[0];
      vTemp[1] = light[1] - intersection4[1];
      vTemp[2] = light[2] - intersection4[2];
      A = vTemp[0] * vTemp[0] + vTemp[1] * vTemp[1] + vTemp[2] * vTemp[2];
      B = 2 * (intersection4[0] + 55) * vTemp[0] + 2 * (intersection4[1] - 100) * vTemp[1] + 2 * (intersection4[2] - 100) * vTemp[2];
      C = (intersection4[0] + 55) * (intersection4[0] + 55) + (intersection4[1] - 100) * (intersection4[1] - 100) + (intersection4[2] - 100) * (intersection4[2] - 100) - 60 * 60;
      D = B * B - 4 * A * C;
      if (D >= 0) {
        if ((-B - sqrt(D)) / (2 * A) > 0 && (-B + sqrt(D)) / (2 * A) > 0) shadow4 = true;
      }
    }

    // check sphere2
    if (choice4 == 0) {
      A = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
      B = 2 * (x - 15) * v[0] + 2 * (y - 60) * v[1] + 2 * (z - 5) * v[2];
      C = (x - 15) * (x - 15) + (y - 60) * (y - 60) + (z - 5) * (z - 5) - 45 * 45;
      D = B * B - 4 * A * C;
      if (D >= 0) choice4 = 2;
      coefficient = (-B - sqrt(D)) / (2 * A);
      if (coefficient <= 0) coefficient = (B - sqrt(D)) / (2 * A);
      intersection4[0] = x + v[0] * coefficient;
      intersection4[1] = y + v[1] * coefficient;
      intersection4[2] = z + v[2] * coefficient;

      // check shadow
      vTemp[0] = light[0] - intersection4[0];
      vTemp[1] = light[1] - intersection4[1];
      vTemp[2] = light[2] - intersection4[2];
      A = vTemp[0] * vTemp[0] + vTemp[1] * vTemp[1] + vTemp[2] * vTemp[2];
      B = 2 * (intersection4[0] - 15) * vTemp[0] + 2 * (intersection4[1] - 60) * vTemp[1] + 2 * (intersection4[2] - 5) * vTemp[2];
      C = (intersection4[0] - 15) * (intersection4[0] - 15) + (intersection4[1] - 60) * (intersection4[1] - 60) + (intersection4[2] - 5) * (intersection4[2] - 5) - 45 * 45;
      D = B * B - 4 * A * C;
      if (D >=0) {
        if ((-B - sqrt(D)) / (2 * A) > 0 && (-B + sqrt(D)) / (2 * A) > 0) shadow4 = true;
      }

      A = vTemp[0] * vTemp[0] + vTemp[1] * vTemp[1] + vTemp[2] * vTemp[2];
      B = 2 * (intersection4[0] + 55) * vTemp[0] + 2 * (intersection4[1] - 100) * vTemp[1] + 2 * (intersection4[2] - 100) * vTemp[2];
      C = (intersection4[0] + 55) * (intersection4[0] + 55) + (intersection4[1] - 100) * (intersection4[1] - 100) + (intersection4[2] - 100) * (intersection4[2] - 100) - 60 * 60;
      D = B * B - 4 * A * C;
      if (D >= 0) {
        if ((-B - sqrt(D)) / (2 * A) > 0 && (-B + sqrt(D)) / (2 * A) > 0) shadow4 = true;
      }
    }


    // check floor
    if (choice4 == 0) {
      D = -y / v[1];
      if (D >= 0) {
        xx = x + v[0] * D;
        yy = y + v[1] * D;
        zz = z + v[2] * D;
        if (xx <= 200 && xx >= -200 && zz <= 250 && zz >= -250) {
          choice4 = 3;
          intersection4[0] = xx;
          intersection4[1] = yy;
          intersection4[2] = zz;
        }
      }

      // check shadow
      vTemp[0] = light[0] - intersection4[0];
      vTemp[1] = light[1] - intersection4[1];
      vTemp[2] = light[2] - intersection4[2];
      A = vTemp[0] * vTemp[0] + vTemp[1] * vTemp[1] + vTemp[2] * vTemp[2];
      B = 2 * (intersection4[0] - 15) * vTemp[0] + 2 * (intersection4[1] - 60) * vTemp[1] + 2 * (intersection4[2] - 5) * vTemp[2];
      C = (intersection4[0] - 15) * (intersection4[0] - 15) + (intersection4[1] - 60) * (intersection4[1] - 60) + (intersection4[2] - 5) * (intersection4[2] - 5) - 45 * 45;
      D = B * B - 4 * A * C;
      if (D >=0) {
        if ((-B - sqrt(D)) / (2 * A) > 0 && (-B + sqrt(D)) / (2 * A) > 0) shadow4 = true;
      }

      A = vTemp[0] * vTemp[0] + vTemp[1] * vTemp[1] + vTemp[2] * vTemp[2];
      B = 2 * (intersection4[0] + 55) * vTemp[0] + 2 * (intersection4[1] - 100) * vTemp[1] + 2 * (intersection4[2] - 100) * vTemp[2];
      C = (intersection4[0] + 55) * (intersection4[0] + 55) + (intersection4[1] - 100) * (intersection4[1] - 100) + (intersection4[2] - 100) * (intersection4[2] - 100) - 60 * 60;
      D = B * B - 4 * A * C;
      if (D >= 0) {
        if ((-B - sqrt(D)) / (2 * A) > 0 && (-B + sqrt(D)) / (2 * A) > 0) shadow4 = true;
      }
    }



    // draw scene
    switch(choice1) {
    case 0:
      color1 = new float[]{69.0, 150.0, 243.0};
      break;
    case 1:
      // direction of incoming light
      S = new float[]{light[0] - intersection1[0], light[1] - intersection1[1], light[2] - intersection1[2]};
      // normalize S
      temp = sqrt(S[0] * S[0] + S[1] * S[1] + S[2] * S[2]);
      S[0] = S[0] / temp;
      S[1] = S[1] / temp;
      S[2] = S[2] / temp;

      // normal vector
      N = new float[]{intersection1[0] + 55.0, intersection1[1] - 100.0, intersection1[2] - 100.0};
      // normalize N
      temp = sqrt(N[0] * N[0] + N[1] * N[1] + N[2] * N[2]);
      N[0] = N[0] / temp;
      N[1] = N[1] / temp;
      N[2] = N[2] / temp;

      // destination (to the camera)
      V = new float[]{eye[0] - intersection1[0], eye[1] - intersection1[1], eye[2] - intersection1[2]};
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
      if (dotTemp >= 0) diffuse = new float[]{150.0 * dotTemp * kd, 150.0 * dotTemp * kd, 150.0 * dotTemp * kd};
      else diffuse = new float[]{0.0, 0.0, 0.0};
      dotTemp = R[0] * V[0] + R[1] * V[1] + R[2] * V[2];
      if (dotTemp >= 0) specDot = pow(dotTemp, specExp);
      else specDot = 0.0;
      specular = new float[]{255.0 * specDot * ks, 255.0 * specDot * ks, 255.0 * specDot * ks};
      if (!shadow1) color1 = new float[]{ambient[0] + diffuse[0] + specular[0], ambient[1] + diffuse[1] + specular[1], ambient[2] + diffuse[2] + specular[2]};
      else color1 = new float[]{ambient[0], ambient[1], ambient[2]};
      break;
    case 2:
      // direction of incoming light
      S = new float[]{light[0] - intersection1[0], light[1] - intersection1[1], light[2] - intersection1[2]};
      // normalize S
      temp = sqrt(S[0] * S[0] + S[1] * S[1] + S[2] * S[2]);
      S[0] = S[0] / temp;
      S[1] = S[1] / temp;
      S[2] = S[2] / temp;

      // normal vector
      N = new float[]{intersection1[0] - 15.0, intersection1[1] - 60.0, intersection1[2] - 5.0};
      // normalize N
      temp = sqrt(N[0] * N[0] + N[1] * N[1] + N[2] * N[2]);
      N[0] = N[0] / temp;
      N[1] = N[1] / temp;
      N[2] = N[2] / temp;

      // destination (to the camera)
      V = new float[]{eye[0] - intersection1[0], eye[1] - intersection1[1], eye[2] - intersection1[2]};
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
      if (dotTemp >= 0) diffuse = new float[]{250.0 * dotTemp * kd, 250.0 * dotTemp * kd, 250.0 * dotTemp * kd};
      else diffuse = new float[]{0.0, 0.0, 0.0};
      dotTemp = R[0] * V[0] + R[1] * V[1] + R[2] * V[2];
      if (dotTemp >= 0) specDot = pow(dotTemp, specExp);
      else specDot = 0.0;
      specular = new float[]{255.0 * specDot * ks, 255.0 * specDot * ks, 255.0 * specDot * ks};
      if (!shadow1) color1 = new float[]{ambient[0] + diffuse[0] + specular[0], ambient[1] + diffuse[1] + specular[1], ambient[2] + diffuse[2] + specular[2]};
      else color1 = new float[]{ambient[0], ambient[1], ambient[2]};
      break;
    case 3:
      // direction of incoming light
      S = new float[]{light[0] - intersection1[0], light[1] - intersection1[1], light[2] - intersection1[2]};
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
      V = new float[]{eye[0] - intersection1[0], eye[1] - intersection1[1], eye[2] - intersection1[2]};
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
      if (dotTemp >= 0) diffuse = new float[]{255.0 * dotTemp * kd, 228.0 * dotTemp * kd, 108.0 * dotTemp * kd};
      else diffuse = new float[]{0.0, 0.0, 0.0};
      dotTemp = R[0] * V[0] + R[1] * V[1] + R[2] * V[2];
      if (dotTemp >= 0) specDot = pow(dotTemp, specExp);
      else specDot = 0.0;
      specular = new float[]{255.0 * specDot * ks, 255.0 * specDot * ks, 255.0 * specDot * ks};
      if (!shadow1) color1 = new float[]{ambient[0] + diffuse[0] + specular[0], ambient[1] + diffuse[1] + specular[1], ambient[2] + diffuse[2] + specular[2]};
      else color1 = new float[]{ambient[0], ambient[1], ambient[2]};
      break;
    }

    switch(choice2) {
    case 0:
      color2 = new float[]{69.0, 150.0, 243.0};
      break;
    case 1:
      // direction of incoming light
      S = new float[]{light[0] - intersection2[0], light[1] - intersection2[1], light[2] - intersection2[2]};
      // normalize S
      temp = sqrt(S[0] * S[0] + S[1] * S[1] + S[2] * S[2]);
      S[0] = S[0] / temp;
      S[1] = S[1] / temp;
      S[2] = S[2] / temp;

      // normal vector
      N = new float[]{intersection2[0] + 55.0, intersection2[1] - 100.0, intersection2[2] - 100.0};
      // normalize N
      temp = sqrt(N[0] * N[0] + N[1] * N[1] + N[2] * N[2]);
      N[0] = N[0] / temp;
      N[1] = N[1] / temp;
      N[2] = N[2] / temp;

      // destination (to the camera)
      V = new float[]{eye[0] - intersection2[0], eye[1] - intersection2[1], eye[2] - intersection2[2]};
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
      if (dotTemp >= 0) diffuse = new float[]{150.0 * dotTemp * kd, 150.0 * dotTemp * kd, 150.0 * dotTemp * kd};
      else diffuse = new float[]{0.0, 0.0, 0.0};
      dotTemp = R[0] * V[0] + R[1] * V[1] + R[2] * V[2];
      if (dotTemp >= 0) specDot = pow(dotTemp, specExp);
      else specDot = 0.0;
      specular = new float[]{255.0 * specDot * ks, 255.0 * specDot * ks, 255.0 * specDot * ks};
      if (!shadow2) color2 = new float[]{ambient[0] + diffuse[0] + specular[0], ambient[1] + diffuse[1] + specular[1], ambient[2] + diffuse[2] + specular[2]};
      else color2 = new float[]{ambient[0], ambient[1], ambient[2]};
      break;
    case 2:
      // direction of incoming light
      S = new float[]{light[0] - intersection2[0], light[1] - intersection2[1], light[2] - intersection2[2]};
      // normalize S
      temp = sqrt(S[0] * S[0] + S[1] * S[1] + S[2] * S[2]);
      S[0] = S[0] / temp;
      S[1] = S[1] / temp;
      S[2] = S[2] / temp;

      // normal vector
      N = new float[]{intersection2[0] - 15.0, intersection2[1] - 60.0, intersection2[2] - 5.0};
      // normalize N
      temp = sqrt(N[0] * N[0] + N[1] * N[1] + N[2] * N[2]);
      N[0] = N[0] / temp;
      N[1] = N[1] / temp;
      N[2] = N[2] / temp;

      // destination (to the camera)
      V = new float[]{eye[0] - intersection2[0], eye[1] - intersection2[1], eye[2] - intersection2[2]};
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
      if (dotTemp >= 0) diffuse = new float[]{250.0 * dotTemp * kd, 250.0 * dotTemp * kd, 250.0 * dotTemp * kd};
      else diffuse = new float[]{0.0, 0.0, 0.0};
      dotTemp = R[0] * V[0] + R[1] * V[1] + R[2] * V[2];
      if (dotTemp >= 0) specDot = pow(dotTemp, specExp);
      else specDot = 0.0;
      specular = new float[]{255.0 * specDot * ks, 255.0 * specDot * ks, 255.0 * specDot * ks};
      if (!shadow2) color2 = new float[]{ambient[0] + diffuse[0] + specular[0], ambient[1] + diffuse[1] + specular[1], ambient[2] + diffuse[2] + specular[2]};
      else color2 = new float[]{ambient[0], ambient[1], ambient[2]};
      break;
    case 3:
      // direction of incoming light
      S = new float[]{light[0] - intersection2[0], light[1] - intersection2[1], light[2] - intersection2[2]};
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
      V = new float[]{eye[0] - intersection2[0], eye[1] - intersection2[1], eye[2] - intersection2[2]};
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
      if (dotTemp >= 0) diffuse = new float[]{255.0 * dotTemp * kd, 228.0 * dotTemp * kd, 108.0 * dotTemp * kd};
      else diffuse = new float[]{0.0, 0.0, 0.0};
      dotTemp = R[0] * V[0] + R[1] * V[1] + R[2] * V[2];
      if (dotTemp >= 0) specDot = pow(dotTemp, specExp);
      else specDot = 0.0;
      specular = new float[]{255.0 * specDot * ks, 255.0 * specDot * ks, 255.0 * specDot * ks};
      if (!shadow2) color2 = new float[]{ambient[0] + diffuse[0] + specular[0], ambient[1] + diffuse[1] + specular[1], ambient[2] + diffuse[2] + specular[2]};
      else color2 = new float[]{ambient[0], ambient[1], ambient[2]};
      break;
    }

    switch(choice3) {
    case 0:
      color3 = new float[]{69.0, 150.0, 243.0};
      break;
    case 1:
      // direction of incoming light
      S = new float[]{light[0] - intersection3[0], light[1] - intersection3[1], light[2] - intersection3[2]};
      // normalize S
      temp = sqrt(S[0] * S[0] + S[1] * S[1] + S[2] * S[2]);
      S[0] = S[0] / temp;
      S[1] = S[1] / temp;
      S[2] = S[2] / temp;

      // normal vector
      N = new float[]{intersection3[0] + 55.0, intersection3[1] - 100.0, intersection3[2] - 100.0};
      // normalize N
      temp = sqrt(N[0] * N[0] + N[1] * N[1] + N[2] * N[2]);
      N[0] = N[0] / temp;
      N[1] = N[1] / temp;
      N[2] = N[2] / temp;

      // destination (to the camera)
      V = new float[]{eye[0] - intersection3[0], eye[1] - intersection3[1], eye[2] - intersection3[2]};
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
      if (dotTemp >= 0) diffuse = new float[]{150.0 * dotTemp * kd, 150.0 * dotTemp * kd, 150.0 * dotTemp * kd};
      else diffuse = new float[]{0.0, 0.0, 0.0};
      dotTemp = R[0] * V[0] + R[1] * V[1] + R[2] * V[2];
      if (dotTemp >= 0) specDot = pow(dotTemp, specExp);
      else specDot = 0.0;
      specular = new float[]{255.0 * specDot * ks, 255.0 * specDot * ks, 255.0 * specDot * ks};
      if (!shadow3) color3 = new float[]{ambient[0] + diffuse[0] + specular[0], ambient[1] + diffuse[1] + specular[1], ambient[2] + diffuse[2] + specular[2]};
      else color3 = new float[]{ambient[0], ambient[1], ambient[2]};
      break;
    case 2:
      // direction of incoming light
      S = new float[]{light[0] - intersection3[0], light[1] - intersection3[1], light[2] - intersection3[2]};
      // normalize S
      temp = sqrt(S[0] * S[0] + S[1] * S[1] + S[2] * S[2]);
      S[0] = S[0] / temp;
      S[1] = S[1] / temp;
      S[2] = S[2] / temp;

      // normal vector
      N = new float[]{intersection3[0] - 15.0, intersection3[1] - 60.0, intersection3[2] - 5.0};
      // normalize N
      temp = sqrt(N[0] * N[0] + N[1] * N[1] + N[2] * N[2]);
      N[0] = N[0] / temp;
      N[1] = N[1] / temp;
      N[2] = N[2] / temp;

      // destination (to the camera)
      V = new float[]{eye[0] - intersection3[0], eye[1] - intersection3[1], eye[2] - intersection3[2]};
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
      if (dotTemp >= 0) diffuse = new float[]{250.0 * dotTemp * kd, 250.0 * dotTemp * kd, 250.0 * dotTemp * kd};
      else diffuse = new float[]{0.0, 0.0, 0.0};
      dotTemp = R[0] * V[0] + R[1] * V[1] + R[2] * V[2];
      if (dotTemp >= 0) specDot = pow(dotTemp, specExp);
      else specDot = 0.0;
      specular = new float[]{255.0 * specDot * ks, 255.0 * specDot * ks, 255.0 * specDot * ks};
      if (!shadow3) color3 = new float[]{ambient[0] + diffuse[0] + specular[0], ambient[1] + diffuse[1] + specular[1], ambient[2] + diffuse[2] + specular[2]};
      else color3 = new float[]{ambient[0], ambient[1], ambient[2]};
      break;
    case 3:
      // direction of incoming light
      S = new float[]{light[0] - intersection3[0], light[1] - intersection3[1], light[2] - intersection3[2]};
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
      V = new float[]{eye[0] - intersection3[0], eye[1] - intersection3[1], eye[2] - intersection3[2]};
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
      if (dotTemp >= 0) diffuse = new float[]{255.0 * dotTemp * kd, 228.0 * dotTemp * kd, 108.0 * dotTemp * kd};
      else diffuse = new float[]{0.0, 0.0, 0.0};
      dotTemp = R[0] * V[0] + R[1] * V[1] + R[2] * V[2];
      if (dotTemp >= 0) specDot = pow(dotTemp, specExp);
      else specDot = 0.0;
      specular = new float[]{255.0 * specDot * ks, 255.0 * specDot * ks, 255.0 * specDot * ks};
      if (!shadow3) color3 = new float[]{ambient[0] + diffuse[0] + specular[0], ambient[1] + diffuse[1] + specular[1], ambient[2] + diffuse[2] + specular[2]};
      else color3 = new float[]{ambient[0], ambient[1], ambient[2]};
      break;
    }

    switch(choice4) {
    case 0:
      color4 = new float[]{69.0, 150.0, 243.0};
      break;
    case 1:
      // direction of incoming light
      S = new float[]{light[0] - intersection4[0], light[1] - intersection4[1], light[2] - intersection4[2]};
      // normalize S
      temp = sqrt(S[0] * S[0] + S[1] * S[1] + S[2] * S[2]);
      S[0] = S[0] / temp;
      S[1] = S[1] / temp;
      S[2] = S[2] / temp;

      // normal vector
      N = new float[]{intersection4[0] + 55.0, intersection4[1] - 100.0, intersection4[2] - 100.0};
      // normalize N
      temp = sqrt(N[0] * N[0] + N[1] * N[1] + N[2] * N[2]);
      N[0] = N[0] / temp;
      N[1] = N[1] / temp;
      N[2] = N[2] / temp;

      // destination (to the camera)
      V = new float[]{eye[0] - intersection4[0], eye[1] - intersection4[1], eye[2] - intersection4[2]};
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
      if (dotTemp >= 0) diffuse = new float[]{150.0 * dotTemp * kd, 150.0 * dotTemp * kd, 150.0 * dotTemp * kd};
      else diffuse = new float[]{0.0, 0.0, 0.0};
      dotTemp = R[0] * V[0] + R[1] * V[1] + R[2] * V[2];
      if (dotTemp >= 0) specDot = pow(dotTemp, specExp);
      else specDot = 0.0;
      specular = new float[]{255.0 * specDot * ks, 255.0 * specDot * ks, 255.0 * specDot * ks};
      if (!shadow4) color4 = new float[]{ambient[0] + diffuse[0] + specular[0], ambient[1] + diffuse[1] + specular[1], ambient[2] + diffuse[2] + specular[2]};
      else color4 = new float[]{ambient[0], ambient[1], ambient[2]};
      break;
    case 2:
      // direction of incoming light
      S = new float[]{light[0] - intersection4[0], light[1] - intersection4[1], light[2] - intersection4[2]};
      // normalize S
      temp = sqrt(S[0] * S[0] + S[1] * S[1] + S[2] * S[2]);
      S[0] = S[0] / temp;
      S[1] = S[1] / temp;
      S[2] = S[2] / temp;

      // normal vector
      N = new float[]{intersection4[0] - 15.0, intersection4[1] - 60.0, intersection4[2] - 5.0};
      // normalize N
      temp = sqrt(N[0] * N[0] + N[1] * N[1] + N[2] * N[2]);
      N[0] = N[0] / temp;
      N[1] = N[1] / temp;
      N[2] = N[2] / temp;

      // destination (to the camera)
      V = new float[]{eye[0] - intersection4[0], eye[1] - intersection4[1], eye[2] - intersection4[2]};
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
      if (dotTemp >= 0) diffuse = new float[]{250.0 * dotTemp * kd, 250.0 * dotTemp * kd, 250.0 * dotTemp * kd};
      else diffuse = new float[]{0.0, 0.0, 0.0};
      dotTemp = R[0] * V[0] + R[1] * V[1] + R[2] * V[2];
      if (dotTemp >= 0) specDot = pow(dotTemp, specExp);
      else specDot = 0.0;
      specular = new float[]{255.0 * specDot * ks, 255.0 * specDot * ks, 255.0 * specDot * ks};
      if (!shadow4) color4 = new float[]{ambient[0] + diffuse[0] + specular[0], ambient[1] + diffuse[1] + specular[1], ambient[2] + diffuse[2] + specular[2]};
      else color4 = new float[]{ambient[0], ambient[1], ambient[2]};
      break;
    case 3:
      // direction of incoming light
      S = new float[]{light[0] - intersection4[0], light[1] - intersection4[1], light[2] - intersection4[2]};
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
      V = new float[]{eye[0] - intersection4[0], eye[1] - intersection4[1], eye[2] - intersection4[2]};
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
      if (dotTemp >= 0) diffuse = new float[]{255.0 * dotTemp * kd, 228.0 * dotTemp * kd, 108.0 * dotTemp * kd};
      else diffuse = new float[]{0.0, 0.0, 0.0};
      dotTemp = R[0] * V[0] + R[1] * V[1] + R[2] * V[2];
      if (dotTemp >= 0) specDot = pow(dotTemp, specExp);
      else specDot = 0.0;
      specular = new float[]{255.0 * specDot * ks, 255.0 * specDot * ks, 255.0 * specDot * ks};
      if (!shadow4) color4 = new float[]{ambient[0] + diffuse[0] + specular[0], ambient[1] + diffuse[1] + specular[1], ambient[2] + diffuse[2] + specular[2]};
      else color4 = new float[]{ambient[0], ambient[1], ambient[2]};
      break;
    }

    c = new float[3];
    c[0] = (color1[0] + color2[0] + color3[0]+ color4[0]) / 4.0;
    c[1] = (color1[1] + color2[1] + color3[1]+ color4[1]) / 4.0;
    c[2] = (color1[2] + color2[2] + color3[2]+ color4[2]) / 4.0;
    set(i, j, color(c[0], c[1], c[2]));
  }
}

//save("Assignment_3_Supersampling.png");
