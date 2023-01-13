size(560, 315);

// finding the near clipping plane (upper left pixel position)
float[] eye = new float[]{-30.0, 150.0, 275.0}; // camera position
float[] lookAt = new float[]{-40.0, 75.0, 100.0}; // camera look at

// for extras (move the camera)
eye = new float[]{300.0, 150.0, 100.0};
lookAt = new float[]{0.0, 0.0, 0.0};

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
/*println(upper_left[0] + " " + upper_left[1] + " " + upper_left[2]);
 println((upper_left[0] - vector_up[0] * 314) + " " + (upper_left[1] - vector_up[1] * 314) + " " + (upper_left[2] - vector_up[2] * 314));
 println(near);*/

float x, y, z; // pixel position
float A, B, C, D; // coeffecients for finding intersections
float[] v = new float[3];
float xx, yy, zz;
int choice;
/* ray tracing
 sphere1: (x+55)^2 + (y-100)^2 + (z-100)^2 = 60^2
 */
for (int i = 0; i < width; i++) {
  for (int j = 0; j < height; j++) {
    choice = 0;
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
    if (D >= 0) choice = 1;

    // check sphere2
    B = 2 * (x - 15) * v[0] + 2 * (y - 60) * v[1] + 2 * (z - 5) * v[2];
    C = (x - 15) * (x - 15) + (y - 60) * (y - 60) + (z - 5) * (z - 5) - 45 * 45;
    D = B * B - 4 * A * C;
    if (D >= 0 && choice == 0) choice = 2;

    // check cylinder
    if (choice == 0) {
      for (int k = 75; k >= 0; k--) {
        D = (k - y) / v[1];
        if (D >= 0) {
          xx = x + v[0] * D;
          yy = y + v[1] * D;
          zz = z + v[2] * D;
          if ((xx - 100) * (xx - 100) + (yy - k) * (yy - k) + (zz + 100) * (zz + 100) <= 30 * 30) {
            choice = 4;
            break;
          }
        }
      }
    }
    
    // check floor
    D = -y / v[1];
    if (D >= 0 && choice == 0) {
      xx = x + v[0] * D;
      yy = y + v[1] * D;
      zz = z + v[2] * D;
      if (xx <= 200 && xx >= -200 && zz <= 250 && zz >= -250) choice = 3;
    }
    
    // draw scene
    switch(choice) {
    case 0:
      set(i, j, color(69, 150, 243));
      break;
    case 1:
      set(i, j, color(150, 150, 150));
      break;
    case 2:
      set(i, j, color(250, 250, 250));
      break;
    case 3:
      set(i, j, color(255, 228, 108));
      break;
    case 4:
      set(i, j, color(128, 0, 0));
      break;
    }
  }
}

//save("Assignment_2_Extras_Michael_Lee.png");
