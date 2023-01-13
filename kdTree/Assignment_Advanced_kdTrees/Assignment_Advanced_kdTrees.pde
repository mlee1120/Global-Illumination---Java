import java.util.HashSet;
boolean no_intersection;
int index_intersected;
float[] intersection_coor = new float[3];
int le = 0, in = 0;
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

//
HashMap<Integer, HashSet<Integer>> table = new HashMap();

// KD tree
KdTree tree;

float zmax, zmin, xmax, xmin, ymax, ymin;


// to update the canvas or not
boolean toDraw;

float buildstart, buildend, renderstart, renderend;

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
  tree = new KdTree();
  tree.root = new Node();
  buildstart = millis();
  buildTree(tree.root, triangles, 0, triangles.size() - 1, "x");
  for (int i = 0; i < triangles.size(); i++) {
    buildLeaf(tree.root, i, 0);
  }
  buildend = millis();
  println("construct: " + (buildend - buildstart));

  for (int i = 0; i < vertices.size(); i++) {
    if (i == 0) {
      xmax = vertices.get(i)[0];
      xmin = xmax;
      ymax = vertices.get(i)[1];
      ymin = ymax;
      zmax = vertices.get(i)[2];
      zmin = zmax;
    } else {
      if (vertices.get(i)[0] > xmax) xmax = vertices.get(i)[0];
      if (vertices.get(i)[0] < xmin) xmin = vertices.get(i)[0];
      if (vertices.get(i)[1] > ymax) ymax = vertices.get(i)[1];
      if (vertices.get(i)[1] < ymin) ymin = vertices.get(i)[1];
      if (vertices.get(i)[2] > zmax) zmax = vertices.get(i)[2];
      if (vertices.get(i)[2] < zmin) zmin = vertices.get(i)[2];
    }
  }
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
    println("render: " + (renderend - renderstart));
    //save("kd_res4.png");
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
  //float[] vTemp = new float[3];
  // coeffecients for finding intersections
  float A, B, C, D;
  float coefficient;


  pixel_position[0] = upper_left[0] + i * vector_right[0] - j * vector_up[0] + dx * vector_right[0] - dy * vector_up[0];
  pixel_position[1] = upper_left[1] + i * vector_right[1] - j * vector_up[1] + dx * vector_right[1] - dy * vector_up[1];
  pixel_position[2] = upper_left[2] + i * vector_right[2] - j * vector_up[2] + dx * vector_right[2] - dy * vector_up[2];
  vector[0] = pixel_position[0] - eye[0];
  vector[1] = pixel_position[1] - eye[1];
  vector[2] = pixel_position[2] - eye[2];
  vector = normalize(vector);

  // three vertices of a triangle
  float[] vertex0, vertex1, vertex2;
  // vector from vertex1 to vertex0 and vector from vertex1 to vertex2
  float[] vector0 = new float[3];
  float[] vector1 = new float[3];
  // normal of the triangle
  float[] norm = new float[3];
  float[] intersection = new float[3];

  // auxiliary parameters for checking if the intersection is inside the triangle
  float uu, vv, a0, b0, c0, a1, b1, c1, a2, b2, c2;
  // auxiliary parameter for finding intersections


  boolean x, y;
  if ((xmax - pixel_position[0]) / vector[0] > (xmin - pixel_position[0]) / vector[0]) x = true;
  else x = false;
  if ((ymax - pixel_position[1]) / vector[1] > (ymin - pixel_position[1]) / vector[1]) y = true;
  else y = false;
  float coefficient1 = (zmax - pixel_position[2]) / vector[2];

  HashSet<Integer> checked = new HashSet();
  HashSet<Integer> now = new HashSet();
  HashSet<Integer> to_check = new HashSet();
  no_intersection = true;
  boolean first = true;
  do {
    to_check.clear();
    now = searchTree(tree.root, new float[]{pixel_position[0] + coefficient1 * vector[0], pixel_position[1] + coefficient1 * vector[1], pixel_position[2] + coefficient1 * vector[2]}, 0);
    for (int k : now) {
      if (!checked.contains(k)) {
        checked.add(k);
        to_check.add(k);
        for (int l = 0; l < 3; l++) {
          int temp = triangles.get(k)[l];
          for (int m : table.get(temp)) {
            if (!checked.contains(m)) to_check.add(m);
          }
        }
      }
    }
    for (int k : to_check) {
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
      intersection[0] = pixel_position[0] + vector[0] * coefficient;
      intersection[1] = pixel_position[1] + vector[1] * coefficient;
      intersection[2] = pixel_position[2] + vector[2] * coefficient;
      a0 = vertex0[0] - vertex2[0];
      b0 = vertex1[0] - vertex2[0];
      c0 = intersection[0] - vertex2[0];
      a1 = vertex0[1] - vertex2[1];
      b1 = vertex1[1] - vertex2[1];
      c1 = intersection[1] - vertex2[1];
      a2 = vertex0[2] - vertex2[2];
      b2 = vertex1[2] - vertex2[2];
      c2 = intersection[2] - vertex2[2];
      if (b0/a0 - b1/a1 == 0)println(546);
      vv = (c0/a0 - c1/a1) / (b0/a0 - b1/a1);
      uu = (c0 - b0 * vv) / a0;
      vv = (c2/a2 - c1/a1) / (b2/a2 - b1/a1);
      uu = (c2 - b2 * vv) / a2;

      if (vv >= 0.0 && vv <= 1.0 && uu >= 0.0 && uu <= 1.0 && uu + vv <= 1.0) {
        if (first) {
          intersection_coor[0] = intersection[0];
          intersection_coor[1] = intersection[1];
          intersection_coor[2] = intersection[2];
          index_intersected = k;
          first = false;
        } else {
          float d1 = sqrt(pow(intersection[0] - eye[0], 2) + pow(intersection[1] - eye[1], 2) + pow(intersection[1] - eye[1], 2));
          float d2 = sqrt(pow(intersection_coor[0] - eye[0], 2) + pow(intersection_coor[1] - eye[1], 2) + pow(intersection_coor[1] - eye[1], 2));
          if (d1 < d2) {
            intersection_coor[0] = intersection[0];
            intersection_coor[1] = intersection[1];
            intersection_coor[2] = intersection[2];
            index_intersected = k;
          }
        }
      }
    }
    if (!first) {
      no_intersection = false;
      break;
    }
    coefficient1 += 0.1;
  }
  while (pixel_position[2] + coefficient1 * vector[2] >= zmin && no_intersection && xy(x, y, pixel_position, coefficient1, vector));
}


void calculate_color() {
  // reflection coefficients
  float ka= 0.45, kd = 0.45, ks = 0.2, specExp = 10.0;
  float[] ambient, diffuse, specular, c = new float[]{0.0, 0.0, 0.0};
  float[] S, N, V, R;
  float dotTemp, specDot;
  float[] vertex0, vertex1, vertex2, vector0 = new float[3], vector1 = new float[3];
  if (no_intersection) colors.add(new float[]{69.0, 150.0, 243.0});
  else {
    vertex0 = vertices.get(triangles.get(index_intersected)[0]);
    vertex1 = vertices.get(triangles.get(index_intersected)[1]);
    vertex2 = vertices.get(triangles.get(index_intersected)[2]);
    vector0[0] = vertex0[0] - vertex1[0];
    vector0[1] = vertex0[1] - vertex1[1];
    vector0[2] = vertex0[2] - vertex1[2];
    vector1[0] = vertex2[0] - vertex1[0];
    vector1[1] = vertex2[1] - vertex1[1];
    vector1[2] = vertex2[2] - vertex1[2];
    vector0 = normalize(vector0);
    vector1 = normalize(vector1);
    N = normalize(crossProduct(vector1, vector0));
    V = new float[]{eye[0] - intersection_coor[0], eye[1] - intersection_coor[1], eye[2] - intersection_coor[2]};
    V = normalize(V);
    if (dotProduct(N, V) < 0) N = normalize(crossProduct(vector0, vector1));

    ambient = new float[]{255.0 * ka, 255.0 * ka, 255.0 * ka};
    c[0] += ambient[0];
    c[1] += ambient[1];
    c[2] += ambient[2];
    S = new float[]{lights.get(0)[0] - intersection_coor[0], lights.get(0)[1] - intersection_coor[1], lights.get(0)[2] - intersection_coor[2]};
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
  //lines = loadStrings("bun_zipper.ply");
  String[] line;
  float scale = 1000.0;
  float y_trans = -60.0;
  float x_trans = -12.0;
  HashSet<String> duplicate = new HashSet();
  for (String s : lines) {
    line = s.split(" ");
    if (line.length == 5) {
      vertices.add(new float[3]);
      vertices.get(vertices.size() - 1)[0] = Float.parseFloat(line[0]) * scale + x_trans;
      vertices.get(vertices.size() - 1)[1] = Float.parseFloat(line[1]) * scale + y_trans;
      vertices.get(vertices.size() - 1)[2] = Float.parseFloat(line[2]) * scale;
    } else {
      if (!duplicate.contains(s)) {
        duplicate.add(s);
        duplicate.add(line[0] + " " + line[1] + " " + line[3] + " " + line[2] + " ");
        duplicate.add(line[0] + " " + line[2] + " " + line[1] + " " + line[3] + " ");
        duplicate.add(line[0] + " " + line[2] + " " + line[3] + " " + line[1] + " ");
        duplicate.add(line[0] + " " + line[3] + " " + line[1] + " " + line[2] + " ");
        duplicate.add(line[0] + " " + line[3] + " " + line[2] + " " + line[1] + " ");
        triangles.add(new int[3]);
        triangles.get(triangles.size() - 1)[0] = Integer.parseInt(line[1]);
        triangles.get(triangles.size() - 1)[1] = Integer.parseInt(line[2]);
        triangles.get(triangles.size() - 1)[2] = Integer.parseInt(line[3]);
      }
    }
  }
  for (int i = 0; i < triangles.size(); i++) {
    for (int j = 0; j < 3; j++) {
      if (!table.containsKey(triangles.get(i)[j])) table.put(triangles.get(i)[j], new HashSet());
      table.get(triangles.get(i)[j]).add(i);
    }
  }
}

ArrayList<int[]> mergeSort(ArrayList<int[]> aList, int start, int end, String xyz) {
  ArrayList<int[]>  sortedList = new ArrayList();
  if (start == end) sortedList.add(aList.get(start));
  else {
    int middle = (start + end) / 2;
    ArrayList<int[]> leftList = mergeSort(aList, start, middle, xyz);
    ArrayList<int[]> rightList = mergeSort(aList, middle + 1, end, xyz);
    int leftIndex = 0, rightIndex = 0;
    int sort;
    if (xyz.equals("x")) sort = 0;
    else if (xyz.equals("y")) sort = 1;
    else sort = 2;
    float leftSum, rightSum;
    while (leftIndex < leftList.size() && rightIndex < rightList.size()) {
      leftSum = 0;
      rightSum = 0;
      leftSum += vertices.get(leftList.get(leftIndex)[0])[sort];
      leftSum += vertices.get(leftList.get(leftIndex)[1])[sort];
      leftSum += vertices.get(leftList.get(leftIndex)[2])[sort];
      rightSum += vertices.get(rightList.get(rightIndex)[0])[sort];
      rightSum += vertices.get(rightList.get(rightIndex)[1])[sort];
      rightSum += vertices.get(rightList.get(rightIndex)[2])[sort];
      if (leftSum / 3.0 < rightSum / 3.0) {
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

void buildTree(Node aNode, ArrayList<int[]> aList, int start, int end, String xyz) {
  if (start == end) aNode.name = "leaf";
  else {
    // sort accordingly
    ArrayList<int[]> tempList = mergeSort(aList, start, end, xyz);
    int middle = (tempList.size() - 1) / 2;
    // checking same priority
    float sum, sum_next;
    int sort;
    while (middle < tempList.size() - 1) {
      sum = 0.0;
      sum_next = 0.0;
      if (xyz.equals("x")) sort = 0;
      else if (xyz.equals("y")) sort = 1;
      else sort = 2;
      sum += vertices.get(tempList.get(middle)[0])[sort];
      sum += vertices.get(tempList.get(middle)[1])[sort];
      sum += vertices.get(tempList.get(middle)[2])[sort];
      sum_next += vertices.get(tempList.get(middle + 1)[0])[sort];
      sum_next += vertices.get(tempList.get(middle + 1)[1])[sort];
      sum_next += vertices.get(tempList.get(middle + 1)[2])[sort];
      if (sum / 3.0 == sum_next / 3.0) {
        middle += 1;
      } else break;
    }


    // add a node
    aNode.indices[0] = tempList.get(middle)[0];
    aNode.indices[1] = tempList.get(middle)[1];
    aNode.indices[2] = tempList.get(middle)[2];

    aNode.name = "internal";

    // not necessary, but good for debugging
    tree.size += 1;

    // handle the left and the right sub-trees
    aNode.leftChild = new Node();
    aNode.rightChild = new Node();

    if (xyz.equals("x")) {
      buildTree(aNode.leftChild, tempList, 0, middle, "y");
      if (middle == tempList.size() - 1) aNode.rightChild.name = "no";
      else buildTree(aNode.rightChild, tempList, middle + 1, tempList.size() - 1, "y");
    } else if (xyz.equals("y")) {
      buildTree(aNode.leftChild, tempList, 0, middle, "z");
      if (middle == tempList.size() - 1) aNode.rightChild.name = "no";
      else buildTree(aNode.rightChild, tempList, middle + 1, tempList.size() - 1, "z");
    } else if (xyz.equals("z")) {
      buildTree(aNode.leftChild, tempList, 0, middle, "x");
      if (middle == tempList.size() - 1) aNode.rightChild.name = "no";
      else buildTree(aNode.rightChild, tempList, middle + 1, tempList.size() - 1, "x");
    }
  }
}

void buildLeaf(Node current, int index, int xyz) {
  if (current.name.equals("internal")) {
    float sumNode = 0.0, sumTri = 0.0;
    for (int i = 0; i < 3; i++) {
      sumNode += vertices.get(current.indices[i])[xyz];
      sumTri += vertices.get(triangles.get(index)[i])[xyz];
    }
    xyz += 1;
    if (xyz == 3) xyz = 0;
    if (sumNode / 3.0 < sumTri / 3.0) buildLeaf(current.rightChild, index, xyz);
    else buildLeaf(current.leftChild, index, xyz);
  } else if (current.name.equals("leaf")) {
    if (current.index.size() != 0) println("!!!");
    current.index.add(index);
  } else println("fuck");
}

HashSet<Integer> searchTree(Node current, float[] ray, int xyz) {
  if (current.name.equals("internal")) {
    float temp = 0.0;
    for (int i = 0; i < 3; i++) {
      temp += vertices.get(current.indices[i])[xyz];
    }
    temp /= 3.0;
    if (ray[xyz] > temp) {
      xyz += 1;
      if (xyz == 3) xyz = 0;
      return searchTree(current.rightChild, ray, xyz);
    } else {
      xyz += 1;
      if (xyz == 3) xyz = 0;
      return searchTree(current.leftChild, ray, xyz);
    }
  } else return current.index;
}

boolean xy(boolean x, boolean y, float[] pp, float co, float[] ve) {
  boolean b1, b2;
  if (x) b1 = pp[0] + co * ve[0] <= xmax;
  else b1 = pp[0] + co * ve[0] >= xmin;
  if (y) b2 = pp[1] + co * ve[1] <= ymax;
  else b2 = pp[1] + co * ve[1] >= ymin;
  return b1 && b2;
}
