import java.util.HashSet;

// light position
ArrayList<float[]> lights;

ArrayList<float[]> lights_i;

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

// this function is executed once for setting up some parameters
void setup() {
  // canvas size
  size(900, 900);
  //size(1280, 720);
  // lights' positions
  lights = new ArrayList();
  lights.add(new float[]{-900, 1200, 750});
  lights.add(new float[]{20, 30, 6});
  lights.add(new float[]{-20, 30, 6});
  lights.add(new float[]{0, 90, -5});
  
  lights_i = new ArrayList();
  //lights_i.add(new float[]{1.0, 1.0, 1.0});
  lights_i.add(new float[]{0.5, 0.5, 0.5});
  lights_i.add(new float[]{1.0, 0.75, 0.1});
  lights_i.add(new float[]{1.0, 0.75, 0.1});
  lights_i.add(new float[]{1.0, 0.75, 0.1});
  
  
  // camera position
  eye = new float[]{0.0, 180.0, 720.0};
  eye = new float[]{0.0, 270.0, 720.0};

  // camera look at
  lookAt = new float[]{0.0, 70.0, 1.0};
  lookAt = new float[]{0.0, 160.0, 1.0};

  // create all primitives
  Primitive temp;
  objects = new ArrayList();
  objects.add(new Primitive("background", new float[]{0.0, 200.0, -1000.0}, null, new float[]{0.0, 0.0, 0.0}, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0));

  objects.add(new Primitive("sphere", new float[]{0.0, 135.0, -60.0}, new float[]{360.0}, new float[]{60.0, 60.0, 60.0}, 0.5, 0.5, 0.3, 30.0, 0.0, 0.9, 1.33));

  temp = new Primitive("cylinder", new float[]{-175.0, 30.0, 0.0}, new float[]{60.0, 15.0}, new float[]{216.0, 216.0, 216.0}, 0.5, 0.6, 0.3, 30.0, 0.0, 0.0, 1.00);
  temp.ifTexture = true;
  temp.texture = loadImage("cylinder_left_1.png");
  objects.add(temp);
  temp = new Primitive("cylinder", new float[]{175.0, 30.0, 0.0}, new float[]{60.0, 15.0}, new float[]{216.0, 216.0, 216.0}, 0.5, 0.6, 0.3, 30.0, 0.0, 0.0, 1.00);
  temp.ifTexture = true;
  temp.texture = loadImage("cylinder_right_1.png");
  objects.add(temp);

  temp = new Primitive("cylinder", new float[]{-35.0, 30.0, 0.0}, new float[]{60.0, 10.0}, new float[]{216.0, 216.0, 216.0}, 0.5, 0.6, 0.3, 30.0, 0.0, 0.0, 1.00);
  temp.ifTexture = true;
  temp.texture = loadImage("cylinder_left_2.png");
  objects.add(temp);
  temp = new Primitive("cylinder", new float[]{35.0, 30.0, 0.0}, new float[]{60.0, 10.0}, new float[]{216.0, 216.0, 216.0}, 0.5, 0.6, 0.3, 30.0, 0.0, 0.0, 1.00);
  temp.ifTexture = true;
  temp.texture = loadImage("cylinder_right_2.png");
  objects.add(temp);

  temp = new Primitive("cylinder", new float[]{-130.0, 25.0, -50.0}, new float[]{50.0, 7.5}, new float[]{216.0, 216.0, 216.0}, 0.5, 0.6, 0.3, 30.0, 0.0, 0.0, 1.00);
  temp.ifTexture = true;
  temp.texture = loadImage("cylinder_left_2.png");
  objects.add(temp);

  temp = new Primitive("cylinder", new float[]{-100.0, 40.0, -100.0}, new float[]{80.0, 10.0}, new float[]{216.0, 216.0, 216.0}, 0.5, 0.6, 0.3, 30.0, 0.0, 0.0, 1.00);
  temp.ifTexture = true;
  temp.texture = loadImage("cylinder_left_3.png");
  objects.add(temp);

  temp = new Primitive("cylinder", new float[]{-100.0, 32.0, 0.0}, new float[]{24.0, 5.0}, new float[]{216.0, 216.0, 216.0}, 0.5, 0.6, 0.3, 30.0, 0.0, 0.0, 1.00);
  temp.ifTexture = true;
  temp.texture = loadImage("cylinder_left_4.png");
  objects.add(temp);
  temp = new Primitive("cylinder", new float[]{100.0, 32.0, 0.0}, new float[]{24.0, 5.0}, new float[]{216.0, 216.0, 216.0}, 0.5, 0.6, 0.3, 30.0, 0.0, 0.0, 1.00);
  temp.ifTexture = true;
  temp.texture = loadImage("cylinder_right_4.png");
  objects.add(temp);

  temp = new Primitive("cylinder", new float[]{-87.5, 37.0, -50.0}, new float[]{74.0, 7.5}, new float[]{216.0, 216.0, 216.0}, 0.5, 0.6, 0.3, 30.0, 0.0, 0.0, 1.00);
  temp.ifTexture = true;
  temp.texture = loadImage("cylinder_left_5.png");
  objects.add(temp);
  temp = new Primitive("cylinder", new float[]{87.5, 37.0, -50.0}, new float[]{74.0, 7.5}, new float[]{216.0, 216.0, 216.0}, 0.5, 0.6, 0.3, 30.0, 0.0, 0.0, 1.00);
  temp.ifTexture = true;
  temp.texture = loadImage("cylinder_right_5.png");
  objects.add(temp);

  temp = new Primitive("cylinder", new float[]{110.0, 40.0, -50.0}, new float[]{80.0, 15}, new float[]{216.0, 216.0, 216.0}, 0.5, 0.6, 0.3, 30.0, 0.0, 0.0, 1.00);
  temp.ifTexture = true;
  temp.texture = loadImage("cylinder_right_6.png");
  objects.add(temp);

  temp = new Primitive("cylinder", new float[]{62.5, 60.0, -65.0}, new float[]{120.0, 10.0}, new float[]{216.0, 216.0, 216.0}, 0.5, 0.6, 0.3, 30.0, 0.0, 0.0, 1.00);
  temp.ifTexture = true;
  temp.texture = loadImage("cylinder_right_7.png");
  objects.add(temp);

  temp = new Primitive("cylinder", new float[]{-20.0, 152.0, -65.0}, new float[]{24.0, 5.0}, new float[]{216.0, 216.0, 216.0}, 0.5, 0.6, 0.3, 30.0, 0.0, 0.0, 1.00);
  temp.ifTexture = true;
  temp.texture = loadImage("cylinder_8.png");
  objects.add(temp);
  temp = new Primitive("cylinder", new float[]{20.0, 152.0, -65.0}, new float[]{24.0, 5.0}, new float[]{216.0, 216.0, 216.0}, 0.5, 0.6, 0.3, 30.0, 0.0, 0.0, 1.00);
  temp.ifTexture = true;
  temp.texture = loadImage("cylinder_8.png");
  objects.add(temp);
  temp = new Primitive("cylinder", new float[]{-30.0, 200.0, -120.0}, new float[]{30.0, 5.0}, new float[]{216.0, 216.0, 216.0}, 0.5, 0.6, 0.3, 30.0, 0.0, 0.0, 1.00);
  temp.ifTexture = true;
  temp.texture = loadImage("cylinder_8.png");
  objects.add(temp);

  temp = new Primitive("cylinder", new float[]{-10.0, 235.0, -120.0}, new float[]{50.0, 10.0}, new float[]{216.0, 216.0, 216.0}, 0.5, 0.6, 0.3, 30.0, 0.0, 0.0, 1.00);
  temp.ifTexture = true;
  temp.texture = loadImage("cylinder_9.png");
  objects.add(temp);

  temp = new Primitive("cylinder", new float[]{-10.0, 160.0, -120.0}, new float[]{100.0, 18.0}, new float[]{216.0, 216.0, 216.0}, 0.5, 0.6, 0.3, 30.0, 0.0, 0.0, 1.00);
  temp.ifTexture = true;
  temp.texture = loadImage("cylinder_10.png");
  objects.add(temp);

  temp = new Primitive("cylinder", new float[]{30.0, 265.0, -100.0}, new float[]{30.0, 10.0}, new float[]{216.0, 216.0, 216.0}, 0.5, 0.6, 0.3, 30.0, 0.0, 0.0, 1.00);
  temp.ifTexture = true;
  temp.texture = loadImage("cylinder_11.png");
  objects.add(temp);

  temp = new Primitive("cylinder", new float[]{30.0, 220.0, -100.0}, new float[]{60.0, 15.0}, new float[]{216.0, 216.0, 216.0}, 0.5, 0.6, 0.3, 30.0, 0.0, 0.0, 1.00);
  temp.ifTexture = true;
  temp.texture = loadImage("cylinder_12.png");
  objects.add(temp);

  temp = new Primitive("cylinder", new float[]{30.0, 95.0, -100.0}, new float[]{190.0, 20.0}, new float[]{216.0, 216.0, 216.0}, 0.5, 0.6, 0.3, 30.0, 0.0, 0.0, 1.00);
  temp.ifTexture = true;
  temp.texture = loadImage("cylinder_13.png");
  objects.add(temp);

  objects.add(new Primitive("cone", new float[]{-130.0, 60.0, -50.0}, new float[]{10.0, 7.5}, new float[]{0.0, 104.0, 185.0}, 0.5, 0.6, 0.15, 30.0, 0.0, 0.0, 1.00));
  objects.add(new Primitive("cone", new float[]{-35.0, 80.0, 0.0}, new float[]{20.0, 10.0}, new float[]{0.0, 104.0, 185.0}, 0.5, 0.6, 0.15, 30.0, 0.0, 0.0, 1.00));
  objects.add(new Primitive("cone", new float[]{-175.0, 84.0, 0.0}, new float[]{24.0, 15.0}, new float[]{0.0, 104.0, 185.0}, 0.5, 0.6, 0.15, 30.0, 0.0, 0.0, 1.00));
  objects.add(new Primitive("cone", new float[]{-100.0, 100.0, -100.0}, new float[]{20.0, 10.0}, new float[]{0.0, 104.0, 185.0}, 0.5, 0.6, 0.15, 30.0, 0.0, 0.0, 1.00));
  objects.add(new Primitive("cone", new float[]{-20.0, 184.0, -65.0}, new float[]{20.0, 5.0}, new float[]{0.0, 104.0, 185.0}, 0.5, 0.6, 0.15, 30.0, 0.0, 0.0, 1.00));
  objects.add(new Primitive("cone", new float[]{-30.0, 235.0, -120.0}, new float[]{20.0, 5.0}, new float[]{0.0, 104.0, 185.0}, 0.5, 0.6, 0.15, 30.0, 0.0, 0.0, 1.00));
  objects.add(new Primitive("cone", new float[]{-10.0, 320.0, -120.0}, new float[]{60.0, 10.0}, new float[]{0.0, 104.0, 185.0}, 0.5, 0.6, 0.15, 30.0, 0.0, 0.0, 1.00));
  objects.add(new Primitive("cone", new float[]{35.0, 80.0, 0.0}, new float[]{20.0, 10.0}, new float[]{0.0, 104.0, 185.0}, 0.5, 0.6, 0.15, 30.0, 0.0, 0.0, 1.00));
  objects.add(new Primitive("cone", new float[]{175.0, 84.0, 0.0}, new float[]{24.0, 15.0}, new float[]{0.0, 104.0, 185.0}, 0.5, 0.6, 0.15, 30.0, 0.0, 0.0, 1.00));
  objects.add(new Primitive("cone", new float[]{110.0, 110.0, -50.0}, new float[]{30.0, 15.0}, new float[]{0.0, 104.0, 185.0}, 0.5, 0.6, 0.15, 30.0, 0.0, 0.0, 1.00));
  objects.add(new Primitive("cone", new float[]{62.5, 160.0, -65.0}, new float[]{40.0, 10.0}, new float[]{0.0, 104.0, 185.0}, 0.5, 0.6, 0.15, 30.0, 0.0, 0.0, 1.00));
  objects.add(new Primitive("cone", new float[]{20.0, 184.0, -65.0}, new float[]{20.0, 5.0}, new float[]{0.0, 104.0, 185.0}, 0.5, 0.6, 0.15, 30.0, 0.0, 0.0, 1.00));
  objects.add(new Primitive("cone", new float[]{30.0, 380.0, -100.0}, new float[]{100.0, 10.0}, new float[]{0.0, 104.0, 185.0}, 0.5, 0.6, 0.15, 30.0, 0.0, 0.0, 1.00));

  temp = new Primitive("plane", new float[]{-102.5, 20.0, 0.0}, new float[]{115.0, 40.0, 0.0, 1.0}, new float[]{216.0, 216.0, 216.0}, 0.5, 0.6, 0.3, 30.0, 0.0, 0.0, 1.0);
  temp.ifTexture = true;
  temp.texture = loadImage("left_wall.png");
  objects.add(temp);

  temp = new Primitive("plane", new float[]{0.0, 0.0, 180.0}, new float[]{1080.0, 0.0, 630.0, 1.0}, new float[]{30.0, 30.0, 30.0}, 0.5, 0.5, 0.3, 30.0, 0.3, 0.0, 1.0);
  temp.ifTexture = true;
  temp.texture = loadImage("ground.png");
  objects.add(temp);

  temp = new Primitive("plane", new float[]{102.5, 20.0, 0.0}, new float[]{115.0, 40.0, 0.0, 1.0}, new float[]{216.0, 216.0, 216.0}, 0.5, 0.6, 0.3, 30.0, 0.0, 0.0, 1.0);
  temp.ifTexture = true;
  temp.texture = loadImage("right_wall.png");
  objects.add(temp);

  temp = new Primitive("plane", new float[]{0.0, 270.0, -135.0}, new float[]{1080.0, 540.0, 0.0, 1.0}, new float[]{216.0, 216.0, 0.0}, 0.5, 0.6, 0.3, 30.0, 0.0, 0.0, 1.0);
  temp.ifTexture = true;
  temp.texture = loadImage("star.png");
  objects.add(temp);
  
  temp = new Primitive("plane", new float[]{0.0, 35.0, -50.0}, new float[]{160.0, 70.0, 0.0, 1.0}, new float[]{216.0, 216.0, 216.0}, 0.5, 0.6, 0.3, 30.0, 0.0, 0.0, 1.0);
  temp.ifTexture = true;
  temp.texture = loadImage("wall_2.png");
  objects.add(temp);
  
  temp = new Primitive("plane", new float[]{0.0, 70.0, -65.0}, new float[]{45.0, 140.0, 0.0, 1.0}, new float[]{216.0, 216.0, 216.0}, 0.5, 0.6, 0.3, 30.0, 0.0, 0.0, 1.0);
  temp.ifTexture = true;
  temp.texture = loadImage("wall_3.png");
  objects.add(temp);

  temp = new Primitive("plane", new float[]{-47.5, 70.0, -65.0}, new float[]{50.0, 70.0, 0.0, 1.0}, new float[]{216.0, 216.0, 216.0}, 0.5, 0.6, 0.3, 30.0, 0.0, 0.0, 1.0);
  temp.ifTexture = true;
  temp.texture = loadImage("wall_4.png");
  objects.add(temp);
  
  temp = new Primitive("plane", new float[]{-72.5, 80.0, -70.0}, new float[]{30.0, 160.0, 0.0, 1.0}, new float[]{216.0, 216.0, 216.0}, 0.5, 0.6, 0.3, 30.0, 0.0, 0.0, 1.0);
  temp.ifTexture = true;
  temp.texture = loadImage("wall_5.png");
  objects.add(temp);
  
  temp = new Primitive("plane", new float[]{37.5, 75.0, -65.0}, new float[]{30.0, 100.0, 0.0, 1.0}, new float[]{216.0, 216.0, 216.0}, 0.5, 0.6, 0.3, 30.0, 0.0, 0.0, 1.0);
  temp.ifTexture = true;
  temp.texture = loadImage("wall_6.png");
  objects.add(temp);
 
  temp = new Primitive("plane", new float[]{-40.0, 100.0, -70.0}, new float[]{35.0, 100.0, 0.0, 1.0}, new float[]{216.0, 216.0, 216.0}, 0.5, 0.6, 0.3, 30.0, 0.0, 0.0, 1.0);
  temp.ifTexture = true;
  temp.texture = loadImage("wall_7.png");
  objects.add(temp);

  temp = new Primitive("plane", new float[]{0.0, 30.0, 0.0}, new float[]{50.0, 60.0, 0.0, 1.0}, new float[]{216.0, 216.0, 216.0}, 0.5, 0.6, 0.3, 30.0, 0.0, 0.0, 1.0);
  temp.ifTexture = true;
  temp.texture = loadImage("door.png");
  objects.add(temp);

  temp = new Primitive("trapezoid", new float[]{0.0, 90.0, 0.0}, new float[]{60.0, 20.0, 50.0}, new float[]{216.0, 216.0, 216.0}, 0.5, 0.6, 0.3, 30.0, 0.0, 0.0, 1.00);
  temp.ifTexture = true;
  temp.texture = loadImage("door_top.png");
  objects.add(temp);

  temp = new Primitive("trapezoid", new float[]{0.0, 170.0, -65.0}, new float[]{60.0, 15.0, 30.0}, new float[]{216.0, 216.0, 216.0}, 0.5, 0.6, 0.3, 30.0, 0.0, 0.0, 1.00);
  temp.ifTexture = true;
  temp.texture = loadImage("wall_top.png");
  objects.add(temp);

  temp = new Primitive("circle", new float[]{0.0, 70.0, 1.0}, new float[]{20.0}, new float[]{230.0, 230.0, 0.0}, 0.5, 0.6, 0.3, 30.0, 0.0, 0.0, 1.00);
  temp.ifTexture = true;
  temp.texture = loadImage("50.png");
  objects.add(temp);
  objects.add(new Primitive("triangle", new float[]{-72.5, 167.7, -75.0}, new float[]{-72.5, 183.0, -85.0, -87.5, 160.0, -70.0, -57.5, 160.0, -70.0}, new float[]{0.0, 104.0, 185.0}, 0.5, 0.6, 0.3, 30.0, 0.0, 0.0, 1.00));


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
}

void draw() {
  if (toDraw) {
    float start = millis();
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
          if (p.type.equals("sphere")) {
            rgb[0] += temp[0];
            rgb[1] += temp[1];
            rgb[2] += temp[2];
          }
        }
        buffer.get(i).set(j, rgb);
        //if (rgb[0] / supersample > max_color) max_color = rgb[0] / supersample;
        //if (rgb[1] / supersample > max_color) max_color = rgb[1] / supersample;
        //if (rgb[2] / supersample > max_color) max_color = rgb[2] / supersample;
      }
    }
    max_color = 255;
    for (int i = 0; i < width; i++) {
      for (int j = 0; j < height; j++) {
        set(i, j, color(buffer.get(i).get(j)[0] / supersample / max_color * 255, buffer.get(i).get(j)[1] / supersample / max_color * 255, buffer.get(i).get(j)[2] / supersample / max_color * 255));
      }
    }
    //save("Final_Night_Glass_15.png");
    toDraw = false;
    println((millis() - start) / 1000.0);
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
  return normalize(new float[]{-v[0] + 2 * n[0] * dotTemp, -v[1] + 2 * n[1] * dotTemp, -v[2] + 2 * n[2] * dotTemp});
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
      coefficient = min((-B + sqrt(D)) / (2 * A), (-B - sqrt(D)) / (2 * A));
      if (coefficient < 0) coefficient = max((-B + sqrt(D)) / (2 * A), (-B - sqrt(D)) / (2 * A));
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
  } else if (p.type.equals("cone")) {
    float theta = atan(p.geometry[1] / p.geometry[0]);
    float[] V = new float[]{0.0, -1.0, 0.0};
    float[] CO = new float[]{origin[0] - p.center[0], origin[1] - p.center[1], origin[2] - p.center[2]};
    A = pow(dotProduct(ray, V), 2) - dotProduct(ray, ray) * pow(cos(theta), 2);
    B = 2 * (dotProduct(ray, V) * dotProduct(CO, V) - dotProduct(CO, ray) * pow(cos(theta), 2));
    C = pow(dotProduct(CO, V), 2) - dotProduct(CO, CO) * pow(cos(theta), 2);
    D = B * B - 4 * A * C;
    if (D >= 0) {
      coefficient = min((-B + sqrt(D)) / (2 * A), (-B - sqrt(D)) / (2 * A));
      if (coefficient < 0) coefficient = max((-B + sqrt(D)) / (2 * A), (-B - sqrt(D)) / (2 * A));
      if (coefficient > 0 && origin[1] + ray[1] * coefficient <= p.center[1] && origin[1] + ray[1] * coefficient >= p.center[1] - p.geometry[0]) {
        intersection[0] = origin[0] + ray[0] * coefficient;
        intersection[1] = origin[1] + ray[1] * coefficient;
        intersection[2] = origin[2] + ray[2] * coefficient;
      }
    }
  } else if (p.type.equals("cylinder")) {
    A = pow(ray[0], 2) + pow(ray[2], 2);
    B = 2 * (ray[0] * (origin[0] - p.center[0]) + ray[2] * (origin[2] - p.center[2]));
    C = pow(origin[0], 2) + pow(p.center[0], 2) - 2 * origin[0] * p.center[0] + pow(origin[2], 2) + pow(p.center[2], 2) - 2 * origin[2] * p.center[2] - pow(p.geometry[1], 2);
    D = B * B - 4 * A * C;
    if (D >= 0) {
      coefficient = min((-B + sqrt(D)) / (2 * A), (-B - sqrt(D)) / (2 * A));
      if (coefficient < 0) coefficient = max((-B + sqrt(D)) / (2 * A), (-B - sqrt(D)) / (2 * A));
      if (coefficient > 0 && origin[1] + ray[1] * coefficient <= p.center[1] + p.geometry[0] / 2 && origin[1] + ray[1] * coefficient >= p.center[1] - p.geometry[0] / 2) {
        intersection[0] = origin[0] + ray[0] * coefficient;
        intersection[1] = origin[1] + ray[1] * coefficient;
        intersection[2] = origin[2] + ray[2] * coefficient;
      }
    }
  } else if (p.type.equals("trapezoid")) {
    // center rectangle
    coefficient = (p.center[2] - origin[2]) / ray[2];
    if (coefficient > 0) {
      float[] temp = new float[3];
      temp[0] = origin[0] + ray[0] * coefficient;
      temp[1] = origin[1] + ray[1] * coefficient;
      temp[2] = origin[2] + ray[2] * coefficient;
      if (temp[0] <= p.geometry[1] / 2 + p.center[0] - (p.geometry[2] - p.geometry[1]) / 2 * (temp[1] - (p.center[1] + p.geometry[0] / 2)) / p.geometry[0] && temp[0] >= -p.geometry[1] / 2 + p.center[0] + (p.geometry[2] - p.geometry[1]) / 2 * (temp[1] - (p.center[1] + p.geometry[0] / 2)) / p.geometry[0] && temp[1] <= p.geometry[0] / 2  + p.center[1] && temp[1] >= -p.geometry[0] / 2 + p.center[1]) {
        intersection[0] = origin[0] + ray[0] * coefficient;
        intersection[1] = origin[1] + ray[1] * coefficient;
        intersection[2] = origin[2] + ray[2] * coefficient;
      }
    }
  } else if (p.type.equals("triangle")) {
    float[] vertex0, vertex1, vertex2;
    // vector from vertex1 to vertex0 and vector from vertex1 to vertex2
    float[] vector0 = new float[3];
    float[] vector1 = new float[3];
    // normal of the triangle
    float[] norm = new float[3];
    // auxiliary parameters for checking if the intersection is inside the triangle
    float uu, vv, a0, b0, c0, a1, b1, c1;
    float[] i_temp = new float[3];

    // left triangle
    vertex0 = new float[]{p.geometry[0], p.geometry[1], p.geometry[2]};
    vertex1 = new float[]{p.geometry[3], p.geometry[4], p.geometry[5]};
    vertex2 = new float[]{p.geometry[6], p.geometry[7], p.geometry[8]};
    vector0[0] = vertex1[0] - vertex0[0];
    vector0[1] = vertex1[1] - vertex0[1];
    vector0[2] = vertex1[2] - vertex0[2];
    vector1[0] = vertex2[0] - vertex0[0];
    vector1[1] = vertex2[1] - vertex0[1];
    vector1[2] = vertex2[2] - vertex0[2];
    vector0 = normalize(vector0);
    vector1 = normalize(vector1);
    norm = normalize(crossProduct(vector1, vector0));
    // formula of the triangle => Ax + By + Cz = D
    A = norm[0];
    B = norm[1];
    C = norm[2];
    D = A * vertex0[0] + B * vertex0[1] + C * vertex0[2];
    coefficient = (D - A * origin[0] - B * origin[1] - C * origin[2]) / (A * ray[0] + B * ray[1] + C * ray[2]);
    if (coefficient > 0) {
      i_temp[0] = origin[0] + ray[0] * coefficient;
      i_temp[1] = origin[1] + ray[1] * coefficient;
      i_temp[2] = origin[2] + ray[2] * coefficient;
      a0 = vertex0[0] - vertex2[0];
      b0 = vertex1[0] - vertex2[0];
      c0 = i_temp[0] - vertex2[0];
      a1 = vertex0[1] - vertex2[1];
      b1 = vertex1[1] - vertex2[1];
      c1 = i_temp[1] - vertex2[1];
      if (a0 == 0) {
        vv = c0 / b0;
        uu = (c1 - b1 * vv) / a1;
      } else if (b1 == 0) {
        uu = c1 / a1;
        vv = (c0 - a0 * uu) / b0;
      } else {
        vv = (c0/a0 - c1/a1) / (b0/a0 - b1/a1);
        uu = (c1 - b1* vv) / a1;
      }
      if (vv >= 0.0 && uu >= 0.0 && uu + vv <= 1.0) {
        intersection[0] = i_temp[0];
        intersection[1] = i_temp[1];
        intersection[2] = i_temp[2];
      }
    }
  } else if (p.type.equals("circle")) {
    // center rectangle
    coefficient = (p.center[2] - origin[2]) / ray[2];
    if (coefficient > 0) {
      float[] temp = new float[3];
      temp[0] = origin[0] + ray[0] * coefficient;
      temp[1] = origin[1] + ray[1] * coefficient;
      temp[2] = origin[2] + ray[2] * coefficient;
      if (calculate_distance(temp, p.center) <= p.geometry[0]) {
        intersection[0] = origin[0] + ray[0] * coefficient;
        intersection[1] = origin[1] + ray[1] * coefficient;
        intersection[2] = origin[2] + ray[2] * coefficient;
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
    ambient = new float[]{p.rgb[0] * p.ka * lights_i.get(0)[0], p.rgb[1] * p.ka * lights_i.get(0)[1], p.rgb[2] * p.ka * lights_i.get(0)[2]};
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
        if (dotTemp >= 0) diffuse = new float[]{p.rgb[0] * dotTemp * p.kd * lights_i.get(i)[0], p.rgb[1] * dotTemp * p.kd * lights_i.get(i)[1], p.rgb[2] * dotTemp * p.kd * lights_i.get(i)[2]};
        else diffuse = new float[]{0.0, 0.0, 0.0};
        dotTemp = dotProduct(R, V);
        if (dotTemp >= 0) specDot = pow(dotTemp, p.specExp);
        else specDot = 0.0;
        specular = new float[]{255.0 * specDot * p.ks * lights_i.get(i)[0], 255.0 * specDot * p.ks * lights_i.get(i)[1], 255.0 * specDot * p.ks * lights_i.get(i)[2]};
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
        final_color[0] += p.kt * refracted_color[0];
        final_color[1] += p.kt * refracted_color[1];
        final_color[2] += p.kt * refracted_color[2];
      }
    }
  } else if (p.type.equals("plane")) {
    if (p.ifTexture) {
      float u, v;
      if (p.geometry[0] == 0) {
        u = (intersection[2] - (p.center[2] - p.geometry[2] / 2)) / p.geometry[2];
        v = (intersection[1] - (p.center[1] - p.geometry[1] / 2)) / p.geometry[1];
      } else if (p.geometry[1] == 0) {
        u = (intersection[0] - (p.center[0] - p.geometry[0] / 2)) / p.geometry[0];
        v = (intersection[2] - (p.center[2] - p.geometry[2] / 2)) / p.geometry[2];
      } else {
        u = (intersection[0] - (p.center[0] - p.geometry[0] / 2)) / p.geometry[0];
        v = 1.0 - (intersection[1] - (p.center[1] - p.geometry[1] / 2)) / p.geometry[1];
      }
      int pixel_x = round(u * (p.texture.width - 1));
      int pixel_y = round(v * (p.texture.height - 1));
      float r = red(p.texture.pixels[pixel_x + pixel_y * p.texture.width]);
      float g = green(p.texture.pixels[pixel_x + pixel_y * p.texture.width]);
      float b = blue(p.texture.pixels[pixel_x + pixel_y * p.texture.width]);
      if (p.center[1] == 270.0 && p.center[2] == -135.0) ambient = new float[]{r, g, b};
      else if (p.kr != 0 && r == 0.0 && g == 0.0 && b == 0.0) ambient = new float[]{18.0 * p.ka * lights_i.get(0)[0], 63.0 * p.ka * lights_i.get(0)[1], 135.0 * p.ka * lights_i.get(0)[2]};
      else ambient = new float[]{r * p.ka * lights_i.get(0)[0], g * p.ka * lights_i.get(0)[1], b * p.ka * lights_i.get(0)[2]};
      final_color[0] += ambient[0];
      final_color[1] += ambient[1];
      final_color[2] += ambient[2];
      V = normalize(new float[]{origin[0] - intersection[0], origin[1] - intersection[1], origin[2] - intersection[2]});
      if (p.geometry[0] == 0.0) N = normalize(new float[]{p.geometry[3], 0.0, 0.0});
      else if (p.geometry[1] == 0.0) N = normalize(new float[]{0.0, p.geometry[3], 0.0});
      else N = normalize(new float[]{0.0, 0.0, p.geometry[3]});
      if (p.kr != 0 && r == 0.0 && g == 0.0 && b == 0.0) {
        float frequency = 54;
        float noise_right = noise((u + 0.01) * frequency, v * frequency);
        float noise_left = noise((u - 0.01) * frequency, v * frequency);
        float noise_top = noise(u * frequency, (v - 0.01) * frequency);
        float noise_bottom = noise(u * frequency, (v + 0.01) * frequency);
        float theta1 = acos(noise_right - noise_left);
        float theta2 = acos(noise_bottom - noise_top);
        float[] vTan = normalize(new float[]{sin(theta1), cos(theta1), 0.0});
        float[] vBi = normalize(new float[]{0.0, cos(theta2), sin(theta2)});
        N = normalize(crossProduct(vBi, vTan));
      }
      for (int i = 0; i < lights.size(); i++) {
        if ((!shadows.get(i) || softshadows.get(i)) && !(p.center[1] == 270.0 && p.center[2] == -135.0)) {
          S = normalize(new float[]{lights.get(i)[0] - intersection[0], lights.get(i)[1] - intersection[1], lights.get(i)[2] - intersection[2]});
          R = reflection(S, N);
          dotTemp = dotProduct(S, N);
          if (dotTemp >= 0) {
            if (p.kr != 0 && r == 0.0 && g == 0.0 && b == 0.0) diffuse = new float[]{18.0 * dotTemp * p.kd * lights_i.get(i)[0], 63.0 * dotTemp * p.kd * lights_i.get(i)[1], 135.0 * dotTemp * p.kd * lights_i.get(i)[2]};
            else diffuse = new float[]{r * dotTemp * p.kd * lights_i.get(i)[0], g * dotTemp * p.kd * lights_i.get(i)[1], b * dotTemp * p.kd * lights_i.get(i)[2]};
          } else diffuse = new float[]{0.0, 0.0, 0.0};
          dotTemp = dotProduct(R, V);
          if (dotTemp >= 0) specDot = pow(dotTemp, p.specExp);
          else specDot = 0.0;
          specular = new float[]{255.0 * specDot * p.ks * lights_i.get(i)[0], 255.0 * specDot * p.ks * lights_i.get(i)[1], 255.0 * specDot * p.ks * lights_i.get(i)[2]};
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
      if (p.kr != 0 && r == 0.0 && g == 0.0 && b == 0.0) {
        float[] reflected_color = reflect(intersection, reflection(V, new float[]{0.0, 1.0, 0.0}));
        final_color[0] += p.kr * reflected_color[0];
        final_color[1] += p.kr * reflected_color[1];
        final_color[2] += p.kr * reflected_color[2];
      }
    } else {
      ambient = new float[]{p.rgb[0] * p.ka * lights_i.get(0)[0], p.rgb[1] * p.ka * lights_i.get(0)[1], p.rgb[2] * p.ka * lights_i.get(0)[2]};
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
          if (dotTemp >= 0) diffuse = new float[]{p.rgb[0] * dotTemp * p.kd * lights_i.get(i)[0], p.rgb[1] * dotTemp * p.kd * lights_i.get(i)[1], p.rgb[2] * dotTemp * p.kd * lights_i.get(i)[2]};
          else diffuse = new float[]{0.0, 0.0, 0.0};
          dotTemp = dotProduct(R, V);
          if (dotTemp >= 0) specDot = pow(dotTemp, p.specExp);
          else specDot = 0.0;
          specular = new float[]{255.0 * specDot * p.ks * lights_i.get(i)[0], 255.0 * specDot * p.ks * lights_i.get(i)[1], 255.0 * specDot * p.ks * lights_i.get(i)[2]};
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
  } else if (p.type.equals("cone")) {
    float v = (intersection[1] - (p.center[1] - p.geometry[0])) / p.geometry[0];
    if (mod(v, 6 / p.geometry[0]) / (6 / p.geometry[0]) > 0.93) ambient = new float[]{60.0 * p.ka * lights_i.get(0)[0], 60.0 * p.ka * lights_i.get(0)[1], 60.0 * p.ka * lights_i.get(0)[2]};
    else ambient = new float[]{p.rgb[0] * p.ka * lights_i.get(0)[0], p.rgb[1] * p.ka * lights_i.get(0)[1], p.rgb[2] * p.ka * lights_i.get(0)[2]};
    final_color[0] += ambient[0];
    final_color[1] += ambient[1];
    final_color[2] += ambient[2];
    if (intersection[0] == p.center[0] && intersection[1] == p.center[1] && intersection[2] == p.center[2]) N = new float[]{0.0, 1.0, 0.0};
    else {
      float[] temp1 = normalize(new float[]{intersection[0] - p.center[0], intersection[1] - p.center[1], intersection[2] - p.center[2]});
      float[] temp2 = new float[]{0.0, -1.0, 0.0};
      temp2 = normalize(crossProduct(temp1, temp2));
      N = normalize(crossProduct(temp1, temp2));
    }
    V = normalize(new float[]{origin[0] - intersection[0], origin[1] - intersection[1], origin[2] - intersection[2]});
    for (int i = 0; i < lights.size(); i++) {
      if (!shadows.get(i) || softshadows.get(i)) {
        S = normalize(new float[]{lights.get(i)[0] - intersection[0], lights.get(i)[1] - intersection[1], lights.get(i)[2] - intersection[2]});
        R = reflection(S, N);
        dotTemp = dotProduct(S, N);
        if (dotTemp >= 0) {
          if (mod(v, 6 / p.geometry[0]) / (6 / p.geometry[0]) > 0.93) diffuse = new float[]{60.0 * dotTemp * p.kd, 60.0 * dotTemp * p.kd, 60.0 * dotTemp * p.kd};
          else diffuse = new float[]{p.rgb[0] * dotTemp * p.kd * lights_i.get(i)[0], p.rgb[1] * dotTemp * p.kd * lights_i.get(i)[1], p.rgb[2] * dotTemp * p.kd * lights_i.get(i)[2]};
        } else diffuse = new float[]{0.0, 0.0, 0.0};
        dotTemp = dotProduct(R, V);
        if (dotTemp >= 0) specDot = pow(dotTemp, p.specExp);
        else specDot = 0.0;
        specular = new float[]{255.0 * specDot * p.ks * lights_i.get(i)[0], 255.0 * specDot * p.ks * lights_i.get(i)[1], 255.0 * specDot * p.ks * lights_i.get(i)[2]};
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
  } else if (p.type.equals("cylinder")) {
    if (p.ifTexture) {
      float u;
      if (intersection[0] == p.center[0]) u = 0.0;
      else {
        float[] v1 = normalize(new float[]{intersection[0] - p.center[0], 0.0, intersection[2] - p.center[2]});
        float[] v2 = new float[]{0.0, 0.0, 1.0};
        float theta = acos(dotProduct(v1, v2));
        if (intersection[0] > p.center[0]) u = 1.0 * theta / 2 / PI;
        else u = 1.0 - 1.0 * theta / 2 / PI;
      }
      float v = 1.0 - (intersection[1] - (p.center[1] - p.geometry[0] / 2)) / p.geometry[0];
      int pixel_x = round(u * (p.texture.width - 1));
      int pixel_y = round(v * (p.texture.height - 1));
      p.rgb[0] = red(p.texture.pixels[pixel_x + pixel_y * p.texture.width]);
      p.rgb[1] = green(p.texture.pixels[pixel_x + pixel_y * p.texture.width]);
      p.rgb[2] = blue(p.texture.pixels[pixel_x + pixel_y * p.texture.width]);
    }
    ambient = new float[]{p.rgb[0] * p.ka * lights_i.get(0)[0], p.rgb[1] * p.ka * lights_i.get(0)[1], p.rgb[2] * p.ka * lights_i.get(0)[2]};
    final_color[0] += ambient[0];
    final_color[1] += ambient[1];
    final_color[2] += ambient[2];
    N = normalize(new float[]{intersection[0] - p.center[0], 0.0, intersection[2] - p.center[2]});
    V = normalize(new float[]{origin[0] - intersection[0], origin[1] - intersection[1], origin[2] - intersection[2]});
    for (int i = 0; i < lights.size(); i++) {
      if (!shadows.get(i) || softshadows.get(i)) {
        S = normalize(new float[]{lights.get(i)[0] - intersection[0], lights.get(i)[1] - intersection[1], lights.get(i)[2] - intersection[2]});
        R = reflection(S, N);
        dotTemp = dotProduct(S, N);
        if (dotTemp >= 0) diffuse = new float[]{p.rgb[0] * dotTemp * p.kd * lights_i.get(i)[0], p.rgb[1] * dotTemp * p.kd * lights_i.get(i)[1], p.rgb[2] * dotTemp * p.kd * lights_i.get(i)[2]};
        else diffuse = new float[]{0.0, 0.0, 0.0};
        dotTemp = dotProduct(R, V);
        if (dotTemp >= 0) specDot = pow(dotTemp, p.specExp);
        else specDot = 0.0;
        specular = new float[]{255.0 * specDot * p.ks * lights_i.get(i)[0], 255.0 * specDot * p.ks * lights_i.get(i)[1], 255.0 * specDot * p.ks * lights_i.get(i)[2]};
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
  } else if (p.type.equals("trapezoid")) {
    float r = 0.0, g = 0.0, b = 0.0;
    if (p.ifTexture) {
      float u = (intersection[0] - (p.center[0] - p.geometry[2] / 2)) / p.geometry[2];
      float v = 1.0 - (intersection[1] - (p.center[1] - p.geometry[0] / 2)) / p.geometry[0];
      int pixel_x = round(u * (p.texture.width - 1));
      int pixel_y = round(v * (p.texture.height - 1));
      r = red(p.texture.pixels[pixel_x + pixel_y * p.texture.width]);
      g = green(p.texture.pixels[pixel_x + pixel_y * p.texture.width]);
      b = blue(p.texture.pixels[pixel_x + pixel_y * p.texture.width]);
      ambient = new float[]{r * p.ka * lights_i.get(0)[0], g * p.ka * lights_i.get(0)[1], b * p.ka * lights_i.get(0)[2]};
    } else ambient = new float[]{p.rgb[0] * p.ka * lights_i.get(0)[0], p.rgb[1] * p.ka * lights_i.get(0)[1], p.rgb[2] * p.ka * lights_i.get(0)[2]};
    final_color[0] += ambient[0];
    final_color[1] += ambient[1];
    final_color[2] += ambient[2];
    N = new float[]{0.0, 0.0, 1.0};
    V = normalize(new float[]{origin[0] - intersection[0], origin[1] - intersection[1], origin[2] - intersection[2]});
    for (int i = 0; i < lights.size(); i++) {
      if (!shadows.get(i) || softshadows.get(i)) {
        S = normalize(new float[]{lights.get(i)[0] - intersection[0], lights.get(i)[1] - intersection[1], lights.get(i)[2] - intersection[2]});
        R = reflection(S, N);
        dotTemp = dotProduct(S, N);
        if (dotTemp >= 0) {
          if (p.ifTexture) diffuse = new float[]{r * dotTemp * p.kd * lights_i.get(i)[0], g * dotTemp * p.kd * lights_i.get(i)[1], b * dotTemp * p.kd * lights_i.get(i)[2]};
          else diffuse = new float[]{p.rgb[0] * dotTemp * p.kd * lights_i.get(i)[0], p.rgb[1] * dotTemp * p.kd * lights_i.get(i)[1], p.rgb[2] * dotTemp * p.kd * lights_i.get(i)[2]};
        } else diffuse = new float[]{0.0, 0.0, 0.0};
        dotTemp = dotProduct(R, V);
        if (dotTemp >= 0) specDot = pow(dotTemp, p.specExp);
        else specDot = 0.0;
        specular = new float[]{255.0 * specDot * p.ks * lights_i.get(i)[0], 255.0 * specDot * p.ks * lights_i.get(i)[1], 255.0 * specDot * p.ks * lights_i.get(i)[2]};
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
  } else if (p.type.equals("circle")) {
    float u = (intersection[0] - (p.center[0] - p.geometry[0])) / (p.geometry[0] * 2);
    float v = 1.0 - (intersection[1] - (p.center[1] - p.geometry[0])) / (p.geometry[0] * 2);
    int pixel_x = round(u * (p.texture.width - 1));
    int pixel_y = round(v * (p.texture.height - 1));
    float r = red(p.texture.pixels[pixel_x + pixel_y * p.texture.width]);
    float g = green(p.texture.pixels[pixel_x + pixel_y * p.texture.width]);
    float b = blue(p.texture.pixels[pixel_x + pixel_y * p.texture.width]);

    ambient = new float[]{p.rgb[0] * p.ka * lights_i.get(0)[0], p.rgb[1] * p.ka * lights_i.get(0)[1], p.rgb[2] * p.ka * lights_i.get(0)[2]};
    final_color[0] += ambient[0];
    final_color[1] += ambient[1];
    final_color[2] += ambient[2];
    N = normalize(new float[]{r / 255, g / 255, b / 255});
    V = normalize(new float[]{origin[0] - intersection[0], origin[1] - intersection[1], origin[2] - intersection[2]});
    for (int i = 0; i < lights.size(); i++) {
      if (!shadows.get(i) || softshadows.get(i)) {
        S = normalize(new float[]{lights.get(i)[0] - intersection[0], lights.get(i)[1] - intersection[1], lights.get(i)[2] - intersection[2]});
        R = reflection(S, N);
        dotTemp = dotProduct(S, N);
        if (dotTemp >= 0) diffuse = new float[]{p.rgb[0] * dotTemp * p.kd * lights_i.get(i)[0], p.rgb[1] * dotTemp * p.kd * lights_i.get(i)[1], p.rgb[2] * dotTemp * p.kd * lights_i.get(i)[2]};
        else diffuse = new float[]{0.0, 0.0, 0.0};
        dotTemp = dotProduct(R, V);
        if (dotTemp >= 0) specDot = pow(dotTemp, p.specExp);
        else specDot = 0.0;
        specular = new float[]{255.0 * specDot * p.ks * lights_i.get(i)[0], 255.0 * specDot * p.ks * lights_i.get(i)[1], 255.0 * specDot * p.ks * lights_i.get(i)[2]};
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
  } else if (p.type.equals("triangle")) {
    float v = (intersection[1] - (p.geometry[7])) / (p.geometry[1] - p.geometry[7]);
    if (mod(v, 6 / (p.geometry[1] - p.geometry[7])) / (6 / (p.geometry[1] - p.geometry[7])) > 0.93) ambient = new float[]{60.0 * p.ka * lights_i.get(0)[0], 60.0 * p.ka * lights_i.get(0)[1], 60.0 * p.ka * lights_i.get(0)[2]};
    else ambient = new float[]{p.rgb[0] * p.ka * lights_i.get(0)[0], p.rgb[1] * p.ka * lights_i.get(0)[1], p.rgb[2] * p.ka * lights_i.get(0)[2]};
    final_color[0] += ambient[0];
    final_color[1] += ambient[1];
    final_color[2] += ambient[2];
    N = normalize(crossProduct(new float[]{p.geometry[3] - p.geometry[0], p.geometry[4] - p.geometry[1], p.geometry[5] - p.geometry[2]}, new float[]{p.geometry[6] - p.geometry[0], p.geometry[7] - p.geometry[1], p.geometry[8] - p.geometry[2]}));
    V = normalize(new float[]{origin[0] - intersection[0], origin[1] - intersection[1], origin[2] - intersection[2]});
    for (int i = 0; i < lights.size(); i++) {
      if (!shadows.get(i) || softshadows.get(i)) {
        S = normalize(new float[]{lights.get(i)[0] - intersection[0], lights.get(i)[1] - intersection[1], lights.get(i)[2] - intersection[2]});
        R = reflection(S, N);
        dotTemp = dotProduct(S, N);
        if (dotTemp >= 0) {
          if (mod(v, 6 / (p.geometry[1] - p.geometry[7])) / (6 / (p.geometry[1] - p.geometry[7])) > 0.93) diffuse = new float[]{60.0 * dotTemp * p.kd * lights_i.get(i)[0], 60.0 * dotTemp * p.kd * lights_i.get(i)[1], 60.0 * dotTemp * p.kd * lights_i.get(i)[2]};
          else diffuse = new float[]{p.rgb[0] * dotTemp * p.kd * lights_i.get(i)[0], p.rgb[1] * dotTemp * p.kd * lights_i.get(i)[1], p.rgb[2] * dotTemp * p.kd * lights_i.get(i)[2]};
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
  }
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
