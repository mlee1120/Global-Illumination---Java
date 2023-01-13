systen: Processing 4
I only documented the code for Final_Day in detail because other versions are similar.

Instructions:
 eye         => camera position
 lookAt      => camera lookat
 supersample => number of samples for supersampling
 softshadow  => opacity of soft shadow
 range       => determine dynamic range ("low"/"mid"/"high")
 size()      => canvas size
 lights      => add multiple light sources (just provide point light's position)
 objects     => add primitives:
                first primitive must be the background: (only needs color)
                  objects.add(new Primitive("background", new float[]{0.0, 200.0, -1000.0}, null, new float[]{0.0, 0.0, 0.0}, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0));

                for adding spheres, here is an example:
                  objects.add(new Primitive("sphere", new float[]{-55.0, 100.0, 100.0}, new float[]{60.0}, new float[]{150.0, 150.0, 150.0}, 0.5, 0.6, 0.3, 30.0, 0.0, 0.6, 1.08));
                  center position; geometry is its radiu; color; ka; kd; ks; specExp; reflection coefficient; reraction coefficient; index of refraction

                for adding cylinders, here is an example:
                  Primitive temp = new Primitive("cylinder", new float[]{175.0, 30.0, 0.0}, new float[]{60.0, 15.0}, new float[]{216.0, 216.0, 216.0}, 0.5, 0.6, 0.3, 30.0, 0.0, 0.0, 1.00);
                  temp.ifTexture = true;
                  temp.texture = loadImage("cylinder_right_1.png");
                  objects.add(temp);
                  center is its center position; geometry is its height and top/base radius; color; ka; kd; ks; specExp; reflection coefficient; reraction coefficient; index of refraction

                for adding cones, here is an example:
                  objects.add(new Primitive("cone", new float[]{-130.0, 60.0, -50.0}, new float[]{10.0, 7.5}, new float[]{0.0, 104.0, 185.0}, 0.5, 0.6, 0.15, 30.0, 0.0, 0.0, 1.00));
                  center is its top position; geometry is its height and base radius; color; ka; kd; ks; specExp; reflection coefficient; reraction coefficient; index of refraction

                for adding planes, here is an example:
                  Primitive temp = new Primitive("plane", new float[]{-102.5, 20.0, 0.0}, new float[]{115.0, 40.0, 0.0, 1.0}, new float[]{216.0, 216.0, 216.0}, 0.5, 0.6, 0.3, 30.0, 0.0, 0.0, 1.0);
                  temp.ifTexture = true;
                  temp.texture = loadImage("left_wall.png");
                  center is its center position; geometry is its sizes in xyz; color; ka; kd; ks; specExp; reflection coefficient; reraction coefficient; index of refraction

                for adding trapezoids, here is an example:
                  Primitive temp = new Primitive("trapezoid", new float[]{0.0, 90.0, 0.0}, new float[]{60.0, 20.0, 50.0}, new float[]{216.0, 216.0, 216.0}, 0.5, 0.6, 0.3, 30.0, 0.0, 0.0, 1.00);
                  temp.ifTexture = true;
                  temp.texture = loadImage("door_top.png");
                  objects.add(temp);
                  center is its center position; geometry is its height, top length, and base length; color; ka; kd; ks; specExp; reflection coefficient; reraction coefficient; index of refraction

                for adding circles, here is an example:
                  temp = new Primitive("circle", new float[]{0.0, 70.0, 1.0}, new float[]{20.0}, new float[]{230.0, 230.0, 0.0}, 0.5, 0.6, 0.3, 30.0, 0.0, 0.0, 1.00);
                  temp.ifTexture = true;
                  // normal map !!!
                  temp.texture = loadImage("50.png");
                  objects.add(temp);
                  center is its center position; geometry is its radius; color; ka; kd; ks; specExp; reflection coefficient; reraction coefficient; index of refraction

                for adding triangles, here is an example:
                  objects.add(new Primitive("triangle", new float[]{-72.5, 167.7, -75.0}, new float[]{-72.5, 183.0, -85.0, -87.5, 160.0, -70.0, -57.5, 160.0, -70.0}, new float[]{0.0, 104.0, 185.0}, 0.5, 0.6, 0.3, 30.0, 0.0, 0.0, 1.00));
                  center is its center of gravity position; geometry is its three vertices' coordinates; color; ka; kd; ks; specExp; reflection coefficient; reraction coefficient; index of refraction