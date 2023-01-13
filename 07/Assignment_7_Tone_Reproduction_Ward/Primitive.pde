/**
 * This file illustates Primitive.pde.
 *
 * @author Michael Lee, ml3406@rit.edu
 */

/**
 * a class represents primitive (object)
 */
class Primitive {
  // the type if this Primitive (plane, sphere, background)
  String type;
  
  // the center of this Primitive
  float[] center;
  
  // the geometry information of this Primitive
  float[] geometry;
  
  // the color of this Primitive
  float[] rgb;
  
  // the reflection and refraction coefficients and indices of this Primitive
  float ka;
  float kd;
  float ks;
  float specExp;
  float kr;
  float kt;
  float n;
  
  /**
   * The constructor initializes all variables.
   */
  public Primitive(String type, float[] center, float[] geometry, float[] rgb, float ka, float kd, float ks, float specExp, float kr, float kt, float n) {
    this.type = type;
    this.center = center;
    this.geometry = geometry;
    this.rgb = rgb;
    this.ka = ka;
    this.kd = kd;
    this.ks = ks;
    this.specExp = specExp;
    this.kr = kr;
    this.kt = kt;
    this.n = n;
  }
}
