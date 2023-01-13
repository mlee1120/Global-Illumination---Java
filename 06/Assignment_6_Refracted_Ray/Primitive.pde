class Primitive {
  String type;
  float[] center;
  float[] geometry;
  float[] rgb;
  float ka;
  float kd;
  float ks;
  float specExp;
  float kr;
  float kt;
  float n;
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
