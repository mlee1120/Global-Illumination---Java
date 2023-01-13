class Primitive{
  String type;
  float[] center;
  float[] geometry;
  float ka;
  float kd;
  float ks;
  float specExp;
  public Primitive(String type, float[] center, float[] geometry, float ka, float kd, float ks, float specExp) {
    this.type = type;
    this.center = center;
    this.geometry = geometry;
    this.ka = ka;
    this.kd = kd;
    this.ks = ks;
    this.specExp = specExp;
  }
}
