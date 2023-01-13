/**
 * This file illustates Quaternion.pde.
 *
 * @author Michael Lee, ml3406@rit.edu
 */

/**
 * a class represents quaternion (for rotation)
 */
class Quaternion {
  /** scalar */
  public float w;

  /** x vector */
  public float x;

  /** y vector */
  public float y;

  /** z vector */
  public float z;

  /**
   * The constructor initializes all variables.
   */
  public Quaternion(float w, float x, float y, float z) {
    this.w = w;
    this.x = x;
    this.y = y;
    this.z = z;
  }

  /**
   * This function normalizes this quaternion.
   */
  public void normalize() {
    float rootSumSquare = sqrt(w * w + x * x + y * y + z * z);
    float difference = 0.0000000001;
    float tempW, tempX, tempY, tempZ;
    do {
      tempW = w / rootSumSquare;
      tempX = x / rootSumSquare;
      tempY = y / rootSumSquare;
      tempZ = z / rootSumSquare;
      rootSumSquare += difference++;
    } while (tempW * tempW + tempX * tempX + tempY * tempY + tempZ * tempZ > 1.0);
    w = tempW;
    x = tempX;
    y = tempY;
    z = tempZ;
  }

  /**
   * This function returns the inverse quaternion of this quaternion.
   */
  public Quaternion inverse() {
    return new Quaternion(w, -x, -y, -z);
  }

  /**
   * This function calculates and returns the production of this quaternion and another one.
   *
   * @param q - the other quaternion
   */
  public Quaternion product(Quaternion q) {
    float tempW, tempX, tempY, tempZ;
    tempW = w * q.w - x * q.x - y * q.y - z * q.z;
    tempX = w * q.x + x * q.w + y * q.z - z * q.y;
    tempY = w * q.y - x * q.z + y * q.w + z * q.x;
    tempZ = w * q.z + x * q.y - y * q.x + z * q.w;
    return new Quaternion(tempW, tempX, tempY, tempZ);
  }
}
