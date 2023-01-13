class Node {
  // indices of three points of the triangle
  int[] indices;
  HashSet<Integer> index;
  String name;
  Node leftChild;
  Node rightChild;
  public Node() {
    indices = new int[]{-1, -1, -1};
    index = new HashSet();
  }
}
