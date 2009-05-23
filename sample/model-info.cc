#include "TriMesh.h"
#include <iostream>
using namespace std;

int main(int argc, char *argv[])
{
  if (argc < 2) {
    cout << "Usage: " << argv[0] << " ply_file" << endl;
  }

  const char *filename = argv[1];

  // Read model from ply file
  TriMesh *m = TriMesh::read(filename);
  if (!m) {
    cerr << "Unable to read model from " << filename << endl;
    return 1;
  }

  cout << "Number of vertices: " << m->vertices.size() << endl;
  cout << "Vertex 0 is at " << m->vertices[0] << endl;

  // Convert triangle strips to faces, if necessary
  m->need_faces();
  cout << "Face 0 has vertices " << m->faces[0][0] << ", "
       << m->faces[0][1] << ", and " << m->faces[0][2] << endl;

  m->need_normals();
  cout << "Vertex 0 has normal " << m->normals[0] << endl;
}
