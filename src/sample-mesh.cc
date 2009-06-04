#include <stdio.h>
#include <stdlib.h>
#include "TriMesh.h"
#include "TriMesh_algo.h"
#include "XForm.h"
#include "GLCamera.h"
#include "ICP.h"
#include <GL/glut.h>
#include <string>
#include <iostream>
#include <set>
#include <map>
#include <limits>
#include "ANN.h"

using namespace std;

int choose_sample(TriMesh* mesh, int o, int n) {
  // we prefer points with larger curvatures.

  float ok1 = mesh->curv1[o];
  float ok2 = mesh->curv2[o];
  float nk1 = mesh->curv1[n];
  float nk2 = mesh->curv2[n];

  if (abs(nk1*nk2) > abs(ok1*ok2)) {
    return n;
  } else {
    return o;
  }
}

void sample_mesh(TriMesh *mesh, float grid_div,
                 map<vector<int>, int>& samples) {
  mesh->need_bbox();
  mesh->need_neighbors();
  mesh->need_adjacentfaces();

  // size in each dimension
  float xs = mesh->bbox.max[0] - mesh->bbox.min[0];
  float ys = mesh->bbox.max[1] - mesh->bbox.min[1];
  float zs = mesh->bbox.max[2] - mesh->bbox.min[2];

  // divide smallest edge into grid_div pieces
  float gs = min(xs, min(ys, zs)) / grid_div;

  printf("xs %f, ys %f, zs %f, gs %f\n", xs, ys, zs, gs);
  printf("grid: %d x %d x %d\n", ceil(xs/gs), ceil(ys/gs), ceil(zs/gs));

  long gx, gy, gz;
  point gc;  // grid center
  int numpruned = 0;
  int numfound = 0;

  for (int i = 0; i < mesh->vertices.size(); i++) {
    // (x, y, z) -> grid coordinates
    gx = (mesh->vertices[i][0] - mesh->bbox.min[0]) / gs;
    gy = (mesh->vertices[i][1] - mesh->bbox.min[1]) / gs;
    gz = (mesh->vertices[i][2] - mesh->bbox.min[2]) / gs;

    gx = int(gx);
    gy = int(gy);
    gz = int(gz);

    float k1 = mesh->curv1[i];
    float k2 = mesh->curv2[i];
    //k1 = k1 / sqrt(1 + k1 * k1);
    //k2 = k2 / sqrt(1 + k2 * k2);
    if (abs(k2/k1) > .75) {
      numpruned++;
      continue;
    }

    vector<int> g;
    g.push_back(gx);
    g.push_back(gy);
    g.push_back(gz);
    if (samples.find(g) == samples.end()) {
      samples[g] = i;
    } else {
      numfound++;
      samples[g] = choose_sample(mesh, samples[g], i);
    }
  }

  printf("numpruned: %d, numfound: %d, size: %d\n", numpruned, numfound,
         samples.size());

}

int main(int argc, char** argv) {
  if (argc < 2) {
    cerr << "Usage: " << argv[0] << " ply_file" << endl;
    return 1;
  }

  const char* filename = argv[1];
  TriMesh *mesh = TriMesh::read(filename);
  if (!mesh) {
    cerr << "Unable to parse " << filename << endl;
  }
  mesh->need_curvatures();
  float smoothsigma = 2.0;
  smoothsigma *= mesh->feature_size();
  cout << "feature size " << mesh->feature_size() << endl;
  //diffuse_curv(mesh, smoothsigma);

  mesh->need_normals();
  mesh->need_tstrips();
  mesh->need_bsphere();

  char* pname = argv[1];
  char* cur = argv[1];
  while(*cur) {
    if (*cur == '/') {
      pname = cur + 1;
    }
    cur++;
  }
  string projname = string(pname);
  projname.replace(projname.length() - 4, 4, "");
  printf("Project name %s\n", projname.c_str());

  map<vector<int>, int> samples;

  sample_mesh(mesh, 72, samples);

  string ofname = projname + "-sample.txt";
  ofstream fout(ofname.c_str());

  for (map<vector<int>, int>::iterator i = samples.begin();
       i != samples.end(); ++i) {
    fout << i->second << endl;
  }


  return 0;

}
