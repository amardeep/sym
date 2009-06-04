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
#include "common.h"

using namespace std;

string projname;


// Signature using euler's angles
vector<float> xform_sig(xform x, float scale) {
  vec angles = EulerAngles(x);
  vector<float> sig;

  sig.push_back(angles[0] * scale);
  sig.push_back(angles[1] * scale);
  sig.push_back(angles[2] * scale);
  sig.push_back(x[12]);
  sig.push_back(x[13]);
  sig.push_back(x[14]);
  return sig;
}

// Signature using complete matrix params
vector<float> xform_sig_old(xform x, float scale) {
  vector<float> sig;

  for (int m = 0; m < 16; ++m) {
    if (m % 4 == 3) continue;
    float o = x[m];
    if (m < 12) o *= scale;
    sig.push_back(0);
  }
  return sig;
}

int main(int argc, char **argv)
{
  if (argc < 2) {
    cerr << "Usage: " << argv[0] << " ply_file sample_file" << endl;
    return 1;
  }

  const char* filename = argv[1];

  TriMesh *mesh = TriMesh::read(filename);
  if (!mesh) {
    cerr << "Unable to parse " << filename << endl;
  }

  mesh->need_normals();

  char* pname = argv[1];
  char* cur = argv[1];
  while(*cur) {
    if (*cur == '/') {
      pname = cur + 1;
    }
    cur++;
  }
  projname = string(pname);
  projname.replace(projname.length() - 4, 4, "");
  printf("Project name %s\n", projname.c_str());

  // create a copy of the face
  TriMesh *copy = new TriMesh();
  copyMesh(mesh, copy);

  // reflect model around yz plane
  scale(copy, -1, 1, 1);
  faceflip(copy);
  copy->need_normals();

  // need curvature
  mesh->need_curvatures();
  copy->need_curvatures();

  float smoothsigma = 3.0;
  smoothsigma *= mesh->feature_size();
  cout << "feature size " << mesh->feature_size() << endl;
  //diffuse_curv(mesh, smoothsigma);
  //diffuse_curv(copy, smoothsigma);

  // read samples
  ifstream fin(argv[2]);
  vector<int> samples;
  int s;
  while (fin >> s) {
    samples.push_back(s);
  }
  printf("Read %d samples\n", samples.size());

  // take 1000 random points to match
  set<int> pdash;
  int index = 0;
  int ss = samples.size();
  while (pdash.size() < 3000) {
    int index = (int)((double)rand() / ((double)RAND_MAX + 1) * ss);
    pdash.insert(samples[index]);
  }
  printf("Selected 1000 random points\n");

  // compute correspondences
  int nv = mesh->vertices.size();
  int nmatches = 0;
  float eps = .00001;

  // output file for correspondences and transformations
  string corrfile = projname + "-corr.txt";
  string transfile = projname + "-trans.txt";
  ofstream fcorr(corrfile.c_str());
  ofstream ftrans(transfile.c_str());

  vector<vector<float> > vtrans;

  ANNkd_tree* kdtree = create_kdtree(mesh);
  float kdradius = 50;
  int kdnn = 30;

  TriMesh *m = mesh;
  for (set<int>::iterator iter = pdash.begin(); iter != pdash.end(); ++iter) {
    int i = *iter;

    // get nearest neighbours
    vector<int> matches;
    nearest_neighbors(kdtree, mesh, i, kdradius, kdnn, &matches);
    int num = matches.size();
    for (int iter2 = 0; iter2 < matches.size(); ++iter2) {
      int j = matches[iter2];
      vec v1 = m->vertices[i];
      vec v2 = m->vertices[j];
      vec v = v1 - v2;
      v = normalize(v);
      float e1 = v DOT (m->normals[i] + m->normals[j]);
      float e2 = v DOT (m->pdir1[i] + m->pdir1[j]);
      float e3 = v DOT (m->pdir2[i] + m->pdir2[j]);
      float e = sqrt(e1*e1 + e2*e2 + e3*e3);
      if (e < 1.5) {
        nmatches++;
        fcorr << i << " " << j << endl;
        xform comp = get_xform(mesh, copy, i, j);
        vector<float> sig = xform_sig(comp, 100);
        vtrans.push_back(sig);
      }
    }
  }

  // write transform signatures to file
  ftrans << vtrans.size() << " " << 6 << endl;
  for (int i = 0; i < vtrans.size(); ++i) {
    for (int m = 0; m < 6; ++m) {
      ftrans << vtrans[i][m] << " ";
    }
    ftrans << endl;
  }

  printf("\nnmatches: %d/%d\n", nmatches, nv * 1000);
}
