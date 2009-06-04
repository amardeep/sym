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

  float smoothsigma = 2.0;
  smoothsigma *= mesh->feature_size();
  cout << "feature size " << mesh->feature_size() << endl;
  diffuse_curv(mesh, smoothsigma);
  diffuse_curv(copy, smoothsigma);

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
  while (pdash.size() < 1000) {
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

  for (set<int>::iterator iter = pdash.begin(); iter != pdash.end(); ++iter) {
    int i = *iter;
    float ak1 = mesh->curv1[i];
    float ak2 = mesh->curv2[i];
    ak1 = ak1 / sqrt(1 + ak1 * ak1);
    ak2 = ak2 / sqrt(1 + ak2 * ak2);
    //for (int j = 0; j < nv; ++j) {
    for (int iter2 = 0; iter2 < samples.size(); ++iter2) {
      int j = samples[iter2];
      if (i == j) continue;
      float bk1 = copy->curv1[j];
      float bk2 = copy->curv2[j];
      bk1 = bk1 / sqrt(1 + bk1 * bk1);
      bk2 = bk2 / sqrt(1 + bk2 * bk2);
      //dist = abs(ak1-bk1) + abs(ak2-bk2);
      float dist = sqrt((ak1-bk1)*(ak1-bk1) + (ak2-bk2)*(ak2-bk2));
      if (dist < eps) {
        nmatches++;
        //printf("m: %f %f - %f %f (%d %d)\n", ak1, ak2, bk1, bk2, i, j);
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
