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

using namespace std;

string projname;

void copyMesh(const TriMesh* src, TriMesh* dest) {
  // For now, we just copy the faces and vertices.
  dest->vertices = src->vertices;
  dest->faces = src->faces;
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
  float eps = .0005;

  // output file for correspondences and transformations
  string corrfile = projname + "-corr.txt";
  string transfile = projname + "-trans.txt";
  ofstream fcorr(corrfile.c_str());
  ofstream ftrans(transfile.c_str());

  for (set<int>::iterator iter = pdash.begin(); iter != pdash.end(); ++iter) {
    int i = *iter;
    float ak1 = mesh->curv1[i];
    float ak2 = mesh->curv2[i];
    //ak1 = ak1 / sqrt(1 + ak1 * ak1);
    //ak2 = ak2 / sqrt(1 + ak2 * ak2);
    for (int j = 0; j < nv; ++j) {
      if (i == j) continue;
      float bk1 = copy->curv1[j];
      float bk2 = copy->curv2[j];
      //bk1 = bk1 / sqrt(1 + bk1 * bk1);
      //bk2 = bk2 / sqrt(1 + bk2 * bk2);
      if (abs(ak1-bk1) < eps && abs(ak2-bk2) < eps) {
        nmatches++;
        //printf("m: %f %f - %f %f (%d %d)\n", ak1, ak2, bk1, bk2, i, j);
        fcorr << i << " " << j << endl;
        vec t = mesh->vertices[i] - copy->vertices[j];
        xform r1(copy->pdir1[j][0], copy->pdir1[j][1], copy->pdir1[j][2], 0,
                 copy->pdir2[j][0], copy->pdir2[j][1], copy->pdir2[j][2], 0,
                 copy->normals[j][0], copy->normals[j][1], copy->normals[j][2], 0,
                 0, 0, 0, 1);
        r1 = inv(r1);
        xform r2(mesh->pdir1[i][0], mesh->pdir1[i][1], mesh->pdir1[i][2], 0,
                 mesh->pdir2[i][0], mesh->pdir2[i][1], mesh->pdir2[i][2], 0,
                 mesh->normals[i][0], mesh->normals[i][1], mesh->normals[i][2], 0,
                 0, 0, 0, 1);
        xform comp = r1*r2*xform::trans(t[0], t[1], t[2]);
        for (int m = 0; m < 16; ++m) {
          if (m%4 == 3) continue;
          float o = comp[m];
          if (m < 12) o *= 100;
          ftrans << o << " ";
        }
        ftrans << endl;
      }
    }
  }
  printf("\nnmatches: %d/%d\n", nmatches, nv * 1000);
  printf("%d %d\n", nmatches, 12);
}
