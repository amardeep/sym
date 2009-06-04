#ifndef SYM_COMMON_H
#define SYM_COMMON_H

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
#include "ANN.h"


void copyMesh(const TriMesh* src, TriMesh* dest) {
  // For now, we just copy the faces and vertices.
  dest->vertices = src->vertices;
  dest->faces = src->faces;
}

// Get transform which transfrom m (mesh) into c (copy)
// i and j are the corresponding vertices in m and c to align.
xform get_xform(TriMesh* m, TriMesh* c, int i, int j) {
  using namespace std;
  xform r1(m->pdir1[i][0], m->pdir1[i][1], m->pdir1[i][2], 0,
           m->pdir2[i][0], m->pdir2[i][1], m->pdir2[i][2], 0,
           m->normals[i][0], m->normals[i][1], m->normals[i][2], 0,
           0, 0, 0, 1);

  xform r2(c->pdir1[j][0], c->pdir1[j][1], c->pdir1[j][2], 0,
           c->pdir2[j][0], c->pdir2[j][1], c->pdir2[j][2], 0,
           c->normals[j][0], c->normals[j][1], c->normals[j][2], 0,
           0, 0, 0, 1);
  xform rot = r2 * inv(r1);

  vec t = c->vertices[j] - rot * m->vertices[i];

  return xform::trans(t[0], t[1], t[2]) * rot;
}


// Compute euler angles from a transformation matrix.
vec EulerAngles(const xform &r){
  float eps = 0.00000001;
  float beta = atan2(-r[2], sqrt(r[0]*r[0]+r[1]*r[1]));
  float alpha = atan2(r[1]/(cos(beta)+eps), r[0]/(cos(beta)+eps));
  float gamma = atan2(r[6]/(cos(beta)+eps), r[10]/(cos(beta)+eps));

  vec angles;
  angles[0] = alpha;
  angles[1] = beta;
  angles[2] = gamma;
  return angles;
}

// Compute xform matrix from euler angles.
xform EulerMatrix(vec angles){
  xform xrot, yrot, zrot;

  zrot[0]  = cos(angles[0]);
  zrot[1]  = sin(angles[0]);
  zrot[4]  = -sin(angles[0]);
  zrot[5]  = cos(angles[0]);
  zrot[10] = 1;

  yrot[0]  = cos(angles[1]);
  yrot[10] = cos(angles[1]);
  yrot[2]  = -sin(angles[1]);
  yrot[8]  = sin(angles[1]);
  yrot[5]  = 1;

  xrot[5]  = cos(angles[2]);
  xrot[6]  = sin(angles[2]);
  xrot[9]  = -sin(angles[2]);
  xrot[10] = cos(angles[2]);
  xrot[0]  = 1;

  return zrot * yrot * xrot;
}

// Read modes from file
void read_modes(const char* filename, float scale, vector<xform>* modes) {
  using namespace std;

  ifstream fin(filename);

  float ignore;
  vec angles;
  vec t;

  while(true) {
    fin >> ignore >> angles[0] >> angles[1] >> angles[2]
        >> t[0] >> t[1] >> t[2];
    if (!fin.good()) break;

    angles /= scale;
    xform mode = EulerMatrix(angles);
    mode = xform::trans(t[0], t[1], t[2]) * mode;
    modes->push_back(mode);
  }
}


// Creates a kd tree from curvature values.
ANNkd_tree* create_kdtree(TriMesh* mesh) {
  int nv = mesh->vertices.size();
  // compute kd tree for curvature values
  ANNpointArray points = annAllocPts(nv, 2);
  // fill point memory
  for (int row = 0; row < nv; ++row) {
    float k1 = mesh->curv1[row];
    float k2 = mesh->curv2[row];
    // NOTE(amardeep): In my experience, uncommenting below resulted
    // in worse experience
    //k1 = k1 / sqrt(1 + k1 * k1);
    //k2 = k2 / sqrt(1 + k2 * k2);
    points[row][0] = k1;
    points[row][1] = k2;
  }
  // create kd-tree
  return new ANNkd_tree(points, nv, 2);
}


// Creates a kd tree from curvature values.
// Restrict to vertices from samples.
ANNkd_tree* create_kdtree(TriMesh* mesh, const vector<int>& samples) {
  int nv = samples.size();
  // compute kd tree for curvature values
  ANNpointArray points = annAllocPts(nv, 2);
  // fill point memory
  for (int iter = 0; iter < nv; ++iter) {
    int row = samples[iter];
    float k1 = mesh->curv1[row];
    float k2 = mesh->curv2[row];
    //k1 = k1 / sqrt(1 + k1 * k1);
    //k2 = k2 / sqrt(1 + k2 * k2);
    points[iter][0] = k1;
    points[iter][1] = k2;
  }
  // create kd-tree
  return new ANNkd_tree(points, nv, 2);
}


void nearest_neighbors(ANNkd_tree* kdtree, TriMesh* mesh, int row,
                       float rad, int num_matches,
                       vector<int>* matches) {
  // copy row into annpoint
  ANNpoint point = annAllocPt(2);
  point[0] = mesh->curv1[row];
  point[1] = mesh->curv2[row];

  // search in kdtree
  ANNidxArray idx = new ANNidx[num_matches];
  ANNdistArray dists = new ANNdist[num_matches];
  kdtree->annkFRSearch(point, rad, num_matches, idx, dists);

  for (int i = 0; i < num_matches; ++i) {
    if (idx[i] == ANN_NULL_IDX) break;
    matches->push_back(idx[i]);
  }

  delete[] idx;
  delete[] dists;
}


#endif
