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


// Get transform which transfrom m (mesh) into c (copy)
// i and j are the corresponding vertices in m and c to align.
xform get_xform(TriMesh* m, TriMesh* c, int i, int j) {
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



#endif
