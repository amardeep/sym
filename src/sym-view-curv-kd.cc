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
#include <limits>
#include "ANN.h"

using namespace std;


// Globals
vector<TriMesh *> meshes;
vector<xform> xforms;
vector<bool> visible;
vector<string> xffilenames;

vector<pair<int, int> > corr;
int corrindex = 0;

TriMesh::BSphere global_bsph;
xform global_xf;
GLCamera camera;

int current_mesh = -1;

bool draw_edges = false;
bool draw_2side = false;
bool draw_shiny = true;
bool draw_lit = true;
bool draw_falsecolor = false;
bool draw_index = false;
bool white_bg = false;


// Make some mesh current
void set_current(int i)
{
	if (i >= 0 && i < meshes.size() && visible[i])
		current_mesh = i;
	else
		current_mesh = -1;
	camera.stopspin();
}


// Change visiblility of a mesh
void toggle_vis(int i)
{
	if (i >= 0 && i < meshes.size())
		visible[i] = !visible[i];
	if (current_mesh == i && !visible[i])
		set_current(-1);
}


// Signal a redraw
void need_redraw()
{
	glutPostRedisplay();
}


// Clear the screen
void cls()
{
	glDisable(GL_DITHER);
	glDisable(GL_BLEND);
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_NORMALIZE);
	glDisable(GL_LIGHTING);
	glDisable(GL_NORMALIZE);
	glDisable(GL_COLOR_MATERIAL);
	if (white_bg)
		glClearColor(1, 1, 1, 0);
	else
		glClearColor(0.08, 0.08, 0.08, 0);
	glClearDepth(1);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

}


// Set up lights and materials
void setup_lighting(int id)
{
	Color c(1.0f);
	if (draw_falsecolor)
		c = Color::hsv(-3.88 * id, 0.6 + 0.2 * sin(0.42 * id), 1.0);
	glColor3fv(c);

	if (!draw_lit || meshes[id]->normals.empty()) {
		glDisable(GL_LIGHTING);
		return;
	}

	GLfloat mat_specular[4] = { 0.18, 0.18, 0.18, 0.18 };
	if (!draw_shiny) {
		mat_specular[0] = mat_specular[1] =
		mat_specular[2] = mat_specular[3] = 0.0f;
	}
	GLfloat mat_shininess[] = { 64 };
	GLfloat global_ambient[] = { 0.02, 0.02, 0.05, 0.05 };
	GLfloat light0_ambient[] = { 0, 0, 0, 0 };
	GLfloat light0_diffuse[] = { 0.85, 0.85, 0.8, 0.85 };
	if (current_mesh >= 0 && id != current_mesh) {
		light0_diffuse[0] *= 0.5f;
		light0_diffuse[1] *= 0.5f;
		light0_diffuse[2] *= 0.5f;
	}
	GLfloat light1_diffuse[] = { -0.01, -0.01, -0.03, -0.03 };
	GLfloat light0_specular[] = { 0.85, 0.85, 0.85, 0.85 };
	glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);
	glLightfv(GL_LIGHT0, GL_AMBIENT, light0_ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light0_specular);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, light1_diffuse);
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, global_ambient);
	glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_FALSE);
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, draw_2side);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_NORMALIZE);
}


// Draw triangle strips.  They are stored as length followed by values.
void draw_tstrips(const TriMesh *themesh)
{
	const int *t = &themesh->tstrips[0];
	const int *end = t + themesh->tstrips.size();
	while (likely(t < end)) {
		int striplen = *t++;
		glDrawElements(GL_TRIANGLE_STRIP, striplen, GL_UNSIGNED_INT, t);
		t += striplen;
	}
}


// Draw the mesh
void draw_mesh(int i)
{
	const TriMesh *themesh = meshes[i];

	glPushMatrix();
	glMultMatrixd(xforms[i]);

	glDepthFunc(GL_LESS);
	glEnable(GL_DEPTH_TEST);

	if (draw_2side) {
		glDisable(GL_CULL_FACE);
	} else {
		glCullFace(GL_BACK);
		glEnable(GL_CULL_FACE);
	}

	// Vertices
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_FLOAT,
			sizeof(themesh->vertices[0]),
			&themesh->vertices[0][0]);

	// Normals
	if (!themesh->normals.empty() && !draw_index) {
		glEnableClientState(GL_NORMAL_ARRAY);
		glNormalPointer(GL_FLOAT,
				sizeof(themesh->normals[0]),
				&themesh->normals[0][0]);
	} else {
		glDisableClientState(GL_NORMAL_ARRAY);
	}

	// Colors
	if (!themesh->colors.empty() && !draw_falsecolor && !draw_index) {
		glEnableClientState(GL_COLOR_ARRAY);
		glColorPointer(3, GL_FLOAT,
			       sizeof(themesh->colors[0]),
			       &themesh->colors[0][0]);
	} else {
		glDisableClientState(GL_COLOR_ARRAY);
	}

	// Main drawing pass
	if (themesh->tstrips.empty()) {
		// No triangles - draw as points
		glPointSize(1);
		glDrawArrays(GL_POINTS, 0, themesh->vertices.size());
		glPopMatrix();
		return;
	}

	if (draw_edges) {
		glPolygonOffset(10.0f, 10.0f);
		glEnable(GL_POLYGON_OFFSET_FILL);
	}

	draw_tstrips(themesh);
	glDisable(GL_POLYGON_OFFSET_FILL);

	// Edge drawing pass
	if (draw_edges) {
		glPolygonMode(GL_FRONT, GL_LINE);
		glDisableClientState(GL_COLOR_ARRAY);
		glDisable(GL_COLOR_MATERIAL);
		GLfloat global_ambient[] = { 0.2, 0.2, 0.2, 1.0 };
		GLfloat light0_diffuse[] = { 0.8, 0.8, 0.8, 0.0 };
		GLfloat light1_diffuse[] = { -0.2, -0.2, -0.2, 0.0 };
		GLfloat light0_specular[] = { 0.0f, 0.0f, 0.0f, 0.0f };
		glLightModelfv(GL_LIGHT_MODEL_AMBIENT, global_ambient);
		glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
		glLightfv(GL_LIGHT1, GL_DIFFUSE, light1_diffuse);
		glLightfv(GL_LIGHT0, GL_SPECULAR, light0_specular);
		GLfloat mat_diffuse[4] = { 0.0f, 0.0f, 1.0f, 1.0f };
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_diffuse);
		glColor3f(0, 0, 1); // Used iff unlit
		draw_tstrips(themesh);
		glPolygonMode(GL_FRONT, GL_FILL);
	}

        const TriMesh* m1 = themesh;
        int j = corrindex;
        float k1 = m1->curv1[corr[j].first];
        float k2 = m1->curv2[corr[j].first];
        printf("curv: %f %f\n", k1, k2);
        if (abs(k2/k1) > .75) cout << "discard" << endl;

        while(true) {
          int a = corr[j].first;
          int b = corr[j].second;
          if (a != corr[corrindex].first) {
            corrindex = j;
            break;
          }

          vec v1 = m1->vertices[a];
          vec v2 = m1->vertices[b];
          vec v = v1 - v2;
          v = normalize(v);
          float e1 = v DOT (m1->normals[a] + m1->normals[b]);
          float e2 = v DOT (m1->pdir1[a] + m1->pdir1[b]);
          float e3 = v DOT (m1->pdir2[b] + m1->pdir2[b]);
          float e = sqrt(e1*e1 + e2*e2 + e3*e3);

          glPushMatrix();
          glColor4f (.5, .5, 1.0, .6);
          glTranslatef(m1->vertices[a][0],
                       m1->vertices[a][1],
                       m1->vertices[a][2]);
          glutSolidSphere(.004, 20, 20);
          //glutSolidSphere(.004, 20, 20);
          //glutSolidSphere(40, 20, 20);
          //cout << m1->vertices[a] << endl;
          glPopMatrix();

          glPushMatrix();
          glColor4f (1, .5, .5, .6);
          //cout << "error: " << e << endl;
          if (e < 1.5) glColor4f (.5, 1, .5, .6);
          glTranslatef(m1->vertices[b][0],
                       m1->vertices[b][1],
                       m1->vertices[b][2]);
          glutSolidSphere(.004, 20, 20);
          glPopMatrix();
          ++j;
          if (j >= corr.size()) j=0;
        }
        glColor3f (1, 1, 1.0);


	glPopMatrix();
}


// Draw the scene
void redraw()
{
  timestamp t = now();
  camera.setupGL(global_xf * global_bsph.center, global_bsph.r);
  glPushMatrix();
  glMultMatrixd(global_xf);
  cls();
  for (int i = 0; i < meshes.size(); i++) {
    if (!visible[i])
      continue;
    setup_lighting(i);
    draw_mesh(i);
  }
  glPopMatrix();
  glutSwapBuffers();
  printf("\r                        \r%.1f msec.", 1000.0f * (now() - t));
  fflush(stdout);
}


// Update global bounding sphere.
void update_bsph()
{
	point boxmin(1e38, 1e38, 1e38);
	point boxmax(-1e38, -1e38, -1e38);
	bool some_vis = false;
	for (int i = 0; i < meshes.size(); i++) {
		if (!visible[i])
			continue;
		some_vis = true;
		point c = xforms[i] * meshes[i]->bsphere.center;
		float r = meshes[i]->bsphere.r;
		for (int j = 0; j < 3; j++) {
			boxmin[j] = min(boxmin[j], c[j]-r);
			boxmax[j] = max(boxmax[j], c[j]+r);
		}
	}
	if (!some_vis)
		return;
	point &gc = global_bsph.center;
	float &gr = global_bsph.r;
	gc = 0.5f * (boxmin + boxmax);
	gr = 0.0f;
	for (int i = 0; i < meshes.size(); i++) {
		if (!visible[i])
			continue;
		point c = xforms[i] * meshes[i]->bsphere.center;
		float r = meshes[i]->bsphere.r;
		gr = max(gr, dist(c, gc) + r);
	}
}


// Set the view...
void resetview()
{
	camera.stopspin();
	for (int i = 0; i < meshes.size(); i++)
		if (!xforms[i].read(xffilenames[i]))
			xforms[i] = xform();
	update_bsph();
	global_xf = xform::trans(0, 0, -5.0f * global_bsph.r) *
		    xform::trans(-global_bsph.center);

	// Special case for 1 mesh
	if (meshes.size() == 1 && xforms[0].read(xffilenames[0])) {
		global_xf = xforms[0];
		xforms[0] = xform();
	}
}


// Save the current image to a PPM file.
// Uses the next available filename matching filenamepattern
void dump_image()
{
	// Find first non-used filename
	const char filenamepattern[] = "img%d.ppm";
	int imgnum = 0;
	FILE *f;
	while (1) {
		char filename[1024];
		sprintf(filename, filenamepattern, imgnum++);
		f = fopen(filename, "rb");
		if (!f) {
			f = fopen(filename, "wb");
			printf("\n\nSaving image %s... ", filename);
			fflush(stdout);
			break;
		}
		fclose(f);
	}

	// Read pixels
	GLint V[4];
	glGetIntegerv(GL_VIEWPORT, V);
	GLint width = V[2], height = V[3];
	char *buf = new char[width*height*3];
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	glReadPixels(V[0], V[1], width, height, GL_RGB, GL_UNSIGNED_BYTE, buf);

	// Flip top-to-bottom
	for (int i = 0; i < height/2; i++) {
		char *row1 = buf + 3 * width * i;
		char *row2 = buf + 3 * width * (height - 1 - i);
		for (int j = 0; j < 3 * width; j++)
			swap(row1[j], row2[j]);
	}

	// Write out file
	fprintf(f, "P6\n#\n%d %d\n255\n", width, height);
	fwrite(buf, width*height*3, 1, f);
	fclose(f);
	delete [] buf;

	printf("Done.\n\n");
}


// Save all transforms
void save_xforms()
{
	if (xforms.size() == 1) {
		printf("Writing %s\n", xffilenames[0].c_str());
		global_xf.write(xffilenames[0]);
		return;
	}
	for (int i = 0; i < xforms.size(); i++) {
		printf("Writing %s\n", xffilenames[i].c_str());
		xforms[i].write(xffilenames[i]);
	}
}


// ICP
void do_icp(int n)
{
	camera.stopspin();

	if (current_mesh < 0 || current_mesh >= meshes.size())
		return;
	if (n < 0 || n >= meshes.size())
		return;
	if (!visible[n] || !visible[current_mesh] || n == current_mesh)
		return;
	ICP(meshes[n], meshes[current_mesh], xforms[n], xforms[current_mesh], 2);
	update_bsph();
	need_redraw();
}


// Handle mouse button and motion events
static unsigned buttonstate = 0;

void doubleclick(int button, int x, int y)
{
	// Render and read back ID reference image
	camera.setupGL(global_xf * global_bsph.center, global_bsph.r);
	glDisable(GL_BLEND);
	glDisable(GL_LIGHTING);
	glClearColor(1,1,1,1);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_DEPTH_TEST);
	draw_index = true;
	glPushMatrix();
	glMultMatrixd(global_xf);
	for (int i = 0; i < meshes.size(); i++) {
		if (!visible[i])
			continue;
		glColor3ub((i >> 16) & 0xff,
			   (i >> 8)  & 0xff,
			    i        & 0xff);
		draw_mesh(i);
	}
	glPopMatrix();
	draw_index = false;
	GLint V[4];
	glGetIntegerv(GL_VIEWPORT, V);
	y = int(V[1] + V[3]) - 1 - y;
	unsigned char pix[3];
	glReadPixels(x, y, 1, 1, GL_RGB, GL_UNSIGNED_BYTE, pix);
	int n = (pix[0] << 16) + (pix[1] << 8) + pix[2];

	if (button == 0 || buttonstate == (1 << 0)) {
		// Double left click - select a mesh
		set_current(n);
	} else if (button == 2 || buttonstate == (1 << 2)) {
		// Double right click - ICP current to clicked-on
		do_icp(n);
	}
}

void mousemotionfunc(int x, int y)
{
	static const Mouse::button physical_to_logical_map[] = {
		Mouse::NONE, Mouse::ROTATE, Mouse::MOVEXY, Mouse::MOVEZ,
		Mouse::MOVEZ, Mouse::MOVEXY, Mouse::MOVEXY, Mouse::MOVEXY,
	};

	Mouse::button b = Mouse::NONE;
	if (buttonstate & (1 << 3))
		b = Mouse::WHEELUP;
	else if (buttonstate & (1 << 4))
		b = Mouse::WHEELDOWN;
	else if (buttonstate & (1 << 30))
		b = Mouse::LIGHT;
	else
		b = physical_to_logical_map[buttonstate & 7];

	if (current_mesh < 0) {
		camera.mouse(x, y, b,
			     global_xf * global_bsph.center, global_bsph.r,
			     global_xf);
	} else {
		xform tmp_xf = global_xf * xforms[current_mesh];
		camera.mouse(x, y, b,
			     tmp_xf * meshes[current_mesh]->bsphere.center,
			     meshes[current_mesh]->bsphere.r,
			     tmp_xf);
		xforms[current_mesh] = inv(global_xf) * tmp_xf;
		update_bsph();
	}
	if (b != Mouse::NONE)
		need_redraw();
}

void mousebuttonfunc(int button, int state, int x, int y)
{
	static timestamp last_click_time;
	static unsigned last_click_buttonstate = 0;
	static float doubleclick_threshold = 0.25f;

	if (glutGetModifiers() & GLUT_ACTIVE_CTRL)
		buttonstate |= (1 << 30);
	else
		buttonstate &= ~(1 << 30);

	if (state == GLUT_DOWN) {
		buttonstate |= (1 << button);
		if (buttonstate == last_click_buttonstate &&
		    now() - last_click_time < doubleclick_threshold) {
			doubleclick(button, x, y);
			last_click_buttonstate = 0;
		} else {
			last_click_time = now();
			last_click_buttonstate = buttonstate;
		}
	} else {
		buttonstate &= ~(1 << button);
	}

	mousemotionfunc(x, y);
}


// Idle callback
void idle()
{
  xform tmp_xf = global_xf;
  if (current_mesh >= 0)
    tmp_xf = global_xf * xforms[current_mesh];

  if (camera.autospin(tmp_xf))
    need_redraw();
  else
    usleep(10000);

  if (current_mesh >= 0) {
    xforms[current_mesh] = inv(global_xf) * tmp_xf;
    update_bsph();
  } else {
    global_xf = tmp_xf;
  }
}


// Keyboard
#define Ctrl (1-'a')
void keyboardfunc(unsigned char key, int x, int y)
{
  switch (key) {
    case ' ':
      if (current_mesh < 0)
        resetview();
      else
        set_current(-1);
      break;
    case '@': // Shift-2
      draw_2side = !draw_2side; break;
    case 'e':
      draw_edges = !draw_edges; break;
    case 'u':
      //corrindex += 10;
      if (corrindex >= corr.size()) corrindex = 0;
      break;
    case 'f':
      draw_falsecolor = !draw_falsecolor; break;
    case 'l':
      draw_lit = !draw_lit; break;
    case 's':
      draw_shiny = !draw_shiny; break;
    case 'w':
      white_bg = !white_bg; break;
    case 'I':
      dump_image(); break;
    case Ctrl+'x':
      save_xforms();
      break;
    case '\033': // Esc
    case Ctrl+'q':
    case 'Q':
    case 'q':
        exit(0);
    default:
      if (key >= '1' && key <= '9') {
        int m = key - '1';
        toggle_vis(m);
      }
  }
  need_redraw();
}


void copyMesh(const TriMesh* src, TriMesh* dest) {
  // For now, we just copy the faces and vertices.
  dest->vertices = src->vertices;
  dest->faces = src->faces;
}

ANNkd_tree* create_kdtree(TriMesh* mesh) {
  int nv = mesh->vertices.size();
  // compute kd tree for curvature values
  ANNpointArray points = annAllocPts(nv, 2);
  // fill point memory
  for (int row = 0; row < nv; ++row) {
    points[row][0] = mesh->curv1[row];
    points[row][1] = mesh->curv2[row];
  }
  // create kd-tree
  return new ANNkd_tree(points, nv, 2);
}


const int nn = 30; // no. of neighbours
void nearest_neighbors(ANNkd_tree* kdtree, TriMesh* mesh, int row, float rad,
                      vector<int>* matches) {
  // copy row into annpoint
  ANNpoint point = annAllocPt(2);
  point[0] = mesh->curv1[row];
  point[1] = mesh->curv2[row];

  // search in kdtree
  ANNidxArray idx = new ANNidx[nn];
  ANNdistArray dists = new ANNdist[nn];
  kdtree->annkFRSearch(point, rad, nn, idx, dists);

  for (int i = 0; i < nn; ++i) {
    if (idx[i] == ANN_NULL_IDX) break;
    matches->push_back(idx[i]);
    cout << idx[i] << endl;
  }

  delete[] idx;
  delete[] dists;
}


xform get_xform(TriMesh* m, TriMesh* c, int i, int j) {
  vec t = m->vertices[i] - c->vertices[j];
  xform r1(c->pdir1[j][0], c->pdir1[j][1], c->pdir1[j][2], 0,
           c->pdir2[j][0], c->pdir2[j][1], c->pdir2[j][2], 0,
           c->normals[j][0], c->normals[j][1], c->normals[j][2], 0,
           0, 0, 0, 1);
  r1 = inv(r1);
  xform r2(m->pdir1[i][0], m->pdir1[i][1], m->pdir1[i][2], 0,
           m->pdir2[i][0], m->pdir2[i][1], m->pdir2[i][2], 0,
           m->normals[i][0], m->normals[i][1], m->normals[i][2], 0,
           0, 0, 0, 1);
  xform comp = r1*r2*xform::trans(t[0], t[1], t[2]);
  return comp;
}


int main(int argc, char **argv)
{
  glutInitWindowSize(512, 512);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
  glutInit(&argc, argv);

  if (argc < 2) {
    cerr << "Usage: " << argv[0] << " ply_file" << endl;
    return 1;
  }

  const char* filename = argv[1];

  TriMesh *mesh = TriMesh::read(filename);
  if (!mesh) {
    cerr << "Unable to parse " << filename << endl;
  }

  point cm = point_center_of_mass(mesh->vertices);
  cout << "com: " << cm << endl;
  //trans(mesh, cm * -1.0f);

  //lmsmooth(mesh, 50);

  mesh->need_curvatures();
  float smoothsigma = 2.0;
  smoothsigma *= mesh->feature_size();
  cout << "feature size " << mesh->feature_size() << endl;
  diffuse_curv(mesh, smoothsigma);

  mesh->need_normals();
  mesh->need_tstrips();
  mesh->need_bsphere();
  meshes.push_back(mesh);


  string xffilename = xfname(filename);
  xffilenames.push_back(xffilename);
  printf("%s %s\n", filename, xffilename.c_str());
  xforms.push_back(xform());
  visible.push_back(true);


  // create a copy of the face
  TriMesh *copy = new TriMesh();
  copyMesh(mesh, copy);
  scale(copy, -1, 1, 1);
  faceflip(copy);
  copy->need_normals();
  copy->need_tstrips();
  copy->need_bsphere();
  //pca_rotate(copy);
  //meshes.push_back(copy);

  string copyfilename = string(filename);
  copyfilename.replace(copyfilename.length() - 4, 4, "-ref.ply");
  xffilename = xfname(copyfilename.c_str());
  xffilenames.push_back(xffilename);
  printf("%s %s\n", filename, xffilename.c_str());
  xforms.push_back(xform());
  visible.push_back(true);

  // print few curvature values
  mesh->need_curvatures();
  copy->need_curvatures();
  //mesh->need_dcurv();
  //copy->need_dcurv();
  diffuse_curv(copy, smoothsigma);
  for (int i = 0; i < mesh->vertices.size(); ++i) {
    if (i > 10) break;

    printf("Curvature %d: %f %f\n", mesh->curv1[i], mesh->curv2[i]);
  }

  // let's see how many points we can prune
  int numprune = 0;
  float gamma = .75;
  float mingc = numeric_limits<float>::infinity();
  float maxgc = 0;
  for (int i = 0; i < mesh->vertices.size(); ++i) {
    float k1 = mesh->curv1[i];
    float k2 = mesh->curv2[i];
    k1 = k1 / sqrt(1 + k1 * k1);
    k2 = k2 / sqrt(1 + k2 * k2);
    float gc = k1 * k2;
    if (abs(gc) < mingc) {
      mingc = abs(gc);
    }
    if (abs(gc) > maxgc) {
      maxgc = abs(gc);
    }
    float ratio = abs(k2/k1);
    if (ratio > gamma) {
      numprune++;
    }
  }
  printf("Num pruned: %d/%d\n", numprune, mesh->vertices.size());
  printf("gc %f %f\n", mingc, maxgc);

  int nv = mesh->vertices.size();

  
  // pick up 1000 random points to match
  set<int> pdash;

  int index = 0;
  while (pdash.size() < 100) {
    int index = (int)((double)rand() / ((double)RAND_MAX + 1) * nv);
    //index++;
    pdash.insert(index);
  }

  //no, lets pick samples
  ifstream fin("sample.txt");
  vector<int> samples;
  pdash.clear();
  while(fin >> index) {
    pdash.insert(index);
    samples.push_back(index);
  }

  //float eps = .0001;
  // create kd tree from curvature values
  // ANNkd_tree* tree = create_kdtree(mesh);
  // set<pair<int, int> > scorr;
  // ofstream fout("o.txt");
  // ofstream fcluster("m.txt");
  // TriMesh *m = mesh;
  // for (set<int>::iterator iter = pdash.begin(); iter != pdash.end(); ++iter) {
  //   int i = *iter;
  //   vector<int> matches;
  //   nearest_neighbors(tree, mesh, i, eps, &matches);
  //   for (int iter2 = 0; iter2 < matches.size(); ++iter2) {
  //     int j = matches[iter2];
  //     if (i < j) {
  //       scorr.insert(make_pair(i, j));
  //     } else {
  //       scorr.insert(make_pair(j, i));
  //     }
  //     corr.push_back(make_pair(i, j));

  //     vec v1 = m->vertices[i];
  //     vec v2 = m->vertices[j];
  //     vec vd = v2 - v1;
  //     float lvd = len(vd);
  //     vec v = normalize(vd);
  //     float e1 = v DOT (m->normals[i] + m->normals[j]);
  //     float e2 = v DOT (m->pdir1[i] + m->pdir1[j]);
  //     float e3 = v DOT (m->pdir2[i] + m->pdir2[j]);
  //     float e = sqrt(e1*e1 + e2*e2 + e3*e3);

  //     if (e > 1.5) continue;
  //     vec c = ((v1 DOT v) + lvd/2) * v;
  //     fout << c[0] << " " << c[1] << " " << c[2] << endl;

  //     xform comp = get_xform(mesh, copy, i , j);
  //     cout << "here" << endl;
  //     for (int m = 0; m < 16; ++m) {
  //       if (m%4 == 3) continue;
  //       float o = comp[m];
  //       //if (m < 12) o *= 100;
  //       fcluster << o << " ";
  //     }
  //     fcluster << endl;

  //   }
  // }
  // printf("corr: %d %d\n", corr.size(), scorr.size());


  // compute correspondences
  int nmatches = 0;

  // for (set<pair<int, int> >::iterator iter = scorr.begin();
  //      iter != scorr.end(); ++iter) {
  //   int i = iter->first;
  //   int j = iter->second;
  //   //printf("%d %d\n", i, j);
  // }

  float eps = .0005;
  int num = 0;
  for (set<int>::iterator iter = pdash.begin(); iter != pdash.end(); ++iter) {
    int i = *iter;
    num++;
    if (num > 200) break;
    float ak1 = mesh->curv1[i];
    float ak2 = mesh->curv2[i];
    ak1 = ak1 / sqrt(1 + ak1 * ak1);
    ak2 = ak2 / sqrt(1 + ak2 * ak2);
    if (abs(ak2/ak1) > .75) continue;
    //for (int j = 0; j < nv; ++j) {
    for (int iter2 = 0; iter2 < samples.size(); ++iter2) {
      int j = samples[iter2];
      float bk1 = mesh->curv1[j];
      float bk2 = mesh->curv2[j];
      bk1 = bk1 / sqrt(1 + bk1 * bk1);
      bk2 = bk2 / sqrt(1 + bk2 * bk2);
      float dist;
      //dist = abs(ak1-bk1) + abs(ak2-bk2);
      dist = sqrt((ak1-bk1)*(ak1-bk1) + (ak2-bk2)*(ak2-bk2));
      //dist = abs(ak1*bk1 - ak2*bk2);
      if (dist < eps) {
        nmatches++;
        //printf("m: %f %f - %f %f (%d %d)\n", ak1, ak2, bk1, bk2, i, j);
        corr.push_back(make_pair(i, j));
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
        //cout << t << endl;
      }
    }
  }
  printf("\nnmatches: %d/%d\n", nmatches, nv * 1000);



  glutCreateWindow(argv[1]);
  glutDisplayFunc(redraw);
  glutMouseFunc(mousebuttonfunc);
  glutMotionFunc(mousemotionfunc);
  glutKeyboardFunc(keyboardfunc);
  glutIdleFunc(idle);

  resetview();

  glutMainLoop();
}
