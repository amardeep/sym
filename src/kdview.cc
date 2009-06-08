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
#include "common.h"
#include <argp.h>

using namespace std;


// Globals
vector<TriMesh *> meshes;
vector<xform> xforms;
vector<bool> visible;
vector<string> xffilenames;

vector<int> pdash;
ANNkd_tree* kdtree;
float kdradius = .1;
int kdnn = 30; // no. of neighbours
int current_sample = 0;
int match_index = -1;  // show which match? -1 means all matches

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

bool samples_only = false;
bool use_kdtree = true;
char* samplesfile = 0;
char* plyfile = 0;

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

void draw_sample(TriMesh* m, int i) {
  vector<int> matches;

  // draw point at i
  float rad = m->feature_size();
  glPushMatrix();
  glColor4f (.5, .5, 1.0, .6);
  glTranslatef(m->vertices[i][0],
               m->vertices[i][1],
               m->vertices[i][2]);
  glutSolidSphere(rad*2, 20, 20);
  glPopMatrix();

  // get nearest neighbours
  nearest_neighbors(kdtree, m, i, kdradius, kdnn, &matches);
  int num = matches.size();
  if (match_index > num) {
    match_index = -1;
  } else if (match_index < -1) {
    match_index = matches.size() - 1;
  }
  for (int iter = 0; iter < matches.size(); ++iter) {
    int j = matches[iter];
    if (i == j) continue;
    if (match_index != -1 && iter != match_index) continue;
    vec v1 = m->vertices[i];
    vec v2 = m->vertices[j];
    vec v = v1 - v2;
    v = normalize(v);
    float e1 = v DOT (m->normals[i] + m->normals[j]);
    float e2 = v DOT (m->pdir1[i] + m->pdir1[j]);
    float e3 = v DOT (m->pdir2[i] + m->pdir2[j]);
    float e = sqrt(e1*e1 + e2*e2 + e3*e3);

    if (iter == match_index) cout << "error: " << e << endl;

    float rad = m->feature_size();
    glPushMatrix();
    glColor4f (1, .5, .5, .6);
    if (e < 1.5) {
      glColor4f (.5, 1, .5, .6);
      rad *= 2;
    }
    glTranslatef(m->vertices[j][0],
                 m->vertices[j][1],
                 m->vertices[j][2]);
    glutSolidSphere(rad, 20, 20);
    glPopMatrix();
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

        draw_sample(meshes[i], pdash[current_sample]);
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
      current_sample += 1;
      if (current_sample >= pdash.size()) current_sample = 0;
      printf("current sample is: %d/%d\n", current_sample, pdash.size());
      break;
    case 'U':
      current_sample -= 1;
      if (current_sample <= 0) current_sample = pdash.size() - 1;
      printf("current sample is: %d/%d\n", current_sample, pdash.size());
      break;
    case 'r':
      kdradius *= 2;
      printf("current radius is: %f\n", kdradius);
      break;
    case 'R':
      kdradius /= 2;
      printf("current radius is: %f\n", kdradius);
      break;
    case 'j':
      match_index += 1;
      break;
    case 'J':
      match_index -= 1;
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


int main(int argc, char **argv)
{
  glutInitWindowSize(512, 512);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
  glutInit(&argc, argv);

  if (argc < 2) {
    cerr << "Usage: " << argv[0] << " ply_file sample_file" << endl;
    return 1;
  }

  const char* filename = argv[1];

  TriMesh *mesh = TriMesh::read(filename);
  if (!mesh) {
    cerr << "Unable to parse " << filename << endl;
  }

  printf("before %d vertices\n", mesh->vertices.size());
  smooth_mesh(mesh, .005);
  printf("after %d vertices\n", mesh->vertices.size());

  mesh->need_curvatures();

  float smoothsigma = 5.0;
  smoothsigma *= mesh->feature_size();
  cout << "feature size " << mesh->feature_size() << endl;
  //diffuse_curv(mesh, smoothsigma);

  mesh->need_curvatures();
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
  //meshes.push_back(copy);

  string copyfilename = string(filename);
  copyfilename.replace(copyfilename.length() - 4, 4, "-ref.ply");
  xffilename = xfname(copyfilename.c_str());
  xffilenames.push_back(xffilename);
  printf("%s %s\n", filename, xffilename.c_str());
  xforms.push_back(xform());
  visible.push_back(true);


  // let's see how many points we can prune
  int numprune = 0;
  float gamma = .75;
  float mingc = numeric_limits<float>::infinity();
  float maxgc = 0;
  for (int i = 0; i < mesh->vertices.size(); ++i) {
    float k1 = mesh->curv1[i];
    float k2 = mesh->curv2[i];
    //k1 = k1 / sqrt(1 + k1 * k1);
    //k2 = k2 / sqrt(1 + k2 * k2);
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

  // pick up 1000 random points to match
  int nv = mesh->vertices.size();
  int index = 0;
  for (int i = 0; i < 100; ++i) {
    int index = (int)((double)rand() / ((double)RAND_MAX + 1) * nv);
    pdash.push_back(index);
  }

  // no, lets pick samples
  if (argc > 2) {
    ifstream fin(argv[2]);
    pdash.clear();
    while(fin >> index) {
      pdash.push_back(index);
    }
    cout << "sample size " << pdash.size() << endl;
  }
  // Uncomment one to just use samples to create kdtree or whole mesh
  //kdtree = create_kdtree(mesh, pdash);
  kdtree = create_kdtree(mesh);

  glutCreateWindow(argv[1]);
  glutDisplayFunc(redraw);
  glutMouseFunc(mousebuttonfunc);
  glutMotionFunc(mousemotionfunc);
  glutKeyboardFunc(keyboardfunc);
  glutIdleFunc(idle);

  resetview();

  glutMainLoop();
}
