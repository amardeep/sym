#include <stdio.h>
#include <stdlib.h>
#include "TriMesh.h"
#include "XForm.h"
#include "GLCamera.h"
#include "ICP.h"
#include <GL/glut.h>
#include <string>
#include "TriMesh_algo.h"
#include "KDtree.h"
#include "lineqn.h"
#include <algorithm>
#include <cmath>
#include <float.h>
#include "KDtree.h"
#include <iostream>
#include <set>

using namespace std;
using std::string;

#ifdef DARWIN
#define isfinite( x ) ( x <= FLT_MAX )
#endif

#define BIGNUM 3.3e33

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

// Quick 'n dirty portable random number generator
static inline float tinyrnd()
{
        static unsigned trand = 0;
        trand = 1664525u * trand + 1013904223u;
        return (float) trand / 4294967296.0f;
}


// Compute a "feature size" for the mesh: computed as 1% of
// the reciprocal of the 10-th percentile curvature
float feature_size(TriMesh *mesh)
{
        mesh->need_curvatures();
        int nv = mesh->curv1.size();
        int nsamp = min(nv, 500);

        vector<float> samples;
        samples.reserve(nsamp * 2);

        for (int i = 0; i < nsamp; i++) {
                // Quick 'n dirty portable random number generator
                static unsigned randq = 0;
                randq = unsigned(1664525) * randq + unsigned(1013904223);

                int ind = int(tinyrnd() * nv);
                samples.push_back(fabs(mesh->curv1[ind]));
                samples.push_back(fabs(mesh->curv2[ind]));
        }
        const float frac = 0.1f;
        const float mult = 0.01f;
        int which = int(frac * samples.size());
        nth_element(samples.begin(), samples.begin() + which, samples.end());
        float f = mult / samples[which];
        if (!isfinite(f)) {
                mesh->need_bsphere();
                f = mesh->bsphere.r;
        }

        fprintf(stderr, "Feature size = %f\n", f);
        return f;
}


// Color based on curvature
void colorbycurv(TriMesh *mesh, const char *scale, const char *smooth)
{
        mesh->need_curvatures();
        float smoothsigma = atof(smooth);
        if (smoothsigma > 0.0f) {
                smoothsigma *= mesh->feature_size();
                diffuse_curv(mesh, smoothsigma);
        }
        float cscale = 10.0f * atof(scale) * feature_size(mesh);
        TriMesh::dprintf("Using scale = %f\n", cscale);
        cscale = sqr(cscale);

        int nv = mesh->vertices.size();
#pragma omp parallel for
        for (int i = 0; i < nv; i++) {
                float H = 0.5f * (mesh->curv1[i] + mesh->curv2[i]);
                float K = mesh->curv1[i] * mesh->curv2[i];
                float h = 4.0f / 3.0f * fabs(atan2(H*H-K,H*H*sgn(H)));
                float s = M_2_PI * atan((2.0f*H*H-K)*cscale);
                mesh->colors[i] = Color::hsv(h,s,1.0);
        }

}


void usage(const char *myname)
{
        fprintf(stderr, "Usage: %s infile...\n", myname);
        exit(1);
}



//Takes in a mesh and an index to a vertex
bool isMinMax(TriMesh *mesh, long idx){
        float mcurv=(mesh->curv1[idx]+mesh->curv2[idx])/2;
        float tempCurv;
        //find the min and max of the neighboring points
        if(mesh-> neighbors[idx].size()>0){
                float minCurv=(mesh->curv1[mesh-> neighbors[idx][0]]+mesh->curv2[mesh-> neighbors[idx][0]])/2;
                float maxCurv=minCurv;
                for(int i=1; i<mesh-> neighbors[idx].size(); i++){
                        tempCurv=(mesh->curv1[mesh-> neighbors[idx][i]]+mesh->curv2[mesh-> neighbors[idx][i]])/2;
                        if(tempCurv<minCurv)
                                minCurv=tempCurv;
                        if(tempCurv>maxCurv)
                                maxCurv=tempCurv;
                }
                if((mcurv<minCurv)||(mcurv>maxCurv))
                        return true;
                else
                        return false;
        }
        return false;

}

int gridDiv=30;
vector<long> sampleMesh(TriMesh *mesh){
        mesh->need_bbox();
        mesh->need_neighbors();
        mesh->need_adjacentfaces();
        float xsize, ysize, zsize;
        xsize = mesh->bbox.max[0] -mesh->bbox.min[0];
        ysize = mesh->bbox.max[1] -mesh->bbox.min[1];
        zsize = mesh->bbox.max[2] -mesh->bbox.min[2];
        float gridSize=1;
        if((xsize<ysize)&&(xsize<zsize)){
                gridSize=xsize/gridDiv;
        }
        if((ysize<xsize)&&(ysize<zsize)){
                gridSize=ysize/gridDiv;
        }
        if((zsize<ysize)&&(zsize<xsize)){
                gridSize=zsize/gridDiv;
        }
        printf("Xsize: %f, Ysize: %f, Zsize: %f, gridSize: %f\n",xsize, ysize, zsize, gridSize);
        printf("Grid, %d by %d by %d\n", (int)(xsize/gridSize+1),(int)(ysize/gridSize+1),(int)(zsize/gridSize+1));
        vector<long> sampleGrid[(int)(xsize/gridSize+1)][(int)(ysize/gridSize+1)][(int)(zsize/gridSize+1)];

        long gridx, gridy, gridz;
        point gridCenter;
        KDtree *kd = new KDtree(mesh->vertices);
        /*for(int i=0; i<(int)(xsize/gridSize+1); i++){
                for(int j=0; j<(int)(ysize/gridSize+1); j++){
                        for(int k=0; k<(int)(zsize/gridSize+1); k++){
                                gridCenter[0]=gridSize*((float)i+.5);
                                gridCenter[1]=gridSize*((float)j+.5);
                                gridCenter[2]=gridSize*((float)k+.5);
                        }
                }
        }*/
        for (int i = 0; i < mesh->vertices.size(); i++) {
                gridx=(mesh->vertices[i][0]-mesh->bbox.min[0])/gridSize;
                gridy=(mesh->vertices[i][1]-mesh->bbox.min[1])/gridSize;
                gridz=(mesh->vertices[i][2]-mesh->bbox.min[2])/gridSize;

                gridCenter[0]=gridSize*((float)gridx+.5);
                gridCenter[1]=gridSize*((float)gridy+.5);
                gridCenter[2]=gridSize*((float)gridz+.5);

                if(sampleGrid[gridx][gridy][gridz].size()==0){
                        sampleGrid[gridx][gridy][gridz].push_back(i); //Once for a min distance placeholder
                        if(isMinMax(mesh, i)){
                                sampleGrid[gridx][gridy][gridz].push_back(i); //because it is a local min/max
                        }
                }

                //calculate distance to center;
                if(dist2(mesh->vertices[i], gridCenter)<dist2(mesh->vertices[sampleGrid[gridx][gridy][gridz].at(0)], gridCenter)){
                        sampleGrid[gridx][gridy][gridz].at(0)=i;
                }

                if(isMinMax(mesh, i)){
                        sampleGrid[gridx][gridy][gridz].push_back(i); //because it is a local min/max
                }
        }

        //Move sample point indices from grid to a single fector
        FILE *fid = fopen("sample.txt","w");
        vector<long> samplePoints;
        for(int i=0; i<(int)(xsize/gridSize+1); i++){
                for(int j=0; j<(int)(ysize/gridSize+1); j++){
                        for(int k=0; k<(int)(zsize/gridSize+1); k++){
                                for(int ind=0; ind<sampleGrid[i][j][k].size(); ind++){
                                        samplePoints.push_back(sampleGrid[i][j][k].at(ind));
                                        fprintf(fid, "%ld\n", sampleGrid[i][j][k].at(ind));
                                }
                        }
                }
        }
        fclose(fid);
        return samplePoints;
}

void expandMesh(TriMesh* original, TriMesh* patch, vector<bool> * pointMask){
        //Check to see if neighbors are calculated
        if(original->neighbors.size()<original->vertices.size())
                original->need_neighbors();

        //Get list of point indicies in patch that need to be checked for their immediate neighbors.
        vector<long> updateList;
        for(int i=0; i<pointMask->size(); i++){
                if(pointMask->at(i))
                        updateList.push_back(i);
        }


        //mesh->vertices[i] = xf * mesh->vertices[i];
        //Grow Patch
        for(int i=0; i<updateList.size(); i++){
                //loop through all neighbors of the current patch point
                for(int j=0; j<original->neighbors.at(updateList[i]).size(); j++){
                        if(!pointMask->at(original->neighbors[updateList[i]].at(j))){ //if is not in the point mask already
                                patch->vertices.push_back(original->vertices.at(original->neighbors[updateList[i]].at(j)));
                                pointMask->at(original->neighbors[updateList[i]].at(j))=true;
                        }
                }
        }
}

//Performs ICP on two patches
float icpPatch(TriMesh* patch1, TriMesh * patch2, xform &xf1, xform &xf2){
        int verbose = 0;
        bool do_scale = false;
        bool do_affine = false;
        return ICP(patch1, patch2, xf1, xf2, verbose, do_scale, do_affine);
}

void clearMesh(TriMesh* patch){
        patch->vertices.clear();
        patch->faces.clear();
        patch->tstrips.clear();
        patch->normals.clear();
        patch->pdir1.clear();
        patch->pdir2.clear();
        patch->curv1.clear();
        patch->curv2.clear();
        patch->neighbors.clear();
        patch->adjacentfaces.clear();
        patch->across_edge.clear();
}

void clearMask(vector<bool>* mask){
        for(int i=0; i<mask->size(); i++){
                mask->at(i)=false;
        }
}

//Merges mask1 and mask2 and stores the result in mask1
void mergeMask(vector<bool>* mask1, vector<bool>* mask2){
        for(int i=0; i<mask1->size(); i++){
                mask1->at(i)=mask1->at(i)||mask2->at(i);
        }
}

void rollBackMesh(TriMesh* mesh, TriMesh* patch, vector<bool>* mask){
        clearMesh(patch);
        for(int i=0; i<mask->size(); i++){
                if(mask->at(i))
                        mesh->vertices.push_back(mesh->vertices.at(i));
        }
}


struct patchPair{
        vector<bool> patch1;
        vector<bool> patch2;
        bool model1;
        bool model2;
        xform T;
        float error;
};

float icpThresh=.001;

//Takes in pairs of vertex indicies of the original mesh that are the correspondences from
//a single cluster and the clustered transform.
vector<patchPair> icpIter(TriMesh* themesh, TriMesh* themesh2, vector<int> &pair1, vector<int> &pair2, vector<bool> &model1, vector<bool> &model2, xform &xf, xform &final){
        int start1, start2;
        TriMesh *patch1 = new TriMesh();
        TriMesh *patch2 = new TriMesh();
        vector<patchPair> matchedPairs;

        //Global mask, set of all points visited by alignment process
        vector<bool> pointMask;
        pointMask.resize(themesh->vertices.size());

        //Local masks, reset if the patches need to be reintilized, used to make sure points are not added the the mask twice.
        vector<bool> pointMask1;
        pointMask1.resize(themesh->vertices.size());
        vector<bool> pointMask2;
        pointMask2.resize(themesh->vertices.size());

        clearMask(&pointMask);
        clearMask(&pointMask1);
        clearMask(&pointMask2);

        //Rollback masks, to undo the last patch expantion
        vector<bool> rollback1;
        vector<bool> rollback2;
        rollback1.resize(themesh->vertices.size());
        rollback2.resize(themesh->vertices.size());
        clearMask(&rollback1);
        clearMask(&rollback2);


        float error, errorPrev;
        xform xf2, xfNew;

        for(int i=0; i<pair1.size(); i++){
        //while(pair1.size()>0){ //Iterate until all point coorespondences are tested
                //Test if either point in the correspondence as already been used in any patch

                if(pointMask.at(pair1[i])||pointMask.at(pair2[i]))
                        continue;
                /*if(pointMask1.at(pair1[i])||pointMask2.at(pair1[i]))
                        continue;
                if(pointMask1.at(pair2[i])||pointMask2.at(pair2[i]))
                        continue;*/

                //clear patches
                clearMesh(patch1);
                clearMesh(patch2);

                //clear masks
                clearMask(&pointMask1);
                clearMask(&pointMask2);


                //initilize the two patches
                if(model1.at(i))
                        patch1->vertices.push_back(themesh2->vertices.at(pair1[i]));
                else
                        patch1->vertices.push_back(themesh->vertices.at(pair1[i]));

                if(model2.at(i))
                        patch2->vertices.push_back(themesh2->vertices.at(pair2[i]));
                else
                        patch2->vertices.push_back(themesh->vertices.at(pair2[i]));

                //init local vertex masks.
                pointMask1.at(pair1[i])=true;
                pointMask2.at(pair2[i])=true;

                rollback1.at(pair1[i])=true;
                rollback2.at(pair2[i])=true;

                //initial growth phase
                if(model1.at(i))
                        expandMesh(themesh2, patch1, &pointMask1);
                else
                        expandMesh(themesh, patch1, &pointMask1);
                if(model2.at(i))
                        expandMesh(themesh2, patch2, &pointMask2);
                else
                        expandMesh(themesh, patch2, &pointMask2);

                printf("Size patch1: %d, size patch2: %d\n", (int) patch1->vertices.size(), (int)patch2->vertices.size());

                error=icpPatch(patch1, patch2, xf, xf2);
                errorPrev=error;
                while(error<icpThresh){
                        //Update Masks

                        mergeMask(&rollback1,&pointMask1);
                        mergeMask(&rollback2,&pointMask2);

                        //expandMesh(themesh, patch1, &pointMask1);
                        //expandMesh(themesh2, patch2, &pointMask2);
                        if(model1.at(i))
                                expandMesh(themesh2, patch1, &pointMask1);
                        else
                                expandMesh(themesh, patch1, &pointMask1);
                        if(model2.at(i))
                                expandMesh(themesh2, patch2, &pointMask2);
                        else
                                expandMesh(themesh, patch2, &pointMask2);

                        errorPrev=error;
                        xfNew=xf2; //Iterate transformation
                        error=icpPatch(patch1, patch2, xfNew, xf2);
                }

                //rollBackMesh(themesh, patch1, &rollback1);
                //rollBackMesh(themesh2, patch2, &rollback1);
                patchPair temp;
                temp.patch1=vector<bool>(rollback1);
                temp.patch2=vector<bool>(rollback2);
                temp.error=errorPrev;
                temp.T=xfNew;
                temp.model1=model1.at(i);
                temp.model2=model2.at(i);
                matchedPairs.push_back(temp);
                mergeMask(&pointMask, &pointMask1);
                mergeMask(&pointMask, &pointMask2);

        }
        return matchedPairs;
}

void copyMesh(const TriMesh* src, TriMesh* dest) {
        // For now, we just copy the faces and vertices.
        dest->vertices = src->vertices;
        dest->faces = src->faces;
}


int main(int argc, char *argv[])
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
        meshes.push_back(copy);

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
        mesh->need_dcurv();
        copy->need_dcurv();
        /*for (int i = 0; i < mesh->vertices.size(); ++i) {
                if (i > 10) break;

                printf("Curvature %d: %f %f\n", i, mesh->curv1[i], mesh->curv2[i]);
        }

        printf("Curvature original, %f, %f.\n",  mesh->curv1[1], mesh->curv2[1]);
        printf("Direction 1 (%f, %f, %f), direction 2: (%f, %f, %f).\n", mesh->pdir1.at(1)[0], mesh->pdir1.at(1)[1], mesh->pdir1.at(1)[2], mesh->pdir2.at(1)[0],mesh->pdir2.at(1)[1],mesh->pdir2.at(1)[2]);
        printf("Curvature copy, %f, %f.\n",  copy->curv1[1], copy->curv2[1]);
        printf("Direction 1 (%f, %f, %f), direction 2: (%f, %f, %f).\n", copy->pdir1.at(1)[0], copy->pdir1.at(1)[1], copy->pdir1.at(1)[2], copy->pdir2.at(1)[0],copy->pdir2.at(1)[1],copy->pdir2.at(1)[2]);


        printf("nose left: %f, %f, nose right: %f, %f.\n", mesh->curv1[29404], mesh->curv2[29404], mesh->curv1[20738], mesh->curv2[20738]);
        */

         // let's see how many points we can prune
        int numprune = 0;
        float gamma = .75;
        float mingc = 1.0/0;
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

        // compute kd tree for curvature values

        // pick up 1000 random points to match
        /*set<int> pdash;
        int nv = mesh->vertices.size();
        while (pdash.size() < 1000) {
                int index = (int)((double)rand() / ((double)RAND_MAX + 1) * nv);
                pdash.insert(index);
        }*/

        int nv = mesh->vertices.size();
        vector<long> samplePoints=sampleMesh(mesh);
        printf("Sample points: %d.\n", samplePoints.size());
        // compute correspondences
        int nmatches = 0;
        float eps = .0001;
        corr.clear();
        //for (set<int>::iterator iter = pdash.begin(); iter != pdash.end(); ++iter) {
        //int i = *iter;

        FILE * fid;
        fid=fopen("cluster.txt","w");
        for(int c=0; c<samplePoints.size(); c++){
                int i = samplePoints.at(c);
                float ak1 = mesh->curv1[i];
                float ak2 = mesh->curv2[i];
                //ak1 = ak1 / sqrt(1 + ak1 * ak1);
                //ak2 = ak2 / sqrt(1 + ak2 * ak2);
                for (int j = 0; j < nv; ++j) {
                        float bk1 = copy->curv1[j];
                        float bk2 = copy->curv2[j];
                        bk1 = bk1 / sqrt(1 + bk1 * bk1);
                        bk2 = bk2 / sqrt(1 + bk2 * bk2);
                        if (abs(ak1-bk1) < eps && abs(ak2-bk2) < eps) {
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

                                //fprintf(fid, "%d %d ", i, j);
                                for(int m=0; m<16; m++){
                                        if((m+1)%4==0)
                                                continue;
                                        fprintf(fid, "%f ", comp[m]);
                                }
                                fprintf(fid, "\n");
                                //cout << t << endl;
                        }
                }
        }
        fclose(fid);
        printf("\nnmatches: %d/%d\n", nmatches, nv * 1000);


        /*if (argc < 2)
                usage(argv[0]);

        for (int i = 1; i < argc; i++) {
                const char *filename = argv[i];
                TriMesh *themesh = TriMesh::read(filename);
                if (!themesh)
                        usage(argv[0]);
                themesh->colors.resize(themesh->vertices.size());
                themesh->need_normals();
                themesh->need_tstrips();
                themesh->need_bsphere();
                colorbycurv(themesh, "1.0", "2.0");
                meshes.push_back(themesh);

                vector<long> samplePoints=sampleMesh(themesh);

                vector<int> pair1, pair2;
                vector<bool> model1, model2;
                model1.push_back(false);
                model2.push_back(true);
                pair1.push_back(100);
                pair2.push_back(1000);
                xform xf = xform();
                xform final = xform();
                vector<patchPair> test = icpIter(themesh,themesh, pair1, pair2, model1, model2, xf, final);

                printf(">>Sampled points: %d\n",(int)samplePoints.size());
                string xffilename = xfname(filename);
                xffilenames.push_back(xffilename);

                xforms.push_back(xform());
                visible.push_back(true);
        }*/

        glutCreateWindow(argv[1]);
        glutDisplayFunc(redraw);
        glutMouseFunc(mousebuttonfunc);
        glutMotionFunc(mousemotionfunc);
        glutKeyboardFunc(keyboardfunc);
        glutIdleFunc(idle);

        resetview();

        glutMainLoop();
}
