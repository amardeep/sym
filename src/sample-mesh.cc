void sample_mesh(TriMesh *mesh, float grid_div) {
  mesh->need_bbox();
  mesh->need_neighbors();
  mesh->need_adjacentfaces();

  // size in each dimension
  float xs = mesh->bbox.max[0] - mesh->bbox.min[0];
  float ys = mesh->bbox.max[1] - mesh->bbox.min[1];
  float zs = mesh->bbox.max[2] - mesh->bbox.min[2];

  // divide smallest edge into grid_div pieces
  float gs = min(xs, ys, zs) / grid_div;

  printf("xs %f, ys %f, zs %f, gs %f\n", xs, ys, zs, gs);
  printf("grid: %d x %d x %d\n", ceil(xs/gs), ceil(ys/gs), ceil(zs/gs));

  vector<long> grid[int(xs/gs)][int(ys/gs)][int(zs/gs)];

  long gx, gy, gz;
  point gc;  // grid center

  for (int i = 0; i < mesh->vertices.size(); i++) {
    // (x, y, z) -> grid coordinates
    gx = (mesh->vertices[i][0] - mesh->bbox.min[0]) / gs;
    gy = (mesh->vertices[i][1] - mesh->bbox.min[1]) / gs;
    gz = (mesh->vertices[i][2] - mesh->bbox.min[2]) / gs;

    gc[0] = gs * (gx + .5);
    gc[1] = gs * (gy + .5);
    gc[2] = gs * (gz + .5);

    if (grid[gx][gy][gz].size() == 0) {
      grid[gx][gy][gz].push_back(i); // Once for a min distance placeholder
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



}
int main(int argc, char** argv) {



  return 0;

}
