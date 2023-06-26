#include <iostream>
#include "src/data.cuh"
#include "src/grids.cuh"
#include <cmath>

int main() {

    Idim3 NodeSize(128, 1, 1);
    Idim3 GridSize(1, 1, 1);
    Idim3 id(0, 0, 0);
    int dim = 1;

    data test(NodeSize, GridSize, id, 0, dim, true);
    test.LinkNeighbor();
    test.assign_bc(0, BC::INFO(BC::CELL_CENTER_PERIODIC), BC::INFO(BC::CELL_CENTER_PERIODIC));

    grids grid(NodeSize, GridSize, id, 0, true, DOMAIN_INFO(-1.0, 1.0));
    grid.LinkNeighbor();
    grid.apply_bc_x(0);

    std::cout << grid.dx << "," << grid.dy << "," << grid.dz;

    for (int i = 0; i < test.NodeSize.x; i++) {
        test.set(cos(M_PI * grid.get(i, 0, 0, 0)), i, 0, 0, 0);
    }
    test.apply_bc_x(0);

    FILE *fp = fopen("test.dat", "w");
    for (int i = -3; i < test.NodeSize.x + 3; i++) {
        fprintf(fp, "%lf %lf %lf %lf\n", grid.get(i, 0, 0, 0), grid.get(i, 0, 0, 1), grid.get(i, 0, 0, 2),
                test.get(i, 0, 0, 0));
    }
    fclose(fp);

}
