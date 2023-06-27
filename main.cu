#include <iostream>
#include "src/variable.cuh"
#include "src/grids.cuh"
#include <cmath>
#include <omp.h>

int main() {

    double oerr1 = 0.0;
    double oerr2 = 0.0;

    std::cout << omp_get_max_threads() << std::endl;

    for (int n = 0; n < 5; n++) {

        Idim3 NodeSize(32 * (int) pow(2, n), 1, 1);
        Idim3 GridSize(1, 1, 1);
        Idim3 id(0, 0, 0);

        grids grid(NodeSize, GridSize, id, 0, true, DOMAIN_INFO(-1.0, 1.0));
        grid.LinkNeighbor();
        grid.apply_bc_x(0);
        grid.ToDevice();

        variable phi(NodeSize, GridSize, id, 0, 1, true, 2, &grid);
        phi.f[0]->LinkNeighbor();
        phi.f[0]->assign_bc(0, BC::INFO(BC::CELL_CENTER_PERIODIC), BC::INFO(BC::CELL_CENTER_PERIODIC));

        for (int i = 0; i < phi.NodeSize.x; i++) {
            phi.f[0]->set(cos(M_PI * grid.get(i, 0, 0, 0)), i, 0, 0, 0);
        }
        phi.f[0]->apply_bc_x(0);
        phi.f[0]->ToDevice();

        phi.get_derivative_SEC(1);
        phi.get_derivative_SEC(2);

        double err1, err2;
        for (int i = 1; i < phi.NodeSize.x - 1; i++) {
            err1 += fabs(phi.df[phi.X]->get(i, 0, 0, 0) + M_PI * sin(M_PI * grid.get(i, 0, 0, 0)));
            err2 += fabs(phi.ddf[phi.XX]->get(i, 0, 0, 0) + M_PI * M_PI * cos(M_PI * grid.get(i, 0, 0, 0)));
        }
        err1 = err1 / phi.NodeSize.x;
        err2 = err2 / phi.NodeSize.x;

        if (n > 0) {
            std::cout << (log(oerr1) - log(err1)) / (log(2.0)) << "," << (log(oerr2) - log(err2)) / (log(2.0))
                      << std::endl;
        }

        oerr1 = err1;
        oerr2 = err2;

    }


}
