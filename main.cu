#include <iostream>
#include "src/variable.cuh"
#include "src/grids.cuh"
#include <cmath>
#include <omp.h>
#include <numbers>

int main() {

    double oerr1 = 0.0;
    double oerr2 = 0.0;

    std::cout << omp_get_max_threads() << std::endl;

    for (int n = 0; n < 10; n++) {

//        Idim3 NodeSize(4 * (int) pow(2, n), 1, 1);
        Idim3 NodeSize(256, 1, 1);
        Idim3 GridSize(1, 1, 1);
        Idim3 id(0, 0, 0);

        grids grid(NodeSize, GridSize, id, 0, true, DOMAIN_INFO(-1.0, 1.0), true);
        grid.LinkNeighbor();
        grid.apply_bc_x(0);
        grid.ToDevice();

        variable phi(NodeSize, GridSize, id, 0, 1, true, 2, &grid);
        phi.f[0]->LinkNeighbor();
        phi.f[0]->assign_bc(0, BC::INFO(BC::CELL_CENTER_PERIODIC), BC::INFO(BC::CELL_CENTER_PERIODIC));

        for (int i = 0; i < phi.NodeSize.x; i++) {
            phi.f[0]->set(cos(std::numbers::pi * grid.get(i, 0, 0, 0)), i, 0, 0, 0);
        }
        phi.f[0]->apply_bc_x(0);
        phi.f[0]->ToDevice();

//        phi.get_derivative_SEC(1, 0, DIM::X);
//        phi.get_derivative_SEC(2, 0, DIM::XX);

        phi.get_derivative_CCD(0, DIM::X);

        double err1, err2;
        err1 = 0.0;
        err2 = 0.0;
        for (int i = 1; i < phi.NodeSize.x - 1; i++) {
            err1 += std::pow(phi.df[0]->get(i, 0, 0, DIM::X) +
                             std::numbers::pi * sin(std::numbers::pi * grid.get(i, 0, 0, DIM::X)), 2.0);
            err2 += std::pow(phi.ddf[0]->get(i, 0, 0, DIM::XX) +
                             std::numbers::pi * std::numbers::pi * cos(std::numbers::pi * grid.get(i, 0, 0, DIM::X)),
                             2.0);
        }
        err1 = std::sqrt(err1 / phi.NodeSize.x);
        err2 = std::sqrt(err2 / phi.NodeSize.x);

        if (n > 0) {
            printf("%d, %.4f, %.5e, %.4f, %.5e\n", 4 * (int) pow(2, n), (log(oerr1) - log(err1)) / (log(2.0)), err1,
                   (log(oerr2) - log(err2)) / (log(2.0)), err2);
        }

        oerr1 = err1;
        oerr2 = err2;

    }


}
