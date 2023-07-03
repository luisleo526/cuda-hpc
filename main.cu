#include <iostream>
#include "src/variable.cuh"
#include "src/grids.cuh"
#include <cmath>
#include <omp.h>
#include <numbers>
#include "src/flux.cuh"

int main() {

    double oerr1 = 0.0;
    double oerr2 = 0.0;

    std::cout << omp_get_max_threads() << std::endl;

    for (int n = 0; n < 10; n++) {

        Idim3 NodeSize(4 * (int) pow(2, n), 1, 1);
//        Idim3 NodeSize(256, 1, 1);
        Idim3 GridSize(1, 1, 1);
        Idim3 id(0, 0, 0);

        grids grid(NodeSize, GridSize, id, 0, CUDA::PINNED_MEM, DOMAIN_INFO(-1.0, 1.0), true);
        grid.LinkNeighbor();
        grid.apply_bc_x(0);
        grid.ToDevice();

        variable phi(NodeSize, GridSize, id, 0, 1, CUDA::PINNED_MEM, VAR::W_SECOND_DERI, &grid);
        phi.f[VAR::SCALAR]->LinkNeighbor();
        phi.f[VAR::SCALAR]->assign_bc(0, BC::INFO(BC::CELL_CENTER_PERIODIC), BC::INFO(BC::CELL_CENTER_PERIODIC));

        variable vel(NodeSize, GridSize, id, 0, 1, CUDA::PINNED_MEM, VAR::WO_DERI, &grid);
        vel.f[VAR::U]->LinkNeighbor();
        vel.f[VAR::U]->assign_bc(0, BC::INFO(BC::CELL_CENTER_PERIODIC), BC::INFO(BC::CELL_CENTER_PERIODIC));

        for (int i = 0; i < phi.NodeSize.x; i++) {
            phi.f[VAR::SCALAR]->set(cos(std::numbers::pi * grid.get(i, 0, 0)), i, 0, 0);
            vel.f[VAR::U]->set(1.0, i, 0, 0);
        }
        phi.f[VAR::SCALAR]->apply_bc_x(0);
        vel.f[VAR::U]->apply_bc_x(0);
        phi.f[VAR::SCALAR]->ToDevice();
        vel.f[VAR::U]->ToDevice();

//        phi.get_derivative_SEC(1, VAR::SCALAR, DIM::X);
//        phi.get_derivative_SEC(2, VAR::SCALAR, DIM::XX);

//        phi.get_derivative_CCD(VAR::SCALAR, DIM::X);
        phi.get_derivative_UCCD(VAR::SCALAR, DIM::X, &phi);

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
