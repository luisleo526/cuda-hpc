#include <iostream>
#include "src/variable.cuh"
#include "src/grids.cuh"
#include "src/flux.cuh"
#include "src/runge_kutta_solver.cuh"
#include <cmath>
#include <omp.h>
#include <numbers>
#include <string>
#include <sstream>

int main() {

    double oerr = 0.0;

    std::cout << omp_get_max_threads() << std::endl;

    for (int n = 0; n < 5; n++) {

        int N = 16 * (int) pow(2, n);
        Idim3 NodeSize(N, 1, 1);
//        Idim3 NodeSize(176, 1, 1);
        Idim3 GridSize(1, 1, 1);
        Idim3 id(0, 0, 0);

        grids grid(NodeSize, GridSize, id, 0, CUDA::UNPINNED_MEM, DOMAIN_INFO(-1.0, 1.0), false);
        grid.LinkNeighbor();
        grid.apply_bc_x(0);
        grid.ToDevice();

        variable phi(NodeSize, GridSize, id, 0, 1, CUDA::UNPINNED_MEM, VAR::W_SECOND_DERI, &grid);
        phi.f[VAR::SCALAR]->LinkNeighbor();
        phi.f[VAR::SCALAR]->assign_bc(0, BC::INFO(BC::CELL_CENTER_PERIODIC), BC::INFO(BC::CELL_CENTER_PERIODIC));

        variable vel(NodeSize, GridSize, id, 0, 1, CUDA::UNPINNED_MEM, VAR::WO_DERI, &grid);
        vel.f[VAR::U]->LinkNeighbor();
        vel.f[VAR::U]->assign_bc(0, BC::INFO(BC::CELL_CENTER_PERIODIC), BC::INFO(BC::CELL_CENTER_PERIODIC));

        runge_kutta_solver tsolver(NodeSize, GridSize, id, 0, CUDA::UNPINNED_MEM);

        for (int i = 0; i < phi.NodeSize.x; i++) {
            phi.f[VAR::SCALAR]->set(sin(std::numbers::pi * grid.get(i, 0, 0)), i, 0, 0);
            vel.f[VAR::U]->set(-1.0, i, 0, 0);
        }
        phi.f[VAR::SCALAR]->apply_bc_x(0);
        vel.f[VAR::U]->apply_bc_x(0);

        double time = 0;
        double dt = 0.01 / 512.0;
        int iter = 0;
        int pltid = 0;

        FILE *fp;
        std::ostringstream oss;
        oss << "solution_" << pltid << ".txt";
        fp = fopen(oss.str().c_str(), "w");
        for (int i = 0; i < phi.NodeSize.x; i++) {
            fprintf(fp, "%.4f, %.4f\n", grid.get(i, 0, 0), phi.f[VAR::SCALAR]->get(i, 0, 0));
        }
        fclose(fp);
        oss.str("");
        oss.clear();
        pltid++;

        while (time < 2.0) {
            tsolver.solve_tvdrk3(dt, &phi, VAR::SCALAR, &vel, &total_derivative);
            phi.f[VAR::SCALAR]->apply_bc_x(0);
            time += dt;
            if (++iter * dt > pltid * 0.5) {
                oss << "solution_" << pltid << ".txt";
                fp = fopen(oss.str().c_str(), "w");
                for (int i = 0; i < phi.NodeSize.x; i++) {
                    fprintf(fp, "%.4f, %.4f\n", grid.get(i, 0, 0), phi.f[VAR::SCALAR]->get(i, 0, 0));
                }
                fclose(fp);
                oss.str("");
                oss.clear();
                pltid++;
            }
        }

        oss << "solution_" << pltid << ".txt";
        fp = fopen(oss.str().c_str(), "w");
        for (int i = 0; i < phi.NodeSize.x; i++) {
            fprintf(fp, "%.4f, %.4f\n", grid.get(i, 0, 0), phi.f[VAR::SCALAR]->get(i, 0, 0));
        }
        fclose(fp);
        oss.str("");
        oss.clear();

        double err = 0.0;
        for (int i = 0; i < phi.NodeSize.x; i++) {
            double exact = sin(std::numbers::pi * (grid.get(i, 0, 0) - time + 2.0));
            err += std::pow(phi.f[VAR::SCALAR]->get(i, 0, 0) - exact, 2);

        }
        err = std::sqrt(err / phi.NodeSize.x);

        if (n > 0) {
            printf("%d, %.5e, %.4f\n", N, err, (log(oerr) - log(err)) / log(2.0));
        } else {
            printf("%d, %.5e\n", N, err);
        }

        oerr = err;

    }


}
