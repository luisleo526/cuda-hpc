//
// Created by leo on 6/29/23.
//

#include "flux.cuh"

void total_derivative(variable *phi, variable *vel, int d) {

    assert(phi->derivative_level > 0);
    assert(phi->spatial_dim == vel->spatial_dim);

    phi->get_derivative_UCCD(d, DIM::X, vel);
    if (phi->spatial_dim > DIM::D1) {
        phi->get_derivative_UCCD(d, DIM::Y, vel);
    }
    if (phi->spatial_dim > DIM::D2) {
        phi->get_derivative_UCCD(d, DIM::Z, vel);
    }

    if (phi->spatial_dim == DIM::D1) {
        for (int k = 0; k < phi->NodeSize.z; k++) {
            for (int j = 0; j < phi->NodeSize.y; j++) {
                for (int i = 0; i < phi->NodeSize.x; i++) {
                    double flux = phi->df[d]->get(i, j, k, DIM::X) * vel->f[VAR::U]->get(i, j, k);
                    phi->flux[d]->set(-flux, i, j, k);
                }
            }
        }
    } else if (phi->spatial_dim == DIM::D2) {
        for (int k = 0; k < phi->NodeSize.z; k++) {
            for (int j = 0; j < phi->NodeSize.y; j++) {
                for (int i = 0; i < phi->NodeSize.x; i++) {
                    double flux = phi->df[d]->get(i, j, k, DIM::X) * vel->f[VAR::U]->get(i, j, k) +
                                  phi->df[d]->get(i, j, k, DIM::Y) * vel->f[VAR::V]->get(i, j, k);
                    phi->flux[d]->set(-flux, i, j, k);
                }
            }
        }
    } else {
        for (int k = 0; k < phi->NodeSize.z; k++) {
            for (int j = 0; j < phi->NodeSize.y; j++) {
                for (int i = 0; i < phi->NodeSize.x; i++) {
                    double flux = phi->df[d]->get(i, j, k, DIM::X) * vel->f[VAR::U]->get(i, j, k) +
                                  phi->df[d]->get(i, j, k, DIM::Y) * vel->f[VAR::V]->get(i, j, k) +
                                  phi->df[d]->get(i, j, k, DIM::Z) * vel->f[VAR::W]->get(i, j, k);
                    phi->flux[d]->set(-flux, i, j, k);
                }
            }
        }
    }
}