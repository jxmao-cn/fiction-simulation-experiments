#include <suitesparse/cholmod.h>

#include <cmath>
#include <iostream>
#include <vector>

int main() {
    cholmod_common Common, *c = &Common;
    if (!cholmod_start(c)) return 1;

    const int n = 5;
    const double k = 2e4, h = 0.01, node_mass = 40.0;
    std::vector<double> diag(n), off(n - 1, -k * h * h);
    for (int i = 0; i < n; ++i) diag[i] = node_mass + 2 * k * h * h;
    diag[0] += k * h * h;
    diag[n - 1] += k * h * h;

    std::vector<int> colp(n + 1), rowi;
    std::vector<double> val;
    for (int j = 0; j < n; ++j) {
        colp[j] = static_cast<int>(rowi.size());
        if (j > 0)      { rowi.push_back(j - 1); val.push_back(off[j - 1]); }
                        { rowi.push_back(j);     val.push_back(diag[j]); }
        if (j < n - 1)  { rowi.push_back(j + 1); val.push_back(off[j]); }
    }
    colp[n] = static_cast<int>(rowi.size());
    const int nnz = static_cast<int>(rowi.size());

    cholmod_sparse* A = cholmod_allocate_sparse(n, n, nnz, true, true, -1, CHOLMOD_REAL, c);
    for (int j = 0; j <= n; ++j) static_cast<int*>(A->p)[j] = colp[j];
    for (int t = 0; t < nnz; ++t) {
        static_cast<int*>(A->i)[t] = rowi[t];
        static_cast<double*>(A->x)[t] = val[t];
    }

    cholmod_dense* b = cholmod_allocate_dense(n, 1, n, CHOLMOD_REAL, c);
    for (int i = 0; i < n; ++i) static_cast<double*>(b->x)[i] = node_mass * 9.81;

    cholmod_factor* L = cholmod_analyze(A, c);
    cholmod_factorize(A, L, c);
    std::cout << "cholmod status = " << c->status << " (0 = OK), L->minor = "
              << static_cast<long>(L->minor) << " / " << n << "\n";
    cholmod_dense* y = cholmod_solve(CHOLMOD_A, L, b, c);
    double* yx = static_cast<double*>(y->x);

    double max_res = 0.0;
    for (int j = 0; j < n; ++j) {
        double Ay = 0.0;
        for (int t = colp[j]; t < colp[j + 1]; ++t) Ay += val[t] * yx[rowi[t]];
        max_res = std::max(max_res, std::abs(static_cast<double*>(b->x)[j] - Ay));
    }
    std::cout.precision(12);
    for (int i = 0; i < n; ++i) std::cout << "y[" << i << "] = " << yx[i] << "\n";
    std::cout << "max |b - A y| = " << max_res << "\n";

    cholmod_free_dense(&y, c);
    cholmod_free_factor(&L, c);
    cholmod_free_dense(&b, c);
    cholmod_free_sparse(&A, c);
    cholmod_finish(c);
    return (max_res < 1e-10 && c->status == CHOLMOD_OK) ? 0 : 1;
}
