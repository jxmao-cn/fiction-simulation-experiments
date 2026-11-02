#include <iostream>
#include <iomanip>
#include <array>
#include <cmath>

constexpr int N = 3;
const double mu = 1.0;

double p(std::array<double, N> x) {
    double r2 = 0;
    for (double v : x) r2 += v * v;
    return -mu / std::sqrt(r2);
}

int main() {
    std::array<double, N> x = {1, 2, 2};
    double r2 = 0;
    for (double v : x) r2 += v * v;
    double r = std::sqrt(r2), r3 = r2 * r, r5 = r2 * r2 * r;
    double h = 1e-4 * r;

    std::cout << std::scientific << std::setprecision(6);
    std::cout << "r = " << r << "  h = " << h << "\n\n梯度检验:\n";
    for (int i = 0; i < N; ++i) {
        auto xp = x, xm = x;
        xp[i] += h;
        xm[i] -= h;
        double num = (p(xp) - p(xm)) / (2 * h);
        double ana = mu * x[i] / r3;
        std::cout << "  dP/dx" << i + 1 << "  解析 " << ana
                  << "  数值 " << num << "  误差 " << std::abs(ana - num) << "\n";
    }

    std::cout << "\n海森检验:\n";
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) {
            auto a = x, b = x, c = x, d = x;
            a[i] += h; a[j] += h;
            b[i] += h; b[j] -= h;
            c[i] -= h; c[j] += h;
            d[i] -= h; d[j] -= h;
            double num = (p(a) - p(b) - p(c) + p(d)) / (4 * h * h);
            double ana = (i == j ? mu / r3 : 0.0) - 3 * mu * x[i] * x[j] / r5;
            std::cout << "  H[" << i + 1 << "][" << j + 1 << "]  解析 " << ana
                      << "  数值 " << num << "  误差 " << std::abs(ana - num) << "\n";
        }

    double tr = 0;
    for (int i = 0; i < N; ++i) tr += mu / r3 - 3 * mu * x[i] * x[i] / r5;
    std::cout << "\ntrace(H) = " << tr << "  （n=3 应为 0 → ∇²p = 0）\n";
}
