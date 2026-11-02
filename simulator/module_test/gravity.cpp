#include <iostream>
#include <iomanip>

int main() {
    const double G = 6.674e-11;

    double m1, m2, r;
    std::cout << "输入 m1(kg) m2(kg) r(m): ";
    if (!(std::cin >> m1 >> m2 >> r) || r <= 0) {
        std::cerr << "无效输入\n";
        return 1;
    }

    double F = G * m1 * m2 / (r * r);
    std::cout << std::scientific << std::setprecision(6)
              << "引力 F = " << F << " N\n";

    double mEarth = 5.972e24, mMoon = 7.348e22, dEM = 3.844e8;
    std::cout << "校验(地月): F = "
              << G * mEarth * mMoon / (dEM * dEM) << " N\n";
    return 0;
}
