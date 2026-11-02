#pragma once
// =============================================================================
// SpringEnergy.h — MassSpring 能量（纯 3D 版）
// =============================================================================
// E = sum 0.5*k*l2*(|d|^2/l2 - 1)^2，边拓扑 + makePSD 用法示例：
// 局部 Hessian 为 3D 边块 6×6，自由度索引 3*node+d
// =============================================================================

#include "energy.h"
#include "utils.h"

#include <array>
#include <vector>

namespace fiction {

class ExampleSpringEnergy final : public IEnergy {
public:
    ExampleSpringEnergy(const std::vector<std::array<int,2>>& edges,
                        const Eigen::VectorXd& l2,
                        const Eigen::VectorXd& k)
        : edges_(edges), l2_(l2), k_(k) {}
    const char* name() const override { return "ExampleSpring"; }
    double hessianWeight(double h) const override { return h*h; }
    double val(const Positions& x) const override {
        double E=0;
        for(size_t i=0;i<edges_.size();++i){
            Vec3 diff = x.row(edges_[i][0]).transpose() - x.row(edges_[i][1]).transpose();
            E += l2_[i]*0.5*k_[i]*std::pow(diff.squaredNorm()/l2_[i]-1.0,2);
        }
        return E;
    }
    void addGrad(const Positions& x, Positions& grad) const override {
        for(size_t i=0;i<edges_.size();++i){
            Vec3 diff = x.row(edges_[i][0]).transpose() - x.row(edges_[i][1]).transpose();
            Vec3 gdiff = 2*k_[i]*(diff.squaredNorm()/l2_[i]-1)*diff;
            grad.row(edges_[i][0]) += gdiff.transpose();
            grad.row(edges_[i][1]) -= gdiff.transpose();
        }
    }
    void addHess(const Positions& x, std::vector<Triplet>& trips) const override {
        for(size_t i=0;i<edges_.size();++i){
            Vec3 diff = x.row(edges_[i][0]).transpose() - x.row(edges_[i][1]).transpose();
            Eigen::Matrix3d Hdiff = 2*k_[i]/l2_[i]*(2*diff*diff.transpose()
                + (diff.squaredNorm()-l2_[i])*Eigen::Matrix3d::Identity());
            Eigen::Matrix<double,6,6> Hl;
            Hl <<  Hdiff, -Hdiff,
                  -Hdiff,  Hdiff;
            Hl = makePSD(Hl); // 3D 边块 6×6
            for(int nI=0;nI<2;++nI) for(int nJ=0;nJ<2;++nJ)
              for(int r=0;r<3;++r) for(int c=0;c<3;++c)
                trips.emplace_back(3*edges_[i][nI]+r, 3*edges_[i][nJ]+c, Hl(nI*3+r,nJ*3+c));
        }
    }
    double initStepSize(const Positions& /*x*/, const Positions& /*p*/) const override {
        return 1.0;
    }
private:
    std::vector<std::array<int,2>> edges_;
    Eigen::VectorXd l2_, k_;
};

} // namespace fiction
