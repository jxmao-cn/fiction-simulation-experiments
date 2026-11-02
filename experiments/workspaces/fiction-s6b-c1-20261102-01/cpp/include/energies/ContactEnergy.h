#pragma once
// =============================================================================
// ContactEnergy.h — IPC Toolkit v1.6.1 网格/平面接触能量（隔离实验）
// =============================================================================
// 仅用于 S6b C1 隔离 workspace：把 toolkit 的候选检测、屏障势与 CCD 接入
// fiction 的 IEnergy / TimeIntegrator。生产默认路径不包含本文件。
//
// 约定：
//   - val/addGrad/addHess 返回未乘 h² 的屏障量，TimeIntegrator 统一乘 h²；
//   - 接触集在每个 Newton 迭代由 preIteration(x) 刷新，线搜索期间冻结；
//   - 力诊断用未缩放的 -∇B，不能用乘过 h² 的梯度；
//   - 几何相邻由 toolkit 宽相排除，体级过滤用碰撞网格顶点 patch。
// =============================================================================

#include "energy.h"

#include <ipc/ipc.hpp>
#include <ipc/barrier/adaptive_stiffness.hpp>
#include <ipc/candidates/candidates.hpp>
#include <ipc/collisions/normal/normal_collisions.hpp>
#include <ipc/potentials/barrier_potential.hpp>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <stdexcept>
#include <unordered_set>
#include <vector>

namespace fiction {

class ContactEnergy : public IEnergy {
public:
    struct Options {
        double dhat = 1.7e-3;       ///< 激活距离 [m]
        double dmin = 0.0;          ///< 最小距离（平面接触必须为 0）
        double kappa = 0.0;         ///< <=0 时在 initKappa 中确定
        bool useAreaWeighting = false;
    };

    ContactEnergy(const Positions& x0,
                  const Eigen::MatrixXi& faces,
                  const Eigen::VectorXi& bodyIdFull,
                  Options opt)
        : opt_(opt)
    {
        if (faces.rows() == 0 || faces.cols() != 3)
            throw std::invalid_argument("ContactEnergy: 需要 [S,3] 表面三角面");
        if (!(opt_.dhat > 0.0) || !(opt_.dmin >= 0.0) || opt_.dmin >= opt_.dhat)
            throw std::invalid_argument("ContactEnergy: 需要 0 <= dmin < dhat");

        edges_ = uniqueEdges(faces);
        const Eigen::MatrixXd X0 = toMatrix(x0);
        cm_ = ipc::CollisionMesh::build_from_full_mesh(X0, edges_, faces);
        cm_.init_adjacencies();
        cm_.init_area_jacobians();

        Eigen::VectorXi patch = Eigen::VectorXi::Zero(cm_.num_vertices());
        if (bodyIdFull.size() == x0.rows()) {
            for (int i = 0; i < cm_.num_vertices(); ++i)
                patch[i] = bodyIdFull[cm_.to_full_vertex_id(i)];
        } else if (bodyIdFull.size() != 0) {
            throw std::invalid_argument("ContactEnergy: bodyIdFull 大小必须为 0 或 N");
        }
        cm_.can_collide = ipc::make_vertex_patches_filter(std::move(patch));

        barrier_ = std::make_shared<ipc::BarrierPotential>(
            opt_.dhat, opt_.kappa > 0.0 ? opt_.kappa : 1.0, false);
        collisions_.set_use_area_weighting(opt_.useAreaWeighting);
        refresh(x0);
    }

    const char* name() const override { return "ContactEnergy"; }

    void addPlane(const Eigen::Vector3d& normal, double offset) {
        cm_.planes.emplace_back(normal.normalized(), offset);
    }

    /// 在当前 x 上重建候选与活跃接触集；必须在每次 Newton 迭代组装前调用。
    void preIteration(const Positions& x) override { refresh(x); }

    void refresh(const Positions& x) {
        const Eigen::MatrixXd V = cm_.vertices(toMatrix(x));
        candidates_.build(cm_, V, 0.5 * (opt_.dhat + opt_.dmin));
        collisions_.build(candidates_, cm_, V, opt_.dhat, opt_.dmin);
        V_ = V;
        minDistanceSq_ = collisions_.compute_minimum_distance(cm_, V_);
    }

    double val(const Positions& x) const override {
        // 线搜索在 xTrial 上求值；碰撞集冻结，但位置必须用当前实参。
        const Eigen::MatrixXd V = cm_.vertices(toMatrix(x));
        return barrier_->operator()(collisions_, cm_, V);
    }

    void addGrad(const Positions& /*x*/, Positions& grad) const override {
        const Eigen::VectorXd gc = barrier_->gradient(collisions_, cm_, V_);
        const Eigen::VectorXd gf = cm_.to_full_dof(gc);
        for (int i = 0; i < grad.rows(); ++i)
            for (int d = 0; d < kDim; ++d)
                grad(i, d) += gf(3 * i + d);
    }

    void addHess(const Positions& /*x*/, std::vector<Triplet>& triplets) const override {
        const Eigen::SparseMatrix<double> Hc = barrier_->hessian(
            collisions_, cm_, V_, ipc::PSDProjectionMethod::CLAMP);
        const Eigen::SparseMatrix<double> Hf = cm_.to_full_dof(Hc);
        for (int k = 0; k < Hf.outerSize(); ++k)
            for (Eigen::SparseMatrix<double>::InnerIterator it(Hf, k); it; ++it)
                triplets.emplace_back(it.row(), it.col(), it.value());
    }

    /// CCD + 非候选保守步长；返回 [0,1]。
    double initStepSize(const Positions& x, const Positions& p) const override {
        const Eigen::MatrixXd V0 = cm_.vertices(toMatrix(x));
        const Eigen::MatrixXd V1 = cm_.vertices(toMatrix(x + p));
        const double aCcd = ipc::compute_collision_free_stepsize(
            cm_, V0, V1, opt_.dmin);
        const double aNc = candidates_.compute_noncandidate_conservative_stepsize(
            cm_, cm_.map_displacements(toMatrix(p)), opt_.dhat);
        lastAlphaCcd_ = aCcd;
        lastAlphaNc_ = aNc;
        lastAlpha_ = std::min({1.0, aCcd, aNc});
        return lastAlpha_;
    }

    /// 用 t=0 的弹性+重力梯度确定一次 κ，之后冻结。
    double initKappa(const Positions& x0,
                     double bboxDiagonal,
                     double averageMass,
                     const Eigen::VectorXd& gradEnergy) {
        refresh(x0);
        const Eigen::VectorXd gb = fullGradient(x0);
        double kappaMax = 0.0;
        const ipc::ClampedLogBarrier barrier;
        const double kappa = ipc::initial_barrier_stiffness(
            bboxDiagonal, barrier, opt_.dhat, averageMass, gradEnergy, gb, kappaMax,
            1e11, opt_.dmin);
        opt_.kappa = kappa;
        barrier_->set_stiffness(kappa);
        kappaMax_ = kappaMax;
        return kappa;
    }

    /// 接触力（未乘 h²），返回 3N 向量。
    Eigen::VectorXd forces(const Positions& x) const {
        const Eigen::MatrixXd V = cm_.vertices(toMatrix(x));
        const Eigen::VectorXd gc = barrier_->gradient(collisions_, cm_, V);
        return -cm_.to_full_dof(gc);
    }

    int contactCount() const { return static_cast<int>(collisions_.size()); }
    int candidateCount() const { return static_cast<int>(candidates_.size()); }
    double minDistance() const {
        return std::isfinite(minDistanceSq_) ? std::sqrt(std::max(0.0, minDistanceSq_))
                                             : std::numeric_limits<double>::infinity();
    }
    double kappa() const { return opt_.kappa; }
    double lastAlphaCcd() const { return lastAlphaCcd_; }
    double lastAlphaNc() const { return lastAlphaNc_; }
    double lastAlpha() const { return lastAlpha_; }
    double kappaMax() const { return kappaMax_; }
    int collisionVertices() const { return static_cast<int>(cm_.num_vertices()); }
    int surfaceEdges() const { return static_cast<int>(edges_.rows()); }
    const ipc::CollisionMesh& collisionMesh() const { return cm_; }

private:
    static Eigen::MatrixXd toMatrix(const Positions& x) {
        return Eigen::MatrixXd(x);
    }

    static Eigen::MatrixXi uniqueEdges(const Eigen::MatrixXi& faces) {
        std::unordered_set<std::uint64_t> seen;
        std::vector<std::array<int, 2>> edges;
        seen.reserve(static_cast<size_t>(faces.rows()) * 3);
        auto key = [](int a, int b) {
            if (a > b) std::swap(a, b);
            return (static_cast<std::uint64_t>(a) << 32) | static_cast<std::uint32_t>(b);
        };
        for (int f = 0; f < faces.rows(); ++f) {
            for (int e = 0; e < 3; ++e) {
                const int a = faces(f, e);
                const int b = faces(f, (e + 1) % 3);
                if (seen.insert(key(a, b)).second) edges.push_back({a, b});
            }
        }
        Eigen::MatrixXi out(static_cast<int>(edges.size()), 2);
        for (int i = 0; i < out.rows(); ++i) {
            out(i, 0) = edges[i][0];
            out(i, 1) = edges[i][1];
        }
        return out;
    }

    Eigen::VectorXd fullGradient(const Positions& x) const {
        const Eigen::MatrixXd V = cm_.vertices(toMatrix(x));
        const Eigen::VectorXd gc = barrier_->gradient(collisions_, cm_, V);
        return cm_.to_full_dof(gc);
    }

    Options opt_;
    Eigen::MatrixXi edges_;
    ipc::CollisionMesh cm_;
    ipc::Candidates candidates_;
    ipc::NormalCollisions collisions_;
    std::shared_ptr<ipc::BarrierPotential> barrier_;
    Eigen::MatrixXd V_;
    double minDistanceSq_ = std::numeric_limits<double>::infinity();
    double kappaMax_ = 0.0;
    mutable double lastAlphaCcd_ = 1.0;
    mutable double lastAlphaNc_ = 1.0;
    mutable double lastAlpha_ = 1.0;
};

} // namespace fiction
