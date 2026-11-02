// =============================================================================
// material_validation.cpp — S0 四面体 Neo-Hookean 局部材料验证
// 独立可执行：material_validation [可选 metrics.csv]
// 任一检查失败、结果非有限或 CSV 写入失败均返回非零。
// =============================================================================

#include "energies/NeoHookeanTetEnergy.h"
#include "scene/mass_spring_builder.h"

#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {
using namespace fiction;
using Mat12 = Eigen::Matrix<double, 12, 12>;

int failures = 0;
void check(bool ok, const std::string& name) {
    std::cout << (ok ? "  [PASS] " : "  [FAIL] ") << name << '\n';
    if (!ok) ++failures;
}
bool finite(double x) { return std::isfinite(x); }

Eigen::MatrixXd denseRaw(const NeoHookeanTetEnergy& e, const Positions& x) {
    std::vector<Triplet> ts;
    e.addHessRaw(x, ts);
    Eigen::MatrixXd h = Eigen::MatrixXd::Zero(3 * x.rows(), 3 * x.rows());
    for (const auto& t : ts) h(t.row(), t.col()) += t.value();
    return h;
}
Eigen::MatrixXd densePsd(const NeoHookeanTetEnergy& e, const Positions& x) {
    std::vector<Triplet> ts;
    e.addHess(x, ts);
    Eigen::MatrixXd h = Eigen::MatrixXd::Zero(3 * x.rows(), 3 * x.rows());
    for (const auto& t : ts) h(t.row(), t.col()) += t.value();
    return h;
}
Positions gradient(const NeoHookeanTetEnergy& e, const Positions& x) {
    Positions g = Positions::Zero(x.rows(), 3);
    e.addGrad(x, g);
    return g;
}
Positions referenceTet() {
    Positions x(4, 3);
    x << 0,0,0, 1,0,0, 0,1,0, 0,0,1;
    return x;
}
const std::vector<std::array<int,4>> oneTet{{0,1,2,3}};

void checkRestAndObjectivity() {
    std::cout << "\n=== S0: rest, translation, objectivity ===\n";
    const Positions rest = referenceTet();
    NeoHookeanTetEnergy e(oneTet, rest, 2e5, 0.32);
    const double er = e.val(rest);
    const Positions gr = gradient(e, rest);
    std::cout << "  rest: energy=" << er << " gradient_norm=" << gr.norm() << '\n';
    check(finite(er) && gr.allFinite() && std::abs(er) < 1e-10 && gr.norm() < 1e-9,
          "参考态能量与弹性梯度为零");
    Positions translated = rest;
    translated.rowwise() += Eigen::RowVector3d(0.37, -1.2, 0.8);
    check(finite(e.val(translated)) && std::abs(e.val(translated)-er) < 1e-10 &&
          (gradient(e, translated)-gr).norm() < 1e-9,
          "只平移当前位形不改变能量与梯度");

    const Eigen::Matrix3d R = Eigen::AngleAxisd(0.73, Eigen::Vector3d(1,2,3).normalized()).toRotationMatrix();
    Positions rotated = rest * R.transpose();
    rotated.rowwise() += Eigen::RowVector3d(0.4,-1.0,2.0);
    // 参考网格固定，只转当前位形；这是检验客观性的关键。
    const double erig = e.val(rotated);
    const Positions grig = gradient(e, rotated);
    const double rigidGrad = grig.norm();
    std::cout << "  fixed-reference rigid rotation: energy=" << erig << " gradient_norm=" << rigidGrad
              << " |dE-Erest|=" << std::abs(erig-er) << '\n';
    check(finite(erig) && grig.allFinite() && finite(rigidGrad) && std::abs(erig) < 1e-9 && rigidGrad < 1e-8,
          "固定参考网格、仅刚体旋转当前位形：零能量/零梯度");

    Eigen::Matrix3d A;
    A << 1.12, 0.18, 0.03, 0.02, 0.94, 0.08, 0.0, 0.04, 1.06;
    Positions def = rest * A.transpose();
    Positions defRot = def * R.transpose();
    defRot.rowwise() += Eigen::RowVector3d(-0.3,0.7,0.2);
    const double e0 = e.val(def), e1 = e.val(defRot);
    const Positions g0 = gradient(e, def), g1 = gradient(e, defRot);
    const double ecov = std::abs(e0-e1) / std::max(1.0,std::abs(e0));
    const double gcov = (g1-g0*R.transpose()).norm() / std::max(1.0,g0.norm());
    std::cout << "  deformed rotation: E=" << e0 << " rotated_E=" << e1 << " relative dE=" << ecov
              << ", gradient covariance=" << gcov << '\n';
    check(finite(e0) && finite(e1) && g0.allFinite() && g1.allFinite() && ecov < 1e-10 && gcov < 1e-9,
          "变形后叠加刚体旋转：能量不变且梯度随旋转协变");
}

struct CsvRow { std::string state; double eps, gradErr, hessErr; bool pass; };
std::vector<CsvRow> rows;
void finiteDifferenceCase(const std::string& label, const Eigen::Matrix3d& F) {
    const Positions rest = referenceTet();
    NeoHookeanTetEnergy e(oneTet, rest, 2e5, 0.32);
    const Positions x = rest * F.transpose();
    const Positions ga = gradient(e,x);
    const Eigen::MatrixXd H = denseRaw(e,x); // 与 FD 比较的必须是未投影精确 Hessian
    const double es = e.val(x);
    const double gscale = std::max(ga.norm(),1e-12);
    const double hscale = std::max(H.norm(),1e-12);
    const std::vector<double> epsilons{1e-4,1e-5,1e-6,1e-7};
    std::vector<bool> pass;
    std::cout << "\n  FD " << label << " J=" << F.determinant() << " E=" << std::setprecision(10) << es << '\n';
    for (double eps : epsilons) {
        Eigen::VectorXd gfd(12);
        Eigen::MatrixXd hfd(12,12);
        for (int k=0;k<12;++k) {
            Positions xp=x,xm=x;
            xp(k/3,k%3)+=eps; xm(k/3,k%3)-=eps;
            gfd[k]=(e.val(xp)-e.val(xm))/(2*eps);
            hfd.col(k)=(Eigen::Map<const Eigen::VectorXd>(gradient(e,xp).data(),12)-
                        Eigen::Map<const Eigen::VectorXd>(gradient(e,xm).data(),12))/(2*eps);
        }
        const double ge=(Eigen::Map<const Eigen::VectorXd>(ga.data(),12)-gfd).norm()/gscale;
        const double he=(H-hfd).norm()/hscale;
        const bool ok=finite(ge)&&finite(he)&&gfd.allFinite()&&hfd.allFinite()&&ge<1e-6&&he<1e-4;
        rows.push_back({label,eps,ge,he,ok}); pass.push_back(ok);
        std::cout << "    eps=" << eps << " grad_rel=" << ge << " raw_hess_rel=" << he << (ok?" PASS":" FAIL") << '\n';
    }
    int consecutive=0,best=0;
    for(bool p:pass) { consecutive=p?consecutive+1:0; best=std::max(best,consecutive); }
    check(F.determinant()>0.2 && F.allFinite() && finite(es),label+" 正体积且计算有限");
    check(best>=2,label+" 梯度<1e-6、原始 Hessian<1e-4 至少连续两档 eps");
}

void checkPsdAndParameterResponse() {
    std::cout << "\n=== S0: PSD and parameter response ===\n";
    const Positions rest=referenceTet();
    const std::vector<std::pair<std::string,Eigen::Matrix3d>> psdStates={
        {"stretch",(Eigen::Matrix3d)(Eigen::DiagonalMatrix<double,3>(1.15,1.0,1.0))},
        {"shear",[] { Eigen::Matrix3d s=Eigen::Matrix3d::Identity(); s(0,1)=0.3; return s; }()},
        {"compression",(Eigen::Matrix3d)(Eigen::DiagonalMatrix<double,3>(0.8,0.8,0.9))},
    };
    bool sawIndefiniteRaw=false;
    for(const auto& [name,Fs]:psdStates) {
        Positions xs=rest*Fs.transpose();
        NeoHookeanTetEnergy es0(oneTet,rest,2e5,0.32);
        const Eigen::MatrixXd Hr=denseRaw(es0,xs);
        const Eigen::MatrixXd Hp=densePsd(es0,xs);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(0.5*(Hp+Hp.transpose()));
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> rawEig(0.5*(Hr+Hr.transpose()));
        const double sym=(Hp-Hp.transpose()).norm()/std::max(1.0,Hp.norm());
        const double mineig=eig.info()==Eigen::Success?eig.eigenvalues().minCoeff():-std::numeric_limits<double>::infinity();
        const double rawMin=rawEig.info()==Eigen::Success?rawEig.eigenvalues().minCoeff():std::numeric_limits<double>::quiet_NaN();
        Eigen::MatrixXd clipped=Eigen::MatrixXd::Zero(12,12);
        if(rawEig.info()==Eigen::Success) clipped=rawEig.eigenvectors()*rawEig.eigenvalues().cwiseMax(0.0).asDiagonal()*rawEig.eigenvectors().transpose();
        const double projectionErr=(Hp-clipped).norm()/std::max(1.0,clipped.norm());
        if(finite(rawMin)&&rawMin < -1e-8*std::max(1.0,Hr.norm())) sawIndefiniteRaw=true;
        std::cout << "  " << name << " raw_min_eig="<<rawMin<<" projected_symmetry="<<sym
                  <<" projected_min_eig="<<mineig<<" spectral_clip_error="<<projectionErr<<'\n';
        check(Hr.allFinite()&&Hp.allFinite()&&finite(sym)&&finite(mineig)&&finite(rawMin)&&finite(projectionErr)&&
              sym<1e-10&&mineig>-1e-8*std::max(1.0,Hp.norm())&&projectionErr<1e-10,
              name+" PSD 投影符合独立谱裁剪，且有限、对称、半正定");
    }
    check(sawIndefiniteRaw,"拉伸/剪切/压缩状态中至少一个原始 Hessian 确实不定，避免恒等投影假通过");

    Eigen::Matrix3d mixed; mixed<<1.1,0.15,0, 0,0.96,0.07, 0.02,0,1.03;
    const Positions x=rest*mixed.transpose();
    NeoHookeanTetEnergy e1(oneTet,rest,1e5,0.30), e2(oneTet,rest,2e5,0.30);
    const double ratio=e2.val(x)/e1.val(x);
    const double gradRatio=gradient(e2,x).norm()/gradient(e1,x).norm();
    const double hessRatio=denseRaw(e2,x).norm()/denseRaw(e1,x).norm();
    std::cout << "  fixed nu=0.30, E scale ratios: energy="<<ratio<<" gradient="<<gradRatio
              <<" raw Hessian="<<hessRatio<<" (expected 2)\n";
    check(finite(ratio)&&finite(gradRatio)&&finite(hessRatio)&&std::abs(ratio-2.0)<1e-10&&
          std::abs(gradRatio-2.0)<1e-10&&std::abs(hessRatio-2.0)<1e-10,
          "固定泊松比，杨氏模量翻倍使能量/梯度/原始 Hessian 翻倍");

    Eigen::Matrix3d hydro=1.04*Eigen::Matrix3d::Identity();
    Positions xh=rest*hydro.transpose();
    NeoHookeanTetEnergy nu1(oneTet,rest,2e5,0.20),nu2(oneTet,rest,2e5,0.40);
    const Positions gn1=gradient(nu1,xh),gn2=gradient(nu2,xh);
    // 沿 x(s)=s*x_ref 的体积加载方向取梯度方向导数；节点力总和会因平移不变性为零。
    const double p1=(gn1.array()*rest.array()).sum();
    const double p2=(gn2.array()*rest.array()).sum();
    auto dEdScale=[&](double nu) {
        const double mu=2e5/(2*(1+nu));
        const double lam=2e5*nu/((1+nu)*(1-2*nu));
        const double s=1.04, J=s*s*s;
        return (1.0/6.0)*3.0*(mu*(s-1.0/s)+lam*std::log(J)/s);
    };
    const double p1ref=dEdScale(0.20),p2ref=dEdScale(0.40);
    std::cout << "  hydrostatic stretch 1.04: nu=.20 E="<<nu1.val(xh)<<" dE/ds="<<p1<<" analytic="<<p1ref
              <<"; nu=.40 E="<<nu2.val(xh)<<" dE/ds="<<p2<<" analytic="<<p2ref<<'\n';
    check(finite(nu1.val(xh))&&finite(nu2.val(xh))&&gn1.allFinite()&&gn2.allFinite()&&
          finite(p1ref)&&finite(p2ref)&&std::abs(p1-p1ref)<1e-8*std::max(1.0,std::abs(p1ref))&&
          std::abs(p2-p2ref)<1e-8*std::max(1.0,std::abs(p2ref))&&p2>p1,
          "固定 E 时 nu=.20/.40 的体积加载导数符合解析值且随 nu 增大");
}

void checkAssemblyAndMass() {
    std::cout << "\n=== S0: two-tet assembly and mass/volume ===\n";
    MshMesh mesh; mesh.x.resize(5,3);
    mesh.x<<0,0,0, 1,0,0, 0,1,0, 0,0,1, 0,0,-1;
    mesh.tets={{0,1,2,3},{0,2,1,4}};
    const auto model=buildMassSpringModel(mesh,1200.0,0.0);
    NeoHookeanTetEnergy all(mesh.tets,mesh.x,2e5,0.3);
    NeoHookeanTetEnergy a({mesh.tets[0]},mesh.x,2e5,0.3),b({mesh.tets[1]},mesh.x,2e5,0.3);
    Eigen::Matrix3d F; F<<1.07,0.08,0, 0.02,0.97,0.03, 0,0.01,1.04;
    Positions x=mesh.x*F.transpose();
    const Positions g=gradient(all,x), gs=gradient(a,x)+gradient(b,x);
    const Eigen::MatrixXd H=denseRaw(all,x), Hs=denseRaw(a,x)+denseRaw(b,x);
    const double ev=std::abs(all.val(x)-a.val(x)-b.val(x));
    const double gv=(g-gs).norm(),hv=(H-Hs).norm();
    std::cout<<"  assembly errors: energy="<<ev<<" grad="<<gv<<" raw Hess="<<hv
             <<" volume="<<model.volume<<" mass="<<model.totalMass<<'\n';
    check(finite(ev)&&finite(gv)&&finite(hv)&&ev<1e-10&&gv<1e-9&&hv<1e-8,
          "全局两单元能量/梯度/原始 Hessian 等于逐单元贡献之和");
    Eigen::VectorXd expectedMass(5); expectedMass<<100,100,100,50,50;
    const double massErr=(model.masses-expectedMass).cwiseAbs().maxCoeff();
    std::cout<<"  nodal masses actual="<<model.masses.transpose()<<" expected="<<expectedMass.transpose()
             <<" max_error="<<massErr<<'\n';
    check(std::abs(model.volume-1.0/3.0)<1e-12&&std::abs(model.totalMass-400.0)<1e-9&&massErr<1e-10,
          "builder 总体积、总质量及逐节点质量符合 rho*V/4 分摊");
    Eigen::VectorXd tx=Eigen::VectorXd::Zero(15);
    for(int i=0;i<5;++i) tx.segment<3>(3*i)=Eigen::Vector3d(0.4,-0.2,0.7);
    check(H.allFinite()&&(H*tx).norm()<1e-8*std::max(1.0,H.norm()),"参考配置附近刚体平移为全局 Hessian 零模态");
}
}

int main(int argc,char** argv) {
    std::cout<<std::scientific<<std::setprecision(6)<<"S0 Neo-Hookean material validation\n";
    checkRestAndObjectivity();
    finiteDifferenceCase("拉伸",(Eigen::Matrix3d)(Eigen::DiagonalMatrix<double,3>(1.15,1.0,1.0)));
    Eigen::Matrix3d shear=Eigen::Matrix3d::Identity(); shear(0,1)=0.3;
    finiteDifferenceCase("剪切",shear);
    finiteDifferenceCase("压缩",(Eigen::Matrix3d)(Eigen::DiagonalMatrix<double,3>(0.8,0.8,0.9)));
    checkPsdAndParameterResponse();
    checkAssemblyAndMass();
    if(argc>2) { std::cerr<<"usage: material_validation [metrics.csv]\n"; return 2; }
    if(argc==2) {
        std::ofstream out(argv[1]);
        if(!out) { std::cerr<<"无法写入 CSV: "<<argv[1]<<'\n'; return 2; }
        out<<"state,epsilon,gradient_relative_error,raw_hessian_relative_error,pass\n";
        for(const auto&r:rows) out<<r.state<<','<<std::setprecision(17)<<r.eps<<','<<r.gradErr<<','<<r.hessErr<<','<<(r.pass?1:0)<<'\n';
        out.flush();
        if(!out.good()) { std::cerr<<"CSV flush/写入失败: "<<argv[1]<<'\n'; return 2; }
    }
    std::cout<<"\nS0 result: "<<(failures?"FAIL":"PASS")<<" ("<<failures<<" failed checks)\n";
    return failures?1:0;
}
