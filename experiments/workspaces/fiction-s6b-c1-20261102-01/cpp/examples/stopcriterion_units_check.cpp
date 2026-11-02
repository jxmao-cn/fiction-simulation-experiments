// 隔离停止判据的量纲和错误早停负例。
#include "time_integrator.h"
#include <cmath>
#include <iostream>
#include <memory>
using namespace fiction;
class Quadratic final : public IEnergy {
public:
 const char* name() const override { return "QuadraticTinyStep"; }
 double val(const Positions& x) const override { return .5*x.squaredNorm(); }
 void addGrad(const Positions& x, Positions& g) const override {g+=x;}
 void addHess(const Positions& x,std::vector<Triplet>& t) const override {
  for(int i=0;i<x.size();++i)t.emplace_back(i,i,1.);
 }
 double initStepSize(const Positions&,const Positions&) const override {return 1e-14;}
};
TimeIntegrator::State run(int n,double h,bool stat,double threshold,int maxIter) {
 TimeIntegrator::Config c;c.h=h;c.staticSolve=stat;c.maxIter=maxIter;c.verbose=false;
 c.tol=1e-12;c.stopMode="combined";c.gradScale=1;c.gradTolRel=threshold;
 TimeIntegrator it(c,Eigen::VectorXd::Ones(n));it.addEnergy(std::make_shared<Quadratic>());
 return it.stepForward(Positions::Ones(n,3),Velocities::Zero(n,3),DBCInfo{});
}
int main(){
 bool ok=true;
 for(double h:{.005,.01,.02}){
  auto s=run(1,h,false,.1,0);
  bool pass=std::abs(s.report.forceResidualInf-1)<1e-12 &&
   std::abs(s.report.freeGradientInf-h*h)<1e-15;
  std::cout<<"dt="<<h<<" IP="<<s.report.freeGradientInf<<" force_N="<<s.report.forceResidualInf<<" pass="<<pass<<'\n';ok &= pass;
 }
 auto strict=run(1,1,true,.1,4),loose=run(1,1,true,2.,4),copies=run(4,1,true,.1,4);
 bool rejected=!strict.report.converged && strict.report.maxIterHit;
 bool accepted=loose.report.converged && loose.report.stopReason=="combined";
 bool invariant=!copies.report.converged && std::abs(copies.report.forceResidualInf-strict.report.forceResidualInf)<1e-12;
 std::cout<<"large_force_tiny_motion_rejected="<<rejected<<" configured_gradTol_effective="<<accepted<<" four_copy_invariant="<<invariant<<'\n';
 return ok && rejected && accepted && invariant ? 0:1;
}
