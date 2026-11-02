// C0：直接调用上游距离、屏障和CCD；不运行物理求解器。
#include <ipc/barrier/barrier.hpp>
#include <ipc/ccd/point_static_plane.hpp>
#include <ipc/distance/point_plane.hpp>
#include <ipc/distance/point_triangle.hpp>
#include <ipc/distance/signed/point_plane.hpp>
#include <cmath>
#include <iomanip>
#include <iostream>
int main() {
    bool pass=true;
    auto check=[&](const char* name,bool ok,double value){
        std::cout<<name<<" value="<<std::setprecision(17)<<value<<" pass="<<ok<<'\n';pass &= ok;
    };
    const double d=.005, dh=.01, q=d*d, qh=dh*dh;
    const double expected=-(q-qh)*(q-qh)*std::log(q/qh);
    double b=ipc::barrier(q,qh);
    check("barrier_manual",std::abs(b-expected)/std::abs(expected)<1e-8,b);
    double g=2*d*ipc::barrier_first_derivative(q,qh);
    double H=4*q*ipc::barrier_second_derivative(q,qh)+2*ipc::barrier_first_derivative(q,qh);
    const double e=1e-6;
    auto B=[&](double z){return ipc::barrier(z*z,qh);};
    auto G=[&](double z){return 2*z*ipc::barrier_first_derivative(z*z,qh);};
    const double gfd=(B(d+e)-B(d-e))/(2*e),hfd=(G(d+e)-G(d-e))/(2*e);
    check("barrier_gradient_fd",std::abs(gfd-g)/std::abs(g)<1e-5,std::abs(gfd-g)/std::abs(g));
    check("barrier_hessian_fd",std::abs(hfd-H)/std::abs(H)<1e-5,std::abs(hfd-H)/std::abs(H));
    Eigen::Vector3d p(0,0,d),minus(0,0,-d),o(0,0,0),n(0,0,1),a(0,0,0),c(1,0,0),f(0,1,0);
    check("distance_is_squared",std::abs(ipc::point_plane_distance(p,o,n)-q)<1e-15,ipc::point_plane_distance(p,o,n));
    check("signed_plane_positive",ipc::point_plane_signed_distance(p,a,c,f)>0,ipc::point_plane_signed_distance(p,a,c,f));
    check("signed_plane_negative",ipc::point_plane_signed_distance(minus,a,c,f)<0,ipc::point_plane_signed_distance(minus,a,c,f));
    // 顶点-三角元：投影落在三角形内部时，平方距离应等于法向偏移的平方。
    Eigen::Vector3d pt(0.2,0.2,d),t0(0,0,0),t1(1,0,0),t2(0,1,0);
    double qt=ipc::point_triangle_distance(pt,t0,t1,t2);
    check("point_triangle_distance_manual",std::abs(qt-q)<1e-15,qt);
    const ipc::Vector12d qtg=ipc::point_triangle_distance_gradient(pt,t0,t1,t2);
    const double bpq=ipc::barrier_first_derivative(qt,qh);
    const Eigen::Vector3d gpt=bpq*qtg.head<3>();
    Eigen::Vector3d gpt_fd;
    for(int i=0;i<3;++i){
        Eigen::Vector3d pp=pt; pp[i]+=e;
        const double b1=ipc::barrier(ipc::point_triangle_distance(pp,t0,t1,t2),qh);
        pp=pt; pp[i]-=e;
        const double b0=ipc::barrier(ipc::point_triangle_distance(pp,t0,t1,t2),qh);
        gpt_fd[i]=(b1-b0)/(2*e);
    }
    check("point_triangle_barrier_gradient_fd",(gpt_fd-gpt).norm()/gpt.norm()<1e-5,(gpt_fd-gpt).norm()/gpt.norm());
    double toi=1;
    ipc::VectorMax3d p0(3),away(3),toward(3),zero(3);
    p0<<0.,0.,.01; away<<0.,0.,.02; toward<<0.,0.,-.01; zero.setZero();
    bool hit=ipc::point_static_plane_ccd(p0,away,o,n,toi);
    check("ccd_away",!hit,toi);
    toi=1;hit=ipc::point_static_plane_ccd(p0,toward,o,n,toi);
    check("ccd_approach",hit&&toi>=0&&toi<.5&&(.01+toi*(-.02))>0,toi);
    toi=1;hit=ipc::point_static_plane_ccd(zero,zero,o,n,toi);
    check("ccd_initial_touch",hit&&toi==0,toi);
    return pass?0:1;
}
