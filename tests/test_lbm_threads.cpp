#include "TestUtil.hpp"
#include "HomeBank.hpp"
#include "lbm/River.hpp"
#include <cstring>
using namespace lbm;
using namespace lbmtest;
static bool bits(double a, double b) { return std::memcmp(&a, &b, sizeof a) == 0; }
int main() {
    for (int geometry = 0; geometry < 3; ++geometry) {
        DomainSpec spec; spec.nx = 33; spec.ny = 25;
        if (geometry == 0) spec.xMin = spec.xMax = spec.yMin = spec.yMax = SideKind::Periodic;
        if (geometry == 1) {
            spec.phi = [](double x, double y) { return std::hypot(x-16., y-12.)-4.2; };
            spec.bodyAt = [](double, double) { return 1; };
        }
        if (geometry == 2) { spec.yMin = SideKind::PressureOutlet; spec.yMax = SideKind::MovingWall; }
        Domain domain(spec);
        for (double kt : {0., 1e-8}) for (int threads : {2, 6}) {
            CollisionConfig cfg = RiverParams::defaultCollision(); cfg.tau = .62; cfg.kT = kt;
            Solver serial(domain, cfg, 1., 1984), threaded(domain, cfg, 1., 1984);
            threaded.setThreadCount(threads);
            std::vector<double> fx(domain.fluidCount()), fy(fx.size());
            for (int n=0; n<domain.fluidCount(); ++n) {
                fx[n]=1e-7*std::sin(n); fy[n]=1e-7*std::cos(n);
                const double u=.001*std::sin(n), v=.001*std::cos(n);
                serial.initializeNode(n,1.,u,v); threaded.initializeNode(n,1.,u,v);
            }
            for (Solver* s : {&serial, &threaded}) {
                s->setForceField(&fx,&fy); s->setUniformForce(1e-8,-2e-8);
                for (size_t li=0; li<domain.links().size(); ++li) s->setWallVelocity(li,0.,-.003);
            }
            bool same=true;
            for(int step=0; step<31; ++step) {
                serial.step(); threaded.step();
                same &= bits(serial.maxTauEffective(),threaded.maxTauEffective());
                for(int n=0;n<domain.fluidCount();++n) {
                    same &= bits(serial.rho()[n],threaded.rho()[n]);
                    same &= bits(serial.ux()[n],threaded.ux()[n]);
                    same &= bits(serial.uy()[n],threaded.uy()[n]);
                    for(int i=0;i<Q;++i) same &= bits(serial.population(n,i),threaded.population(n,i));
                }
                for(int d=0;d<2;++d) same &= bits(serial.bodyForce(1)[d],threaded.bodyForce(1)[d]);
                if(step==15) threaded.setThreadCount(threads==2 ? 6 : 2);
            }
            check(same,"thread parity geometry="+std::to_string(geometry)+" thermal="+std::to_string(kt)+" threads="+std::to_string(threads));
        }
    }
    RiverParams p; p.rows=112; p.nestedBaseRows=112; p.tracers.enabled=false;
    RiverSimulation coarse(homeBankSegments(),1440,900,p);
    for(int r:{2,4}) {
        p.rows=112*r; RiverSimulation fine(homeBankSegments(),1440,900,p);
        check(fine.units().nx==r*coarse.units().nx && bits(fine.units().x0,coarse.units().x0),"nested dimensions and origin");
        checkNear(fine.units().h*r,coarse.units().h,0.,"nested spacing");
        checkNear((fine.units().tau-.5)/r,coarse.units().tau-.5,1e-16,"physical viscosity scaling");
    }
    return finish();
}
