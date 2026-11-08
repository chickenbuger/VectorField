#include "../../opgnel_project_smoke/SharpInterfaceSimulation.h"
#include <chrono>
#include <cstdint>
#include <functional>
#include <iomanip>
#include <iostream>
#include <string>

using namespace fieldlab;
namespace {
constexpr double Pi = 3.14159265358979323846;
constexpr double Dt = 1.0 / 120.0;
struct Review {
    int checks=0, failed=0;
    void Check(bool pass,const std::string& name) {
        ++checks;
        if(!pass) { ++failed; std::cout<<"FAIL "<<name<<'\n'; }
    }
};
Vec2 Center(int i,int j,double h) { return {-2+(i-.5)*h,-2+(j-.5)*h}; }
double Beta(const SharpInterfaceSimulation& sim,int component,int i,int j) {
    const double a=sim.Phi(i,j);
    const double b=sim.Phi(i+(component==0),j+(component==1));
    const double ra=a<0?sim.LiquidDensity():sim.GasDensity();
    const double rb=b<0?sim.LiquidDensity():sim.GasDensity();
    if((a<0)==(b<0)) return 1/ra;
    const double theta=std::abs(a)/(std::abs(a)+std::abs(b));
    return 1/(theta*ra+(1-theta)*rb);
}
double Divergence(const SharpInterfaceSimulation& sim) {
    double squares=0;
    for(int j=1;j<=sim.Size();++j)for(int i=1;i<=sim.Size();++i) {
        double d=(sim.FaceVelocity(0,i,j)-sim.FaceVelocity(0,i-1,j)
            +sim.FaceVelocity(1,i,j)-sim.FaceVelocity(1,i,j-1))/sim.CellWidth();
        squares+=d*d;
    }
    return std::sqrt(squares/(sim.Size()*sim.Size()));
}
double Energy(const SharpInterfaceSimulation& sim) {
    double result=0,h=sim.CellWidth();
    for(int c=0;c<2;++c)for(int j=1;j<=sim.Size()-(c==1);++j)
        for(int i=1;i<=sim.Size()-(c==0);++i) {
            double u=sim.FaceVelocity(c,i,j);
            result+=.5*h*h/Beta(sim,c,i,j)*u*u;
        }
    return result;
}
double SpeedMaximum(const SharpInterfaceSimulation& sim) {
    double result=0;
    for(int c=0;c<2;++c)for(int j=0;j<=sim.Size();++j)for(int i=0;i<=sim.Size();++i)
        result=(std::max)(result,std::abs(sim.FaceVelocity(c,i,j)));
    return result;
}
bool WallsClosed(const SharpInterfaceSimulation& sim) {
    for(int i=1;i<=sim.Size();++i) {
        if(sim.FaceVelocity(0,0,i)!=0||sim.FaceVelocity(0,sim.Size(),i)!=0
            ||sim.FaceVelocity(1,i,0)!=0||sim.FaceVelocity(1,i,sim.Size())!=0)return false;
    }
    return true;
}
double Area(const SharpInterfaceSimulation& sim) {
    double sum=0;
    for(int j=1;j<=sim.Size();++j)for(int i=1;i<=sim.Size();++i)sum+=sim.LiquidFraction(i,j);
    return sum*sim.CellWidth()*sim.CellWidth();
}
double PressureError(const SharpInterfaceSimulation& sim,const std::function<double(Vec2)>& exact,double dt) {
    double offset=0;
    for(int j=1;j<=sim.Size();++j)for(int i=1;i<=sim.Size();++i)
        offset+=sim.Pressure(i,j)*dt-exact(Center(i,j,sim.CellWidth()));
    offset/=sim.Size()*sim.Size();
    double result=0;
    for(int j=1;j<=sim.Size();++j)for(int i=1;i<=sim.Size();++i)
        result=(std::max)(result,std::abs(sim.Pressure(i,j)*dt-exact(Center(i,j,sim.CellWidth()))-offset));
    return result;
}
std::uint64_t rng=91249811;
double Random() {rng=rng*6364136223846793005ULL+1442695040888963407ULL;return double(rng>>11)/9007199254740992.0;}
}

int main() {
    Review review;
    std::cout<<std::setprecision(17);
    int randomCases=0,crossingFaces=0;
    double maxBetaRelativeError=0,maxThetaError=0,maxRelativeEnergyIncrease=0,maxProjectionDivRatio=0,maxTransmissionFluxError=0;
    double maxReportedDivDifference=0,maxManufacturedVelocity=0,maxManufacturedPressureError=0;
    for(int n:{8,16,32,48})for(double density:{1.2,20.,80.,1000.})for(double angle:{0.,.5410520681182421,Pi/2}) {
        SharpInterfaceSimulation sim(n);
        sim.SetGasDensity(density);sim.SetSubcellInterface(true);
        sim.SetPhi([=](Vec2 p){return std::cos(angle)*p.x+std::sin(angle)*p.y+.137;});
        sim.SetVelocity([](Vec2){return Vec2{};});
        for(int c=0;c<2;++c)for(int j=1;j<=n-(c==1);++j)for(int i=1;i<=n-(c==0);++i) {
            const double expected=Beta(sim,c,i,j),actual=sim.FaceBeta(c,i,j);
            maxBetaRelativeError=(std::max)(maxBetaRelativeError,std::abs(actual-expected)/expected);
            const double a=sim.Phi(i,j),b=sim.Phi(i+(c==0),j+(c==1));
            if((a<0)!=(b<0)) {
                ++crossingFaces;
                maxThetaError=(std::max)(maxThetaError,std::abs(sim.FaceTheta(c,i,j)-std::abs(a)/(std::abs(a)+std::abs(b))));
                const double pa=(a<0?1000:density)*a,pb=(b<0?1000:density)*b;
                const double actualFlux=actual*(pb-pa)/sim.CellWidth();
                const double exactFlux=c==0?std::cos(angle):std::sin(angle);
                maxTransmissionFluxError=(std::max)(maxTransmissionFluxError,std::abs(actualFlux-exactFlux));
            }
            sim.AddFaceVelocity(c,i,j,2*Random()-1);
        }
        const double beforeEnergy=Energy(sim),beforeDiv=Divergence(sim);
        sim.ProjectOnly(Dt);
        const double afterEnergy=Energy(sim),afterDiv=Divergence(sim);
        const double divRatio=afterDiv/(beforeDiv+1e-30);
        maxProjectionDivRatio=(std::max)(maxProjectionDivRatio,divRatio);
        maxRelativeEnergyIncrease=(std::max)(maxRelativeEnergyIncrease,(afterEnergy-beforeEnergy)/(1+beforeEnergy));
        maxReportedDivDifference=(std::max)(maxReportedDivDifference,std::abs(afterDiv-sim.divergenceAfterProjection));
        review.Check(sim.IsFinite()&&WallsClosed(sim)&&divRatio<2e-7&&afterEnergy<=beforeEnergy*(1+1e-10),"random variable-density projection");
        ++randomCases;
    }
    review.Check(crossingFaces>0&&maxBetaRelativeError<5e-13&&maxThetaError<5e-13,"independent crossing fraction and density resistance");
    review.Check(maxTransmissionFluxError<2e-12,"continuous piecewise-linear pressure transmits the exact normal flux");
    review.Check(maxReportedDivDifference<1e-12,"reported divergence matches independently summed MAC arrays");
    std::cout<<"Independent random projections: "<<randomCases<<"; checked crossing faces: "<<crossingFaces<<'\n'
        <<"Maximum relative beta error: "<<maxBetaRelativeError<<"; theta absolute error: "<<maxThetaError<<'\n'
        <<"Maximum piecewise-linear continuous-pressure transmission flux error: "<<maxTransmissionFluxError<<'\n'
        <<"Maximum post/pre divergence ratio: "<<maxProjectionDivRatio<<"; relative kinetic energy increase: "<<maxRelativeEnergyIncrease<<'\n'
        <<"Maximum reported versus actual divergence RMS difference: "<<maxReportedDivDifference<<'\n';

    SharpInterfaceSimulation endpoint(16);
    endpoint.SetPhi([](Vec2 p){return -.125-p.x;});
    review.Check(std::abs(endpoint.FaceTheta(0,8,8))<1e-14&&std::abs(endpoint.FaceBeta(0,8,8)-.001)<1e-14,
        "zero level set at left cell center has theta zero and full water segment");
    endpoint.SetPhi([](Vec2 p){return p.x+.125;});
    review.Check(std::abs(endpoint.FaceTheta(0,7,8)-1)<1e-14&&std::abs(endpoint.FaceBeta(0,7,8)-.001)<1e-14,
        "zero level set at right cell center has theta one and full water segment");
    std::cout<<"Cell-center interface endpoint cases: 2.\n";

    for(int n:{16,32,48})for(double density:{1.2,80.}) {
        SharpInterfaceSimulation sim(n);sim.SetGasDensity(density);sim.SetSubcellInterface(true);
        sim.SetPhi([](Vec2 p){return p.x+.41*p.y+.137;});
        sim.SetVelocity([](Vec2){return Vec2{};});
        const auto q=[](Vec2 p){return .03*std::cos(Pi*(p.x+2)/4)*std::cos(2*Pi*(p.y+2)/4);};
        for(int c=0;c<2;++c)for(int j=1;j<=n-(c==1);++j)for(int i=1;i<=n-(c==0);++i) {
            Vec2 a=Center(i,j,sim.CellWidth()),b=Center(i+(c==0),j+(c==1),sim.CellWidth());
            sim.AddFaceVelocity(c,i,j,Beta(sim,c,i,j)*(q(b)-q(a))/sim.CellWidth());
        }
        sim.ProjectOnly(Dt);
        const double speed=SpeedMaximum(sim),pError=PressureError(sim,q,Dt);
        maxManufacturedVelocity=(std::max)(maxManufacturedVelocity,speed);
        maxManufacturedPressureError=(std::max)(maxManufacturedPressureError,pError);
        review.Check(speed<2e-9&&pError<2e-7,"manufactured pressure impulse recovered up to gauge");
    }
    std::cout<<"Manufactured pressure cases: 6; max postprojection speed: "<<maxManufacturedVelocity
        <<"; max pressure-impulse gauge error: "<<maxManufacturedPressureError<<'\n';

    double maxHydrostaticSpeed=0,maxHydrostaticPressureError=0;
    for(int n:{16,32,48})for(double density:{1.2,80.}) {
        SharpInterfaceSimulation sim(n);sim.SetGasDensity(density);sim.SetSubcellInterface(true);
        constexpr double surface=-.137,g=-9.81;
        sim.SetPhi([=](Vec2 p){return p.y-surface;});
        sim.SetVelocity([=](Vec2){return Vec2{0,g*Dt};});
        sim.ProjectOnly(Dt);
        const auto hydro=[=](Vec2 p){return Dt*g*(p.y<surface?1000.:density)*(p.y-surface);};
        const double speed=SpeedMaximum(sim),pError=PressureError(sim,hydro,Dt);
        maxHydrostaticSpeed=(std::max)(maxHydrostaticSpeed,speed);
        maxHydrostaticPressureError=(std::max)(maxHydrostaticPressureError,pError);
        review.Check(speed<2e-7&&pError<2e-5,"subcell flat-interface hydrostatic balance");
    }
    std::cout<<"Hydrostatic cases: 6; max speed: "<<maxHydrostaticSpeed
        <<"; max pressure-impulse gauge error: "<<maxHydrostaticPressureError<<'\n';

    for(int n:{16,32,48}) {
        double exactError=0,staircaseError=0;
        for(bool subcell:{false,true}) {
            SharpInterfaceSimulation sim(n);sim.SetGasDensity(1.2);sim.SetSubcellInterface(subcell);
            constexpr double surface=-.137,g=-9.81;
            sim.SetPhi([=](Vec2 p){return p.y-surface;});
            sim.SetVelocity([=](Vec2){return Vec2{0,g*Dt};});
            sim.ProjectOnly(Dt);
            const auto hydro=[=](Vec2 p){return Dt*g*(p.y<surface?1000.:1.2)*(p.y-surface);};
            const double error=PressureError(sim,hydro,Dt);
            if(subcell)exactError=error;else staircaseError=error;
        }
        review.Check(exactError<staircaseError*.01,"subcell hydrostatic pressure more accurate than midpoint density");
        std::cout<<"Matched hydrostatic pressure n="<<n<<" rho_g=1.2: midpoint error="<<staircaseError
            <<"; subcell error="<<exactError<<'\n';
    }

    for(double density:{80.,1.2}) {
        SharpInterfaceSimulation sim(48);sim.SetGasDensity(density);
        const double initialArea=Area(sim);
        bool finite=true,closed=true,bounded=true,converged=true;
        double maximumDiv=0,maximumAreaRelativeError=0,maximumSpeed=0,maximumRawAreaError=0,maximumAreaShift=0;
        auto started=std::chrono::steady_clock::now();
        for(int step=0;step<600;++step) {
            sim.Step(Dt);finite&=sim.IsFinite();closed&=WallsClosed(sim);converged&=sim.pressureConverged;
            maximumDiv=(std::max)(maximumDiv,Divergence(sim));
            maximumAreaRelativeError=(std::max)(maximumAreaRelativeError,std::abs(Area(sim)-initialArea)/(initialArea+1e-30));
            maximumSpeed=(std::max)(maximumSpeed,SpeedMaximum(sim));
            maximumRawAreaError=(std::max)(maximumRawAreaError,std::abs(sim.RawAreaError()));
            maximumAreaShift=(std::max)(maximumAreaShift,std::abs(sim.AreaCorrectionShift()));
            for(int j=1;j<=sim.Size();++j)for(int i=1;i<=sim.Size();++i) {
                const double fraction=sim.LiquidFraction(i,j);
                bounded&=std::isfinite(fraction)&&fraction>=0&&fraction<=1;
            }
        }
        const double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
        review.Check(finite&&closed&&bounded&&converged&&maximumDiv<1e-7&&maximumAreaRelativeError<1e-6&&std::abs(sim.Time()-5)<1e-9,"default 48x48 full evolution at density "+std::to_string(density));
        std::cout<<"Full 48x48 evolution rho_g="<<density<<": 600 x 1/120 s; elapsed "<<elapsed<<" s; finite="<<finite
            <<" closed_walls="<<closed<<" bounded_fraction="<<bounded<<" all_final_substep_pcg_converged="<<converged<<"; max divergence="<<maximumDiv
            <<"; max relative corrected reconstructed water-area drift="<<maximumAreaRelativeError<<"; max face speed="<<maximumSpeed
            <<"; max absolute pre-correction relative area error="<<maximumRawAreaError<<"; max absolute correction shift="<<maximumAreaShift<<'\n';
    }
    std::cout<<"Independent interface review: "<<review.checks<<" checks, "<<review.failed<<" failures.\n";
    return review.failed?1:0;
}
