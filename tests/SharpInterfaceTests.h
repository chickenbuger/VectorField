#ifndef FIELDLAB_SHARP_INTERFACE_TESTS_H
#define FIELDLAB_SHARP_INTERFACE_TESTS_H
#include "../opgnel_project_smoke/SharpInterfaceSimulation.h"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <string>

namespace fieldlab {
namespace sharp_interface_test_detail {
inline double ActualDivergence(const SharpInterfaceSimulation& a) {
    double sum=0;const double h=a.CellWidth();
    for(int j=1;j<=a.Size();++j)for(int i=1;i<=a.Size();++i){const double d=(a.FaceVelocity(0,i,j)-a.FaceVelocity(0,i-1,j)+a.FaceVelocity(1,i,j)-a.FaceVelocity(1,i,j-1))/h;sum+=d*d;}
    return std::sqrt(sum/(a.Size()*a.Size()));
}
inline bool Closed(const SharpInterfaceSimulation& a) {
    for(int i=1;i<=a.Size();++i)if(a.FaceVelocity(0,0,i)!=0||a.FaceVelocity(0,a.Size(),i)!=0||a.FaceVelocity(1,i,0)!=0||a.FaceVelocity(1,i,a.Size())!=0)return false;
    return true;
}
inline double ActualEnergy(const SharpInterfaceSimulation& a) {
    double result=0;const double h=a.CellWidth();
    for(int c=0;c<2;++c)for(int j=1;j<=a.Size()-(c==1);++j)for(int i=1;i<=a.Size()-(c==0);++i){const double v=a.FaceVelocity(c,i,j);result+=0.5*h*h*v*v/a.FaceBeta(c,i,j);}
    return result;
}
inline double PressureError(const SharpInterfaceSimulation& a,const std::function<double(Vec2)>& exact) {
    double gauge=0;const int n=a.Size();const double h=a.CellWidth();
    for(int j=1;j<=n;++j)for(int i=1;i<=n;++i)gauge+=a.Pressure(i,j)-exact({-2+(i-0.5)*h,-2+(j-0.5)*h});
    gauge/=n*n;double e=0;
    for(int j=1;j<=n;++j)for(int i=1;i<=n;++i){const double d=a.Pressure(i,j)-exact({-2+(i-0.5)*h,-2+(j-0.5)*h})-gauge;e+=d*d;}
    return std::sqrt(e/(n*n));
}
inline void Hydrostatic(SharpInterfaceSimulation& a,double y0,double dt=1.0/120.0) {
    a.SetPhi([=](Vec2 p){return p.y-y0;});a.SetJetEnabled(false);
    a.SetVelocity([=](Vec2){return Vec2{0,-9.81*dt};});a.ProjectOnly(dt);
}
}

inline int RunSharpInterfaceTests() {
    using namespace sharp_interface_test_detail;
    int passed=0,failed=0;
    const auto check=[&](const std::string& name,bool ok){std::cout<<(ok?"PASS ":"FAIL ")<<name<<'\n';ok?++passed:++failed;};
    const double dt=1.0/120.0,pi=3.14159265358979323846;
    SharpInterfaceSimulation a(24);
    check("sharp interface defaults are finite and use separate material densities",a.IsFinite()&&a.GasDensity()==80&&a.LiquidDensity()==1000);
    check("sharp initial water area lies inside shared domain",a.WaterArea()>0&&a.WaterArea()<16);
    check("shared MAC domain has closed normal walls",Closed(a));
    bool fractionBounded=true;
    for(int j=1;j<=a.Size();++j)for(int i=1;i<=a.Size();++i)fractionBounded=fractionBounded&&a.LiquidFraction(i,j)>=0&&a.LiquidFraction(i,j)<=1;
    check("material display fractions stay in zero to one",fractionBounded);
    a.SetPhi([](Vec2 p){return p.y+0.137;});
    const int row=static_cast<int>(std::floor((-0.137+2)/a.CellWidth()+0.5));
    const double left=a.Phi(12,row),right=a.Phi(12,row+1),theta=std::abs(left)/(std::abs(left)+std::abs(right));
    const double beta=1/(theta*1000+(1-theta)*80);
    check("SetPhi preserves a planar signed distance without reinitialization",std::abs(a.Phi(7,11)-(-2+10.5*a.CellWidth()+0.137))<1e-14);
    check("subcell crossing fraction matches linear zero location",std::abs(a.FaceTheta(1,12,row)-theta)<1e-14&&theta>0&&theta<1);
    check("interface face uses integrated density resistance",std::abs(a.FaceBeta(1,12,row)-beta)<1e-14);
    check("pure water face retains inverse water density",std::abs(a.FaceBeta(0,5,3)-0.001)<1e-14);
    check("pure gas face retains inverse gas density",std::abs(a.FaceBeta(0,5,22)-0.0125)<1e-14);
    a.SetSubcellInterface(false);
    check("midpoint comparison uses one half density lengths",std::abs(a.FaceBeta(1,12,row)-2.0/1080)<1e-14);
    a.SetSubcellInterface(true);
    const double qLeft=1000*left,qRight=80*right;
    const double sharpFlux=a.FaceBeta(1,12,row)*(qRight-qLeft)/a.CellWidth();
    a.SetSubcellInterface(false);const double midpointFlux=a.FaceBeta(1,12,row)*(qRight-qLeft)/a.CellWidth();
    check("subcell coefficient reproduces continuous piecewise-linear interface flux",std::abs(sharpFlux-1)<1e-12);
    check("subcell flux is more accurate than midpoint on off-center interface",std::abs(sharpFlux-1)<std::abs(midpointFlux-1)*1e-6&&std::abs(midpointFlux-1)>0.01);
    std::cout<<std::setprecision(12)<<"Sharp planar flux: exact=1, subcell="<<sharpFlux<<", midpoint="<<midpointFlux<<"; theta="<<theta<<".\n";

    SharpInterfaceSimulation endpoint(16);const double center=-2+7.5*endpoint.CellWidth();
    endpoint.SetPhi([=](Vec2 p){return center-p.y;});
    check("interface exactly on a cell center reports endpoint fraction zero",endpoint.FaceTheta(1,5,8)==0&&std::abs(endpoint.FaceBeta(1,5,8)-0.001)<1e-14);
    endpoint.SetPhi([=](Vec2 p){return p.y-center;});
    check("reverse cell-center endpoint reports fraction one",endpoint.FaceTheta(1,5,7)==1&&std::abs(endpoint.FaceBeta(1,5,7)-0.001)<1e-14);

    a.SetSubcellInterface(true);a.SetVelocity([](Vec2){return Vec2{};});a.ProjectOnly(dt);
    check("zero velocity produces zero pressure and no spurious flow",a.MaxSpeed()<1e-14&&a.pressureIterations==0&&a.pressureConverged);

    double maximumHydroSharp=0,maximumHydroMid=0,maximumHydroSpeed=0;
    bool hydroFine=true,hydroImproved=true;
    for(double gas:{1.2,80.0,1000.0}) {
        SharpInterfaceSimulation sharp(32),midpoint(32);sharp.SetGasDensity(gas);midpoint.SetGasDensity(gas);midpoint.SetSubcellInterface(false);
        Hydrostatic(sharp,-0.137);Hydrostatic(midpoint,-0.137);
        const auto exact=[=](Vec2 p){const double y=p.y+0.137;return -9.81*(y<0?1000:gas)*y;};
        const double e=PressureError(sharp,exact),m=PressureError(midpoint,exact);
        maximumHydroSharp=(std::max)(maximumHydroSharp,e);maximumHydroMid=(std::max)(maximumHydroMid,m);maximumHydroSpeed=(std::max)(maximumHydroSpeed,sharp.MaxSpeed());
        hydroFine=hydroFine&&e<1e-5&&sharp.MaxSpeed()<1e-8&&sharp.pressureConverged;
        if(gas<1000)hydroImproved=hydroImproved&&e<m*1e-5&&m>1e-2;
    }
    check("hydrostatic pressure matches analytic two-density column",hydroFine);
    check("subcell hydrostatic pressure is more accurate than midpoint",hydroImproved);
    std::cout<<"Hydrostatic 32x32 gas density 1.2/80/1000: max sharp pressure RMS="<<maximumHydroSharp<<", max midpoint RMS="<<maximumHydroMid<<", max residual speed="<<maximumHydroSpeed<<".\n";

    SharpInterfaceSimulation manufactured(24);manufactured.SetGasDensity(1.2);
    manufactured.SetPhi([](Vec2 p){return p.x+0.37*p.y-0.11;});manufactured.SetVelocity([](Vec2){return Vec2{};});
    const auto q=[=](int i,int j){return 0.4*std::cos(pi*(i-0.5)/manufactured.Size())*std::cos(2*pi*(j-0.5)/manufactured.Size());};
    for(int j=1;j<=manufactured.Size();++j)for(int i=1;i<manufactured.Size();++i)manufactured.AddFaceVelocity(0,i,j,manufactured.FaceBeta(0,i,j)*(q(i+1,j)-q(i,j))/manufactured.CellWidth());
    for(int j=1;j<manufactured.Size();++j)for(int i=1;i<=manufactured.Size();++i)manufactured.AddFaceVelocity(1,i,j,manufactured.FaceBeta(1,i,j)*(q(i,j+1)-q(i,j))/manufactured.CellWidth());
    manufactured.ProjectOnly(dt);
    const double manufacturedError=PressureError(manufactured,[=](Vec2 p){return 0.4/dt*std::cos(pi*(p.x+2)/4)*std::cos(2*pi*(p.y+2)/4);});
    check("variable-density manufactured pressure is recovered up to gauge",manufacturedError<1e-5);
    check("manufactured gradient velocity projects to zero",manufactured.MaxSpeed()<1e-8);
    check("explicit pressure residual confirms PCG convergence",manufactured.pressureConverged&&manufactured.pressureResidual<1e-8);
    std::cout<<"Manufactured slanted-interface pressure RMS="<<manufacturedError<<", explicit relative residual="<<manufactured.pressureResidual<<"; iterations="<<manufactured.pressureIterations<<".\n";

    SharpInterfaceSimulation random(24);random.SetGasDensity(1.2);random.SetPhi([](Vec2 p){return std::hypot(p.x+0.13,p.y-0.17)-0.91;});
    random.SetVelocity([](Vec2 p){return Vec2{std::sin(3.7*p.x+1.4*p.y),std::cos(2.2*p.y-0.8*p.x)};});
    const double energyBefore=ActualEnergy(random),divBefore=ActualDivergence(random);random.ProjectOnly(dt);const double energyAfter=ActualEnergy(random),divAfter=ActualDivergence(random);
    check("projection removes independently summed divergence",divAfter<divBefore*1e-7);
    check("projection preserves prescribed closed wall velocities",Closed(random));
    check("density-weighted actual-face kinetic energy cannot increase",energyAfter<=energyBefore*(1+1e-10));
    check("public divergence agrees with independent MAC sum",std::abs(divAfter-random.DivergenceRms())<1e-13);
    check("public energy agrees with independently weighted face sum",std::abs(energyAfter-random.KineticEnergy())<1e-10);
    std::cout<<"Independent face projection: divergence "<<divBefore<<" -> "<<divAfter<<"; kinetic energy "<<energyBefore<<" -> "<<energyAfter<<".\n";

    SharpInterfaceSimulation evolution(48);const double initialArea=evolution.InitialWaterArea();bool finite=true,walls=true,converged=true;double maxArea=0,maxRaw=0,maxCfl=0,maxDiv=0;
    const auto start=std::chrono::steady_clock::now();
    for(int step=0;step<600;++step){if(step==240)evolution.Stir({-0.8,-0.6},{0.08,0.11});evolution.Step(dt);finite=finite&&evolution.IsFinite();walls=walls&&Closed(evolution);converged=converged&&evolution.pressureConverged;maxArea=(std::max)(maxArea,std::abs(evolution.RelativeAreaError()));maxRaw=(std::max)(maxRaw,std::abs(evolution.RawAreaError()));maxCfl=(std::max)(maxCfl,evolution.lastCfl);maxDiv=(std::max)(maxDiv,ActualDivergence(evolution));}
    const double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    check("default 48x48 complete 600-step evolution stays finite",finite);
    check("complete evolution keeps both phases inside closed domain",walls&&evolution.WaterArea()>0&&evolution.WaterArea()<16);
    check("all default pressure solves meet explicit residual tolerance",converged);
    check("global level-set area correction meets documented reference area",maxArea<1e-8&&std::abs(evolution.InitialWaterArea()-initialArea)<1e-12);
    check("default adaptive transport remains within estimated CFL target",maxCfl<=0.6500001);
    check("fixed-step interface simulation reaches five seconds",std::abs(evolution.Time()-5)<1e-9);
    std::cout<<"Sharp 48x48 600-step scenario: "<<elapsed<<" s; max corrected area relative error="<<maxArea<<", max raw reinit/advection area error="<<maxRaw<<", max estimated CFL="<<maxCfl<<", max divergence="<<maxDiv<<"; initial area="<<initialArea<<".\n";
    SharpInterfaceSimulation air(32);air.SetGasDensity(1.2);bool airFinite=true,airConverged=true;
    for(int step=0;step<240;++step){air.Step(dt);airFinite=airFinite&&air.IsFinite();airConverged=airConverged&&air.pressureConverged;}
    check("room-air-like density ratio stays finite in complete evolution",airFinite&&airConverged&&std::abs(air.RelativeAreaError())<1e-8);
    SharpInterfaceSimulation off(16);off.SetProjectionEnabled(false);off.Step(dt);
    check("projection-off comparison exposes uncorrected divergence",off.DivergenceRms()>1e-4&&off.pressureIterations==0&&!off.ProjectionEnabled());
    SharpInterfaceSimulation uniform(16);uniform.SetGasDensity(1000);bool sameDensity=true;
    for(int j=1;j<=16;++j)for(int i=1;i<16;++i)sameDensity=sameDensity&&std::abs(uniform.FaceBeta(0,i,j)-0.001)<1e-14;
    check("equal material densities reduce to constant-coefficient pressure",sameDensity);
    const double time=uniform.Time();uniform.Step(-1);uniform.Step(std::numeric_limits<double>::infinity());
    check("invalid time increments leave interface state unchanged",uniform.Time()==time&&uniform.IsFinite());
    bool rejected=false;try{uniform.SetGasDensity(0);}catch(const std::invalid_argument&){rejected=true;}
    check("invalid density is rejected before modifying state",rejected&&uniform.GasDensity()==1000);
    rejected=false;try{uniform.SetPhi([](Vec2){return std::numeric_limits<double>::quiet_NaN();});}catch(const std::invalid_argument&){rejected=true;}
    check("invalid level-set injection is atomic",rejected&&uniform.IsFinite());
    SharpInterfaceSimulation pure(16);pure.SetGravity({0,0});pure.SetJetEnabled(false);pure.SetPhi([](Vec2){return -1.0;});pure.Step(dt);
    check("pure water retains full area without artificial offset",pure.WaterArea()==16&&pure.AreaCorrectionShift()==0&&pure.IsFinite());
    pure.SetPhi([](Vec2){return 1.0;});pure.Step(dt);
    check("pure gas retains empty water area without artificial offset",pure.WaterArea()==0&&pure.AreaCorrectionShift()==0&&pure.IsFinite());
    SharpInterfaceSimulation alternating(16);const double cell=alternating.CellWidth();
    alternating.SetVelocity([=](Vec2 p){return Vec2{4*std::cos(pi*(p.x+2)/cell),0};});
    check("CFL bound sees alternating face speeds hidden by cell averaging",alternating.FaceSpeedBound()>3.999999&&std::abs(alternating.U(8,8))<1e-10);
    // Scratch reuse must not carry an old pressure correction into equilibrium.
    SharpInterfaceSimulation reused(16);
    reused.SetVelocity([](Vec2 p){return Vec2{std::sin(2.3*p.x),std::cos(1.7*p.y)};});
    reused.ProjectOnly(dt);reused.SetVelocity([](Vec2){return Vec2{};});reused.ProjectOnly(dt);
    double restPressure=0;
    for(int j=1;j<=16;++j)for(int i=1;i<=16;++i)restPressure=(std::max)(restPressure,std::abs(reused.Pressure(i,j)));
    check("reused pressure workspace clears a previous nonzero correction",reused.IsFinite()&&reused.pressureConverged&&reused.MaxSpeed()<1e-13&&restPressure<1e-13);
    SharpInterfaceSimulation contrast(16);contrast.SetLiquidDensity(5000);contrast.SetGasDensity(0.1);
    contrast.SetPhi([](Vec2 p){return std::hypot(p.x+0.19,p.y-0.21)-0.83;});
    contrast.SetVelocity([](Vec2 p){return Vec2{std::sin(2.1*p.x+0.4*p.y),std::cos(1.9*p.y-0.2*p.x)};});
    const double contrastBefore=ActualDivergence(contrast),contrastEnergy=ActualEnergy(contrast);contrast.ProjectOnly(dt);
    check("MIC pressure solve handles the supported extreme density contrast",contrast.IsFinite()&&contrast.pressureConverged&&Closed(contrast)&&ActualDivergence(contrast)<contrastBefore*1e-7&&ActualEnergy(contrast)<=contrastEnergy*(1+1e-10));
    std::cout<<"SharpInterface tests: "<<passed<<" passed, "<<failed<<" failed, "<<passed+failed<<" total.\n";
    return failed?1:0;
}
} // namespace fieldlab
#endif
