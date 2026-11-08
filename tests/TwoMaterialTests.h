#ifndef FIELDLAB_TWO_MATERIAL_TESTS_H
#define FIELDLAB_TWO_MATERIAL_TESTS_H

#include "../opgnel_project_smoke/TwoMaterialSimulation.h"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <string>

namespace fieldlab {
namespace two_material_test_detail {
inline bool Near(double a, double b, double tolerance = 1e-9) {
    return std::isfinite(a) && std::abs(a-b) <= tolerance;
}
inline bool NearVector(Vec2 a, Vec2 b, double tolerance = 1e-9) {
    return Near(a.x,b.x,tolerance) && Near(a.y,b.y,tolerance);
}
inline Vec2 ActualWaterMomentum(const TwoMaterialSimulation& sim) {
    Vec2 value{};
    for(const auto& p:sim.Liquid().Particles()) value=value+p.velocity*sim.Liquid().ParticleMass();
    return value;
}
inline Vec2 ActualGasMomentum(const TwoMaterialSimulation& sim) {
    Vec2 value{};
    const auto& gas=sim.Gas(); const int n=gas.Size();
    const double mass=sim.GasDensity()*gas.CellWidth()*gas.CellWidth();
    for(int j=1;j<=n;++j)for(int i=1;i<n;++i)value.x+=mass*gas.MacFaceVelocity(0,i,j);
    for(int j=1;j<n;++j)for(int i=1;i<=n;++i)value.y+=mass*gas.MacFaceVelocity(1,i,j);
    return value;
}
inline double GasDifference(const TwoMaterialSimulation& a,const TwoMaterialSimulation& b) {
    double sum=0;
    for(int j=1;j<=a.Gas().Size();++j)for(int i=1;i<=a.Gas().Size();++i){
        const double du=a.Gas().U(i,j)-b.Gas().U(i,j),dv=a.Gas().V(i,j)-b.Gas().V(i,j);
        sum+=du*du+dv*dv;
    }
    return std::sqrt(sum/(a.Gas().Size()*a.Gas().Size()));
}
inline double PositionDifference(const TwoMaterialSimulation& a,const TwoMaterialSimulation& b){
    double sum=0;
    if(a.Liquid().Count()!=b.Liquid().Count())return std::numeric_limits<double>::infinity();
    for(std::size_t i=0;i<a.Liquid().Count();++i){
        const Vec2 d=a.Liquid().Particles()[i].position-b.Liquid().Particles()[i].position;
        sum+=d.x*d.x+d.y*d.y;
    }
    return a.Liquid().Count()?std::sqrt(sum/a.Liquid().Count()):0;
}
inline bool LiquidInside(const TwoMaterialSimulation& sim){
    const double lo=-2+sim.Liquid().ParticleRadius()-1e-8,hi=2-sim.Liquid().ParticleRadius()+1e-8;
    for(const auto&p:sim.Liquid().Particles())if(p.position.x<lo||p.position.x>hi||p.position.y<lo||p.position.y>hi)return false;
    return true;
}
inline bool WallsClosed(const TwoMaterialSimulation& sim){
    for(int i=0;i<=sim.Gas().Size();++i){
        const double q=double(i)/sim.Gas().Size();
        if(std::abs(sim.Gas().SampleVelocity({0,q}).x)>1e-10||std::abs(sim.Gas().SampleVelocity({1,q}).x)>1e-10
            ||std::abs(sim.Gas().SampleVelocity({q,0}).y)>1e-10||std::abs(sim.Gas().SampleVelocity({q,1}).y)>1e-10)return false;
    }
    return true;
}
inline TwoMaterialSimulation SingleParticle(int n=24,Vec2 p={0.13,-0.19},Vec2 v={0,0}){
    TwoMaterialSimulation sim(n);sim.SetJetEnabled(false);sim.Liquid().SetGravity({0,0});
    sim.Liquid().SetParticleState({{p,v}});sim.RebuildOccupancy();return sim;
}
}

inline int RunTwoMaterialTests(){
    using namespace two_material_test_detail;
    int passed=0,failed=0;
    const auto check=[&](const std::string& name,bool result){
        std::cout<<(result?"PASS ":"FAIL ")<<name<<'\n';result?++passed:++failed;
    };
    std::cout<<std::setprecision(12);
    TwoMaterialSimulation initial(24);
    check("coupled default uses 400 liquid particles",initial.Liquid().Count()==400);
    check("default educational density ratio is 12.5",Near(initial.LiquidDensity()/initial.GasDensity(),12.5));
    check("coupled initial state finite and particles inside",initial.IsFinite()&&LiquidInside(initial));
    bool fractionBounds=true;double phiMax=0;
    for(int j=1;j<=initial.Gas().Size();++j)for(int i=1;i<=initial.Gas().Size();++i){
        const double phi=initial.Occupancy(i,j);fractionBounds&=phi>=0&&phi<=1;phiMax=(std::max)(phiMax,phi);
    }
    check("poly6 occupancy is bounded and resolves a liquid interior",fractionBounds&&phiMax>0.85);

    auto equilibrium=SingleParticle();equilibrium.Gas().SetPresetVelocity(FieldType::Uniform,1.0);
    equilibrium.Liquid().SetParticleState({{{0.13,-0.19},{1.0,0.25}}});
    const double equilibriumEnergy=equilibrium.TotalKineticEnergy();equilibrium.ExchangeOnly(0.3);
    check("matching gas and water velocity produces no impulse",equilibrium.couplingWaterImpulse.Length()<1e-12&&equilibrium.couplingGasImpulse.Length()<1e-12);
    check("equilibrium exchange does not change kinetic energy",Near(equilibrium.TotalKineticEnergy(),equilibriumEnergy,1e-10));

    auto exact=SingleParticle(16,{0,-0.125});exact.SetDragRate(300);
    exact.Gas().AddMacFaceVelocity(0,8,8,2.0);
    const double mp=exact.Liquid().ParticleMass(),mg=exact.GasDensity()*0.25*0.25;
    const double common=2.0*mg/(mp+mg);
    exact.ExchangeOnly(1000.0);
    check("single-face stiff drag reaches exact common velocity",Near(exact.Liquid().Particles()[0].velocity.x,common,1e-12)
        &&Near(exact.Gas().MacFaceVelocity(0,8,8),common,1e-12));
    check("single-face impulse conserves independently measured momentum",NearVector(ActualWaterMomentum(exact)+ActualGasMomentum(exact),{2.0*mg,0},1e-10));

    auto random=SingleParticle(24);
    std::vector<WaterParticle> particles;
    unsigned int randomState=1297;
    const auto number=[&](){randomState=1664525u*randomState+1013904223u;return double(randomState)/4294967295.0;};
    for(int k=0;k<90;++k)particles.push_back({{-1.8+3.6*number(),-1.8+3.6*number()},{-3+6*number(),-3+6*number()}});
    random.Liquid().SetParticleState(particles);
    for(int component=0;component<2;++component)for(int j=1;j<=random.Gas().Size();++j)for(int i=1;i<=random.Gas().Size();++i)
        random.Gas().AddMacFaceVelocity(component,i,j,-4+8*number());
    const Vec2 rwBefore=ActualWaterMomentum(random),rgBefore=ActualGasMomentum(random);
    const double randomEnergy=random.TotalKineticEnergy();random.ExchangeOnly(0.2);
    const Vec2 rwDelta=ActualWaterMomentum(random)-rwBefore,rgDelta=ActualGasMomentum(random)-rgBefore;
    const double randomResidual=(rwDelta+rgDelta).Length();
    check("nonuniform multi-particle exchange changes both media",rwDelta.Length()>1e-5&&rgDelta.Length()>1e-5);
    check("random-state actual-array exchange conserves momentum",randomResidual<1e-8);
    check("reported impulses equal independent measured state changes",NearVector(random.couplingWaterImpulse,rwDelta,1e-11)&&NearVector(random.couplingGasImpulse,rgDelta,1e-11));
    check("two-way drag cannot increase kinetic energy",random.TotalKineticEnergy()<=randomEnergy+1e-8);
    std::cout<<"Coupling random-state exchange residual: "<<randomResidual<<"; energy "<<randomEnergy<<" -> "<<random.TotalKineticEnergy()<<".\n";

    auto walls=SingleParticle(32);walls.Liquid().SetParticleState({{{-1.967,-1.967},{2,-1}},{{1.967,1.967},{-1,3}}});
    walls.Gas().SetPresetVelocity(FieldType::Rotation,0.8);
    bool normalized=true;
    for(const auto& p:walls.Liquid().Particles())for(int c=0;c<2;++c){
        const auto s=walls.Gas().InteriorVelocityStencil((p.position+Vec2{2,2})/4.0,c);double sum=0;
        for(const auto& f:s){sum+=f.weight;normalized&=walls.Gas().IsInteriorFace(c,f.i,f.j);}
        normalized&=!s.empty()&&Near(sum,1.0);
    }
    check("wall-truncated MAC stencils exclude wall DOFs and normalize",normalized);
    const Vec2 wallMomentum=ActualWaterMomentum(walls)+ActualGasMomentum(walls);const double wallEnergy=walls.TotalKineticEnergy();
    walls.ExchangeOnly(0.5);
    check("near-wall exchange preserves actual momentum",NearVector(ActualWaterMomentum(walls)+ActualGasMomentum(walls),wallMomentum,1e-8));
    check("near-wall exchange dissipates energy and retains closed walls",walls.TotalKineticEnergy()<=wallEnergy+1e-8&&WallsClosed(walls));

    auto off=SingleParticle();off.Gas().SetPresetVelocity(FieldType::Uniform,2);off.SetCouplingMode(CouplingMode::Off);
    const Vec2 offP=off.TotalMomentum();const double offE=off.TotalKineticEnergy();off.ExchangeOnly(1);
    check("off mode leaves both velocities unchanged",NearVector(off.TotalMomentum(),offP)&&Near(off.TotalKineticEnergy(),offE));
    check("off mode reports zero exchanged impulse",off.couplingWaterImpulse.Length()==0&&off.couplingGasImpulse.Length()==0);
    auto one=SingleParticle();one.Gas().SetPresetVelocity(FieldType::Uniform,2);one.SetCouplingMode(CouplingMode::OneWay);
    const Vec2 oneGas=ActualGasMomentum(one);one.ExchangeOnly(0.1);
    check("one-way mode transfers gas velocity to the liquid",one.Liquid().Particles()[0].velocity.Length()>0.1);
    check("one-way reservoir intentionally has no gas reaction",NearVector(oneGas,ActualGasMomentum(one))&&one.couplingGasImpulse.Length()==0);
    check("one-way diagnostic exposes missing momentum reaction",one.couplingResidual.Length()>0.01);

    auto weak=SingleParticle();weak.Gas().SetPresetVelocity(FieldType::Uniform,1);
    auto strong=weak;weak.SetGasDensity(1.2);weak.ExchangeOnly(0.1);strong.ExchangeOnly(0.1);
    check("physical-like high density ratio exchange remains finite and conservative",weak.IsFinite()&&weak.couplingResidual.Length()<1e-8);
    check("smaller gas inertia gives a smaller liquid impulse",weak.couplingWaterImpulse.Length()<strong.couplingWaterImpulse.Length());
    auto zero=SingleParticle();zero.Gas().SetPresetVelocity(FieldType::Uniform,1);zero.SetDragRate(0);const Vec2 zeroP=zero.TotalMomentum();zero.ExchangeOnly(0.1);
    check("zero drag rate disables exchange",NearVector(zero.TotalMomentum(),zeroP)&&zero.couplingWaterImpulse.Length()==0);
    auto timeOne=SingleParticle();timeOne.Gas().SetPresetVelocity(FieldType::Uniform,1);auto timeTwo=timeOne;
    timeOne.ExchangeOnly(0.02);timeTwo.ExchangeOnly(0.01);timeTwo.ExchangeOnly(0.01);
    check("one particle exponential exchange composes over two half steps",NearVector(timeOne.Liquid().Particles()[0].velocity,timeTwo.Liquid().Particles()[0].velocity,1e-10)
        &&GasDifference(timeOne,timeTwo)<1e-10);
    auto invalid=SingleParticle();invalid.Gas().SetPresetVelocity(FieldType::Uniform,1);const Vec2 invalidP=invalid.TotalMomentum();
    invalid.ExchangeOnly(0);invalid.ExchangeOnly(-1);invalid.ExchangeOnly(std::numeric_limits<double>::quiet_NaN());
    check("invalid exchange durations do not change actual state",NearVector(invalid.TotalMomentum(),invalidP)&&invalid.IsFinite());
    bool densityRejected=false,dragRejected=false;
    try{invalid.SetGasDensity(0);}catch(const std::invalid_argument&){densityRejected=true;}
    try{invalid.SetDragRate(-1);}catch(const std::invalid_argument&){dragRejected=true;}
    check("invalid material density and drag are rejected",densityRejected&&dragRejected);
    auto empty=SingleParticle();empty.Liquid().SetParticleState({});empty.RebuildOccupancy();empty.ExchangeOnly(0.1);
    check("empty liquid has zero occupancy and exchange",Near(empty.Occupancy(10,10),0)&&empty.couplingResidual.Length()==0&&empty.IsFinite());

    TwoMaterialSimulation marker(24);marker.SetJetEnabled(false);marker.SetCouplingMode(CouplingMode::Off);marker.Gas().FillDye(1);marker.Step(1.0/120.0);
    bool denseExcluded=true,airRetained=false;int denseCells=0;
    for(int j=1;j<=marker.Gas().Size();++j)for(int i=1;i<=marker.Gas().Size();++i){
        if(marker.Occupancy(i,j)>0.85){denseExcluded&=Near(marker.Gas().Dye(i,j),0);++denseCells;}
        if(marker.Occupancy(i,j)<0.01&&marker.Gas().Dye(i,j)>0.9)airRetained=true;
    }
    check("liquid occupancy removes tracer from dense liquid cells",denseCells>0&&denseExcluded&&marker.removedDyeTracer>0);
    check("dye exclusion retains tracer in the gas region",airRetained);

    TwoMaterialSimulation decoupled(24),oneWay(24),twoWay(24);
    decoupled.SetCouplingMode(CouplingMode::Off);oneWay.SetCouplingMode(CouplingMode::OneWay);
    for(int k=0;k<360;++k){decoupled.Step(1.0/120);oneWay.Step(1.0/120);twoWay.Step(1.0/120);}
    const double offOneGas=GasDifference(decoupled,oneWay),offTwoGas=GasDifference(decoupled,twoWay);
    const double offOneWater=PositionDifference(decoupled,oneWay),offTwoWater=PositionDifference(decoupled,twoWay);
    check("one-way and off have identical gas velocity evolution",offOneGas<1e-12);
    check("coupling to driven gas changes the liquid trajectory",offOneWater>0.001);
    check("two-way coupling changes water trajectories",offTwoWater>0.001);
    check("moving liquid changes the gas velocity in two-way mode",offTwoGas>0.001);
    check("three full evolution modes remain finite",decoupled.IsFinite()&&oneWay.IsFinite()&&twoWay.IsFinite());
    std::cout<<"Coupled mode comparison at 3s (24x24): water position RMS off/one="<<offOneWater
        <<", off/two="<<offTwoWater<<"; gas velocity RMS off/one="<<offOneGas<<", off/two="<<offTwoGas<<".\n";
    TwoMaterialSimulation noJet(24);noJet.SetCouplingMode(CouplingMode::OneWay);noJet.SetJetEnabled(false);
    for(int k=0;k<360;++k)noJet.Step(1.0/120);
    const double jetEffect=PositionDifference(oneWay,noJet);
    check("matched jet on and off isolates the source effect on water",jetEffect>0.001);
    std::cout<<"Matched OneWay gas jet on/off water position RMS at 3s (24x24): "<<jetEffect<<".\n";

    TwoMaterialSimulation longRun(24);const double initialMass=longRun.LiquidMass();const auto start=std::chrono::steady_clock::now();
    bool allFinite=true,allInside=true,closed=true;double maxResidual=0;
    for(int k=0;k<600;++k){
        if(k==240)longRun.Liquid().Stir({0.4,-1.3},{0.12,0.08});
        longRun.Step(1.0/120);allFinite&=longRun.IsFinite();allInside&=LiquidInside(longRun);closed&=WallsClosed(longRun);
        maxResidual=(std::max)(maxResidual,longRun.couplingResidual.Length());
    }
    check("600 coupled steps remain finite",allFinite);
    check("600 coupled steps preserve water particle mass and count",longRun.Liquid().Count()==400&&Near(longRun.LiquidMass(),initialMass,1e-10));
    check("600 coupled steps keep particles inside and gas walls closed",allInside&&closed);
    check("two-way full-run exchange residual remains below 1e-7",maxResidual<1e-7);
    check("fixed-step coupled time reaches exactly five seconds",Near(longRun.Time(),5.0,1e-10));
    std::cout<<"Coupled 600-step scenario (24x24): "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()
        <<" s; max exchange residual "<<maxResidual<<"; water mass "<<longRun.LiquidMass()
        <<"; final gas divergence "<<longRun.Gas().DivergenceRms()<<".\n";
    TwoMaterialSimulation physical(16);physical.SetGasDensity(1.2);bool physicalFinite=true;
    for(int k=0;k<240;++k){physical.Step(1.0/120);physicalFinite&=physical.IsFinite();}
    check("water to room-air-like density ratio stays stable in full evolution",physicalFinite&&physical.maxCouplingResidual<1e-7);
    std::cout<<"TwoMaterial tests: "<<passed<<" passed, "<<failed<<" failed, "<<(passed+failed)<<" total.\n";
    return failed?1:0;
}
}
#endif
