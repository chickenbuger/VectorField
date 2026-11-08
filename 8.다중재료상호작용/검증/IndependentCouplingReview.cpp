#include "../../opgnel_project_smoke/TwoMaterialSimulation.h"
#include <iomanip>
#include <iostream>
#include <cstdint>
#include <chrono>

using namespace fieldlab;
struct Totals { Vec2 gas, liquid; double energy=0; };
Totals ReadActual(const TwoMaterialSimulation& sim) {
    Totals result;
    const auto& gas=sim.Gas();
    const double mf=sim.GasDensity()*gas.CellWidth()*gas.CellWidth();
    for(int component=0;component<2;++component) {
        for(int j=1;j<=gas.Size();++j) for(int i=1;i<=gas.Size();++i) {
            if(!gas.IsInteriorFace(component,i,j))continue;
            double v=gas.MacFaceVelocity(component,i,j);
            if(component==0)result.gas.x+=mf*v;else result.gas.y+=mf*v;
            result.energy+=0.5*mf*v*v;
        }
    }
    const double mp=sim.Liquid().ParticleMass();
    for(const auto& particle:sim.Liquid().Particles()) {
        result.liquid=result.liquid+particle.velocity*mp;
        result.energy+=0.5*mp*(particle.velocity.x*particle.velocity.x+particle.velocity.y*particle.velocity.y);
    }
    return result;
}
std::uint64_t rng=9471337;
double Random() {rng=rng*6364136223846793005ULL+1442695040888963407ULL;return static_cast<double>(rng>>11)/9007199254740992.0;}
void Fill(TwoMaterialSimulation& sim) {
    auto& gas=sim.Gas();
    for(int c=0;c<2;++c)for(int j=1;j<=gas.Size();++j)for(int i=1;i<=gas.Size();++i)
        if(gas.IsInteriorFace(c,i,j))gas.AddMacFaceVelocity(c,i,j,4*(Random()-0.5));
    gas.RefreshBoundaries();
    std::vector<WaterParticle> particles;
    for(int j=0;j<5;++j)for(int i=0;i<5;++i)
        particles.push_back({{-1.6+i*0.71,-1.5+j*0.68},{4*(Random()-0.5),4*(Random()-0.5)}});
    const double wall=2-sim.Liquid().ParticleRadius();
    for(Vec2 p:std::vector<Vec2>{{-wall,-wall},{wall,-wall},{-wall,wall},{wall,wall},{-wall,0},{wall,0},{0,-wall},{0,wall}})
        particles.push_back({p,{4*(Random()-0.5),4*(Random()-0.5)}});
    sim.Liquid().SetParticleState(particles);
}
int main() {
    int scenarios=0,failed=0;
    double maxRelativeMomentum=0,maxEnergyIncrease=0;
    for(int resolution:{8,12,24,48})for(double density:{0.1,1.2,80.,1000.})for(double dt:{1e-5,1./120,0.2,10.}) {
        TwoMaterialSimulation sim(resolution);sim.SetJetEnabled(false);sim.SetGasDensity(density);Fill(sim);
        Totals before=ReadActual(sim);sim.ExchangeOnly(dt);Totals after=ReadActual(sim);
        Vec2 actualResidual=after.gas+after.liquid-before.gas-before.liquid;
        double scale=1+before.gas.Length()+before.liquid.Length();
        double relative=actualResidual.Length()/scale;
        double energyIncrease=(after.energy-before.energy)/(1+before.energy);
        maxRelativeMomentum=(std::max)(maxRelativeMomentum,relative);
        maxEnergyIncrease=(std::max)(maxEnergyIncrease,energyIncrease);
        bool okay=relative<2e-12 && energyIncrease<2e-12 && sim.IsFinite()
            &&(sim.couplingResidual-actualResidual).Length()<2e-9*scale;
        if(!okay){++failed;std::cout<<"FAIL n="<<resolution<<" rho="<<density<<" dt="<<dt<<" p="<<relative<<" dE="<<energyIncrease<<'\n';}
        ++scenarios;
    }
    TwoMaterialSimulation equilibrium(24);
    equilibrium.Gas().Clear();
    for(int c=0;c<2;++c)for(int j=1;j<=24;++j)for(int i=1;i<=24;++i)
        if(equilibrium.Gas().IsInteriorFace(c,i,j))equilibrium.Gas().AddMacFaceVelocity(c,i,j,c==0?0.4:-0.3);
    equilibrium.Gas().RefreshBoundaries();
    equilibrium.Liquid().SetParticleState({{{-0.3,0.3},{0.4,-0.3}},{{0.4,-0.4},{0.4,-0.3}}});
    auto equilibriumBefore=ReadActual(equilibrium);equilibrium.ExchangeOnly(10);auto equilibriumAfter=ReadActual(equilibrium);
    if((equilibriumBefore.liquid-equilibriumAfter.liquid).Length()>1e-11||std::abs(equilibriumBefore.energy-equilibriumAfter.energy)>1e-10)++failed;
    for(CouplingMode mode:{CouplingMode::Off,CouplingMode::OneWay,CouplingMode::TwoWay}) {
        TwoMaterialSimulation sim(24);sim.SetCouplingMode(mode);Fill(sim);auto before=ReadActual(sim);sim.ExchangeOnly(.1);auto after=ReadActual(sim);
        if(mode!=CouplingMode::TwoWay&&(after.gas-before.gas).Length()!=0)++failed;
        if(mode==CouplingMode::Off&&(after.liquid-before.liquid).Length()!=0)++failed;
        if(mode!=CouplingMode::Off&&(after.liquid-before.liquid).Length()<1e-6)++failed;
    }
    std::cout<<std::setprecision(17)<<"Independent review scenarios: "<<scenarios<<" randomized near-wall cases + equilibrium + 3 modes\n"
        <<"Maximum relative actual-state momentum residual: "<<maxRelativeMomentum<<"\nMaximum relative kinetic energy increase: "<<maxEnergyIncrease<<'\n'
        <<"Independent exchange review failures: "<<failed<<'\n';

    TwoMaterialSimulation production(48);
    const double initialMass=production.Liquid().Count()*production.Liquid().ParticleMass();
    bool finite=true,inside=true,closed=true,massPreserved=true;
    double maxStepRelativeResidual=0,maxStepAbsoluteResidual=0;
    const auto started=std::chrono::steady_clock::now();
    for(int step=0;step<600;++step) {
        production.Step(1.0/120.0);
        const auto actual=ReadActual(production);
        finite &= production.IsFinite()&&std::isfinite(actual.energy)&&std::isfinite(actual.gas.x)&&std::isfinite(actual.gas.y);
        const double limit=2-production.Liquid().ParticleRadius()+1e-10;
        for(const auto& p:production.Liquid().Particles())inside &= std::abs(p.position.x)<=limit&&std::abs(p.position.y)<=limit;
        for(int i=0;i<48;++i) {
            double q=(i+0.5)/48.0;
            closed &= std::abs(production.Gas().SampleVelocity({0,q}).x)<1e-10&&std::abs(production.Gas().SampleVelocity({1,q}).x)<1e-10
                &&std::abs(production.Gas().SampleVelocity({q,0}).y)<1e-10&&std::abs(production.Gas().SampleVelocity({q,1}).y)<1e-10;
        }
        massPreserved &= production.Liquid().Count()==400&&std::abs(production.Liquid().Count()*production.Liquid().ParticleMass()-initialMass)<1e-10;
        double absolute=production.couplingResidual.Length();
        maxStepAbsoluteResidual=(std::max)(maxStepAbsoluteResidual,absolute);
        maxStepRelativeResidual=(std::max)(maxStepRelativeResidual,absolute/(1+production.couplingWaterImpulse.Length()+production.couplingGasImpulse.Length()));
    }
    bool productionPassed=finite&&inside&&closed&&massPreserved&&maxStepRelativeResidual<1e-10&&std::abs(production.Time()-5)<1e-10;
    if(!productionPassed)++failed;
    std::cout<<"Default 48x48 full evolution: 600 x 1/120 s; elapsed "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()<<" s\n"
        <<"Default full-run finite="<<finite<<" inside="<<inside<<" closed_walls="<<closed<<" fixed_water_mass="<<massPreserved<<'\n'
        <<"Default full-run maximum exchange absolute residual: "<<maxStepAbsoluteResidual<<"; relative: "<<maxStepRelativeResidual<<'\n';
    TwoMaterialSimulation noCoupling(48),withCoupling(48);noCoupling.SetCouplingMode(CouplingMode::Off);
    for(int step=0;step<120;++step){noCoupling.Step(1.0/120);withCoupling.Step(1.0/120);}
    double waterDifference=0,gasDifference=0;int gasFaceCount=0;
    for(std::size_t i=0;i<noCoupling.Liquid().Count();++i){Vec2 d=noCoupling.Liquid().Particles()[i].position-withCoupling.Liquid().Particles()[i].position;waterDifference+=d.x*d.x+d.y*d.y;}
    for(int c=0;c<2;++c)for(int j=1;j<=48;++j)for(int i=1;i<=48;++i)if(noCoupling.Gas().IsInteriorFace(c,i,j)){
        double d=noCoupling.Gas().MacFaceVelocity(c,i,j)-withCoupling.Gas().MacFaceVelocity(c,i,j);gasDifference+=d*d;++gasFaceCount;
    }
    waterDifference=std::sqrt(waterDifference/noCoupling.Liquid().Count());gasDifference=std::sqrt(gasDifference/gasFaceCount);
    if(!(waterDifference>1e-4&&gasDifference>1e-4&&noCoupling.IsFinite()&&withCoupling.IsFinite()))++failed;
    std::cout<<"Default 48x48 Off/TwoWay at 1 s: particle position RMS "<<waterDifference<<"; MAC-face gas velocity RMS "<<gasDifference<<'\n'
        <<"Independent review failures: "<<failed<<'\n';
    return failed?1:0;
}
