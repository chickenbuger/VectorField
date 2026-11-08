#include "SharpInterfaceBaseline.h"
#include "SharpInterfaceSimulation.h"
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace fieldlab;
namespace {
using Clock=std::chrono::steady_clock;
constexpr double Dt=1.0/120.0;
constexpr int N=48,Steps=600,Trials=3;
struct State { std::vector<double> phi,pressure,dye,u,v; double area=0; };
struct Run {
    std::string version;double density=80;int trial=0;
    double totalMs=0,medianMs=0,p95Ms=0,maxMs=0;
    double maxDivergence=0,maxArea=0,maxRawArea=0,maxCfl=0;
    std::size_t substeps=0,iterations=0;
    bool finite=true,walls=true,fractions=true,converged=true;
    State state;
};
double Percentile(std::vector<double> samples,double p) {
    std::sort(samples.begin(),samples.end());
    return samples[static_cast<std::size_t>(std::ceil(p*samples.size()))-1];
}
template<class Sim> double ActualDivergence(const Sim& sim) {
    double sum=0;for(int j=1;j<=N;++j)for(int i=1;i<=N;++i){
        const double d=(sim.FaceVelocity(0,i,j)-sim.FaceVelocity(0,i-1,j)+sim.FaceVelocity(1,i,j)-sim.FaceVelocity(1,i,j-1))/sim.CellWidth();sum+=d*d;
    }return std::sqrt(sum/(N*N));
}
template<class Sim> Run Execute(const std::string& version,double density,int trial,std::ofstream& raw) {
    Sim sim(N);sim.SetGasDensity(density);
    Run out;out.version=version;out.density=density;out.trial=trial;
    std::vector<double> stepTimes;stepTimes.reserve(Steps);
    for(int step=0;step<Steps;++step) {
        const auto start=Clock::now();sim.Step(Dt);
        const double ms=std::chrono::duration<double,std::milli>(Clock::now()-start).count();
        stepTimes.push_back(ms);out.totalMs+=ms;out.substeps+=sim.lastSubsteps;out.iterations+=sim.pressureIterations;
        // Validation is outside the timed Step interval and uses actual faces.
        const double div=ActualDivergence(sim);
        out.maxDivergence=(std::max)(out.maxDivergence,div);
        out.maxArea=(std::max)(out.maxArea,std::abs(sim.RelativeAreaError()));
        out.maxRawArea=(std::max)(out.maxRawArea,std::abs(sim.RawAreaError()));
        out.maxCfl=(std::max)(out.maxCfl,sim.lastCfl);
        out.finite=out.finite&&sim.IsFinite();out.converged=out.converged&&sim.pressureConverged;
        for(int k=1;k<=N;++k)out.walls=out.walls&&sim.FaceVelocity(0,0,k)==0&&sim.FaceVelocity(0,N,k)==0&&sim.FaceVelocity(1,k,0)==0&&sim.FaceVelocity(1,k,N)==0;
        for(int j=1;j<=N;++j)for(int i=1;i<=N;++i)out.fractions=out.fractions&&sim.LiquidFraction(i,j)>=0&&sim.LiquidFraction(i,j)<=1;
        raw<<version<<','<<density<<','<<trial<<','<<step+1<<','<<std::setprecision(17)<<ms<<','<<sim.lastSubsteps<<','<<sim.pressureIterations<<','<<sim.pressureResidual<<','<<div<<','<<sim.RelativeAreaError()<<','<<sim.RawAreaError()<<','<<sim.lastCfl<<'\n';
    }
    out.medianMs=Percentile(stepTimes,.5);out.p95Ms=Percentile(stepTimes,.95);out.maxMs=Percentile(stepTimes,1);
    for(int j=1;j<=N;++j)for(int i=1;i<=N;++i){out.state.phi.push_back(sim.Phi(i,j));out.state.pressure.push_back(sim.Pressure(i,j));out.state.dye.push_back(sim.Dye(i,j));}
    for(int j=1;j<=N;++j)for(int i=0;i<=N;++i)out.state.u.push_back(sim.FaceVelocity(0,i,j));
    for(int j=0;j<=N;++j)for(int i=1;i<=N;++i)out.state.v.push_back(sim.FaceVelocity(1,i,j));
    out.state.area=sim.WaterArea();
    out.finite=out.finite&&std::abs(sim.Time()-5)<1e-9;
    return out;
}
struct Difference {double rms=0,max=0,relativeRms=0;};
Difference Compare(const std::vector<double>& a,const std::vector<double>& b,bool gauge=false) {
    Difference d;double mean=0,scale=0;
    if(gauge){for(std::size_t k=0;k<a.size();++k)mean+=a[k]-b[k];mean/=a.size();}
    for(std::size_t k=0;k<a.size();++k){const double e=a[k]-b[k]-mean;d.rms+=e*e;d.max=(std::max)(d.max,std::abs(e));scale+=a[k]*a[k];}
    d.rms=std::sqrt(d.rms/a.size());d.relativeRms=d.rms/(std::max)(1.0,std::sqrt(scale/a.size()));return d;
}
void WriteSummary(std::ofstream& csv,const Run& r) {
    csv<<r.version<<','<<r.density<<','<<r.trial<<','<<N<<','<<Steps<<','<<std::setprecision(17)<<r.totalMs<<','<<r.medianMs<<','<<r.p95Ms<<','<<r.maxMs<<','<<r.substeps<<','<<r.iterations<<','<<r.maxDivergence<<','<<r.maxArea<<','<<r.maxRawArea<<','<<r.maxCfl<<','<<r.finite<<','<<r.walls<<','<<r.fractions<<','<<r.converged<<'\n';
    std::cout<<r.version<<" rho="<<r.density<<" trial="<<r.trial<<": total Step CPU wall="<<r.totalMs<<" ms; step p50="<<r.medianMs<<", p95="<<r.p95Ms<<", max="<<r.maxMs<<" ms; internal substeps="<<r.substeps<<", last-substep PCG iterations sum="<<r.iterations<<".\n"<<std::flush;
}
template<class Sim> void WarmUp() {Sim sim(24);for(int step=0;step<20;++step)sim.Step(Dt);}
int ProjectionCheck(double density) {
    BaselineSharpInterfaceSimulation base(32);SharpInterfaceSimulation next(32);
    base.SetGasDensity(density);next.SetGasDensity(density);
    const auto phi=[](Vec2 p){return p.x+.37*p.y-.11;};
    const auto velocity=[](Vec2 p){return Vec2{std::sin(3.7*p.x+1.4*p.y),std::cos(2.2*p.y-.8*p.x)};};
    base.SetPhi(phi);next.SetPhi(phi);base.SetVelocity(velocity);next.SetVelocity(velocity);
    base.ProjectOnly(Dt);next.ProjectOnly(Dt);
    double maxVelocity=0,pmean=0,pRms=0;
    for(int j=1;j<=32;++j)for(int i=1;i<=32;++i)pmean+=base.Pressure(i,j)-next.Pressure(i,j);pmean/=32*32;
    for(int j=1;j<=32;++j)for(int i=1;i<=32;++i){const double e=base.Pressure(i,j)-next.Pressure(i,j)-pmean;pRms+=e*e;}
    for(int j=1;j<=32;++j)for(int i=0;i<=32;++i)maxVelocity=(std::max)(maxVelocity,std::abs(base.FaceVelocity(0,i,j)-next.FaceVelocity(0,i,j)));
    for(int j=0;j<=32;++j)for(int i=1;i<=32;++i)maxVelocity=(std::max)(maxVelocity,std::abs(base.FaceVelocity(1,i,j)-next.FaceVelocity(1,i,j)));
    pRms=std::sqrt(pRms/(32*32));
    std::cout<<"Fixed-geometry paired projection rho="<<density<<": MAC-face velocity max difference="<<maxVelocity<<"; gauge-adjusted pressure RMS="<<pRms<<"; relative residual baseline="<<base.pressureResidual<<", optimized="<<next.pressureResidual<<".\n";
    return maxVelocity<1e-8&&pRms<1e-4&&next.pressureConverged?0:1;
}
}
int main() {
    std::cout<<std::setprecision(12);
    std::ofstream raw("performance_steps.csv"),summary("performance_summary.csv"),differences("paired_state_differences.csv");
    if(!raw||!summary||!differences)return 2;
    raw<<"version,gas_density,trial,fixed_step,step_ms,internal_substeps,last_substep_pcg_iterations,last_substep_relative_residual,actual_divergence_rms,corrected_relative_area_error,raw_relative_area_error,estimated_cfl\n";
    summary<<"version,gas_density,trial,grid,fixed_steps,total_step_ms,step_p50_ms,step_p95_ms,step_max_ms,total_internal_substeps,last_substep_pcg_iterations_sum,max_actual_divergence_rms,max_corrected_relative_area_error,max_raw_relative_area_error,max_estimated_cfl,finite_and_time,closed_walls,bounded_fraction,all_final_substeps_pcg_converged\n";
    differences<<"gas_density,trial,field,rms_difference,max_difference,relative_rms_difference\n";
    WarmUp<BaselineSharpInterfaceSimulation>();WarmUp<SharpInterfaceSimulation>();
    int failures=ProjectionCheck(80)+ProjectionCheck(1.2);
    std::vector<Run> runs;
    for(double density:{80.0,1.2})for(int trial=1;trial<=Trials;++trial){
        Run base,next;
        if(trial%2){base=Execute<BaselineSharpInterfaceSimulation>("baseline_75cf82b",density,trial,raw);WriteSummary(summary,base);next=Execute<SharpInterfaceSimulation>("optimized",density,trial,raw);WriteSummary(summary,next);}
        else{next=Execute<SharpInterfaceSimulation>("optimized",density,trial,raw);WriteSummary(summary,next);base=Execute<BaselineSharpInterfaceSimulation>("baseline_75cf82b",density,trial,raw);WriteSummary(summary,base);}
        for(const Run* r:{&base,&next})if(!(r->finite&&r->walls&&r->fractions&&r->converged&&r->maxDivergence<1e-7&&r->maxArea<1e-8&&r->maxCfl<.6500001))++failures;
        const auto phi=Compare(base.state.phi,next.state.phi),pressure=Compare(base.state.pressure,next.state.pressure,true),dye=Compare(base.state.dye,next.state.dye),u=Compare(base.state.u,next.state.u),v=Compare(base.state.v,next.state.v);
        const Difference ds[]={phi,pressure,dye,u,v};const char* names[]={"phi","pressure_mean_gauge","dye","u_faces","v_faces"};
        for(int k=0;k<5;++k)differences<<std::setprecision(17)<<density<<','<<trial<<','<<names[k]<<','<<ds[k].rms<<','<<ds[k].max<<','<<ds[k].relativeRms<<'\n';
        // Before measurement: tolerate <=0.1% relative velocity/pressure RMS,
        // phi RMS <=0.1% of a cell, and corrected area difference <=1e-8.
        if(phi.rms>1e-3*(4.0/N)||pressure.relativeRms>1e-3||u.relativeRms>1e-3||v.relativeRms>1e-3||std::abs(base.state.area-next.state.area)>1e-8)++failures;
        std::cout<<"Paired final state rho="<<density<<" trial="<<trial<<": phi RMS="<<phi.rms<<", pressure RMS="<<pressure.rms<<", u/v face RMS="<<u.rms<<'/'<<v.rms<<", area difference="<<std::abs(base.state.area-next.state.area)<<".\n";
        runs.push_back(base);runs.push_back(next);
    }
    for(double density:{80.0,1.2}){
        std::vector<double> baseTimes,nextTimes;
        for(const auto& r:runs)if(r.density==density)(r.version=="optimized"?nextTimes:baseTimes).push_back(r.totalMs);
        const double b=Percentile(baseTimes,.5),n=Percentile(nextTimes,.5);
        std::cout<<"Median total Step time rho="<<density<<": baseline="<<b<<" ms; optimized="<<n<<" ms; speedup="<<b/n<<"x; reduction="<<100*(1-n/b)<<"%.\n";
    }
    for(const auto& r:runs)if(r.trial>1){
        const Run* first=nullptr;for(const auto& candidate:runs)if(candidate.trial==1&&candidate.density==r.density&&candidate.version==r.version){first=&candidate;break;}
        if(!first){++failures;continue;}
        const double maximum=(std::max)({Compare(first->state.phi,r.state.phi).max,Compare(first->state.pressure,r.state.pressure).max,Compare(first->state.dye,r.state.dye).max,Compare(first->state.u,r.state.u).max,Compare(first->state.v,r.state.v).max});
        std::cout<<"Repeat determinism "<<r.version<<" rho="<<r.density<<" trial="<<r.trial<<": maximum final-array difference="<<maximum<<".\n";
        if(maximum>1e-12)++failures;
    }
    std::cout<<"Performance/paired-state validation failures: "<<failures<<".\n";
    return failures?1:0;
}
