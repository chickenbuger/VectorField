#include "SharpInterfaceProfile.h"
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>

int main() {
    fieldlab::ProfiledSharpInterfaceSimulation sim(48);
    sim.SetGasDensity(80.0);
    const auto start=std::chrono::steady_clock::now();
    for(int step=0;step<600;++step)sim.Step(1.0/120.0);
    const double total=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count();
    const auto& p=sim.profile;
    const char* names[]={"snapshots","phi_advection","reinitialization","area_correction","face_coefficients","velocity_advection_and_sources","pressure_projection","dye_transport"};
    const double values[]={p.snapshots,p.phi,p.reinitialize,p.area,p.coefficients,p.velocity,p.pressure,p.dye};
    std::ofstream csv("baseline_phase_profile.csv");
    csv<<"baseline_commit,grid,gas_density,fixed_steps,internal_substeps,phase,milliseconds,percent_of_total\n";
    std::cout<<std::setprecision(12)<<"Baseline 75cf82b instrumented 48x48 rho_g=80: "<<total<<" ms, "<<p.substeps<<" internal substeps.\n";
    for(int i=0;i<8;++i){
        csv<<std::setprecision(17)<<"75cf82b,48,80,600,"<<p.substeps<<','<<names[i]<<','<<values[i]<<','<<100*values[i]/total<<'\n';
        std::cout<<names[i]<<": "<<values[i]<<" ms ("<<100*values[i]/total<<"%).\n";
    }
    std::cout<<"finite="<<sim.IsFinite()<<" time="<<sim.Time()<<" last_pcg_converged="<<sim.pressureConverged<<'\n';
    return sim.IsFinite()&&sim.pressureConverged?0:1;
}
