#ifndef FIELDLAB_SHARP_INTERFACE_SIMULATION_H
#define FIELDLAB_SHARP_INTERFACE_SIMULATION_H

#include "FieldMath.h"
#include <functional>
#include <limits>
#include <vector>

namespace fieldlab {

// Stage 9: one MAC velocity and a signed-distance material interface. Pressure
// is continuous (surface tension = 0); only its normal derivative changes.
// The level set is not a conservative VOF mass/flux representation.
class SharpInterfaceSimulation {
    struct Workspace {
        std::vector<double> b,q,r,z,d,ad,diag,oldU,oldV,oldPhi,oldDye,distance,icInverse,icLeft,icBottom,forward;
        std::vector<unsigned char> fixed;
        std::vector<std::size_t> areaBand;
        explicit Workspace(std::size_t count)
            : b(count),q(count),r(count),z(count),d(count),ad(count),diag(count),
              oldU(count),oldV(count),oldPhi(count),oldDye(count),distance(count),
              icInverse(count),icLeft(count),icBottom(count),forward(count),fixed(count) {
            areaBand.reserve(count);
        }
    };
public:
    explicit SharpInterfaceSimulation(int n=48)
        : n_((std::max)(8,(std::min)(128,n))), s_(n_+2), h_(4.0/n_),
          phi_(Count()), u_(Count()), v_(Count()), dye_(Count()), pressure_(Count()),
          bx_(Count()), by_(Count()), work_(Count()) { Reset(); }

    int Size() const { return n_; }
    double CellWidth() const { return h_; }
    double Time() const { return time_; }
    double FixedTimeStep() const { return 1.0/120.0; }
    double GasDensity() const { return gasDensity_; }
    double LiquidDensity() const { return liquidDensity_; }
    bool SubcellInterface() const { return subcell_; }
    bool ProjectionEnabled() const { return projection_; }
    bool JetEnabled() const { return jet_; }
    double Phi(int i,int j) const { return phi_[I(Cell(i),Cell(j))]; }
    double LiquidFraction(int i,int j) const { return Fraction(Phi(i,j)); }
    double Pressure(int i,int j) const { return pressure_[I(Cell(i),Cell(j))]; }
    double Dye(int i,int j) const { return dye_[I(Cell(i),Cell(j))]; }
    double U(int i,int j) const { i=Cell(i);j=Cell(j);return 0.5*(u_[I(i-1,j)]+u_[I(i,j)]); }
    double V(int i,int j) const { i=Cell(i);j=Cell(j);return 0.5*(v_[I(i,j-1)]+v_[I(i,j)]); }
    Vec2 Gravity() const { return gravity_; }
    double FaceVelocity(int c,int i,int j) const {
        if(c==0 && i>=0 && i<=n_ && j>=1 && j<=n_)return u_[I(i,j)];
        if(c==1 && i>=1 && i<=n_ && j>=0 && j<=n_)return v_[I(i,j)];
        return 0;
    }
    double FaceBeta(int c,int i,int j) const {
        return Interior(c,i,j)?(c==0?bx_:by_)[I(i,j)]:0.0;
    }
    double FaceTheta(int c,int i,int j) const {
        if(!Interior(c,i,j))return 0;
        const double a=Phi(i,j),b=c==0?Phi(i+1,j):Phi(i,j+1);
        return (a<0)!=(b<0)?std::abs(a)/(std::abs(a)+std::abs(b)):0.5;
    }
    void AddFaceVelocity(int c,int i,int j,double increment) {
        if(Interior(c,i,j) && std::isfinite(increment)) {
            double& value=(c==0?u_:v_)[I(i,j)];
            if(std::isfinite(value+increment))value+=increment;
        }
    }
    void SetGasDensity(double density) { RequirePositive(density);gasDensity_=Clamp(density,0.1,5000.0);BuildCoefficients(); }
    void SetLiquidDensity(double density) { RequirePositive(density);liquidDensity_=Clamp(density,0.1,5000.0);BuildCoefficients(); }
    void SetSubcellInterface(bool enabled) { subcell_=enabled;BuildCoefficients(); }
    void SetProjectionEnabled(bool enabled) { projection_=enabled; }
    void SetJetEnabled(bool enabled) { jet_=enabled; }
    void SetGravity(Vec2 gravity) { if(!Finite(gravity))throw std::invalid_argument("Gravity must be finite.");gravity_=gravity; }
    void SetPhi(const std::function<double(Vec2)>& field) {
        auto next=phi_;
        for(int j=1;j<=n_;++j)for(int i=1;i<=n_;++i) {
            const double value=field(Center(i,j));
            if(!std::isfinite(value))throw std::invalid_argument("Level set must be finite.");
            next[I(i,j)]=value;
        }
        phi_=next;ScalarBoundary(phi_);BuildCoefficients();
        initialArea_=WaterArea();rawAreaError_=areaShift_=0;
    }
    void SetVelocity(const std::function<Vec2(Vec2)>& field) {
        auto nextU=u_,nextV=v_;
        for(int j=1;j<=n_;++j)for(int i=1;i<n_;++i) {
            const Vec2 value=field({-2+i*h_,-2+(j-0.5)*h_});
            if(!Finite(value))throw std::invalid_argument("Velocity must be finite.");nextU[I(i,j)]=value.x;
        }
        for(int j=1;j<n_;++j)for(int i=1;i<=n_;++i) {
            const Vec2 value=field({-2+(i-0.5)*h_,-2+j*h_});
            if(!Finite(value))throw std::invalid_argument("Velocity must be finite.");nextV[I(i,j)]=value.y;
        }
        u_=nextU;v_=nextV;VelocityBoundary();
    }

    double divergenceBeforeProjection=0,divergenceAfterProjection=0;
    double pressureResidual=0;
    int pressureIterations=0;
    int lastSubsteps=1;
    double lastCfl=0;
    bool pressureConverged=true;

    void Reset() {
        time_=accumulator_=0;gasDensity_=80;liquidDensity_=1000;
        subcell_=projection_=jet_=true;gravity_={0,-9.81};
        for(auto* a:{&u_,&v_,&dye_,&pressure_})std::fill(a->begin(),a->end(),0.0);
        // A liquid dam; signed distance to the box [-1.75,-0.20] x [-2,0.65].
        SetPhi([](Vec2 p){
            const Vec2 q{std::abs(p.x+0.975)-0.775,std::abs(p.y+0.675)-1.325};
            return Vec2{(std::max)(0.0,q.x),(std::max)(0.0,q.y)}.Length()
                +(std::min)(0.0,(std::max)(q.x,q.y));
        });
        divergenceBeforeProjection=divergenceAfterProjection=pressureResidual=0;
        pressureIterations=0;pressureConverged=true;lastSubsteps=1;lastCfl=0;
    }
    Vec2 SampleVelocity(Vec2 normalized) const { return SampleVelocityFrom(u_,v_,normalized); }
    double SampleDye(Vec2 normalized) const {
        if(!Finite(normalized))return 0;
        return Bilinear(dye_,normalized.x*n_+0.5,normalized.y*n_+0.5);
    }
    void AddDye(double x,double y,double amount,double radius=0.06) {
        if(!ValidSource(x,y,amount,radius)||amount<=0)return;
        radius=(std::max)(radius,0.75/n_);
        for(int j=1;j<=n_;++j)for(int i=1;i<=n_;++i)if(Phi(i,j)>=0)
            dye_[I(i,j)]=Clamp(dye_[I(i,j)]+amount*Splat((i-0.5)/n_,(j-0.5)/n_,x,y,radius),0,10);
        ScalarBoundary(dye_);
    }
    // Source increments act in gas; the pressure solve transmits their effect.
    void AddVelocity(double x,double y,double du,double dv,double radius=0.08) {
        if(!ValidSource(x,y,du,radius)||!std::isfinite(dv))return;
        radius=(std::max)(radius,0.75/n_);
        for(int j=1;j<=n_;++j)for(int i=1;i<n_;++i)if(Phi(i,j)>=0&&Phi(i+1,j)>=0)
            u_[I(i,j)]+=Clamp(du,-10,10)*Splat(double(i)/n_,(j-0.5)/n_,x,y,radius);
        for(int j=1;j<n_;++j)for(int i=1;i<=n_;++i)if(Phi(i,j)>=0&&Phi(i,j+1)>=0)
            v_[I(i,j)]+=Clamp(dv,-10,10)*Splat((i-0.5)/n_,double(j)/n_,x,y,radius);
        VelocityBoundary();
    }
    void Stir(Vec2 world,Vec2 delta) {
        if(!Finite(world)||!Finite(delta))return;
        if(delta.Length()>0.2)delta=delta*(0.2/delta.Length());
        for(int j=1;j<=n_;++j)for(int i=1;i<n_;++i)if(Phi(i,j)<0||Phi(i+1,j)<0)
            u_[I(i,j)]+=25*delta.x*Splat(-2+i*h_,-2+(j-0.5)*h_,world.x,world.y,0.4);
        for(int j=1;j<n_;++j)for(int i=1;i<=n_;++i)if(Phi(i,j)<0||Phi(i,j+1)<0)
            v_[I(i,j)]+=25*delta.y*Splat(-2+(i-0.5)*h_,-2+j*h_,world.x,world.y,0.4);
    }
    double DivergenceRms() const {
        double value=0;for(int j=1;j<=n_;++j)for(int i=1;i<=n_;++i){const double d=Div(i,j);value+=d*d;}
        return std::sqrt(value/(n_*n_));
    }
    double KineticEnergy() const {
        double value=0;
        for(int j=1;j<=n_;++j)for(int i=1;i<n_;++i)value+=0.5*h_*h_*u_[I(i,j)]*u_[I(i,j)]/bx_[I(i,j)];
        for(int j=1;j<n_;++j)for(int i=1;i<=n_;++i)value+=0.5*h_*h_*v_[I(i,j)]*v_[I(i,j)]/by_[I(i,j)];
        return value;
    }
    double WaterArea() const { return AreaAtOffset(0); }
    double InitialWaterArea() const { return initialArea_; }
    double RelativeAreaError() const { return (WaterArea()-initialArea_)/(std::max)(initialArea_,1e-12); }
    double RawAreaError() const { return rawAreaError_; }
    double AreaCorrectionShift() const { return areaShift_; }
    double MaxSpeed() const {
        double value=0;for(int j=1;j<=n_;++j)for(int i=1;i<=n_;++i)value=(std::max)(value,Vec2{U(i,j),V(i,j)}.Length());return value;
    }
    double FaceSpeedBound() const {
        double x=0,y=0;
        for(int j=1;j<=n_;++j)for(int i=1;i<n_;++i)x=(std::max)(x,std::abs(u_[I(i,j)]));
        for(int j=1;j<n_;++j)for(int i=1;i<=n_;++i)y=(std::max)(y,std::abs(v_[I(i,j)]));
        // A componentwise interpolation and 2D transport bound. Cell-center
        // averaging alone could hide alternating-sign face velocities.
        return x+y;
    }

    // Fixed geometry projection is exposed for independent manufactured tests.
    void ProjectOnly(double dt=1.0/120.0) {
        if(!(dt>0)||!std::isfinite(dt))return;
        VelocityBoundary();BuildCoefficients();divergenceBeforeProjection=DivergenceRms();
        auto& b=work_.b;auto& q=work_.q;auto& r=work_.r;auto& z=work_.z;
        auto& d=work_.d;auto& ad=work_.ad;auto& diag=work_.diag;
        std::fill(q.begin(),q.end(),0.0);
        for(int j=1;j<=n_;++j)for(int i=1;i<=n_;++i) {
            const auto k=I(i,j);b[k]=-Div(i,j);
            diag[k]=((i>1?bx_[I(i-1,j)]:0)+(i<n_?bx_[I(i,j)]:0)
                +(j>1?by_[I(i,j-1)]:0)+(j<n_?by_[I(i,j)]:0))/(h_*h_);
        }
        ZeroMean(b);r=b;const double bnorm=std::sqrt(Dot(b,b));
        const double tolerance=(std::max)(1e-11,1e-10*bnorm);
        BuildModifiedIncompleteCholesky(diag);
        ApplyPreconditioner(r,z);
        ZeroMean(z);d=z;double rz=Dot(r,z);pressureIterations=0;
        double norm=bnorm;
        for(int iteration=0;iteration<1000 && norm>tolerance;++iteration) {
            Apply(d,ad);const double denominator=Dot(d,ad);
            if(!(denominator>0)||!std::isfinite(denominator)||!std::isfinite(rz))break;
            const double alpha=rz/denominator;
            for(int j=1;j<=n_;++j)for(int i=1;i<=n_;++i){const auto k=I(i,j);q[k]+=alpha*d[k];r[k]-=alpha*ad[k];}
            ZeroMean(r);norm=std::sqrt(Dot(r,r));pressureIterations=iteration+1;
            if(norm<=tolerance)break;
            ApplyPreconditioner(r,z);
            ZeroMean(z);const double next=Dot(r,z),beta=next/rz;rz=next;
            for(int j=1;j<=n_;++j)for(int i=1;i<=n_;++i)d[I(i,j)]=z[I(i,j)]+beta*d[I(i,j)];
            ZeroMean(d);
        }
        ZeroMean(q);
        // Report an explicit A*q-b residual, not just the CG recurrence estimate.
        Apply(q,ad);double actual=0;
        for(int j=1;j<=n_;++j)for(int i=1;i<=n_;++i){const double e=ad[I(i,j)]-b[I(i,j)];actual+=e*e;pressure_[I(i,j)]=q[I(i,j)]/dt;}
        pressureResidual=std::sqrt(actual)/(std::max)(bnorm,1e-12);
        pressureConverged=std::sqrt(actual)<=(std::max)(2e-10,5e-10*bnorm);
        for(int j=1;j<=n_;++j)for(int i=1;i<n_;++i)u_[I(i,j)]-=bx_[I(i,j)]*(q[I(i+1,j)]-q[I(i,j)])/h_;
        for(int j=1;j<n_;++j)for(int i=1;i<=n_;++i)v_[I(i,j)]-=by_[I(i,j)]*(q[I(i,j+1)]-q[I(i,j)])/h_;
        VelocityBoundary();ScalarBoundary(pressure_);divergenceAfterProjection=DivergenceRms();
    }
    void Step(double dt) {
        if(!(dt>0)||!std::isfinite(dt))return;
        accumulator_+=(std::min)(0.05,dt);const double fixed=FixedTimeStep();int frames=0;
        while(accumulator_+1e-12>=fixed&&frames<6) {
            const double speed=FaceSpeedBound()+(std::abs(gravity_.x)+std::abs(gravity_.y)+(jet_?29.0:0.0))*fixed;
            lastSubsteps=(std::max)(1,(std::min)(16,static_cast<int>(std::ceil(speed*fixed/(0.65*h_)))));
            const double subdt=fixed/lastSubsteps;lastCfl=speed*subdt/h_;
            for(int k=0;k<lastSubsteps;++k)Substep(subdt);
            time_+=fixed;accumulator_=(std::max)(0.0,accumulator_-fixed);++frames;
        }
    }
    bool IsFinite() const {
        if(!std::isfinite(time_)||!std::isfinite(accumulator_)||!std::isfinite(pressureResidual)||!std::isfinite(WaterArea()))return false;
        for(const auto* a:{&phi_,&u_,&v_,&dye_,&pressure_,&bx_,&by_})for(double value:*a)if(!std::isfinite(value))return false;
        return true;
    }

private:
    int n_,s_;double h_,time_=0,accumulator_=0,gasDensity_=80,liquidDensity_=1000;
    bool subcell_=true,projection_=true,jet_=true;Vec2 gravity_{0,-9.81};
    double initialArea_=0,rawAreaError_=0,areaShift_=0;
    std::vector<double> phi_,u_,v_,dye_,pressure_,bx_,by_;
    Workspace work_;
    std::size_t Count() const { return static_cast<std::size_t>(s_*s_); }
    std::size_t I(int i,int j) const { return static_cast<std::size_t>(j*s_+i); }
    int Cell(int i) const { return (std::max)(1,(std::min)(n_,i)); }
    Vec2 Center(int i,int j) const { return {-2+(i-0.5)*h_,-2+(j-0.5)*h_}; }
    static double Clamp(double x,double a,double b) { return (std::max)(a,(std::min)(b,x)); }
    static bool Finite(Vec2 p) { return std::isfinite(p.x)&&std::isfinite(p.y); }
    static void RequirePositive(double x) { if(!(x>0)||!std::isfinite(x))throw std::invalid_argument("Density must be finite and positive."); }
    bool Interior(int c,int i,int j) const { return c==0?i>=1&&i<n_&&j>=1&&j<=n_:c==1&&i>=1&&i<=n_&&j>=1&&j<n_; }
    static bool ValidSource(double x,double y,double value,double r) { return std::isfinite(x)&&std::isfinite(y)&&std::isfinite(value)&&std::isfinite(r)&&r>0; }
    static double Splat(double px,double py,double x,double y,double radius) { const double t=(std::max)(0.0,1-((px-x)*(px-x)+(py-y)*(py-y))/(radius*radius));return t*t; }
    double Density(double phi) const { return phi<0?liquidDensity_:gasDensity_; }
    double Beta(double left,double right) const {
        if((left<0)==(right<0))return 1.0/Density(left);
        const double fraction=subcell_?std::abs(left)/(std::abs(left)+std::abs(right)):0.5;
        return 1.0/(fraction*Density(left)+(1-fraction)*Density(right));
    }
    void BuildCoefficients() {
        std::fill(bx_.begin(),bx_.end(),0.0);std::fill(by_.begin(),by_.end(),0.0);
        for(int j=1;j<=n_;++j)for(int i=1;i<n_;++i)bx_[I(i,j)]=Beta(Phi(i,j),Phi(i+1,j));
        for(int j=1;j<n_;++j)for(int i=1;i<=n_;++i)by_[I(i,j)]=Beta(Phi(i,j),Phi(i,j+1));
    }
    void ScalarBoundary(std::vector<double>& a) const {
        for(int k=1;k<=n_;++k){a[I(0,k)]=a[I(1,k)];a[I(n_+1,k)]=a[I(n_,k)];a[I(k,0)]=a[I(k,1)];a[I(k,n_+1)]=a[I(k,n_)];}
        a[I(0,0)]=a[I(1,1)];a[I(0,n_+1)]=a[I(1,n_)];a[I(n_+1,0)]=a[I(n_,1)];a[I(n_+1,n_+1)]=a[I(n_,n_)];
    }
    void VelocityBoundary() {
        for(int j=0;j<=n_+1;++j)u_[I(0,j)]=u_[I(n_,j)]=u_[I(n_+1,j)]=0;
        for(int i=1;i<n_;++i){u_[I(i,0)]=u_[I(i,1)];u_[I(i,n_+1)]=u_[I(i,n_)];}
        for(int i=0;i<=n_+1;++i)v_[I(i,0)]=v_[I(i,n_)]=v_[I(i,n_+1)]=0;
        for(int j=1;j<n_;++j){v_[I(0,j)]=v_[I(1,j)];v_[I(n_+1,j)]=v_[I(n_,j)];}
    }
    double Bilinear(const std::vector<double>& a,double x,double y) const {
        x=Clamp(x,0,double(n_+1));y=Clamp(y,0,double(n_+1));
        const int i=static_cast<int>(x),j=static_cast<int>(y),r=(std::min)(n_+1,i+1),t=(std::min)(n_+1,j+1);const double sx=x-i,sy=y-j;
        return (1-sy)*((1-sx)*a[I(i,j)]+sx*a[I(r,j)])+sy*((1-sx)*a[I(i,t)]+sx*a[I(r,t)]);
    }
    Vec2 SampleVelocityFrom(const std::vector<double>& u,const std::vector<double>& v,Vec2 p) const {
        if(!Finite(p))return {};p.x=Clamp(p.x,0,1);p.y=Clamp(p.y,0,1);
        return {Bilinear(u,p.x*n_,p.y*n_+0.5),Bilinear(v,p.x*n_+0.5,p.y*n_)};
    }
    double Div(int i,int j) const { return (u_[I(i,j)]-u_[I(i-1,j)]+v_[I(i,j)]-v_[I(i,j-1)])/h_; }
    void ZeroMean(std::vector<double>& a) const {
        double mean=0;for(int j=1;j<=n_;++j)for(int i=1;i<=n_;++i)mean+=a[I(i,j)];mean/=n_*n_;
        for(int j=1;j<=n_;++j)for(int i=1;i<=n_;++i)a[I(i,j)]-=mean;
    }
    double Dot(const std::vector<double>& a,const std::vector<double>& b) const {
        double value=0;for(int j=1;j<=n_;++j)for(int i=1;i<=n_;++i)value+=a[I(i,j)]*b[I(i,j)];return value;
    }
    void Apply(const std::vector<double>& x,std::vector<double>& out) const {
        for(int j=1;j<=n_;++j)for(int i=1;i<=n_;++i) {
            const auto k=I(i,j);double value=0;
            if(i>1)value+=bx_[I(i-1,j)]*(x[k]-x[I(i-1,j)]);
            if(i<n_)value+=bx_[I(i,j)]*(x[k]-x[I(i+1,j)]);
            if(j>1)value+=by_[I(i,j-1)]*(x[k]-x[I(i,j-1)]);
            if(j<n_)value+=by_[I(i,j)]*(x[k]-x[I(i,j+1)]);
            out[k]=value/(h_*h_);
        }
    }
    // Modified IC(0) keeps the original lower-neighbor sparsity. It approximates the
    // pressure matrix only inside PCG; A, the convergence tolerance and the
    // final explicitly recomputed residual remain unchanged.
    void BuildModifiedIncompleteCholesky(const std::vector<double>& diagonal) {
        auto& inverse=work_.icInverse;auto& left=work_.icLeft;auto& bottom=work_.icBottom;
        const double inverseH2=1.0/(h_*h_);
        for(int j=1;j<=n_;++j)for(int i=1;i<=n_;++i) {
            const auto k=I(i,j);
            left[k]=i>1?-bx_[I(i-1,j)]*inverseH2*inverse[I(i-1,j)]:0.0;
            bottom[k]=j>1?-by_[I(i,j-1)]*inverseH2*inverse[I(i,j-1)]:0.0;
            const double droppedLeft=i>1&&j<n_?bx_[I(i-1,j)]*by_[I(i-1,j)]*inverseH2*inverseH2*inverse[I(i-1,j)]*inverse[I(i-1,j)]:0.0;
            const double droppedBottom=j>1&&i<n_?by_[I(i,j-1)]*bx_[I(i,j-1)]*inverseH2*inverseH2*inverse[I(i,j-1)]*inverse[I(i,j-1)]:0.0;
            double reduced=diagonal[k]-left[k]*left[k]-bottom[k]*bottom[k]-0.97*(droppedLeft+droppedBottom);
            // A diagonal fallback protects the auxiliary factor at extreme
            // density contrasts and the Neumann nullspace; it does not modify A.
            if(reduced<0.25*diagonal[k])reduced=diagonal[k];
            inverse[k]=1.0/std::sqrt(reduced);
        }
    }
    void ApplyPreconditioner(const std::vector<double>& residual,std::vector<double>& result) {
        const auto& inverse=work_.icInverse;const auto& left=work_.icLeft;const auto& bottom=work_.icBottom;
        auto& forward=work_.forward;
        for(int j=1;j<=n_;++j)for(int i=1;i<=n_;++i) {
            const auto k=I(i,j);double value=residual[k];
            if(i>1)value-=left[k]*forward[I(i-1,j)];
            if(j>1)value-=bottom[k]*forward[I(i,j-1)];
            forward[k]=value*inverse[k];
        }
        for(int j=n_;j>=1;--j)for(int i=n_;i>=1;--i) {
            const auto k=I(i,j);double value=forward[k];
            if(i<n_)value-=left[I(i+1,j)]*result[I(i+1,j)];
            if(j<n_)value-=bottom[I(i,j+1)]*result[I(i,j+1)];
            result[k]=value*inverse[k];
        }
    }
    double Fraction(double phi) const {
        const double eps=h_;if(phi<=-eps)return 1;if(phi>=eps)return 0;
        const double pi=3.14159265358979323846;
        return 0.5-phi/(2*eps)-std::sin(pi*phi/eps)/(2*pi);
    }
    double AreaAtOffset(double shift) const {
        double area=0;for(int j=1;j<=n_;++j)for(int i=1;i<=n_;++i)area+=Fraction(phi_[I(i,j)]+shift);return area*h_*h_;
    }
    void Reinitialize() {
        auto& d=work_.distance;auto& fixed=work_.fixed;
        std::fill(d.begin(),d.end(),8.0);std::fill(fixed.begin(),fixed.end(),static_cast<unsigned char>(0));
        for(int j=1;j<=n_;++j)for(int i=1;i<=n_;++i) {
            const auto k=I(i,j);const double p=phi_[k];
            for(const Vec2 offset:{Vec2{-1,0},Vec2{1,0},Vec2{0,-1},Vec2{0,1}}) {
                const int a=i+static_cast<int>(offset.x),b=j+static_cast<int>(offset.y);if(a<1||a>n_||b<1||b>n_)continue;
                const double neighbor=phi_[I(a,b)];
                if((p<0)!=(neighbor<0)) { d[k]=(std::min)(d[k],h_*std::abs(p)/(std::abs(p)+std::abs(neighbor)));fixed[k]=1; }
            }
            if(p==0){d[k]=0;fixed[k]=1;}
        }
        for(int iteration=0;iteration<3;++iteration)for(int sweep=0;sweep<4;++sweep) {
            const int sx=(sweep&1)?-1:1,sy=(sweep&2)?-1:1;
            for(int jj=0;jj<n_;++jj)for(int ii=0;ii<n_;++ii) {
                const int i=sx>0?ii+1:n_-ii,j=sy>0?jj+1:n_-jj;const auto k=I(i,j);if(fixed[k])continue;
                const double a=(std::min)(i>1?d[I(i-1,j)]:8.0,i<n_?d[I(i+1,j)]:8.0);
                const double b=(std::min)(j>1?d[I(i,j-1)]:8.0,j<n_?d[I(i,j+1)]:8.0);
                const double delta=std::abs(a-b);const double next=delta>=h_?(std::min)(a,b)+h_:0.5*(a+b+std::sqrt(2*h_*h_-delta*delta));
                d[k]=(std::min)(d[k],next);
            }
        }
        for(int j=1;j<=n_;++j)for(int i=1;i<=n_;++i)phi_[I(i,j)]=(phi_[I(i,j)]<0?-1:1)*d[I(i,j)];
        ScalarBoundary(phi_);
    }
    void CorrectArea() {
        rawAreaError_=(WaterArea()-initialArea_)/(std::max)(initialArea_,1e-12);
        areaShift_=0;
        if(std::abs(rawAreaError_)<1e-12)return;
        // For every permitted shift, |phi| >= 1.5*h has a constant fraction.
        // Evaluate only the remaining band in the unchanged 40 bisections.
        auto& band=work_.areaBand;band.clear();double fullWater=0;
        for(int j=1;j<=n_;++j)for(int i=1;i<=n_;++i) {
            const auto k=I(i,j);const double p=phi_[k];
            if(p<=-1.5*h_)fullWater+=1.0;
            else if(p<1.5*h_)band.push_back(k);
        }
        const auto bandArea=[&](double shift) {
            double area=fullWater;for(const auto k:band)area+=Fraction(phi_[k]+shift);
            return area*h_*h_;
        };
        double lo=-0.5*h_,hi=0.5*h_;
        for(int iteration=0;iteration<40;++iteration){const double mid=0.5*(lo+hi);if(bandArea(mid)>initialArea_)lo=mid;else hi=mid;}
        areaShift_=0.5*(lo+hi);
        for(int j=1;j<=n_;++j)for(int i=1;i<=n_;++i)phi_[I(i,j)]+=areaShift_;
        ScalarBoundary(phi_);
    }
    void Substep(double dt) {
        auto& oldU=work_.oldU;auto& oldV=work_.oldV;auto& oldPhi=work_.oldPhi;auto& oldDye=work_.oldDye;
        oldU=u_;oldV=v_;oldPhi=phi_;oldDye=dye_;
        // Advect material with the previous projected velocity (RK2 backtrace).
        for(int j=1;j<=n_;++j)for(int i=1;i<=n_;++i) {
            const Vec2 p{(i-0.5)/n_,(j-0.5)/n_};
            const Vec2 mid=p-SampleVelocityFrom(oldU,oldV,p)*(0.5*dt/4.0);
            const Vec2 back=p-SampleVelocityFrom(oldU,oldV,mid)*(dt/4.0);
            phi_[I(i,j)]=Bilinear(oldPhi,back.x*n_+0.5,back.y*n_+0.5);
        }
        ScalarBoundary(phi_);Reinitialize();CorrectArea();
        for(int j=1;j<=n_;++j)for(int i=1;i<n_;++i) {
            const Vec2 p{double(i)/n_,(j-0.5)/n_};const Vec2 back=p-SampleVelocityFrom(oldU,oldV,p)*(dt/4.0);
            u_[I(i,j)]=SampleVelocityFrom(oldU,oldV,back).x+dt*gravity_.x;
        }
        for(int j=1;j<n_;++j)for(int i=1;i<=n_;++i) {
            const Vec2 p{(i-0.5)/n_,double(j)/n_};const Vec2 back=p-SampleVelocityFrom(oldU,oldV,p)*(dt/4.0);
            v_[I(i,j)]=SampleVelocityFrom(oldU,oldV,back).y+dt*gravity_.y;
        }
        VelocityBoundary();
        if(jet_){AddVelocity(0.84,0.37,-25*dt,4*dt,0.11);AddDye(0.84,0.37,3*dt,0.09);}
        if(projection_)ProjectOnly(dt);
        else {BuildCoefficients();std::fill(pressure_.begin(),pressure_.end(),0.0);divergenceBeforeProjection=divergenceAfterProjection=DivergenceRms();pressureResidual=0;pressureIterations=0;pressureConverged=false;}
        // Passive gas tracer, with material clipping separate from pressure.
        for(int j=1;j<=n_;++j)for(int i=1;i<=n_;++i) {
            const Vec2 p{(i-0.5)/n_,(j-0.5)/n_};const Vec2 back=p-SampleVelocity(p)*(dt/4.0);
            const double transported=Bilinear(oldDye,back.x*n_+0.5,back.y*n_+0.5);
            const double injected=jet_?3*dt*Splat(p.x,p.y,0.84,0.37,0.09):0;
            dye_[I(i,j)]=Phi(i,j)<0?0:Clamp(transported+injected,0,10);
        }
        ScalarBoundary(dye_);
    }
};

} // namespace fieldlab
#endif
