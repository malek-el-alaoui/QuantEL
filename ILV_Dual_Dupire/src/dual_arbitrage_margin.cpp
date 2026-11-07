#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

static constexpr double S0 = 100.0;
static constexpr double SIG0 = 0.20;
static constexpr double SIG_REF = 0.35;
static constexpr double KAPPA = 0.03;
static constexpr double TAU = 0.09;
static constexpr double T0 = 0.50;
static constexpr double MAX_SHOCK = 0.0125;
static constexpr double SAFETY = 0.25;
static constexpr double VAR_FLOOR = 0.0025;
static constexpr double VAR_CAP = 1.0;

static const std::vector<double> MATURITIES = {0.25, 0.50, 0.75, 1.00};

static double normalCdf(double x) {
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

static double bsCall(double K, double T, double sigma) {
    if (T <= 1e-14) return std::max(S0 - K, 0.0);
    const double sw = std::max(sigma * std::sqrt(T), 1e-14);
    const double d1 = (std::log(S0 / K) + 0.5 * sw * sw) / sw;
    const double d2 = d1 - sw;
    return S0 * normalCdf(d1) - K * normalCdf(d2);
}

static double bumpShape(double K, double T) {
    const double x = std::log(K / S0);
    return std::exp(-0.5 * (x / KAPPA) * (x / KAPPA))
         * std::exp(-0.5 * ((T - T0) / TAU) * ((T - T0) / TAU));
}

static double impliedVol(double K, double T, double shock) {
    return SIG0 + shock * bumpShape(K, T);
}

static double marketCall(double K, double T, double shock) {
    return bsCall(K, T, impliedVol(K, T, shock));
}

static double referenceCall(double K, double T) {
    return bsCall(K, T, SIG_REF);
}

struct Derivs { double ckk, ct; };

static Derivs derivativesMarket(double K, double T, double shock, double hk=0.05, double ht=1e-4) {
    const double c0 = marketCall(K, T, shock);
    const double ckk = (marketCall(K + hk, T, shock) - 2.0 * c0 + marketCall(K - hk, T, shock)) / (hk * hk);
    double ct;
    if (T > ht) ct = (marketCall(K, T + ht, shock) - marketCall(K, T - ht, shock)) / (2.0 * ht);
    else ct = (marketCall(K, T + ht, shock) - c0) / ht;
    return {ckk, ct};
}

static Derivs derivativesRef(double K, double T, double hk=0.05, double ht=1e-4) {
    const double c0 = referenceCall(K, T);
    const double ckk = (referenceCall(K + hk, T) - 2.0 * c0 + referenceCall(K - hk, T)) / (hk * hk);
    double ct;
    if (T > ht) ct = (referenceCall(K, T + ht) - referenceCall(K, T - ht)) / (2.0 * ht);
    else ct = (referenceCall(K, T + ht) - c0) / ht;
    return {ckk, ct};
}

static std::tuple<double,double,double> betaRequired(double shock) {
    double bg = 0.0, bt = 0.0;
    const int nK = 240, nT = 87;
    for (int jt = 0; jt <= nT; ++jt) {
        const double T = 0.08 + (0.95 - 0.08) * jt / nT;
        for (int ik = 0; ik <= nK; ++ik) {
            const double K = 70.0 + (130.0 - 70.0) * ik / nK;
            const auto c = derivativesMarket(K, T, shock);
            const auto r = derivativesRef(K, T);
            if (r.ckk > 1e-8) bg = std::max(bg, std::max(-c.ckk / r.ckk, 0.0));
            if (r.ct > 1e-6) bt = std::max(bt, std::max(-c.ct / r.ct, 0.0));
        }
    }
    return {std::max(bg,bt), bg, bt};
}

static std::pair<double,double> robustBeta() {
    double bstar = 0.0;
    for (int i = 0; i <= 25; ++i) {
        const double a = MAX_SHOCK * i / 25.0;
        auto [b,bg,bt] = betaRequired(a);
        bstar = std::max(bstar, b);
    }
    return {bstar, (1.0 + SAFETY) * bstar};
}

static std::vector<double> linspace(double a, double b, int n) {
    std::vector<double> x(n);
    for (int i = 0; i < n; ++i) x[i] = a + (b-a) * i / (n-1.0);
    return x;
}

static double interp(const std::vector<double>& x, const std::vector<double>& y, double z) {
    if (z <= x.front()) return y.front();
    if (z >= x.back()) return y.back();
    auto it = std::lower_bound(x.begin(), x.end(), z);
    size_t j = static_cast<size_t>(it - x.begin());
    size_t i = j - 1;
    double a = (z - x[i]) / (x[j] - x[i]);
    return y[i] + a * (y[j] - y[i]);
}

struct LVSlice {
    std::vector<double> plus;
    std::vector<double> classic;
    std::vector<double> raw;
};

static const std::vector<double> KLV = linspace(55.0,145.0,361);

static LVSlice localVarianceAt(double T, double shock, double beta) {
    LVSlice out;
    out.plus.resize(KLV.size());
    out.classic.resize(KLV.size());
    out.raw.resize(KLV.size(), std::numeric_limits<double>::quiet_NaN());
    for (size_t i=0;i<KLV.size();++i) {
        double K=KLV[i];
        auto c=derivativesMarket(K,T,shock);
        auto r=derivativesRef(K,T);
        double den=c.ckk+beta*r.ckk;
        double num=c.ct+beta*r.ct;
        double vp=SIG_REF*SIG_REF;
        if (r.ckk>1e-9 && den>1e-10 && num>1e-10)
            vp=2.0*num/(K*K*den);
        vp=std::min(VAR_CAP,std::max(VAR_FLOOR,vp));
        out.plus[i]=vp;

        bool valid=(c.ckk>1e-10 && c.ct>1e-10);
        double vc=VAR_FLOOR;
        if (valid) {
            double vr=2.0*c.ct/(K*K*c.ckk);
            out.raw[i]=vr;
            vc=std::min(VAR_CAP,std::max(VAR_FLOOR,vr));
        }
        if (K<65.0 || K>135.0) vc=SIG0*SIG0;
        out.classic[i]=vc;
    }
    return out;
}

struct LVGrid {
    std::vector<std::vector<double>> plus;
    std::vector<std::vector<double>> classic;
};

static LVGrid precomputeLV(double shock, double beta, int nsteps) {
    LVGrid g;
    g.plus.resize(nsteps);
    g.classic.resize(nsteps);
    for (int j=1;j<=nsteps;++j) {
        double T=double(j)/nsteps;
        auto s=localVarianceAt(T,shock,beta);
        g.plus[j-1]=std::move(s.plus);
        g.classic[j-1]=std::move(s.classic);
    }
    return g;
}

static std::vector<double> solveTridiagonal(
    const std::vector<double>& lower,
    const std::vector<double>& diag,
    const std::vector<double>& upper,
    const std::vector<double>& rhs)
{
    int n=static_cast<int>(diag.size());
    std::vector<double> c=upper,d=rhs,b=diag;
    for(int i=1;i<n;++i){
        double m=lower[i-1]/b[i-1];
        b[i]-=m*c[i-1];
        d[i]-=m*d[i-1];
    }
    std::vector<double> x(n);
    x[n-1]=d[n-1]/b[n-1];
    for(int i=n-2;i>=0;--i) x[i]=(d[i]-c[i]*x[i+1])/b[i];
    return x;
}

static std::map<double,std::vector<double>> forwardPDE(const std::vector<std::vector<double>>& lv) {
    const auto KPDE=linspace(0.0,200.0,801);
    const double dK=KPDE[1]-KPDE[0];
    const int nsteps=static_cast<int>(lv.size());
    const double dt=1.0/nsteps;
    std::vector<double> C(KPDE.size());
    for(size_t i=0;i<KPDE.size();++i) C[i]=std::max(S0-KPDE[i],0.0);
    std::map<int,double> target;
    for(double T:MATURITIES) target[int(std::round(T*nsteps))]=T;
    std::map<double,std::vector<double>> saved;

    int nint=static_cast<int>(KPDE.size())-2;
    for(int j=1;j<=nsteps;++j){
        std::vector<double> lower(nint-1),diag(nint),upper(nint-1),rhs(nint);
        std::vector<double> lam(nint);
        for(int i=0;i<nint;++i){
            double K=KPDE[i+1];
            double v=interp(KLV,lv[j-1],K);
            lam[i]=0.5*v*K*K*dt/(dK*dK);
            diag[i]=1.0+2.0*lam[i];
            rhs[i]=C[i+1];
        }
        for(int i=0;i<nint-1;++i){ upper[i]=-lam[i]; lower[i]=-lam[i+1]; }
        rhs.front() += lam.front()*S0;
        auto x=solveTridiagonal(lower,diag,upper,rhs);
        C.front()=S0; C.back()=0.0;
        for(int i=0;i<nint;++i) C[i+1]=x[i];
        auto it=target.find(j);
        if(it!=target.end()) saved[it->second]=C;
    }
    return saved;
}

static double autocallPDE(const std::vector<std::vector<double>>& lv) {
    const auto SGRID=linspace(0.0,250.0,1001);
    const double dS=SGRID[1]-SGRID[0];
    const int nsteps=static_cast<int>(lv.size());
    const double dt=1.0/nsteps;
    std::vector<double> V(SGRID.size());
    for(size_t i=0;i<SGRID.size();++i){
        double S=SGRID[i];
        if(S>=100.0) V[i]=108.0;
        else if(S>=70.0) V[i]=100.0;
        else V[i]=100.0*S/S0;
    }
    std::map<int,int> obs={{nsteps/4,1},{nsteps/2,2},{3*nsteps/4,3}};
    int nint=static_cast<int>(SGRID.size())-2;

    for(int j=nsteps;j>=1;--j){
        std::vector<double> lower(nint-1),diag(nint),upper(nint-1),rhs(nint),lam(nint);
        for(int i=0;i<nint;++i){
            double S=SGRID[i+1];
            double v=interp(KLV,lv[j-1],S);
            lam[i]=0.5*v*S*S*dt/(dS*dS);
            diag[i]=1.0+2.0*lam[i];
            rhs[i]=V[i+1];
        }
        for(int i=0;i<nint-1;++i){ upper[i]=-lam[i]; lower[i]=-lam[i+1]; }
        double tprev=double(j-1)/nsteps;
        int nextI = (tprev<=0.25+1e-12?1:(tprev<=0.50+1e-12?2:(tprev<=0.75+1e-12?3:4)));
        double high=100.0*(1.0+0.02*nextI);
        rhs.back() += lam.back()*high;
        auto x=solveTridiagonal(lower,diag,upper,rhs);
        V.front()=0.0; V.back()=high;
        for(int i=0;i<nint;++i) V[i+1]=x[i];
        int step=j-1;
        auto it=obs.find(step);
        if(it!=obs.end()){
            double redeem=100.0*(1.0+0.02*it->second);
            for(size_t i=0;i<SGRID.size();++i) if(SGRID[i]>=100.0) V[i]=redeem;
        }
    }
    return interp(SGRID,V,S0);
}

static int quoteButterflyCount(double shock, double& worstGamma) {
    worstGamma=std::numeric_limits<double>::infinity();
    int count=0;
    const double h=1.0;
    for(double T:MATURITIES){
        for(double K=81.0;K<=119.0+1e-12;K+=1.0){
            double g=(marketCall(K-h,T,shock)-2.0*marketCall(K,T,shock)+marketCall(K+h,T,shock))/(h*h);
            worstGamma=std::min(worstGamma,g);
            if(g<0.0) ++count;
        }
    }
    return count;
}

static std::pair<double,double> autocallPrices(double shock,double beta,double vminus) {
    auto g=precomputeLV(shock,beta,252);
    double pplus=autocallPDE(g.plus);
    double pdual=(1.0+beta)*pplus-beta*vminus;
    double pclassic=autocallPDE(g.classic);
    return {pclassic,pdual};
}

int main(){
    try{
        auto [bstar,beta]=robustBeta();
        std::cout<<std::setprecision(12);
        std::cout<<"beta_star="<<bstar<<"\n";
        std::cout<<"beta_margin="<<beta<<"\n";

        double worstGamma;
        int butterflies=quoteButterflyCount(MAX_SHOCK,worstGamma);
        std::cout<<"butterfly_violations="<<butterflies<<" worst_gamma="<<worstGamma<<"\n";

        // Vanilla forward-PDE comparison at the arbitrable 1.25vp surface.
        auto lv=precomputeLV(MAX_SHOCK,beta,1000);
        std::vector<std::vector<double>> minus(1000,std::vector<double>(KLV.size(),SIG_REF*SIG_REF));
        auto cp=forwardPDE(lv.plus);
        auto cm=forwardPDE(minus);
        auto cc=forwardPDE(lv.classic);
        const auto KPDE=linspace(0.0,200.0,801);
        long double ssed=0.0L,ssec=0.0L; int n=0;
        double maxd=0.0,maxc=0.0;
        for(double T:MATURITIES){
            for(double K=80.0;K<=120.0+1e-12;K+=1.0){
                double m=marketCall(K,T,MAX_SHOCK);
                double dual=(1.0+beta)*interp(KPDE,cp[T],K)-beta*interp(KPDE,cm[T],K);
                double clas=interp(KPDE,cc[T],K);
                double ed=dual-m, ec=clas-m;
                ssed+=ed*ed; ssec+=ec*ec; ++n;
                maxd=std::max(maxd,std::abs(ed)); maxc=std::max(maxc,std::abs(ec));
            }
        }
        std::cout<<"dual_rmse="<<std::sqrt(double(ssed/n))<<" classic_rmse="<<std::sqrt(double(ssec/n))<<"\n";
        std::cout<<"dual_max_abs="<<maxd<<" classic_max_abs="<<maxc<<"\n";

        // Autocall vega step-size test near first convexity loss.
        std::vector<std::vector<double>> minus252(252,std::vector<double>(KLV.size(),SIG_REF*SIG_REF));
        double vminus=autocallPDE(minus252);
        const double center=0.0085;
        const std::vector<double> hs={0.0001,0.0002,0.00025,0.0005,0.0010};
        std::ofstream out("autocall_vega_step_stability_cpp.csv");
        out<<"half_bump,classical_vega_per_1vp,dual_vega_per_1vp\n";
        double minVC=1e100,maxVC=-1e100,minVD=1e100,maxVD=-1e100;
        for(double h:hs){
            auto pm=autocallPrices(center-h,beta,vminus);
            auto pp=autocallPrices(center+h,beta,vminus);
            double vc=0.01*(pp.first-pm.first)/(2.0*h);
            double vd=0.01*(pp.second-pm.second)/(2.0*h);
            out<<std::setprecision(15)<<h<<","<<vc<<","<<vd<<"\n";
            minVC=std::min(minVC,vc); maxVC=std::max(maxVC,vc);
            minVD=std::min(minVD,vd); maxVD=std::max(maxVD,vd);
            std::cout<<"h="<<h<<" classic_vega_1vp="<<vc<<" dual_vega_1vp="<<vd<<"\n";
        }
        std::cout<<"classic_vega_range="<<(maxVC-minVC)<<" dual_vega_range="<<(maxVD-minVD)<<"\n";
        return 0;
    }catch(const std::exception& e){
        std::cerr<<"ERROR: "<<e.what()<<"\n";
        return 1;
    }
}
