#include <algorithm>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace fs = std::filesystem;

struct Quote { double T, K, price; };

struct Params {
    double beta{};
    std::vector<double> wPlus;
    std::vector<double> wMinus;
};

struct ScenarioDef {
    std::string name;
    std::string label;
    double betaTrue;
    std::vector<double> wPlusTrue;
    std::vector<double> wMinusTrue;
    double kMin;
    double kMax;
    double kStep;
    double betaGuess;
    double plusVarianceRateGuess;
    double minusVarianceRateGuess;
};

struct SurfaceGrid {
    std::vector<double> T;
    std::vector<double> K;
    std::vector<std::vector<double>> C; // C[j][i]
};

struct ClassicNode {
    double dCdt{};
    double gamma{};
    double rawVariance{std::numeric_limits<double>::quiet_NaN()};
    double regularizedVariance{};
    bool gammaInvalid{};
    bool varianceInvalid{};
    bool clipped{};
};

struct ClassicGrid {
    SurfaceGrid surface;
    std::vector<std::vector<ClassicNode>> node;
};

struct ClassicSolution {
    std::map<double, std::vector<double>> pricesAtMaturity;
};

static constexpr double S0 = 5000.0;
static const std::vector<double> MATURITIES = {0.10, 0.25, 0.50, 1.00};
static constexpr double CLASSIC_VAR_FLOOR = 1e-6;
static constexpr double CLASSIC_VAR_CAP   = 1.0;   // 100% instantaneous vol squared cap.
static constexpr double CLASSIC_PDE_DK    = 25.0;
static constexpr double CLASSIC_PDE_KMAX  = 10000.0;
static constexpr double CLASSIC_PDE_DT    = 0.0005;

static double normalCdf(double x) {
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

static double normalPdf(double x) {
    static const double invSqrt2Pi = 1.0 / std::sqrt(2.0 * M_PI);
    return invSqrt2Pi * std::exp(-0.5 * x * x);
}

// Forward-neutralized demo: r=q=0. Black-Scholes call parameterized by total variance w(T).
static double callFromTotalVariance(double S, double K, double w) {
    if (w <= 1e-14) return std::max(S - K, 0.0);
    const double sw = std::sqrt(w);
    const double d1 = (std::log(S / K) + 0.5 * w) / sw;
    const double d2 = d1 - sw;
    return S * normalCdf(d1) - K * normalCdf(d2);
}

// With r=q=0, C_KK equals the terminal density for an ordinary positive branch.
static double strikeGammaFromTotalVariance(double S, double K, double w) {
    if (w <= 1e-14) return 0.0;
    const double sw = std::sqrt(w);
    const double d2 = (std::log(S / K) - 0.5 * w) / sw;
    return normalPdf(d2) / (K * sw);
}

static size_t maturityIndex(double T) {
    for (size_t j = 0; j < MATURITIES.size(); ++j)
        if (std::abs(T - MATURITIES[j]) < 1e-10) return j;
    throw std::runtime_error("Maturity not found in calibration nodes");
}

static double instantVariance(double T, const std::vector<double>& w) {
    if (T <= MATURITIES.front()) return w.front() / MATURITIES.front();
    for (size_t j = 1; j < MATURITIES.size(); ++j) {
        if (T <= MATURITIES[j] + 1e-14)
            return (w[j] - w[j - 1]) / (MATURITIES[j] - MATURITIES[j - 1]);
    }
    const size_t n = MATURITIES.size();
    return (w[n - 1] - w[n - 2]) / (MATURITIES[n - 1] - MATURITIES[n - 2]);
}

static double modelCall(double K, double T, const Params& p) {
    const size_t j = maturityIndex(T);
    const double alpha = 1.0 + p.beta;
    return alpha * callFromTotalVariance(S0, K, p.wPlus[j])
         - p.beta * callFromTotalVariance(S0, K, p.wMinus[j]);
}

static Params decode(const std::vector<double>& x) {
    if (x.size() != 9) throw std::runtime_error("Expected 9 transformed parameters");
    Params p;
    p.beta = std::exp(x[0]);
    p.wPlus.resize(4);
    p.wMinus.resize(4);
    double acc = 0.0;
    for (int i = 0; i < 4; ++i) {
        acc += std::exp(x[1 + i]);
        p.wPlus[i] = acc;
    }
    acc = 0.0;
    for (int i = 0; i < 4; ++i) {
        acc += std::exp(x[5 + i]);
        p.wMinus[i] = acc;
    }
    return p;
}

static double objective(const std::vector<double>& x, const std::vector<Quote>& quotes) {
    for (double z : x) if (!std::isfinite(z) || z > 4.0 || z < -25.0) return 1e30;
    Params p = decode(x);
    if (!(p.beta > 0.0 && p.beta < 3.0)) return 1e30;
    if (p.wPlus.back() > 1.0 || p.wMinus.back() > 1.0) return 1e30;
    long double sum = 0.0L;
    for (const auto& q : quotes) {
        const double pred = modelCall(q.K, q.T, p);
        const double scale = std::max(10.0, std::abs(q.price));
        const double e = (pred - q.price) / scale;
        sum += static_cast<long double>(e) * e;
    }
    return static_cast<double>(sum / quotes.size());
}

static std::vector<double> residualVector(const std::vector<double>& x, const std::vector<Quote>& quotes) {
    Params p = decode(x);
    std::vector<double> r; r.reserve(quotes.size());
    for (const auto& q : quotes) {
        const double pred = modelCall(q.K, q.T, p);
        const double scale = std::max(10.0, std::abs(q.price));
        r.push_back((pred - q.price) / scale);
    }
    return r;
}

static bool solveLinearSystem(std::vector<std::vector<double>> A, std::vector<double> b,
                              std::vector<double>& x) {
    const int n = static_cast<int>(b.size());
    for (int i = 0; i < n; ++i) {
        int piv = i;
        for (int k = i + 1; k < n; ++k)
            if (std::abs(A[k][i]) > std::abs(A[piv][i])) piv = k;
        if (std::abs(A[piv][i]) < 1e-18) return false;
        if (piv != i) { std::swap(A[piv], A[i]); std::swap(b[piv], b[i]); }
        const double d = A[i][i];
        for (int j = i; j < n; ++j) A[i][j] /= d;
        b[i] /= d;
        for (int k = 0; k < n; ++k) {
            if (k == i) continue;
            const double f = A[k][i];
            if (f == 0.0) continue;
            for (int j = i; j < n; ++j) A[k][j] -= f * A[i][j];
            b[k] -= f * b[i];
        }
    }
    x = b;
    return true;
}

static std::vector<double> levenbergMarquardt(std::vector<double> x,
                                              const std::vector<Quote>& quotes,
                                              int maxIter = 250) {
    const int n = static_cast<int>(x.size());
    const int m = static_cast<int>(quotes.size());
    double lambda = 1e-4;
    double f = objective(x, quotes);
    for (int iter = 0; iter < maxIter; ++iter) {
        const auto r = residualVector(x, quotes);
        std::vector<std::vector<double>> J(m, std::vector<double>(n, 0.0));
        for (int j = 0; j < n; ++j) {
            const double h = 1e-5 * std::max(1.0, std::abs(x[j]));
            auto xp = x, xm = x;
            xp[j] += h; xm[j] -= h;
            const auto rp = residualVector(xp, quotes);
            const auto rm = residualVector(xm, quotes);
            for (int i = 0; i < m; ++i) J[i][j] = (rp[i] - rm[i]) / (2.0 * h);
        }
        std::vector<std::vector<double>> A(n, std::vector<double>(n, 0.0));
        std::vector<double> g(n, 0.0);
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                g[j] += J[i][j] * r[i];
                for (int k = 0; k < n; ++k) A[j][k] += J[i][j] * J[i][k];
            }
        }
        for (int j = 0; j < n; ++j) {
            A[j][j] += lambda * std::max(1e-12, A[j][j]);
            g[j] = -g[j];
        }
        std::vector<double> dx;
        if (!solveLinearSystem(A, g, dx)) { lambda *= 10.0; continue; }
        auto trial = x;
        double normDx = 0.0;
        for (int j = 0; j < n; ++j) {
            trial[j] += dx[j];
            normDx = std::max(normDx, std::abs(dx[j]));
        }
        const double ft = objective(trial, quotes);
        if (ft < f) {
            x = trial; f = ft; lambda = std::max(1e-12, lambda / 3.0);
            if (normDx < 1e-11 || f < 1e-24) break;
        } else lambda = std::min(1e12, lambda * 10.0);
    }
    return x;
}

struct Vertex { std::vector<double> x; double f; };

static std::vector<double> nelderMead(const std::vector<double>& x0,
                                      const std::vector<Quote>& quotes,
                                      int maxIter = 6000) {
    const int n = static_cast<int>(x0.size());
    std::vector<Vertex> simplex(n + 1);
    simplex[0] = {x0, objective(x0, quotes)};
    for (int i = 0; i < n; ++i) {
        auto x = x0; x[i] += (i == 0 ? 0.15 : 0.12);
        simplex[i + 1] = {x, objective(x, quotes)};
    }
    const double alpha = 1.0, gamma = 1.0 + 2.0 / n;
    const double rho = 0.75 - 1.0 / (2.0 * n), sigma = 1.0 - 1.0 / n;
    for (int iter = 0; iter < maxIter; ++iter) {
        std::sort(simplex.begin(), simplex.end(), [](const Vertex& a, const Vertex& b){return a.f < b.f;});
        double fSpread = std::abs(simplex.back().f - simplex.front().f), xSpread = 0.0;
        for (int i = 1; i <= n; ++i)
            for (int d = 0; d < n; ++d)
                xSpread = std::max(xSpread, std::abs(simplex[i].x[d] - simplex[0].x[d]));
        if (fSpread < 1e-22 && xSpread < 1e-10) break;
        std::vector<double> centroid(n, 0.0);
        for (int i = 0; i < n; ++i) for (int d = 0; d < n; ++d) centroid[d] += simplex[i].x[d];
        for (double& z : centroid) z /= n;
        std::vector<double> xr(n);
        for (int d = 0; d < n; ++d) xr[d] = centroid[d] + alpha * (centroid[d] - simplex[n].x[d]);
        const double fr = objective(xr, quotes);
        if (fr < simplex[0].f) {
            std::vector<double> xe(n);
            for (int d = 0; d < n; ++d) xe[d] = centroid[d] + gamma * (xr[d] - centroid[d]);
            const double fe = objective(xe, quotes);
            simplex[n] = (fe < fr) ? Vertex{xe, fe} : Vertex{xr, fr};
        } else if (fr < simplex[n - 1].f) simplex[n] = {xr, fr};
        else {
            const bool outside = fr < simplex[n].f;
            std::vector<double> xc(n);
            if (outside) for (int d = 0; d < n; ++d) xc[d] = centroid[d] + rho * (xr[d] - centroid[d]);
            else for (int d = 0; d < n; ++d) xc[d] = centroid[d] + rho * (simplex[n].x[d] - centroid[d]);
            const double fc = objective(xc, quotes);
            if (fc < (outside ? fr : simplex[n].f)) simplex[n] = {xc, fc};
            else {
                for (int i = 1; i <= n; ++i) {
                    for (int d = 0; d < n; ++d)
                        simplex[i].x[d] = simplex[0].x[d] + sigma * (simplex[i].x[d] - simplex[0].x[d]);
                    simplex[i].f = objective(simplex[i].x, quotes);
                }
            }
        }
    }
    std::sort(simplex.begin(), simplex.end(), [](const Vertex& a, const Vertex& b){return a.f < b.f;});
    return simplex.front().x;
}

static void ensureDirs() {
    fs::create_directories("data");
    fs::create_directories("output/regular");
    fs::create_directories("output/violent");
    fs::create_directories("plots/regular");
    fs::create_directories("plots/violent");
    fs::create_directories("paper");
}

static void generateSyntheticMarket(const ScenarioDef& s) {
    const std::string path = "data/eurostoxx_" + s.name + "_market.csv";
    const double alpha = 1.0 + s.betaTrue;
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot write " + path);
    out << "maturity,strike,market_price\n" << std::setprecision(15);
    for (size_t j = 0; j < MATURITIES.size(); ++j) {
        for (double K = s.kMin; K <= s.kMax + 1e-9; K += s.kStep) {
            const double c = alpha * callFromTotalVariance(S0, K, s.wPlusTrue[j])
                           - s.betaTrue * callFromTotalVariance(S0, K, s.wMinusTrue[j]);
            out << MATURITIES[j] << ',' << K << ',' << c << '\n';
        }
    }
}

static std::vector<Quote> loadMarket(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot read market CSV: " + path);
    std::string line; std::getline(in, line);
    std::vector<Quote> q;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line); std::string a,b,c;
        std::getline(ss,a,','); std::getline(ss,b,','); std::getline(ss,c,',');
        q.push_back({std::stod(a),std::stod(b),std::stod(c)});
    }
    return q;
}

static SurfaceGrid makeSurface(const std::vector<Quote>& quotes) {
    SurfaceGrid s; s.T = MATURITIES;
    std::map<double,double> first;
    for (const auto& q : quotes) if (std::abs(q.T - MATURITIES.front()) < 1e-12) first[q.K] = q.price;
    for (const auto& kv : first) s.K.push_back(kv.first);
    s.C.assign(s.T.size(), std::vector<double>(s.K.size(), 0.0));
    std::map<std::pair<double,double>,double> m;
    for (const auto& q : quotes) m[{q.T,q.K}] = q.price;
    for (size_t j=0;j<s.T.size();++j) for (size_t i=0;i<s.K.size();++i)
        s.C[j][i] = m.at({s.T[j],s.K[i]});
    return s;
}

static Params calibrateExtended(const ScenarioDef& s, const std::vector<Quote>& quotes) {
    const std::vector<double> dt = {0.10,0.15,0.25,0.50};
    std::vector<double> x0(9);
    x0[0] = std::log(s.betaGuess);
    for (int i=0;i<4;++i) x0[1+i] = std::log(s.plusVarianceRateGuess * dt[i]);
    for (int i=0;i<4;++i) x0[5+i] = std::log(s.minusVarianceRateGuess * dt[i]);
    auto x = nelderMead(x0, quotes);
    x = levenbergMarquardt(x, quotes);
    return decode(x);
}

static ClassicGrid computeClassicGrid(const SurfaceGrid& s) {
    ClassicGrid g; g.surface = s;
    const size_t nt=s.T.size(), nk=s.K.size();
    g.node.assign(nt, std::vector<ClassicNode>(nk));
    for (size_t j=0;j<nt;++j) {
        for (size_t i=0;i<nk;++i) {
            ClassicNode n;
            const double prevC = (j==0) ? std::max(S0-s.K[i],0.0) : s.C[j-1][i];
            const double prevT = (j==0) ? 0.0 : s.T[j-1];
            n.dCdt = (s.C[j][i]-prevC)/(s.T[j]-prevT);
            g.node[j][i] = n;
        }
        // centered gamma in the interior; copy nearest interior at boundaries.
        for (size_t i=1;i+1<nk;++i) {
            const double h1=s.K[i]-s.K[i-1], h2=s.K[i+1]-s.K[i];
            double gamma=std::numeric_limits<double>::quiet_NaN();
            if (std::abs(h1-h2)<1e-10)
                gamma=(s.C[j][i-1]-2.0*s.C[j][i]+s.C[j][i+1])/(h1*h1);
            g.node[j][i].gamma=gamma;
        }
        g.node[j][0].gamma=g.node[j][1].gamma;
        g.node[j][nk-1].gamma=g.node[j][nk-2].gamma;
        for (size_t i=0;i<nk;++i) {
            auto& n=g.node[j][i];
            n.gammaInvalid = (!std::isfinite(n.gamma) || n.gamma <= 0.0);
            if (std::isfinite(n.gamma) && std::abs(n.gamma)>1e-16) {
                n.rawVariance=2.0*n.dCdt/(s.K[i]*s.K[i]*n.gamma);
            }
            n.varianceInvalid = (!std::isfinite(n.rawVariance) || n.rawVariance <= 0.0);
            double v=n.rawVariance;
            if (!std::isfinite(v)) v=CLASSIC_VAR_FLOOR;
            const double reg=std::clamp(v, CLASSIC_VAR_FLOOR, CLASSIC_VAR_CAP);
            n.regularizedVariance=reg;
            n.clipped = (!std::isfinite(n.rawVariance) || std::abs(reg-n.rawVariance)>1e-14);
        }
    }
    return g;
}

static double interp1(const std::vector<double>& x, const std::vector<double>& y, double z) {
    if (z<=x.front()) return y.front();
    if (z>=x.back()) return y.back();
    auto it=std::upper_bound(x.begin(),x.end(),z);
    size_t i=static_cast<size_t>(it-x.begin()-1);
    double a=(z-x[i])/(x[i+1]-x[i]);
    return y[i]+a*(y[i+1]-y[i]);
}

static double classicVarianceAt(const ClassicGrid& g, double t, double K) {
    size_t j=0;
    while (j+1<g.surface.T.size() && t>g.surface.T[j]+1e-14) ++j;
    std::vector<double> v(g.surface.K.size());
    for (size_t i=0;i<v.size();++i) v[i]=g.node[j][i].regularizedVariance;
    return interp1(g.surface.K,v,K);
}

static void solveTridiagonal(const std::vector<double>& a, const std::vector<double>& b,
                             const std::vector<double>& c, const std::vector<double>& d,
                             std::vector<double>& x) {
    const size_t n=d.size();
    std::vector<double> cp(n,0.0), dp(n,0.0);
    cp[0]=c[0]/b[0]; dp[0]=d[0]/b[0];
    for (size_t i=1;i<n;++i) {
        const double den=b[i]-a[i]*cp[i-1];
        cp[i]=(i+1<n)?c[i]/den:0.0;
        dp[i]=(d[i]-a[i]*dp[i-1])/den;
    }
    x.assign(n,0.0); x[n-1]=dp[n-1];
    for (size_t ii=n-1;ii-->0;) x[ii]=dp[ii]-cp[ii]*x[ii+1];
}

static ClassicSolution solveClassicForwardPDE(const ClassicGrid& g) {
    const int nK=static_cast<int>(std::round(CLASSIC_PDE_KMAX/CLASSIC_PDE_DK))+1;
    std::vector<double> K(nK), C(nK), Cnew(nK);
    for (int i=0;i<nK;++i) { K[i]=i*CLASSIC_PDE_DK; C[i]=std::max(S0-K[i],0.0); }
    ClassicSolution sol;
    const int steps=static_cast<int>(std::round(MATURITIES.back()/CLASSIC_PDE_DT));
    size_t nextMat=0;
    for (int step=1;step<=steps;++step) {
        const double t=step*CLASSIC_PDE_DT;
        const int n=nK-2;
        std::vector<double> lo(n,0.0),di(n,0.0),up(n,0.0),rhs(n,0.0),x;
        for (int m=0;m<n;++m) {
            const int i=m+1;
            const double v=classicVarianceAt(g,t,K[i]);
            const double lambda=CLASSIC_PDE_DT*0.5*K[i]*K[i]*v/(CLASSIC_PDE_DK*CLASSIC_PDE_DK);
            lo[m]=-lambda; di[m]=1.0+2.0*lambda; up[m]=-lambda; rhs[m]=C[i];
        }
        // Dirichlet boundaries: C(0,t)=S0 and C(Kmax,t)=0.
        rhs[0]-=lo[0]*S0; lo[0]=0.0;
        up[n-1]=0.0;
        solveTridiagonal(lo,di,up,rhs,x);
        Cnew[0]=S0; Cnew[nK-1]=0.0;
        for (int m=0;m<n;++m) Cnew[m+1]=x[m];
        C.swap(Cnew);
        if (nextMat<MATURITIES.size() && std::abs(t-MATURITIES[nextMat])<CLASSIC_PDE_DT*0.25) {
            std::vector<double> vals;
            vals.reserve(g.surface.K.size());
            for (double kq:g.surface.K) {
                const int idx=static_cast<int>(std::round(kq/CLASSIC_PDE_DK));
                vals.push_back(C[idx]);
            }
            sol.pricesAtMaturity[MATURITIES[nextMat]]=vals;
            ++nextMat;
        }
    }
    return sol;
}

static void writeClassicDiagnostics(const std::string& dir, const ClassicGrid& g) {
    std::ofstream out(dir+"/classic_dupire_local_variance.csv");
    out << "maturity,strike,market_price,dCdt,gamma_fd,raw_dupire_variance,regularized_variance,gamma_nonpositive,variance_invalid,clipped\n";
    out << std::setprecision(15);
    for (size_t j=0;j<g.surface.T.size();++j) for (size_t i=0;i<g.surface.K.size();++i) {
        const auto& n=g.node[j][i];
        out << g.surface.T[j] << ',' << g.surface.K[i] << ',' << g.surface.C[j][i] << ','
            << n.dCdt << ',' << n.gamma << ',' << n.rawVariance << ',' << n.regularizedVariance << ','
            << n.gammaInvalid << ',' << n.varianceInvalid << ',' << n.clipped << '\n';
    }
}

static void writeExtendedComplex(const std::string& dir, const SurfaceGrid& s, const Params& p) {
    const double alpha=1.0+p.beta;
    std::ofstream out(dir+"/extended_complex_dupire.csv");
    out << "maturity,strike,market_price,C_plus,C_minus,gamma_market,gamma_envelope,v_plus,v_minus,sigma_ext2_real,sigma_ext2_imag,formal_naive_variance\n";
    out << std::setprecision(15);
    for (size_t j=0;j<s.T.size();++j) {
        const double vp=instantVariance(s.T[j],p.wPlus), vm=instantVariance(s.T[j],p.wMinus);
        for (double K:s.K) {
            const double cp=callFromTotalVariance(S0,K,p.wPlus[j]);
            const double cm=callFromTotalVariance(S0,K,p.wMinus[j]);
            const double gp=strikeGammaFromTotalVariance(S0,K,p.wPlus[j]);
            const double gm=strikeGammaFromTotalVariance(S0,K,p.wMinus[j]);
            const double P=alpha*gp, N=p.beta*gm;
            const double gammaC=P-N, gammaH=P+N;
            const double denom=P*P+N*N;
            const double re=denom>0?(P*P*vp+N*N*vm)/denom:0.0;
            const double im=denom>0?(P*N*(vm-vp))/denom:0.0;
            const double formal=(std::abs(gammaC)>1e-16)?(P*vp-N*vm)/gammaC:std::numeric_limits<double>::quiet_NaN();
            const double c=alpha*cp-p.beta*cm;
            out << s.T[j] << ',' << K << ',' << c << ',' << cp << ',' << cm << ',' << gammaC << ','
                << gammaH << ',' << vp << ',' << vm << ',' << re << ',' << im << ',' << formal << '\n';
        }
    }
}

static void writeComparison(const std::string& dir, const SurfaceGrid& s,
                            const Params& p, const ClassicSolution& classic) {
    std::ofstream out(dir+"/model_comparison.csv");
    out << "maturity,strike,market_price,extended_price,classic_price,extended_error,classic_error,extended_abs_error,classic_abs_error\n";
    out << std::setprecision(15);
    for (size_t j=0;j<s.T.size();++j) {
        const auto it=classic.pricesAtMaturity.find(s.T[j]);
        if (it==classic.pricesAtMaturity.end()) throw std::runtime_error("Missing classic PDE maturity");
        for (size_t i=0;i<s.K.size();++i) {
            const double market=s.C[j][i];
            const double ext=modelCall(s.K[i],s.T[j],p);
            const double cls=it->second[i];
            out << s.T[j] << ',' << s.K[i] << ',' << market << ',' << ext << ',' << cls << ','
                << ext-market << ',' << cls-market << ',' << std::abs(ext-market) << ',' << std::abs(cls-market) << '\n';
        }
    }
}

static void writeCalibratedParams(const std::string& dir, const ScenarioDef& s, const Params& p) {
    std::ofstream out(dir+"/extended_calibrated_parameters.csv");
    out << "parameter,index,maturity,true_value,calibrated_value\n" << std::setprecision(15);
    out << "beta,0,0," << s.betaTrue << ',' << p.beta << '\n';
    for (size_t j=0;j<MATURITIES.size();++j) {
        out << "w_plus,"<<j<<','<<MATURITIES[j]<<','<<s.wPlusTrue[j]<<','<<p.wPlus[j]<<'\n';
        out << "w_minus,"<<j<<','<<MATURITIES[j]<<','<<s.wMinusTrue[j]<<','<<p.wMinus[j]<<'\n';
    }
}

static void writeSummary(const std::string& dir, const ScenarioDef& def, const SurfaceGrid& s,
                         const Params& p, const ClassicGrid& cg, const ClassicSolution& classic) {
    long double sseE=0.0L,sseC=0.0L; double maxE=0.0,maxC=0.0;
    int negativeGamma=0, invalidVar=0, clipped=0; double worstGamma=0.0, minRaw=1e300, maxRaw=-1e300;
    int calendarViol=0; double worstCal=0.0;
    for (size_t j=0;j<s.T.size();++j) {
        for (size_t i=0;i<s.K.size();++i) {
            const double market=s.C[j][i], ext=modelCall(s.K[i],s.T[j],p);
            const double cls=classic.pricesAtMaturity.at(s.T[j])[i];
            const double e1=ext-market,e2=cls-market;
            sseE+=(long double)e1*e1; sseC+=(long double)e2*e2;
            maxE=std::max(maxE,std::abs(e1)); maxC=std::max(maxC,std::abs(e2));
            const auto& n=cg.node[j][i];
            if (n.gammaInvalid) ++negativeGamma;
            if (n.varianceInvalid) ++invalidVar;
            if (n.clipped) ++clipped;
            if (std::isfinite(n.gamma)) worstGamma=std::min(worstGamma,n.gamma);
            if (std::isfinite(n.rawVariance)) { minRaw=std::min(minRaw,n.rawVariance); maxRaw=std::max(maxRaw,n.rawVariance); }
            if (j>0) {
                const double d=s.C[j][i]-s.C[j-1][i];
                if (d<0) {++calendarViol; worstCal=std::min(worstCal,d);}
            }
        }
    }
    const double n=static_cast<double>(s.T.size()*s.K.size());
    std::ofstream out(dir+"/validation_summary.txt");
    out << std::setprecision(15);
    out << def.label << "\nSpot: "<<S0<<"\nr=q=0 forward-neutralized demonstration\n";
    out << "Quotes: "<<static_cast<long long>(n)<<"\n";
    out << "Nonpositive/invalid market gamma nodes used by classic Dupire: "<<negativeGamma<<"\n";
    out << "Worst discrete market gamma: "<<worstGamma<<"\n";
    out << "Calendar violations: "<<calendarViol<<"\nWorst calendar price difference: "<<worstCal<<"\n";
    out << "Raw classic Dupire invalid/nonpositive variance nodes: "<<invalidVar<<"\n";
    out << "Classic variance nodes clipped/floored/capped: "<<clipped<<"\n";
    out << "Finite raw classic variance range: ["<<minRaw<<", "<<maxRaw<<"]\n";
    out << "Extended calibrated beta: "<<p.beta<<"\n";
    out << "Extended max abs pricing error: "<<maxE<<"\n";
    out << "Extended RMSE: "<<std::sqrt((double)(sseE/n))<<"\n";
    out << "Classic admissible diffusion max abs pricing error: "<<maxC<<"\n";
    out << "Classic admissible diffusion RMSE: "<<std::sqrt((double)(sseC/n))<<"\n";
    out << "Classic regularization: variance clipped to ["<<CLASSIC_VAR_FLOOR<<", "<<CLASSIC_VAR_CAP<<"]\n";
    out << "Classic forward PDE: fully implicit, dK="<<CLASSIC_PDE_DK<<", dt="<<CLASSIC_PDE_DT<<", Kmax="<<CLASSIC_PDE_KMAX<<"\n";
}

static void runScenario(const ScenarioDef& def) {
    std::cout << "\n=== "<<def.label<<" ===\n";
    generateSyntheticMarket(def);
    const std::string dataPath="data/eurostoxx_"+def.name+"_market.csv";
    const std::string outDir="output/"+def.name;
    const auto quotes=loadMarket(dataPath);
    const auto surface=makeSurface(quotes);
    std::cout << "quotes="<<quotes.size()<<"; calibrating extended model...\n";
    const Params p=calibrateExtended(def,quotes);
    std::cout << "extended beta="<<std::setprecision(12)<<p.beta<<"\n";
    const auto classicGrid=computeClassicGrid(surface);
    std::cout << "solving admissible classic Dupire forward PDE...\n";
    const auto classic=solveClassicForwardPDE(classicGrid);
    writeClassicDiagnostics(outDir,classicGrid);
    writeExtendedComplex(outDir,surface,p);
    writeComparison(outDir,surface,p,classic);
    writeCalibratedParams(outDir,def,p);
    writeSummary(outDir,def,surface,p,classicGrid,classic);
}

int main() {
    try {
        ensureDirs();
        const ScenarioDef regular{
            "regular",
            "Synthetic EuroStoxx-like arbitrageable surface (moderate violations)",
            0.14459377107413454,
            {0.00412421,0.00997241,0.02206662,0.03357289},
            {0.00491697,0.01632676,0.04630139,0.09493636},
            2500.0,7000.0,100.0,
            0.20,0.20*0.20,0.28*0.28
        };
        const ScenarioDef violent{
            "violent",
            "Synthetic EuroStoxx-like surface with violent butterfly arbitrage",
            1.0,
            {0.016,0.040,0.080,0.160},
            {0.0002,0.0005,0.0010,0.0020},
            3000.0,7000.0,100.0,
            0.70,0.35*0.35,0.07*0.07
        };
        runScenario(regular);
        runScenario(violent);
        std::cout << "\nAll scenarios completed.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "ERROR: "<<e.what()<<"\n";
        return 1;
    }
}
