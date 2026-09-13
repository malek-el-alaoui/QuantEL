#include "market_data.hpp"
#include "lognormal_mixture.hpp"

#include <cstddef>
#include <cctype>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace lnmix {
namespace {

std::vector<std::string> splitCsv(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    bool quoted = false;
    for (std::size_t i = 0; i < line.size(); ++i) {
        const char ch = line[i];
        if (ch == '"') {
            if (quoted && i + 1 < line.size() && line[i + 1] == '"') {
                cur.push_back('"'); ++i;
            } else quoted = !quoted;
        } else if (ch == ',' && !quoted) {
            out.push_back(cur); cur.clear();
        } else cur.push_back(ch);
    }
    out.push_back(cur);
    return out;
}

std::string trim(std::string s) {
    auto notSpace = [](unsigned char c) { return !std::isspace(c); };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), notSpace));
    s.erase(std::find_if(s.rbegin(), s.rend(), notSpace).base(), s.end());
    return s;
}

double numberOrNaN(const std::string& s) {
    const auto x = trim(s);
    if (x.empty() || x == "NA" || x == "NaN" || x == "nan")
        return std::numeric_limits<double>::quiet_NaN();
    return std::stod(x);
}

std::unordered_map<std::string, std::size_t> headerMap(const std::vector<std::string>& h) {
    std::unordered_map<std::string, std::size_t> m;
    for (std::size_t i = 0; i < h.size(); ++i) m[trim(h[i])] = i;
    return m;
}

std::string field(const std::vector<std::string>& r,
                  const std::unordered_map<std::string, std::size_t>& h,
                  const std::string& name,
                  const std::string& fallback = "") {
    auto it = h.find(name);
    if (it == h.end() || it->second >= r.size()) return fallback;
    return trim(r[it->second]);
}

MarketSlice& getSlice(std::map<double, MarketSlice>& slices,
                      double t,
                      const std::string& instrument,
                      const std::string& asOf,
                      const std::string& maturity,
                      double forward) {
    std::pair<std::map<double, MarketSlice>::iterator, bool> result =
        slices.emplace(t, MarketSlice{});
    std::map<double, MarketSlice>::iterator it = result.first;
    const bool inserted = result.second;
    if (inserted) {
        it->second.instrument = instrument;
        it->second.asOf = asOf;
        it->second.maturity = maturity;
        it->second.t = t;
        it->second.forward = forward;
    }
    return it->second;
}

} // namespace

MarketSurface loadStrikeSurfaceCsv(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open data file: " + path);

    std::string line;
    if (!std::getline(in, line)) throw std::runtime_error("Empty data file: " + path);
    const auto h = headerMap(splitCsv(line));
    std::map<double, MarketSlice> slices;
    MarketSurface out;

    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto r = splitCsv(line);
        const auto instrument = field(r, h, "instrument");
        const auto asOf = field(r, h, "asof");
        const auto maturity = field(r, h, "maturity");
        const double t = std::stod(field(r, h, "ttm"));
        const double forward = std::stod(field(r, h, "forward"));
        const double strike = std::stod(field(r, h, "strike"));
        const double marketIv = std::stod(field(r, h, "market_iv"));
        const double sourcePrice = numberOrNaN(field(r, h, "source_call_mid", "NA"));
        const double sourceIv = numberOrNaN(field(r, h, "source_iv", "NA"));
        const auto pillar = field(r, h, "pillar", std::to_string(strike));
        const double k = strike / forward;

        auto& s = getSlice(slices, t, instrument, asOf, maturity, forward);
        s.points.push_back({k, strike, marketIv, sourcePrice, sourceIv, pillar});
        if (out.instrument.empty()) { out.instrument = instrument; out.asOf = asOf; }
    }

    for (std::map<double, MarketSlice>::iterator it = slices.begin(); it != slices.end(); ++it) {
        MarketSlice& s = it->second;
        std::sort(s.points.begin(), s.points.end(),
                  [](const MarketPoint& a, const MarketPoint& b) { return a.k < b.k; });
        out.slices.push_back(std::move(s));
    }
    return out;
}

MarketSurface loadFxDeltaSurfaceCsv(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open data file: " + path);

    std::string line;
    if (!std::getline(in, line)) throw std::runtime_error("Empty data file: " + path);
    const auto h = headerMap(splitCsv(line));
    std::map<double, MarketSlice> slices;
    MarketSurface out;

    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto r = splitCsv(line);
        const auto instrument = field(r, h, "instrument");
        const auto asOf = field(r, h, "asof");
        const auto maturity = field(r, h, "maturity");
        const double t = std::stod(field(r, h, "ttm"));
        const double deltaPct = std::stod(field(r, h, "delta"));
        const auto optionType = field(r, h, "option_type");
        const double vol = std::stod(field(r, h, "market_iv"));
        const auto pillar = field(r, h, "pillar");
        const double totalVol = vol * std::sqrt(t);

        // Demonstration conversion uses forward delta. For a call, Delta_F=N(d1).
        // For a put with absolute delta d, N(d1)=1-d. ATM is set to K/F=1.
        double k = 1.0;
        if (optionType == "Call") {
            const double d1 = LognormalMixtureCalibrator::invNormCdf(deltaPct / 100.0);
            k = std::exp(0.5 * totalVol * totalVol - d1 * totalVol);
        } else if (optionType == "Put") {
            const double d1 = LognormalMixtureCalibrator::invNormCdf(1.0 - deltaPct / 100.0);
            k = std::exp(0.5 * totalVol * totalVol - d1 * totalVol);
        }

        auto& s = getSlice(slices, t, instrument, asOf, maturity, 1.0);
        s.points.push_back({k, k, vol,
                            std::numeric_limits<double>::quiet_NaN(), vol, pillar});
        if (out.instrument.empty()) { out.instrument = instrument; out.asOf = asOf; }
    }

    for (std::map<double, MarketSlice>::iterator it = slices.begin(); it != slices.end(); ++it) {
        MarketSlice& s = it->second;
        std::sort(s.points.begin(), s.points.end(),
                  [](const MarketPoint& a, const MarketPoint& b) { return a.k < b.k; });
        out.slices.push_back(std::move(s));
    }
    return out;
}

} // namespace lnmix
