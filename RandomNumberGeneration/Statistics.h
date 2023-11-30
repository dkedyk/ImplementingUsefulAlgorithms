#ifndef IGMDK_STATISTICS_H
#define IGMDK_STATISTICS_H
#include "Random.h"
#include "../Sorting/Sort.h"
#include "../Utils/Bits.h"
#include "../Utils/Vector.h"
#include "../Heaps/Heap.h"
#include "../NumericalMethods/Matrix.h"
#include "../NumericalMethods/NumericalMethods.h"
#include <cstdlib>
#include <cmath>
namespace igmdk{

template<typename CDF> double invertCDF(CDF const& c, double u,
    double guess = 0, double step0 = 1, double prec = 0.0001)
{
    assert(u > 0 && u < 1);
    auto f = [c, u](double x){return u - c(x);};
    pair<double, double> bracket = findInterval0(f, guess, step0, 100);
    return solveFor0(f, bracket.first, bracket.second, prec).first;
}

double approxNormalCDF(double x)
{
    assert(isfinite(x));
    double uHalf = erf(abs(x)/sqrt(2))/2;
    return 0.5 + (x >= 0 ? uHalf : -uHalf);
}
double approxNormal2SidedConf(double x){return 2 * approxNormalCDF(x) - 1;}

double find2SidedConfZ(double conf)
{
    assert(conf > 0 && conf < 1);
    return invertCDF([](double x){return approxNormalCDF(x);}, 0.5 + conf/2);
}

struct HoefFunctor
{
    double ave, a;
    int n;
    double getTerm(double t, double d)const
    {//limit is 0 safely approached
        if(d < 0) assert(d <= t);
        double result = log(1 + d/t) * t;
        return isfinite(result) ? result : 0;
    }
public:
    double operator()(double d)const
    {
        assert(d >= 0 && d <= 1);
        return (getTerm(1 - ave, -d) + getTerm(ave, d)) * n - log(a);
    }
    static pair<double, double> conf(double ave, int n,
        double confidence = 0.95)
    {
        assert(ave >= 0 && ave <= 1 && confidence > 0 && confidence < 1 &&
            n > 0);
        HoefFunctor f = {ave, (1 - confidence)/2, n};
        double e = numeric_limits<double>::epsilon(), upperD = ave < 1 - e ?
            solveFor0(f, 0, 1 - ave - e).first : 0;
        f.ave = 1 - ave;
        double lowerD = ave > e ? solveFor0(f, 0,
            ave - e).first : 0;
        return make_pair(ave - lowerD, ave + upperD);
    }
};

double approxTCDF(double t, int v)
{//Hill3 method, always < 10-4?
    assert(v > 0);
    if(v == 1) return 0.5 + atan(t)/PI();//Cauchy
    if(v == 2) return 0.5 + t/2/sqrt(2 + t * t);//also exact
    double a = v - 0.5, b = 48 * a * a, w = sqrt(a * log(1 + t * t/v)), w27[6];
    w27[0] = w * w;
    for(int i = 1; i < 6; ++i) w27[i] = w * w27[i - 1];
    double z = w + (w27[0] + 3) * w/b - (4 * w27[5] + 33 * w27[3] +
        240 * w27[1] + 855 * w)/(10 * b *(b + 0.8 * w27[2] + 100));
    return approxNormalCDF(t > 0 ? z : -z);
}
double approxT2SidedConf(double x, int v){return 2 * approxTCDF(x, v) - 1;}
double find2SidedConfT(double conf, int v)
{
    assert(conf > 0 && conf < 1 && v > 0);
    return invertCDF([v](double x){return approxTCDF(x, v);}, 0.5 + conf/2);
}
pair<double, double> getTConf(IncrementalStatistics const&s, double a = 0.05)
{
    assert(s.n > 1 && a > 0 && a < 1);
    double ste = s.getStandardErrorSummary().stddev() *
        find2SidedConfT(1 - a, s.n - 1);
    return make_pair(s.getMean() - ste, s.getMean() + ste);
}
pair<double, double> getTConf(Vector<double> const& data, double a = 0.05)
{
    IncrementalStatistics s;
    for(int i = 0; i < data.getSize(); ++i) s.addValue(data[i]);
    return getTConf(s, a);
}

bool confIncludes(pair<double, double> const& interval, double value)
    {return interval.first <= value && value <= interval.second;}

template<typename FUNCTION> struct SpeedTester
{
    FUNCTION f;
    SpeedTester(FUNCTION const& theFunction = FUNCTION()): f(theFunction){}
    double operator()()const
    {
        int now = clock();
        f();
        return (clock() - now) * 1.0/CLOCKS_PER_SEC;
    }
};

bool normalTestAreEqual(double m1, double v1, double m2, double v2,
    double z = 3)
{
    NormalSummary diff = NormalSummary(m1, v1) - NormalSummary(m2, v2);
    return abs(diff.mean) <= z * diff.stddev();
}

bool signTestAreEqual(double winCount1, double winCount2, double z = 2)
    {return abs(winCount1 - winCount2)/sqrt(winCount1 + winCount2) <= z;}
pair<double, double> countWins(Vector<pair<double, double> > const& data)
{
    int n1 = 0, n2 = 0;
    for(int i = 0; i < data.getSize(); ++i)
    {
        if(data[i].first < data[i].second) ++n1;
        else if(data[i].first > data[i].second) ++n2;
        else{n1 += 0.5; n2 += 0.5;}
    }
    return make_pair(n1, n2);
}
bool signTestPairs(Vector<pair<double, double> > const& data, double z = 2)
{
    pair<double, double> wins = countWins(data);
    return signTestAreEqual(wins.first, wins.second, z);
}

Vector<double> convertToRanks(Vector<double> a)
{//create index array, sort it, and convert indices into ranks
    int n = a.getSize();
    Vector<int> indices(n);
    for(int i = 0; i < n; ++i) indices[i] = i;
    IndexComparator<double> c(a.getArray());
    quickSort(indices.getArray(), 0, n - 1, c);
    for(int i = 0; i < n; ++i)
    {//rank lookahead to scan for ties, then change a entries
        int j = i;
        while(i + 1 < n && c.isEqual(indices[j], indices[i + 1])) ++i;
        double rank = (i + j)/2.0 + 1;
        for(; j <= i; ++j) a[indices[j]] = rank;
    }
    return a;
}

double trimmedMean(Vector<double> data, double c = 0.2,
    bool isSorted = false)
{
    int n = data.getSize(), trim = c * n;
    assert(n > 0 && c >= 0 && c < 0.5);
    if(!isSorted)
    {
        quickSelect(data.getArray(), n, n - trim - 1);
        quickSelect(data.getArray(), n - trim - 1, trim);
    }
    double sum = 0;
    for(int i = trim; i < n - trim; ++i) sum += data[i];
    return sum/(n - 2 * trim);
}
double trimmedMeanStandardError(Vector<double> data, double c = 0.2,
    bool isSorted = false)
{
    int n = data.getSize(), trim = c * n;
    assert(n > 0 && c >= 0 && c < 0.5);
    if(!isSorted)
    {
        quickSelect(data.getArray(), n, n - trim - 1);
        quickSelect(data.getArray(), n - trim - 1, trim);
    }//Windsorise tails
    for(int i = 0; i < trim; ++i) data[i] = data[trim];
    for(int i = n - trim; i < n; ++i) data[i] = data[n - trim - 1];
    IncrementalStatistics s;//calc regular se of values
    for(int i = 0; i < n; ++i) s.addValue(data[i]);
    return s.getStandardErrorSummary().stddev()/(1 - 2 * c);
}
pair<double, double> trimmedMeanConf(Vector<double> const& data, double z = 2)
{
    double tm = trimmedMean(data), ste = trimmedMeanStandardError(data) * z;
    return make_pair(tm - ste, tm + ste);
}

double quantile(Vector<double> data, double q, bool isSorted = false)
{
    assert(data.getSize() > 0);
    if(q < 0) return -numeric_limits<double>::infinity();
    else if(q > 1) return numeric_limits<double>::infinity();
    int n = data.getSize(), u = q * n, l = u - 1;
    if(u == n) u = l;//check corner cases
    else if(u == 0 || double(u) != q * n) l = u;//and border values
    if(!isSorted)
    {
        quickSelect(data.getArray(), n, u);
        if(l != u) quickSelect(data.getArray(), u, l);
    }
    return (data[l] + data[u])/2;
}
double median(Vector<double> const& data, bool isSorted = false)
    {return quantile(data, 0.5, isSorted);}
pair<double, double> quantileConf(Vector<double> const& data, double q = 0.5,
    bool isSorted = false, double z = 2)
{
    double d = z * sqrt(q * (1 - q)/data.getSize());
    return make_pair(quantile(data, q - d, isSorted),
        quantile(data, q + d, isSorted));
}

pair<double, double> normal2SampleDiff(double mean1, double ste1,
    double mean2, double ste2, double z)
{//difference of approximately normal-based confidence intervals
    NormalSummary n1(mean1, ste1 * ste1), n2(mean2, ste2 * ste2),
        diff = n1 - n2;
    double ste = diff.stddev() * z;
    return make_pair(diff.mean - ste, diff.mean + ste);
}
pair<double, double> normalConfDiff(double mean1, pair<double, double> const&
    conf1, double mean2, pair<double, double> const& conf2, double z)
{//difference of approximately normal-based confidence intervals
    return normal2SampleDiff(mean1, (conf1.second - conf1.first)/2/z,
        mean2, (conf2.second - conf2.first)/2/z, z);
}
pair<double, double> median2SampleDiffConf(Vector<double> const& samples1,
    Vector<double> const& samples2, double z = 2)
{
    return normalConfDiff(median(samples1), quantileConf(samples1, 0.5,false,
        z), median(samples2), quantileConf(samples2, 0.5, false, z), z);
}

pair<double, double> wilsonScoreInterval(double p, int n, double z = 2)
{
    assert(p >= 0 && p <= 1 && n > 0 && z >= 0);
    double temp = z * z/n, delta = z * sqrt((p * (1 - p) + temp/4)/n);
    return make_pair((p + temp/2 - delta)/(1 + temp),
        (p + temp/2 + delta)/(1 + temp));
}

pair<double, double> trimmedMean2SampleDiffConf(Vector<double> const&
    samples1, Vector<double> const& samples2, double z = 2)
{
    return normal2SampleDiff(trimmedMean(samples1),
        trimmedMeanStandardError(samples1), trimmedMean(samples2),
        trimmedMeanStandardError(samples2), z);
}


pair<double, double> medianMADN(Vector<double> data)
{//scaled to match normal stdev
    assert(data.getSize() > 2);
    double med = median(data);
    for(int i = 0; i < data.getSize(); ++i) data[i] = abs(data[i] - med);
    return make_pair(med, 1.4826 * median(data));
}

}//end namespace
#endif
