#include "Random.h"
#include "Statistics.h"
#include "PermutationTests.h"
#include "Bootstrap.h"
#include "DistributionTests.h"
#include "../NumericalMethods/NumericalMethods.h"
#include "../NumericalMethods/Matrix.h"
#include "../Utils/Debug.h"
#include "../ExternalMemoryAlgorithms/CSV.h"
#include "testCommon.h"
using namespace igmdk;

struct BernFunctor
{
    double m, t;
    double operator()(double x)const
    {
        return (m + sqrt(x * (1 - x) * 2 * t) + t/3) - x;
    }
};
pair<double, double> numericalBernBounds(IncrementalStatistics const & si,
    double confidence = 0.95, int nFactor = 1)
{
    int n = si.n;
    assert(n > 0 && confidence > 0 && confidence < 1);
    double p = 1 - confidence, t = log(2/p)/n;
    BernFunctor b = {si.getMean(), t};
    double result = solveFor0(b, 0, 1).first;
    DEBUG(result);
    BernFunctor b2 = {1 - si.getMean(), t};
    double result2 = solveFor0(b2, 0, 1).first;
    DEBUG(result2);
    return make_pair(1 - result2, result);
}

void testNormalEval()
{
    DEBUG(approxNormal2SidedConf(2));
    DEBUG(find2SidedConfZ(0.95));
    DEBUG(find2SidedConfZBonf(100, 0.95));
}

void testTrimmedMean()
{
    Vector<double> data;
    for(int i = 0; i < 10; ++i) data.append(i);
    data[0] = -100;
    data[1] = -200;
    data[8] = 300;
    data[9] = 400;
    double tm = trimmedMean(data);
    DEBUG(tm);
}

void testQuantiles()
{
    Vector<double> data;
    for(int i = 0; i < 4; ++i) data.append(i + 1);
    DEBUG(quantile(data, -0.01));
    DEBUG(quantile(data, 0));
    DEBUG(quantile(data, 0.01));
    DEBUG(quantile(data, 0.24));
    DEBUG(quantile(data, 0.25));
    DEBUG(quantile(data, 0.26));
    DEBUG(quantile(data, 0.49));
    DEBUG(quantile(data, 0.5));
    DEBUG(quantile(data, 0.51));
    DEBUG(quantile(data, 0.74));
    DEBUG(quantile(data, 0.75));
    DEBUG(quantile(data, 0.76));
    DEBUG(quantile(data, 0.99));
    DEBUG(quantile(data, 1));
    DEBUG(quantile(data, 1.01));
    for(int i = 4; i < 100; ++i) data.append(i + 1);
    pair<double, double> conf = quantileConf(data, 0.5);
    DEBUG(conf.first);
    DEBUG(median(data));
    DEBUG(conf.second);
    IncrementalStatistics s;
    for(int i = 1; i < 100; ++i) s.addValue(i + 1);
    DEBUG(s.getMean());
    DEBUG(s.getStandardErrorSummary().error95());
    DEBUG(trimmedMean(data));
    DEBUG(2 * trimmedMeanStandardError(data));

    conf = quantileConf(data, 0.01);
    DEBUG(conf.first);
    DEBUG(quantile(data, 0.01));
    DEBUG(conf.second);

    conf = quantileConf(data, 0.99);
    DEBUG(conf.first);
    DEBUG(quantile(data, 0.99));
    DEBUG(conf.second);
}

double approxTCDF2(double t, int v)
{//Gleason method, worst error 2.9e-3 for v = 3 and t = 0.9
    assert(v > 0);
    if(v == 1) return 0.5 + atan(t)/PI();//Cauchy
    if(v == 2) return 0.5 + t/2/sqrt(2 + t * t);//also exact
    double z = sqrt(log(1 + t * t/v)/(v - 1.5 - 0.1/v + 0.5825/v/v)) * (v - 1);
    return approxNormalCDF(t > 0 ? z : -z);
}
void testT()
{//fyi table has limited precision so wont go beyond 5-6 decimals
    DEBUG(approxTCDF(1.96, 3) - 0.975);
    DEBUG(approxTCDF(1.96, 10) - 0.975);
    DEBUG(approxTCDF(1.96, 30) - 0.975);
    DEBUG(approxTCDF(1.96, 100) - 0.975);
    DEBUG(approxTCDF(1.96, 1000) - 0.975);
    DEBUG(approxTCDF(-1.96, 1000) - 0.025);
    DEBUG(approxNormalCDF(-1.96));
    DEBUG(approxNormalCDF(1.96));
    DEBUG(find2SidedConfT(0.95, 3));
    DEBUG(find2SidedConfT(0.95, 10));
    DEBUG(find2SidedConfT(0.95, 30));
    DEBUG(find2SidedConfT(0.95, 100));
    DEBUG(find2SidedConfT(0.95, 1000));

    DEBUG(approxTCDF(12.706, 1) - 0.975);
    DEBUG(approxTCDF2(12.706, 1) - 0.975);
    DEBUG(approxTCDF(4.303, 2) - 0.975);
    DEBUG(approxTCDF2(4.303, 2) - 0.975);
    DEBUG(approxTCDF(0.9, 3) - 0.7828);//worst case for Gleason < 0.0025
    DEBUG(approxTCDF2(0.9, 3) - 0.7828);
    DEBUG(approxTCDF(3.182, 3) - 0.975);
    DEBUG(approxTCDF2(3.182, 3) - 0.975);
    DEBUG(approxTCDF(2.776, 4) - 0.975);
    DEBUG(approxTCDF2(2.776, 4) - 0.975);
    DEBUG(approxTCDF(2.571, 5) - 0.975);
    DEBUG(approxTCDF2(2.571, 5) - 0.975);
    DEBUG(approxTCDF(2.447, 6) - 0.975);
    DEBUG(approxTCDF2(2.447, 6) - 0.975);
    DEBUG(approxTCDF(2.365, 7) - 0.975);
    DEBUG(approxTCDF2(2.365, 7) - 0.975);
    DEBUG(approxTCDF(2.306, 8) - 0.975);
    DEBUG(approxTCDF2(2.306, 8) - 0.975);
    DEBUG(approxTCDF(2.262, 9) - 0.975);
    DEBUG(approxTCDF2(2.262, 9) - 0.975);
    DEBUG(approxTCDF(2.228, 10) - 0.975);
    DEBUG(approxTCDF2(2.228, 10) - 0.975);
}

template<typename SAMPLER>
void testPairedSign(int t = 10000, int n = 100)
{
    IncrementalStatistics s;
    SAMPLER sa;
    for(int j = 0; j < t; ++j)
    {
        Vector<pair<double, double> > samples;
        for(int i = 0; i < n; ++i) samples.append(make_pair(sa(), sa()));
        s.addValue(!signTestPairs(samples));
    }
    DEBUG(s.getMean());
    pair<double, double> conf = wilsonScoreInterval(s.getMean(), s.n);
    DEBUG(conf.first);
    DEBUG(conf.second);
    //ok for this to fail?
    //move uncertainly in % estimate to skepticism?
    assert(conf.first <= 0.05);
}
//ABOUT T CONF:
//unclear when do use this--for unknown small sample fixed data should use
//trimmed mean most likely or bootstrap BCA but useless estimate anyway
//Belle rule says 12 enough for conf but unclear!?
pair<double, double> tTestPairsP(Vector<pair<double, double> > const& data)
{
    IncrementalStatistics s;
    for(int i = 0; i < data.getSize(); ++i)
    {
        s.addValue(data[i].first - data[i].second);
    }
    return getTConf(s);
}
template<typename SAMPLER>
void testPairedT(int t = 10000, int n = 100)
{
    IncrementalStatistics s;
    SAMPLER sa;
    for(int j = 0; j < t; ++j)
    {
        Vector<pair<double, double> > samples;
        for(int i = 0; i < n; ++i) samples.append(make_pair(sa(), sa()));
        s.addValue(!confIncludes(tTestPairsP(samples), 0));
    }
    DEBUG(s.getMean());
    pair<double, double> conf = wilsonScoreInterval(s.getMean(), s.n);
    DEBUG(conf.first);
    DEBUG(conf.second);
    //ok for this to fail?
    //move uncertainly in % estimate to skepticism?
    assert(conf.first <= 0.05);//2.58 too small for uniformly most powerful
}

template<typename LOCATION_F> pair<double, double> permutationPairedConf(
    Vector<pair<double, double> > const& data, double a = 0.05, int b = 10000)
{
    PairedTestLocationPermuter<LOCATION_F> p = {data};
    return permutationConf(p, a, b);
}

double trimmedTestPairsP(Vector<pair<double, double> > const& data)
{
    Vector<double> diffs;
    for(int i = 0; i < data.getSize(); ++i)
    {
        diffs.append(data[i].first - data[i].second);
    }

    double tm = trimmedMean(diffs), ste = trimmedMeanStandardError(diffs);
    double z = tm/ste;
    return 1 - approxT2SidedConf(z, data.getSize());
    //return 1 - approxTCDF(z, data.getSize());
}
template<typename SAMPLER>
void testPairedTrimmed(int t = 10000, int n = 100)
{
    IncrementalStatistics s;
    SAMPLER sa;
    for(int j = 0; j < t; ++j)
    {
        Vector<pair<double, double> > samples;
        for(int i = 0; i < n; ++i) samples.append(make_pair(sa(), sa()));
        s.addValue(trimmedTestPairsP(samples) < 0.05);
    }
    DEBUG(s.getMean());
    pair<double, double> conf = wilsonScoreInterval(s.getMean(), s.n);
    DEBUG(conf.first);
    DEBUG(conf.second);
    //ok for this to fail?
    //move uncertainly in % estimate to skepticism?
    assert(conf.first <= 0.05);//2.58 too small for uniformly most powerful
}


template<typename SAMPLER>
void testPairedPermConf(int t = 10000, int n = 100)
{
    IncrementalStatistics s;
    SAMPLER sa;
    for(int j = 0; j < t; ++j)
    {
        Vector<pair<double, double> > samples;
        for(int i = 0; i < n; ++i) samples.append(make_pair(sa(), sa()));
        s.addValue(!confIncludes(permutationPairedConf<meaner>(samples), 0));
    }
    DEBUG(s.getMean());
    pair<double, double> conf = wilsonScoreInterval(s.getMean(), s.n);
    DEBUG(conf.first);
    DEBUG(conf.second);
    //ok for this to fail?
    //move uncertainly in % estimate to skepticism?
    assert(conf.first <= 0.05);//2.58 too small for uniformly most powerful
}
double permutationPairedTest(Vector<pair<double, double> > const& data,
    int b = 10000)
{
    PairedTestLocationPermuter<meaner> p = {data};
    return permutationTest(p, b);
}
template<typename SAMPLER>
void testPairedPerm(int t = 1000, int n = 100)
{
    IncrementalStatistics s;
    SAMPLER sa;
    for(int j = 0; j < t; ++j)
    {
        Vector<pair<double, double> > samples;
        for(int i = 0; i < n; ++i) samples.append(make_pair(sa(), sa()));
        s.addValue(permutationPairedTest(samples) < 0.05);
    }
    DEBUG(s.getMean());
    pair<double, double> conf = wilsonScoreInterval(s.getMean(), s.n);
    DEBUG(conf.first);
    DEBUG(conf.second);
    //ok for this to fail?
    //move uncertainly in % estimate to skepticism?
    assert(conf.first <= 0.05);
}
void testPairedTests()
{
    DEBUG("testPairedSign<NSampler>");
    testPairedSign<NSampler>();
    DEBUG("testPairedT<NSampler>");
    testPairedT<NSampler>();
    DEBUG("testPairedPermConf<NSampler>");
    testPairedPermConf<NSampler>();
    DEBUG("testPairedTrimmed<NSampler>");
    testPairedTrimmed<NSampler>();
    DEBUG("testPairedPerm<NSampler>");
    testPairedPerm<NSampler>();
}

pair<double, double> normal2SampleConf(Vector<double> const& samples1,
    Vector<double> const& samples2, double z = 2)
{
    IncrementalStatistics s1, s2;
    for(int i = 0; i < samples1.getSize(); ++i) s1.addValue(samples1[i]);
    for(int i = 0; i < samples2.getSize(); ++i) s2.addValue(samples2[i]);
    NormalSummary diff = s1.getStandardErrorSummary() -
        s2.getStandardErrorSummary();
    double ste = diff.stddev() * z;
    return make_pair(diff.mean - ste, diff.mean + ste);
}
template<typename SAMPLER>
void test2SampleDiffNormal(double value = 0, int t = 10000, int n = 100)
{
    DEBUG("n");
    IncrementalStatistics s, s2;
    SAMPLER sa;
    for(int j = 0; j < t; ++j)
    {
        Vector<double> samples, samples2;
        for(int i = 0; i < n; ++i) samples.append(sa());
        for(int i = 0; i < n; ++i) samples2.append(sa());
        pair<double, double> c = normal2SampleConf(samples, samples2);
        s2.addValue(c.second - c.first);
        s.addValue(!confIncludes(normal2SampleConf(samples, samples2), value));
    }
    DEBUG(s2.getMean());
    DEBUG(s.getMean());
    pair<double, double> conf = wilsonScoreInterval(s.getMean(), s.n);
    DEBUG(conf.first);
    DEBUG(conf.second);
    //ok for this to fail?
    //move uncertainly in % estimate to skepticism?
    assert(conf.first <= 0.05);
}


template<typename SAMPLER>
void test2SampleDiffTrimmed(double value = 0, double a = 0.05, int t = 10000, int n = 100)
{
    DEBUG("tm");
    IncrementalStatistics s, s2;
    SAMPLER sa;
    for(int j = 0; j < t; ++j)
    {
        Vector<double> samples, samples2;
        for(int i = 0; i < n; ++i) samples.append(sa());
        for(int i = 0; i < n; ++i) samples2.append(sa());
        pair<double, double> c = trimmedMean2SampleDiffConf(samples, samples2);
        s2.addValue(c.second - c.first);
        s.addValue(!confIncludes(c, value));
    }
    DEBUG(s2.getMean());
    DEBUG(s.getMean());
    pair<double, double> conf = wilsonScoreInterval(s.getMean(), s.n);
    DEBUG(conf.first);
    DEBUG(conf.second);
    //ok for this to fail?
    //move uncertainly in % estimate to skepticism?
    assert(conf.first <= a);
}

bool confIntervalIntersect(pair<double, double> const& i1,
    pair<double, double> const& i2)
{//separated if endpoints bounded from each other
    return !(i1.second < i2.first || i2.second < i1.first);
}
pair<double, double> median2SampleConf(Vector<double> const& samples1,
    Vector<double> const& samples2, double z = 2)
{
    pair<double, double> i1 = quantileConf(samples1),
        i2 = quantileConf(samples2);
    return make_pair(i1.first - i2.second, i2.second - i1.first);
}

template<typename SAMPLER>
void test2SampleDiffMedian(double value = 0, int t = 10000, int n = 100)
{
    DEBUG("med");
    IncrementalStatistics s, s2;
    SAMPLER sa;
    for(int j = 0; j < t; ++j)
    {
        Vector<double> samples, samples2;
        for(int i = 0; i < n; ++i) samples.append(sa());
        for(int i = 0; i < n; ++i) samples2.append(sa());
        pair<double, double> c = median2SampleConf(samples, samples2);
        s2.addValue(c.second - c.first);
        s.addValue(!confIncludes(median2SampleConf(samples, samples2), value));
    }
    DEBUG(s2.getMean());
    DEBUG(s.getMean());
    pair<double, double> conf = wilsonScoreInterval(s.getMean(), s.n);
    DEBUG(conf.first);
    DEBUG(conf.second);
    //ok for this to fail?
    //move uncertainly in % estimate to skepticism?
    assert(conf.first <= 0.05);
}

template<typename SAMPLER>
void test2SampleDiffMedianApprox(double value = 0, int t = 10000, int n = 100)
{
    DEBUG("med approx");
    IncrementalStatistics s, s2;
    SAMPLER sa;
    for(int j = 0; j < t; ++j)
    {
        Vector<double> samples, samples2;
        for(int i = 0; i < n; ++i) samples.append(sa());
        for(int i = 0; i < n; ++i) samples2.append(sa());
        pair<double, double> c = median2SampleDiffConf(samples, samples2);
        s2.addValue(c.second - c.first);
        s.addValue(!confIncludes(c, value));
    }
    DEBUG(s2.getMean());
    DEBUG(s.getMean());
    pair<double, double> conf = wilsonScoreInterval(s.getMean(), s.n);
    DEBUG(conf.first);
    DEBUG(conf.second);
    //ok for this to fail?
    //move uncertainly in % estimate to skepticism?
    assert(conf.first <= 0.06);
}

void test2SampleDiff()
{
    DEBUG("uniform");
    test2SampleDiffNormal<USampler>();
    test2SampleDiffMedian<USampler>();
    test2SampleDiffMedianApprox<USampler>();
    test2SampleDiffTrimmed<USampler>();
    DEBUG("normal");
    test2SampleDiffNormal<NSampler>();
    test2SampleDiffMedian<NSampler>();
    test2SampleDiffMedianApprox<NSampler>();
    test2SampleDiffTrimmed<NSampler>();
    DEBUG("expo");
    test2SampleDiffNormal<ExpoSampler>();
    test2SampleDiffMedian<ExpoSampler>();
    test2SampleDiffMedianApprox<ExpoSampler>();
    test2SampleDiffTrimmed<ExpoSampler>();
    DEBUG("Cauchy");
    test2SampleDiffNormal<CauchySampler>();
    test2SampleDiffMedian<CauchySampler>();
    //med approx fails if diff left skew with right skew?
    test2SampleDiffMedianApprox<CauchySampler>();
    //for Caucny trimmed mean diff sensitive to last value? just not as much as mean which is even more sensitive!?
    test2SampleDiffTrimmed<CauchySampler>(0, 0.11); //fails 0.05

    //DO BOOTSRAP ALSO EXTEND IT TO TAKE VECTOR OF SAMPLES?
}

template<typename F>
pair<double, double> bootSampleConf(Vector<double> const& data)
{
    BasicBooter<F> booter(data);
    return bootstrapMixedInterval(booter);
}
template<typename F, typename SAMPLER>
void testSampleBoot(double value = 0, double al=0.05,int t = 10000, int n = 100)
{
    DEBUG("boot");
    IncrementalStatistics s, s2;
    SAMPLER sa;
    for(int j = 0; j < t; ++j)
    {
        Vector<double> samples, samples2;
        for(int i = 0; i < n; ++i) samples.append(sa());
        pair<double, double> c = bootSampleConf<F>(samples);
        s2.addValue(c.second - c.first);
        s.addValue(!confIncludes(c, value));
    }
    DEBUG(s2.getMean());
    DEBUG(s.getMean());
    pair<double, double> conf = wilsonScoreInterval(s.getMean(), s.n);
    DEBUG(conf.first);
    DEBUG(conf.second);
    //ok for this to fail?
    //move uncertainly in % estimate to skepticism?
    assert(conf.first <= al);
}

pair<double, double> normalSampleConf(Vector<double> const& samples1,
    double z = 2)
{
    IncrementalStatistics s1;
    for(int i = 0; i < samples1.getSize(); ++i) s1.addValue(samples1[i]);
    NormalSummary diff = s1.getStandardErrorSummary();
    double ste = diff.stddev() * z;
    return make_pair(diff.mean - ste, diff.mean + ste);
}
template<typename SAMPLER>
void testSampleNormal(double value = 0, double al=0.05,int t = 10000, int n = 100)
{
    DEBUG("n");
    IncrementalStatistics s, s2;
    SAMPLER sa;
    for(int j = 0; j < t; ++j)
    {
        Vector<double> samples, samples2;
        for(int i = 0; i < n; ++i) samples.append(sa());
        pair<double, double> c = normalSampleConf(samples);
        s2.addValue(c.second - c.first);
        s.addValue(!confIncludes(c, value));
    }
    DEBUG(s2.getMean());
    DEBUG(s.getMean());
    pair<double, double> conf = wilsonScoreInterval(s.getMean(), s.n);
    DEBUG(conf.first);
    DEBUG(conf.second);
    //ok for this to fail?
    //move uncertainly in % estimate to skepticism?
    assert(conf.first <= al);
}

template<typename SAMPLER>
void testSampleTrimmed(double value = 0, double al=0.05,int t = 10000, int n = 100)
{
    DEBUG("tm");
    IncrementalStatistics s, s2;
    SAMPLER sa;
    for(int j = 0; j < t; ++j)
    {
        Vector<double> samples, samples2;
        for(int i = 0; i < n; ++i) samples.append(sa());
        pair<double, double> c = trimmedMeanConf(samples);
        s2.addValue(c.second - c.first);
        s.addValue(!confIncludes(c, value));
    }
    DEBUG(s2.getMean());
    DEBUG(s.getMean());
    pair<double, double> conf = wilsonScoreInterval(s.getMean(), s.n);
    DEBUG(conf.first);
    DEBUG(conf.second);
    //ok for this to fail?
    //move uncertainly in % estimate to skepticism?
    assert(conf.first <= al);
}

template<typename SAMPLER>
void testSampleMedian(double value = 0, double al=0.05,int t = 10000, int n = 100)
{
    DEBUG("med");
    IncrementalStatistics s, s2;
    SAMPLER sa;
    for(int j = 0; j < t; ++j)
    {
        Vector<double> samples, samples2;
        for(int i = 0; i < n; ++i) samples.append(sa());
        pair<double, double> c = quantileConf(samples);
        s2.addValue(c.second - c.first);
        s.addValue(!confIncludes(c, value));
    }
    DEBUG(s2.getMean());
    DEBUG(s.getMean());
    pair<double, double> conf = wilsonScoreInterval(s.getMean(), s.n);
    DEBUG(conf.first);
    DEBUG(conf.second);
    //ok for this to fail?
    //move uncertainly in % estimate to skepticism?
    assert(conf.first <= al);
}

void test1SampleConf()
{
    DEBUG("uniform");
    testSampleNormal<USampler>(0.5);
    testSampleBoot<meaner, USampler>(0.5, 0.05, 100);
    testSampleMedian<USampler>(0.5, 0.06);
    testSampleBoot<meder, USampler>(0.5, 0.15, 100);
    testSampleTrimmed<USampler>(0.5);
    testSampleBoot<trimmer, USampler>(0.5, 0.05, 100);
    DEBUG("normal");
    testSampleNormal<NSampler>();
    testSampleBoot<meaner, NSampler>(0, 0.05, 100);
    testSampleMedian<NSampler>(0, 0.06);
    testSampleBoot<meder, NSampler>(0, 0.15, 100);
    testSampleTrimmed<NSampler>();
    testSampleBoot<trimmer, NSampler>(0, 0.05, 100);
    /*DEBUG("expo");
    testSampleNormal<ExpoSampler>();
    testSampleBoot<meaner, ExpoSampler>();
    testSampleMedian<ExpoSampler>();
    testSampleBoot<meder, ExpoSampler>();
    testSampleTrimmed<ExpoSampler>();
    testSampleBoot<trimmer, ExpoSampler>();*/
    DEBUG("Cauchy");
    testSampleNormal<CauchySampler>();
    testSampleBoot<meaner, CauchySampler>(0, 0.05, 100);
    testSampleMedian<CauchySampler>();
    testSampleBoot<meder, CauchySampler>(0, 0.05, 100);
    testSampleTrimmed<CauchySampler>(0, 0.1);
    testSampleBoot<trimmer, CauchySampler>(0, 0.05, 100);

    //DO BOOTSRAP ALSO EXTEND IT TO TAKE VECTOR OF SAMPLES?
}


void testTExact()
{
    int n = 30;
    meaner f;
    pair<double, pair<double, double> > target = testBootstrapHelper<NSampler, double>(f, n);
    DEBUG("Exact Simulation");
    DEBUG(target.first);
    DEBUG(target.second.first);
    DEBUG(target.second.second);
    NSampler s;
    Vector<double> samples(n);
    for(int j = 0; j < n; ++j) samples[j] = s();
    double stat = f(samples);
    double left = target.first + stat - target.second.second;
    double right = target.first + stat - target.second.first;
    DEBUG(left);
    DEBUG(right);
    pair<double, double> tconf = getTConf(samples);
    DEBUG(tconf.first);
    DEBUG(tconf.second);
}

double confError(double stat, pair<double, double> conf)
{
    return max(stat - conf.first, conf.second - stat);
}



pair<double, double> getNormalConf(IncrementalStatistics const&s, double a = 0.05)
{
    assert(s.n > 1 && a > 0 && a < 1);
    double ste = s.getStandardErrorSummary().stddev() *
        find2SidedConfZ(1 - a);
    return make_pair(s.getMean() - ste, s.getMean() + ste);
}
pair<double, double> getNormalConf(Vector<double> const& data, double a = 0.05)
{
    IncrementalStatistics s;
    for(int i = 0; i < data.getSize(); ++i) s.addValue(data[i]);
    return getNormalConf(s, a);
}

template<typename SAMPLER>
void testMeanConfs(Vector<Vector<string> >& matrix)
{
    int n = 30;
    typedef double X;
    typedef meaner F;
    F f;
    pair<double, pair<double, double> > target = testBootstrapHelper<SAMPLER, X>(f, n);
    DEBUG("Exact Simulation");
    DEBUG(target.first);
    DEBUG(target.second.first);
    DEBUG(target.second.second);
    SAMPLER s;
    BootstrapSimulationResult p("T"), bp("Normal"), bc("Permutation"),
        bt("Bootstrap-t-capped");
    for(int i = 0; i < 10000; ++i)
    {
        Vector<X> samples(n);
        for(int j = 0; j < n; ++j) samples[j] = s();
        double stat = f(samples);
        double left = target.first + stat - target.second.second;
        double right = target.first + stat - target.second.first;
        pair<double, double> exact(left, right);
        BasicBooter<F, X> booter(samples);
        p.addValue(getTConf(samples), target.first, exact, stat);
        bp.addValue(getNormalConf(samples), target.first, exact, stat);
        //booter.resample = booter.data;//reset
        //bm.addValue(bootstrapMixedInterval(booter), target.first, exact, stat);
        booter.resample = booter.data;//reset
        bc.addValue(permutationLocationConf<F>(samples), target.first, exact, stat);
        booter.resample = booter.data;//reset
        bt.addValue(bootstrapTIntervalCapped(booter), target.first, exact, stat);
    }
    p.print(matrix);
    bp.print(matrix);
    bc.print(matrix);
    //bmr.print(matrix);
    //bm.print(matrix);
    bt.print(matrix);
    //br.print(matrix);
}

void testMeanConfsDriver()
{
    Vector<Vector<string> > matrix;
    DEBUG("Normal Mean");
    matrix.append(Vector<string>());
    matrix.lastItem().append("Normal Mean");
    testMeanConfs<NSampler>(matrix);
    DEBUG("LogNormal Mean");
    matrix.append(Vector<string>());
    matrix.lastItem().append("LogNormal Mean");
    testMeanConfs<LogNormalSampler>(matrix);
    DEBUG("2-Normal Mix Mean");
    matrix.append(Vector<string>());
    matrix.lastItem().append("2-Normal Mix Mean");
    testMeanConfs<NMixSampler>(matrix);

    makeConfsReport("Mean", matrix);
}


template<typename SAMPLER>
void testMedianConfs(Vector<Vector<string> >& matrix)
{
    int n = 30;
    typedef double X;
    typedef meder F;
    F f;
    pair<double, pair<double, double> > target = testBootstrapHelper<SAMPLER, X>(f, n);
    DEBUG("Exact Simulation");
    DEBUG(target.first);
    DEBUG(target.second.first);
    DEBUG(target.second.second);
    SAMPLER s;
    BootstrapSimulationResult p("Nonp"), bc("Permutation"),
        bt("Bootstrap-t-capped");
    for(int i = 0; i < 10000; ++i)
    {
        Vector<X> samples(n);
        for(int j = 0; j < n; ++j) samples[j] = s();
        double stat = f(samples);
        double left = target.first + stat - target.second.second;
        double right = target.first + stat - target.second.first;
        pair<double, double> exact(left, right);
        BasicBooter<F, X> booter(samples);
        p.addValue(quantileConf(samples), target.first, exact, stat);

        //bm.addValue(bootstrapMixedInterval(booter), target.first, exact, stat);
        booter.resample = booter.data;//reset
        bc.addValue(permutationLocationConf<F>(samples), target.first, exact, stat);
        booter.resample = booter.data;//reset
        bt.addValue(bootstrapTIntervalCapped(booter), target.first, exact, stat);
    }
    p.print(matrix);
    bc.print(matrix);
    //bmr.print(matrix);
    //bm.print(matrix);
    bt.print(matrix);
    //br.print(matrix);
}

void testMedianConfsDriver()
{
    Vector<Vector<string> > matrix;
    DEBUG("Normal Median");
    matrix.append(Vector<string>());
    matrix.lastItem().append("Normal Median");
    testMedianConfs<NSampler>(matrix);
    DEBUG("Levy Median");
    matrix.append(Vector<string>());
    matrix.lastItem().append("Levy Median");
    testMedianConfs<LevySampler>(matrix);
    DEBUG("2-Normal Mix Median");
    matrix.append(Vector<string>());
    matrix.lastItem().append("2-Normal Mix Median");
    testMedianConfs<NMixSampler>(matrix);

    makeConfsReport("Median", matrix);
}


template<typename SAMPLER>
void testTrimmedMeanConfs(Vector<Vector<string> >& matrix)
{
    int n = 30;
    typedef double X;
    typedef trimmer F;
    F f;
    pair<double, pair<double, double> > target = testBootstrapHelper<SAMPLER, X>(f, n);
    DEBUG("Exact Simulation");
    DEBUG(target.first);
    DEBUG(target.second.first);
    DEBUG(target.second.second);
    SAMPLER s;
    BootstrapSimulationResult p("Normal"), bc("Permutation"),
        bt("Bootstrap-t-capped");
    for(int i = 0; i < 10000; ++i)
    {
        Vector<X> samples(n);
        for(int j = 0; j < n; ++j) samples[j] = s();
        double stat = f(samples);
        double left = target.first + stat - target.second.second;
        double right = target.first + stat - target.second.first;
        pair<double, double> exact(left, right);
        BasicBooter<F, X> booter(samples);
        p.addValue(trimmedMeanConf(samples), target.first, exact, stat);

        //bm.addValue(bootstrapMixedInterval(booter), target.first, exact, stat);
        booter.resample = booter.data;//reset
        bc.addValue(permutationLocationConf<F>(samples), target.first, exact, stat);
        booter.resample = booter.data;//reset
        bt.addValue(bootstrapTIntervalCapped(booter), target.first, exact, stat);
    }
    p.print(matrix);
    bc.print(matrix);
    //bmr.print(matrix);
    //bm.print(matrix);
    bt.print(matrix);
    //br.print(matrix);
}

void testTrimmedMeanConfsDriver()
{
    Vector<Vector<string> > matrix;
    DEBUG("Normal TrimmedMean");
    matrix.append(Vector<string>());
    matrix.lastItem().append("Normal TrimmedMean");
    testTrimmedMeanConfs<NSampler>(matrix);
    DEBUG("LogNormal TrimmedMean");
    matrix.append(Vector<string>());
    matrix.lastItem().append("LogNormal TrimmedMean");
    testTrimmedMeanConfs<LogNormalSampler>(matrix);
    DEBUG("Levy TrimmedMean");
    matrix.append(Vector<string>());
    matrix.lastItem().append("Levy TrimmedMean");
    testTrimmedMeanConfs<LevySampler>(matrix);
    DEBUG("2-Normal Mix TrimmedMean");
    matrix.append(Vector<string>());
    matrix.lastItem().append("2-Normal Mix TrimmedMean");
    testTrimmedMeanConfs<NMixSampler>(matrix);

    makeConfsReport("TrimmedMean", matrix);
}

//test with multisample
double medianDiff(Vector<Vector<double> > const& data)
{
    assert(data.getSize() == 2);
    return median(data[0]) - median(data[1]);
}

pair<double, double> medianIQR(Vector<double> const& data)
{//scaled to match normal stdev //is IRQ 33% efficient relative to std?
    assert(data.getSize() > 2);
    quickSort(data.getArray(), data.getSize());
    double q025 = quantile(data, 0.25, true),
        q075 = quantile(data, 0.75, true);
    return make_pair(quantile(data, 0.5, true), 0.7413 * (q075 - q025));
}

pair<double, double> trimmedMeanTMAD(Vector<double> data)
{//scaled to match normal stdev
    assert(data.getSize() > 2);
    double tm = trimmedMean(data);
    for(int i = 0; i < data.getSize(); ++i) data[i] = abs(data[i] - tm);
    return make_pair(tm, trimmedMean(data));
}

pair<double, double> HoefdingConf(double ave, int n, double confidence = 0.95)
{//Hoeffding inequality
    assert(n > 1 && confidence > 0 && confidence < 1);
    double a = (1 - confidence)/2, d = sqrt(log(1/a)/2/n);
    return make_pair(max(0.0, ave - d), min(1.0, ave + d));
}

template<typename SAMPLER>
void testProportionConfs(Vector<Vector<string> >& matrix)
{
    int n = 30;
    typedef double X;
    typedef meaner F;
    F f;
    pair<double, pair<double, double> > target = testBootstrapHelper<SAMPLER, X>(f, n);
    DEBUG("Exact Simulation");
    DEBUG(target.first);
    DEBUG(target.second.first);
    DEBUG(target.second.second);
    SAMPLER s;
    BootstrapSimulationResult p("Wilson"), bp("Normal"), be("Hoeffding"), bho("Hoef_Ori");
    for(int i = 0; i < 10000000; ++i)
    {
        Vector<X> samples(n);
        for(int j = 0; j < n; ++j) samples[j] = s();
        double stat = f(samples);
        double left = target.first + stat - target.second.second;
        double right = target.first + stat - target.second.first;
        pair<double, double> exact(left, right);

        p.addValue(wilsonScoreInterval(stat, n), target.first, exact, stat);
        bho.addValue(HoefFunctor::conf(stat, n), target.first, exact, stat);
        bp.addValue(getNormalConf(samples), target.first, exact, stat);
        be.addValue(HoefdingConf(stat, n), target.first, exact, stat);
    }
    p.print(matrix);
    bp.print(matrix);
    be.print(matrix);
    bho.print(matrix);
}
struct BerSampler05
{
    double operator()()const{return GlobalRNG().bernoulli(0.5);}
};
struct BerSampler01
{
    double operator()()const{return GlobalRNG().bernoulli(0.1);}
};
struct BerSampler001
{
    double operator()()const{return GlobalRNG().bernoulli(0.01);}
};
void testProportionConfsDriver()
{
    Vector<Vector<string> > matrix;
    DEBUG("Bernoulli05 Mean");
    matrix.append(Vector<string>());
    matrix.lastItem().append("Bernoulli05 Mean");
    testProportionConfs<BerSampler05>(matrix);
    DEBUG("Bernoulli01 Mean");
    matrix.append(Vector<string>());
    matrix.lastItem().append("Bernoulli01 Mean");
    testProportionConfs<BerSampler01>(matrix);
    DEBUG("Bernoulli001 Mean");
    matrix.append(Vector<string>());
    matrix.lastItem().append("Bernoulli001 Mean");
    testProportionConfs<BerSampler001>(matrix);

    makeConfsReport("Proportion", matrix);
}

void testHoefOriConf()
{
    double ave = 0.5;
    int n = 100;
    DEBUG("HOEF");
    double h = sqrt(log(1/0.025)/2/n);
    DEBUG(ave - h);
    DEBUG(ave + h);
    pair<double, double> conf = HoefFunctor::conf(ave, n);
    DEBUG("HORI");
    DEBUG(conf.first);
    DEBUG(conf.second);
    DEBUG("WILSON");
    conf = wilsonScoreInterval(ave, n);
    DEBUG(conf.first);
    DEBUG(conf.second);
}

int main(int argc, char *argv[])
{
    testT();
    return 0;
    testProportionConfsDriver();
    return 0;
    testHoefOriConf();
    return 0;
    testTrimmedMeanConfsDriver();
    return 0;
    testMedianConfsDriver();
    return 0;
    testMeanConfsDriver();
    return 0;
    testPairedTests();
    return 0;
    test1SampleConf();
    return 0;
    test2SampleDiff();
    return 0;
    testTExact();
    return 0;
    testQuantiles();
    return 0;
    testTrimmedMean();
    return 0;
    testNormalEval();
    return 0;
    DEBUG(approxNormal2SidedConf(2));
    DEBUG(signTestAreEqual(7, 10, 0.5));
    return 0;
}
