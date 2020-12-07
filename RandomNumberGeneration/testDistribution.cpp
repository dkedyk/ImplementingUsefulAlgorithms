#include "DistributionTests.h"
#include "Statistics.h"
#include "testCommon.h"
#include "MCMC.h"
using namespace igmdk;

void testChi2Eval()
{//should be 0.95 for all
    DEBUG(evaluateChiSquaredCdf(3.84, 1));
    DEBUG(evaluateChiSquaredCdf(11.1, 5));
    DEBUG(evaluateChiSquaredCdf(18.3, 10));
    DEBUG(evaluateChiSquaredCdf(31.4, 20));
    DEBUG(evaluateChiSquaredCdf(124, 100));
}

void testChiSquared()
{
    int t = 10000, n = 128;
    assert(isPowerOfTwo(n));
    IncrementalStatistics s;
    for(int j = 0; j < t; ++j)
    {
        Vector<double> means;
        //geometic05
        for(int count = 64; count >= 8; count /= 2)
        {
            means.append(count);
        }
        Vector<int> counts(means.getSize());
        for(int i = 0; i < n; ++i)
        {//power not super good - 0.55 also seems to pass the tests with 128
            //int var = GlobalRNG().geometric(0.55);
            int var = GlobalRNG().geometric05();
            if(var - 1 < counts.getSize()) ++counts[var - 1];//geom starts from 1
        }
        s.addValue(chiSquaredP(counts, means) <= 0.05);
    }
    DEBUG(s.getMean());
    pair<double, double> conf = wilsonScoreInterval(s.getMean(), s.n);
    DEBUG(conf.first);
    DEBUG(conf.second);
    //ok for this to fail?
    //move uncertainly in % estimate to skepticism?
    assert(conf.first <= 0.05);
}

void testChiSquaredUniform10()
{
    int t = 10000, n = 100;
    IncrementalStatistics s;
    for(int j = 0; j < t; ++j)
    {
        Vector<double> means(10, n/10);
        Vector<int> counts(means.getSize());
        for(int i = 0; i < n; ++i) ++counts[GlobalRNG().mod(counts.getSize())];
        s.addValue(chiSquaredP(counts, means) <= 0.05);
    }
    DEBUG(s.getMean());
    pair<double, double> conf = wilsonScoreInterval(s.getMean(), s.n);
    DEBUG(conf.first);
    DEBUG(conf.second);
    //ok for this to fail?
    //move uncertainly in % estimate to skepticism?
    assert(conf.first <= 0.05);
}

template<typename CDF, typename SAMPLER>
void testDKWProportion(int t = 10000, int n = 100)
{
    IncrementalStatistics s;
    CDF c;
    SAMPLER sa;
    for(int j = 0; j < t; ++j)
    {
        Vector<double> samples;
        for(int i = 0; i < n; ++i) samples.append(sa());
        s.addValue(DKWPValue(samples, c) <= 0.05);
    }
    DEBUG(s.getMean());
    pair<double, double> conf = wilsonScoreInterval(s.getMean(), s.n);
    DEBUG(conf.first);
    DEBUG(conf.second);
    //ok for this to fail?
    //move uncertainly in % estimate to skepticism?
    assert(conf.first <= 0.05);
}

struct MCMCN01Sampler
{
    struct N01PDF
    {
        double operator()(double x)const{return exp(-x * x/2);}
    };
    mutable GridRWM<N01PDF> s;
    double operator()()const{return s.sample();}
};

struct Uni01CDF
{
    double operator()(double x)const
    {
        if(x <= 0) return 0;
        else if(x >= 1) return 1;
        else return x;
    }
};

void testDistroMatch()
{//do this for most distros as unit test?
    DEBUG("testChiSquared()");
    testChiSquared();
    DEBUG("testChiSquaredUniform10()");
    testChiSquaredUniform10();
    DEBUG("testDKWProportion<Uni01CDF, USampler>()");
    testDKWProportion<Uni01CDF, USampler>();
    DEBUG("testDKWProportion<Norm01CDF, NSampler>()");
    testDKWProportion<Norm01CDF, NSampler>();
    //ggggggggggggggets under 10% o not perfect sampler
    //DEBUG("testDKWProportion<Norm01CDF, MCMCN01Sampler>()");
    //testDKWProportion<Norm01CDF, MCMCN01Sampler>();
    DEBUG("testDKWProportion<T1001CDF, T10Sampler>()");
    testDKWProportion<T1001CDF, T10Sampler>();
    DEBUG("testDKWProportion<Chi10CDF, Chi10Sampler>()");
    testDKWProportion<Chi10CDF, Chi10Sampler>();
}

void testKSEvals()
{
    double corr = exp(-4) - exp(-9);
    DEBUG(corr/exp(-1));
    DEBUG(exp(-16));
    double corr2 = exp(-16) - exp(-25);
}


void testKS()
{
    Vector<double> a, b;
    /*a.append(-0.15);
    a.append(8.60);
    a.append(5.00);
    a.append(3.71);
    a.append(4.29);
    a.append(7.74);
    a.append(2.48);
    a.append(3.25);
    a.append(-1.15);
    a.append(8.38);

    b.append(2.55);
    b.append(12.07);
    b.append(0.46);
    b.append(0.35);
    b.append(2.69);
    b.append(-0.94);
    b.append(1.73);
    b.append(0.73);
    b.append(-0.35);
    b.append(-0.37);*/

    a.append(-2.18);
    a.append(-1.79);
    a.append(-1.66);
    a.append(-0.65);
    a.append(-0.05);
    a.append(0.54);
    a.append(0.85);
    a.append(1.69);

    b.append(-1.91);
    b.append(-1.22);
    b.append(-0.96);
    b.append(-0.72);
    b.append(0.14);
    b.append(0.82);
    b.append(1.45);
    b.append(1.86);


    DEBUG(KS2SamplePValue(a, b));
}

int main(int argc, char *argv[])
{
    testDistroMatch();
    return 0;
    testKSEvals();
    return 0;
    testChi2Eval();
    return 0;
    testKS();
    return 0;
}
