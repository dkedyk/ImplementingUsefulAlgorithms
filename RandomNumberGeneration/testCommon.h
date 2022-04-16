#ifndef IGMDK_TEST_COMMON_H
#define IGMDK_TEST_COMMON_H
#include "Statistics.h"
#include "Correlation.h"
#include "DistributionTests.h"
#include "../ExternalMemoryAlgorithms/CSV.h"
namespace igmdk{

struct meder
{
    double operator()(Vector<double> observations)const
        {return median(observations);}
};

struct trimmer
{
    double operator()(Vector<double> observations)const
        {return trimmedMean(observations);}
};

struct varer
{
    double operator()(Vector<double> observations)const
    {
        IncrementalStatistics s;
        for(int i = 0; i < observations.getSize(); ++i)
            s.addValue(observations[i]);
        return s.stdev();
    }
};

struct meaner
{
    double operator()(Vector<double> const& observations)const
    {
        double result = observations[0];
        for(int i = 1; i < observations.getSize(); ++i)
        {
            result += observations[i];
        }
        return result / observations.getSize();
    }
};

struct maxer
{
    double operator()(Vector<double> const& observations)const
    {//assume uniform
        int n = observations.getSize();
        double result = observations[0];
        for(int i = 1; i < n; ++i)
        {
            if(observations[i] > result) result = observations[i];
        }
        return result * (n - 1)/n;
    }
};

struct correr
{
    double operator()(Vector<pair<double, double> > const& observations)const
    {
        //DEBUG(observations.getSize());
        /*for(int i = 0; i < observations.getSize(); ++i)
        {
            DEBUG(observations[i].first);
            DEBUG(observations[i].second);
        }*/
        //DEBUG(PearsonCorrelation(observations));
       // system("PAUSE");
        return PearsonCorrelation(observations);
    }
};

struct correrS
{
    double operator()(Vector<pair<double, double> > const& observations)const
    {
        //DEBUG(observations.getSize());
        /*for(int i = 0; i < observations.getSize(); ++i)
        {
            DEBUG(observations[i].first);
            DEBUG(observations[i].second);
        }*/
        //DEBUG(PearsonCorrelation(observations));
       // system("PAUSE");
        return SpearmanCorrelation(observations);
    }
};

struct funer
{
    double operator()(Vector<double> const& observations)const
    {
        double result = observations[0];
        for(int i = 1; i < observations.getSize(); ++i)
        {
            result += observations[i];
        }
        result /= observations.getSize();
        return result < 0 ? result : result * 2;
    }
};

struct USampler
{
    double operator()()const{return GlobalRNG().uniform01();}
};

struct Norm01CDF
{
    double operator()(double x)const
    {
        return approxNormalCDF(x);
    }
};
struct NSampler
{
    double operator()()const{return GlobalRNG().normal01();}
};
struct T1001CDF
{
    double operator()(double x)const
    {
        return approxTCDF(x, 10);
    }
};
struct T10Sampler
{
    double operator()()const{return GlobalRNG().t(10);}
};
struct Chi10CDF
{
    double operator()(double x)const
    {
        assert(x >= 0);
        return evaluateChiSquaredCdf(x, 10);
    }
};
struct Chi10Sampler
{
    double operator()()const{return GlobalRNG().chiSquared(10);}
};

struct ExpoSampler
{
    double operator()()const{return GlobalRNG().exponential(1);}
};
struct CauchySampler
{
    double operator()()const{return GlobalRNG().cauchy(0, 1);}
};

struct PoissonSampler
{
    double operator()()const{return GlobalRNG().poisson(10);}
};
struct LogNormalSampler
{
    double operator()()const{return GlobalRNG().logNormal(0, 1);}
};
struct LevySampler
{
    double operator()()const{return GlobalRNG().Levy();}
};
struct NMixSampler
{
    double operator()()const{return GlobalRNG().mod(5) == 0 ?
        GlobalRNG().normal01() : GlobalRNG().normal(1, 10);}
};
struct NErrorPairSampler
{
    pair<double, double> operator()()const
    {
        double x = GlobalRNG().uniform01();
        double y = x - 0.1 * x * x + GlobalRNG().normal(0, 0.1);
        return make_pair(x, y);
    }
};
struct CubicPairSampler
{
    pair<double, double> operator()()const
    {
        double x = GlobalRNG().uniform01();
        double y = x * x * x;
        return make_pair(x, y);
    }
};

template<typename SAMPLER, typename X, typename F>
pair<double, pair<double, double> > testBootstrapHelper(F const& f, int n,
    int nSim = 10000000)
{
    SAMPLER s;
    Vector<X> samples2(nSim);
    for(int i = 0; i < nSim; ++i) samples2[i] = s();
    double stat = f(samples2);
    Vector<double> samples(nSim);
    for(int i = 0; i < nSim; ++i)
    {
        Vector<X> samplesN(n);
        for(int j = 0; j < n; ++j) samplesN[j] = s();
        samples[i] = f(samplesN);
    }
    pair<double, double> quantiles = make_pair(quantile(samples, 0.025),
        quantile(samples, 0.975));
    return make_pair(stat, quantiles);
}

double chebyshevLoss(double a, double length){return length/sqrt(1/a);}
double expLoss(double a, double aTarget, double length)
    {return length * exp(max(0.0, a/aTarget - 1));}
double rootMeanSquare(double al, double ar)
    {return sqrt((al * al + ar * ar)/2);}
double chebyshevLossAsymetric(double al, double ar, double length)
{
    return chebyshevLoss(rootMeanSquare(al, ar), length);
}
struct BootstrapSimulationResult
{
    string name;
    IncrementalStatistics cl, cr, l, ll, lr;
    BootstrapSimulationResult(string const& theName): name(theName){}
    void addValue(pair<double, double> const& conf, double target,
        pair<double, double> const& exact, double stat)
    {
        cl.addValue(target < conf.first);
        cr.addValue(target > conf.second);
        l.addValue((conf.second - conf.first)/(exact.second - exact.first));
        ll.addValue(max(0.0, target - conf.first)/(exact.second - exact.first));
        lr.addValue(max(0.0, conf.second - target)/(exact.second - exact.first));
        //sum of above will be penalized length, which is interesting concept in itself
    }
    void print(Vector<Vector<string> >& matrix)const
    {
        double aTarget = 0.5;
        DEBUG(name);
        matrix.lastItem().append(name);
        double c = cl.getMean() + cr.getMean();
        DEBUG(c);
        matrix.lastItem().append(toStringDouble(c));
        DEBUG(cl.getMean());
        matrix.lastItem().append(toStringDouble(cl.getMean()));
        DEBUG(cr.getMean());
        matrix.lastItem().append(toStringDouble(cr.getMean()));
        DEBUG(l.getMean() - 1);
        matrix.lastItem().append(toStringDouble(l.getMean() - 1));
        double rms = rootMeanSquare(cl.getMean(), cr.getMean());
        /*DEBUG(2 * rms);
        matrix.lastItem().append(toStringDouble(2 * rms));*/
        /*DEBUG(chebyshevLoss(2 * rms, l.getMean()));
        matrix.lastItem().append(toStringDouble(
            chebyshevLoss(2 * rms, l.getMean())));
        DEBUG(linearLoss(2 * rms, aTarget, l.getMean()));
        matrix.lastItem().append(toStringDouble(
            linearLoss(2 * rms, aTarget, l.getMean())));*/
        double chebL = chebyshevLoss(2 * cl.getMean(), ll.getMean()) +
            chebyshevLoss(2 * cr.getMean(), lr.getMean());
        DEBUG(chebL);
        matrix.lastItem().append(toStringDouble(chebL));
        double expL = expLoss(2 * cl.getMean(), aTarget, ll.getMean()) +
            expLoss(2 * cr.getMean(), aTarget, lr.getMean());
        DEBUG(expL);
        matrix.lastItem().append(toStringDouble(expL));
    }
};

void makeConfsReport(string const& filenamePrefix,
    Vector<Vector<string> > const& matrix)
{
    int reportNumber = time(0);
    string filename = filenamePrefix + "Confs" + to_string(reportNumber) +
        ".csv";
    createCSV(matrix, filename.c_str());
    Vector<string> names;
    names.append("A-level");
    names.append("Left A-level");
    names.append("Right A-level");
    names.append("Percent Longer");
    names.append("Asymmetric Cheb Loss");
    names.append("Asymmetric Exp Loss");
    createAugmentedCSVFiles(matrix, names, filename);
}

}//end namespace
#endif

