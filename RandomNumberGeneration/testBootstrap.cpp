#include "Bootstrap.h"
#include "../NumericalMethods/NumericalMethods.h"
#include "../Utils/Debug.h"
#include "../ExternalMemoryAlgorithms/CSV.h"
#include "testCommon.h"
using namespace igmdk;

void testBootstrap()
{//book examples
    {
        DEBUG("mean of uniform 1000");
        Vector<double> data;
        meaner med;
        for(int i = 0; i < 1000; ++i) data.append(GlobalRNG().uniform01());
        double estimate = med(data);
        DEBUG(estimate);
        IncrementalStatistics s;
        for(int i = 0; i < data.getSize(); ++i) s.addValue(data[i]);
        DEBUG(s.getStandardErrorSummary().error95());

        BasicBooter<meaner> booter(data);
        pair<double, double> conf = bootstrapNormalInterval(booter);
        DEBUG((conf.second - conf.first)/2);

    }
    {
        DEBUG("median of uniform 1000");
        Vector<double> data;
        meder med;
        for(int i = 0; i < 1000; ++i) data.append(GlobalRNG().uniform01());
        double estimate = med(data);
        DEBUG(estimate);
        pair<double, double> medianConf = quantileConf(data);
        DEBUG(medianConf.first - estimate);
        DEBUG(medianConf.second - estimate);

        BasicBooter<meder> booter(data);
        pair<double, double> conf = bootstrapNormalInterval(booter);
        DEBUG((conf.second - conf.first)/2);
    }
}

struct BootstrapResult
{
    double fValue, bias, std, iFa, iF1ma;
    double biasFactor(){return bias/std;}
    pair<double, double> normalInterval(double z = 2)const
        {return make_pair(fValue - z * std, fValue + z * std);}
    pair<double, double> normalBiasAdjustedInterval(double z = 2)const
    {
        pair<double, double> result = normalInterval(z);
        result.first -= bias;
        result.second -= bias;
        return result;
    }
    //not presented, for testing only
    pair<double, double> percentileInterval()const
        {return make_pair(iFa, iF1ma);}
    pair<double, double> pivotalInterval()const
        {return make_pair(2 * fValue - iF1ma, 2 * fValue - iFa);}
};
template<typename BOOTER> BootstrapResult bootstrap(BOOTER& booter,
    int b = 3000, double confidence = 0.95)
{
    assert(b > 2);
    int tailSize = b * (1 - confidence)/2;
    if(tailSize < 1) tailSize = 1;
    if(tailSize > b/2 - 1) tailSize = b/2 - 1;
    //max heap to work with the largest value
    Heap<double, ReverseComparator<double> > left;
    Heap<double> right;//min heap for the smallest value
    double q = booter.eval();
    IncrementalStatistics s;
    for(int i = 0; i < b; ++i)
    {
        booter.boot();//resample
        double value = booter.eval();
        s.addValue(value);
        if(left.getSize() < tailSize)
        {//heaps not full
            left.insert(value);
            right.insert(value);
        }
        else
        {//heaps full - replace if more extreme
            if(value < left.getMin()) left.changeKey(0, value);
            else if(value > right.getMin()) right.changeKey(0, value);
        }
    }
    BootstrapResult r = {q, s.getMean() - q, s.stdev(),
        left.getMin(), right.getMin()};
    return r;
}

template<typename BOOTER>
double findBCAAccelerationConstant(BOOTER& booter)
{//for analytical calculators can specialize as well
    int n = booter.data.getSize();
    Vector<double> values(n);
    IncrementalStatistics s;
    booter.resample = booter.data;
    booter.resample.removeLast();
    for(int i = 0; i < n; ++i)
    {
        if(i != n - 1) booter.resample[i] = booter.data[n - 1];
        values[i] = booter.f(booter.resample);
        if(i != n - 1) booter.resample[i] = booter.data[i];
        s.addValue(values[i]);
    }//restore resample
    booter.resample.append(booter.data[n - 1]);
    double cubeSum = 0, squareSum = 0;
    for(int i = 0; i < n; ++i)
    {
        double diff = (s.getMean() - values[i]);
        cubeSum += diff * diff * diff;
        squareSum += diff * diff;
    }
    double result = cubeSum/6/pow(squareSum, 1.5);
    return isfinite(result) ? result : 0;//check for 0 in denominator
}
pair<double, double> bootstrapBCIntervalHelper(Vector<double> const& values,
    double q, int b = 3000, double confidence = 0.95, double acc = 0)
{
    int nSmaller = 0;
    for(int i = 0; i < b; ++i) if(values[i] < q) ++nSmaller;
    double tail = (1 - confidence)/2, tails[2], p0 = (nSmaller + 1) * 1.0/
        (b + 2), z0 = invertCDF([](double x){return approxNormalCDF(x);}, p0);
    for(int i = 0; i < 2; ++i)
    {
        double temp = z0 + invertCDF([](double x){return approxNormalCDF(x);},
            i == 0 ? tail : 1 - tail);
        tails[i] = approxNormalCDF(z0 + temp/(1 - acc * temp));
    }
    int left = max<int>(0, tails[0] * b),
        right = min<int>(b - 1, tails[1] * b);
    return make_pair(values[left], values[right]);
}
template<typename BOOTER> pair<double, double> bootstrapBCInterval(
    BOOTER& booter, int b = 3000, double confidence = 0.95, double acc = 0)
{//computes BC interval when acc = 0
    assert(b > 2 && isfinite(acc));
    double q = booter.eval();
    int nSmaller = 0;
    Vector<double> values(b);
    for(int i = 0; i < b; ++i)
    {
        booter.boot();//resample
        values[i] = booter.eval();
        if(values[i] < q) ++nSmaller;
    }
    quickSort(values.getArray(), b);//safeguard estimate
    return bootstrapBCIntervalHelper(values, q, b, confidence, acc);
}
template<typename BOOTER> pair<double, double> bootstrapBCaInterval(
    BOOTER& booter, int b = 3000, double confidence = 0.95)
{
    return bootstrapBCInterval(booter, b, confidence,
        findBCAAccelerationConstant(booter));
}

template<typename BOOTER> pair<double, double> bootstrapTInterval(
    BOOTER& booter, int b = 1000, int b2 = 50, double confidence = 0.95)
{
    assert(b > 2);
    int tailSize = b * (1 - confidence)/2;
    if(tailSize < 1) tailSize = 1;
    if(tailSize > b/2 - 1) tailSize = b/2 - 1;
    //max heap to work with the largest value
    Heap<double, ReverseComparator<double> > left;
    Heap<double> right;//min heap for the smallest value
    double q = booter.eval();
    IncrementalStatistics s;
    for(int i = 0; i < b; ++i)
    {
        booter.boot();//resample
        BOOTER booter2 = booter.cloneForNested();
        double valStat = booter.eval(), stdeInner = bootstrapStde(booter2, b2),
            value = (valStat - q)/stdeInner;
        if(!isfinite(value)) continue; //possible division by 0
        s.addValue(valStat);
        if(left.getSize() < tailSize)
        {//heaps not full
            left.insert(value);
            right.insert(value);
        }
        else
        {//heaps full - replace if more extreme
            if(value < left.getMin()) left.changeKey(0, value);
            else if(value > right.getMin()) right.changeKey(0, value);
        }
    }//beware - get standard deviation here not standard error
    //protect against 0 standard deviation just in case
    double stde = s.n > 2 ? s.stdev() : 0, leftPoint = right.getSize() > 0 ?
        q - stde * right.getMin() : q,//UNINTUITIVE FIX SIGN BELOW
        rightPoint = left.getSize() > 0 ? q - stde * left.getMin() : q;
    return make_pair(leftPoint, rightPoint);
}

class BTDCFunctor
{
    Vector<IncrementalStatistics> sets;
    bool findLeft;
    double q, aTarget;
    pair<double, double> makeInterval(int i, double a)const
    {
        assert(0 <= i && i < sets.getSize());
        double m = sets[i].getMean(), delta = sets[i].stdev() * sqrt(1/a);
        return make_pair(m - delta, m + delta);
    }
public:
    BTDCFunctor(double theQ, double a): q(theQ), aTarget(a), findLeft(true){}
    void addSet(){sets.append(IncrementalStatistics());}
    void addValue(double value)
    {
        assert(sets.getSize() > 0);
        sets.lastItem().addValue(value);
    }
    void flipSide(){findLeft = !findLeft;}
    double operator()(double a)const
    {//returns < 0 if under target
        int missCount = 0, b = sets.getSize();
        for(int i = 0; i < b; ++i)
        {//make sure that don't confuse miss left and miss right
            pair<double, double> conf = makeInterval(i, a);
            missCount += (findLeft ? q < conf.first : conf.second < q);
        }//each side misses half
        return aTarget/2 - missCount * 1.0/b;
    }
    double findA()const
    {
        double left = aTarget, right = 1, yLeft = this->operator()(left);
        return haveDifferentSign(yLeft, this->operator()(right)) ?
            solveFor0(*this, left, right).first :
            yLeft < 0 ? left : right;
    }
};
template<typename BOOTER> pair<double, double>
    bootstrapDoubleChebyshevInterval(BOOTER& booter, int b = 1000,
    int b2 = 50, double confidence = 0.95)
{
    assert(b > 2 && confidence > 0 && confidence < 1);
    double q = booter.eval();
    BTDCFunctor f(q, 1 - confidence);
    IncrementalStatistics s;
    for(int i = 0; i < b; ++i)
    {
        booter.boot();//resample
        s.addValue(booter.eval());
        BOOTER booter2 = booter.cloneForNested();
        f.addSet();
        for(int j = 0; j < b; ++j)
        {
            booter2.boot();//resample nested
            f.addValue(booter2.eval());
        }
    }
    double stdev = s.stdev(), aLeft = f.findA();
    f.flipSide();//default is left, next find right
    return make_pair(q - stdev * sqrt(1/aLeft), q + stdev * sqrt(1/f.findA()));
}

template<typename SAMPLER, typename F, typename X>
void testBootstrapConfs(Vector<Vector<string> >& matrix)
{
    int n = 5;
    F f;
    pair<double, pair<double, double> > target = testBootstrapHelper<SAMPLER, X>(f, n);
    DEBUG("Exact Simulation");
    DEBUG(target.first);
    DEBUG(target.second.first);
    DEBUG(target.second.second);
    SAMPLER s;
    BootstrapSimulationResult p("Pivotal"), bp("Percentile"), bn("Normal"), bmr("Chebyshev Double"),
        bm("Mixed"), bc("BC"), bca("BCa"), bt("Bootstrap-t"), btc("Bootstrap-t-capped"), br("Percentile Double");
    for(int i = 0; i < 30000; ++i)
    {
        Vector<X> samples(n);
        for(int j = 0; j < n; ++j) samples[j] = s();
        double stat = f(samples);
        double left = target.first + stat - target.second.second;
        double right = target.first + stat - target.second.first;
        pair<double, double> exact(left, right);
        BasicBooter<F, X> booter(samples);
        BootstrapResult r = bootstrap(booter);
        p.addValue(r.pivotalInterval(), target.first, exact, stat);
        bp.addValue(r.percentileInterval(), target.first, exact, stat);
        booter.resample = booter.data;//reset
        bn.addValue(bootstrapNormalInterval(booter), target.first, exact, stat);
        booter.resample = booter.data;//reset
        bm.addValue(bootstrapMixedInterval(booter), target.first, exact, stat);
        booter.resample = booter.data;//reset
        bc.addValue(bootstrapBCInterval(booter), target.first, exact, stat);
        booter.resample = booter.data;//reset
        bca.addValue(bootstrapBCaInterval(booter), target.first, exact, stat);
        booter.resample = booter.data;//reset
        bt.addValue(bootstrapTInterval(booter), target.first, exact, stat);
        booter.resample = booter.data;//reset
        btc.addValue(bootstrapTIntervalCapped(booter), target.first, exact, stat);
        booter.resample = booter.data;//reset
        br.addValue(bootstrapDoublePercentileInterval(booter), target.first, exact, stat);
        booter.resample = booter.data;//reset
        bmr.addValue(bootstrapDoubleChebyshevInterval(booter), target.first, exact, stat);
    }
    p.print(matrix);
    bp.print(matrix);
    bn.print(matrix);
    bm.print(matrix);
    bc.print(matrix);
    bca.print(matrix);
    bt.print(matrix);
    btc.print(matrix);
    br.print(matrix);
    bmr.print(matrix);
}

void testBootstrapConfsDriver()
{
    Vector<Vector<string> > matrix;
    DEBUG("Normal Mean");
    matrix.append(Vector<string>());
    matrix.lastItem().append("Normal Mean");
    testBootstrapConfs<NSampler, meaner, double>(matrix);
    DEBUG("Lognormal Mean");
    matrix.append(Vector<string>());
    matrix.lastItem().append("Lognormal Mean");
    testBootstrapConfs<LogNormalSampler, meaner, double>(matrix);
    DEBUG("Levy Median");
    matrix.append(Vector<string>());
    matrix.lastItem().append("Levy Median");
    testBootstrapConfs<LevySampler, meder, double>(matrix);
    DEBUG("Normal Stdev");
    matrix.append(Vector<string>());
    matrix.lastItem().append("Normal Stdev");
    testBootstrapConfs<NSampler, varer, double>(matrix);
    DEBUG("2-Normal Mix Mean");
    matrix.append(Vector<string>());
    matrix.lastItem().append("2-Normal Mix Mean");
    testBootstrapConfs<NMixSampler, meaner, double>(matrix);
    DEBUG("Normal Error Corr");
    matrix.append(Vector<string>());
    matrix.lastItem().append("Normal Error Corr");
    testBootstrapConfs<NErrorPairSampler, correr, pair<double, double> >(matrix);
    DEBUG("Normal Error Spearman");
    matrix.append(Vector<string>());
    matrix.lastItem().append("Normal Error Spearman");
    testBootstrapConfs<NErrorPairSampler, correrS, pair<double, double> >(matrix);
    /*DEBUG("Normal Broken Line(Mean)");//SAME AS NORMAL NO USE
    matrix.append(Vector<string>());
    matrix.lastItem().append("Normal Broken Line(Mean)");
    testBootstrapConfs<NSampler, funer, double>(matrix);*/
    /*DEBUG("Uniform Max");
    matrix.append(Vector<string>());
    matrix.lastItem().append("Uniform Max");
    testBootstrapConfs<USampler, maxer, double>(matrix);*/
    /*DEBUG("Levy Median Diff");
    matrix.append(Vector<string>());
    matrix.lastItem().append("Levy Median Diff");
    testBootstrapConfs<USampler, maxer, double>(matrix);*/
    /*DEBUG("Poisson10 Mean");
    matrix.append(Vector<string>());
    matrix.lastItem().append("Poisson10 Mean");
    testBootstrapConfs<PoissonSampler, meaner, double>(matrix);*/

    makeConfsReport("Bootstrap", matrix);
}

int main(int argc, char *argv[])
{
    testBootstrapConfsDriver();
    return 0;
    testBootstrap();
    return 0;
}
