#ifndef IGMDK_INTEGRATION_H
#define IGMDK_INTEGRATION_H
#include <cmath>
#include <iomanip>
#include "../Utils/Bits.h"
#include "../Utils/Vector.h"
#include "../Utils/Utils.h"
#include "../Sorting/Sort.h"
#include "../Utils/Queue.h"
#include "../Heaps/Heap.h"
#include "../RandomNumberGeneration/Random.h"
#include "NumericalCommon.h"
#include "Interpolation.h"
namespace igmdk{

template<typename FUNCTION> pair<double, double> integrateCC(FUNCTION const& f,
    double a, double b, int maxEvals = 5000, int minEvals = 17)
{
    ScaledChebAB c(adaptiveChebEstimate(ScaledFunctionM11<FUNCTION>(a, b, f),
        maxEvals, minEvals), a, b);
    return c.integrate();
}

class IntervalGaussLobattoKronrod
{
    double a, b, y[7];
    void generateX(double x[7])const
    {
        x[0] = -1; x[1] = -sqrt(2.0/3); x[2] = -1/sqrt(5); x[3] = 0;
        x[4] = -x[2]; x[5] = -x[1]; x[6] = -x[0];
    }
    template<typename FUNCTION> void initHelper(FUNCTION const& f, double fa,
        double fb)
    {
        assert(a < b);
        y[0] = fa;
        y[6] = fb;
        double x[7];
        generateX(x);
        ScaledFunctionM11<FUNCTION> fM11(a, b, f);
        for(int i = 1; i < 6; ++i) y[i] = fM11(x[i]);
    }
    double integrateGL()const
    {
        double w[7] = {1.0/6, 0, 5.0/6, 0};
        w[4] = w[2]; w[5] = w[1]; w[6] = w[0];
        double result = 0;
        for(int i = 0; i < 7; ++i) result += w[i] * y[i];
        return result * (b - a)/2;
    }
    template<typename FUNCTION> IntervalGaussLobattoKronrod(FUNCTION const& f,
        double theA, double theB, double fa, double fb): a(theA), b(theB)
        {initHelper(f, fa, fb);}
public:
    template<typename FUNCTION> IntervalGaussLobattoKronrod(FUNCTION const& f,
        double theA, double theB): a(theA), b(theB)
        {initHelper(f, f(a), f(b));}
    double integrate()const
    {
        double w[7] = {11.0/210, 72.0/245, 125.0/294, 16.0/35};
        w[4] = w[2]; w[5] = w[1]; w[6] = w[0];
        double result = 0;
        for(int i = 0; i < 7; ++i) result += w[i] * y[i];
        return result * (b - a)/2;//need to scale back
    }
    double length()const{return b - a;}
    double error()const{return abs(integrate() - integrateGL());}
    template<typename FUNCTION>
    Vector<IntervalGaussLobattoKronrod> split(FUNCTION const& f)const
    {
        Vector<IntervalGaussLobattoKronrod> result;
        double x[7];
        generateX(x);
        for(int i = 0; i < 7; ++i) x[i] = a + (x[i] - -1) * (b - a)/2;
        for(int i = 0; i < 6; ++i) result.append(
            IntervalGaussLobattoKronrod(f, x[i], x[i + 1], y[i], y[i + 1]));
        return result;
    }
    static int initEvals(){return 7;}
    static int splitEvals(){return 30;}
};

/*
Want the most robust integrator that does best precision in allowed num of
evals upto traget precision
So use abs precision not sum that cancels out in other rules and dont strive
for super small num of evals -- want precision instead to each small
subinterval must have small enough length
*/

template<typename INTEGRATION_INTERVAL, typename FUNCTION> pair<double, double>
    integrateAdaptiveHeap(FUNCTION const& f, double a, double b, double
    eRelAbs = highPrecEps, int maxEvals = 1000000, int minEvals = -1)
{
    typedef INTEGRATION_INTERVAL II;
    if(minEvals == -1) minEvals = sqrt(maxEvals);
    assert(a < b && maxEvals >= minEvals && minEvals >= II::initEvals());
    II i0(f, a, b);
    double result = i0.integrate(), totalError = i0.error();
    AdaptiveIntevalComparator<II> ic = {(b - a)/minEvals};
    Heap<II, AdaptiveIntevalComparator<II> > h(ic);
    h.insert(i0);
    for(int usedEvals = II::initEvals();
        usedEvals < minEvals || (usedEvals + II::splitEvals() <= maxEvals &&
        !isEEqual(result, result + totalError, eRelAbs));
        usedEvals += II::splitEvals())
    {
        II next = h.deleteMin();
        Vector<II> division = next.split(f);
        result -= next.integrate();
        totalError -= next.error();
        for(int i = 0; i < division.getSize(); ++i)
        {
            h.insert(division[i]);
            result += division[i].integrate();
            totalError += division[i].error();
        }
    }
    return make_pair(result, totalError);
}
template<typename FUNCTION> struct SingularityWrapper
{
    FUNCTION f;
    mutable int sCount;
    SingularityWrapper(FUNCTION const& theF = FUNCTION()): f(theF), sCount(0){}
    double operator()(double x)const
    {
        double y = f(x);
        if(isfinite(y)) return y;
        else
        {
            ++sCount;
            return 0;
        }
    }
};

template<typename FUNCTION> pair<double, double> integrateHybrid(
    FUNCTION const& f, double a, double b, double eRelAbs = highPrecEps,
    int maxEvals = 1000000, int minEvals = -1)
{
    int CCEvals = min(maxEvals/2, 1000);
    pair<double, double> resultCC = integrateCC(f, a, b, CCEvals);
    if(isEEqual(resultCC.first, resultCC.first + resultCC.second, eRelAbs))
        return resultCC;
    pair<double, double> resultGLK = integrateAdaptiveHeap<
        IntervalGaussLobattoKronrod>(f, a, b, eRelAbs, maxEvals - CCEvals,
        minEvals);
    return resultCC.second < resultGLK.second ? resultCC : resultGLK;
}

pair<double, double> integrateFromData(Vector<pair<double, double> > xyPairs)
{
    assert(xyPairs.getSize() >= 3);
    quickSort(xyPairs.getArray(), 0, xyPairs.getSize() - 1,
        PairFirstComparator<double, double>());
    double result[2] = {0, 0}, last = 0;
    for(int j = 2; j >= 1; --j)
        for(int i = 0; i + j < xyPairs.getSize(); i += j) result[j - 1] +=
            last = (xyPairs[i + j].first - xyPairs[i].first) *
            (xyPairs[i + j].second + xyPairs[i].second)/2;
    if(xyPairs.getSize() % 2 == 0) result[1] += last;
    return make_pair(result[0], abs(result[0] -  result[1]));
}

struct Cheb1DIntegrator
{//ensure getting error estimate
    template<typename FUNCTION> pair<double, double> operator()(
        FUNCTION const& f, double a, double b, int maxEvals)
        {return integrateCC(f, a, b, maxEvals, min(17, (maxEvals - 1)/2 + 1));}
};
template<typename FUNCTION, typename INTEGRATOR1D = Cheb1DIntegrator>
class RecursiveIntegralFunction
{
    FUNCTION f;
    mutable Vector<double> xBound;
    typedef Vector<pair<double, double> > BOX;
    BOX box;
    mutable Vector<pair<double, long long> >* errors;//copy-proof
    int maxEvalsPerDim;//default for low D power of 2 + 1
    pair<double, double> integrateHelper()const
    {
        INTEGRATOR1D i1D;
        int i = xBound.getSize();
        return i1D(*this, box[i].first, box[i].second, maxEvalsPerDim);
    }
public:
    RecursiveIntegralFunction(BOX const& theBox, int theMaxEvalsPerDim = 33,
        FUNCTION const& theF = FUNCTION()): f(theF), box(theBox), errors(0),
        maxEvalsPerDim(theMaxEvalsPerDim){}
    pair<double, double> integrate()const//the main function
    {
        errors = new Vector<pair<double, long long> >(box.getSize());
        pair<double, double> result = integrateHelper();
        double error = 0;
        for(int i = box.getSize() - 1; i >= 0; --i) error = (box[i].second -
            box[i].first) * (error + (*errors)[i].first/(*errors)[i].second);
        result.second += error;
        delete errors;
        return result;
    }
    double operator()(double x)const//called by INTEGRATOR1D
    {
        xBound.append(x);
        bool recurse = xBound.getSize() < box.getSize();
        pair<double, double> result = recurse ? integrateHelper() :
            make_pair(f(xBound), numeric_limits<double>::epsilon());
        if(!recurse) result.second *= max(1.0, abs(result.first));
        xBound.removeLast();
        int i = xBound.getSize();
        (*errors)[i].first += result.second;
        ++(*errors)[i].second;
        return result.first;
    }
};

double boxVolume(Vector<pair<double, double> > const& box)
{
    double result = 1;
    for(int i = 0; i < box.getSize(); ++i)
        result *= box[i].second - box[i].first;
    return result;
}
struct InsideTrue{
    bool operator()(Vector<double> const& dummy)const{return true;}};
template<typename TEST, typename FUNCTION> pair<double, double>
    MonteCarloIntegrate(Vector<pair<double, double> > const& box, int n,
    TEST const& isInside = TEST(), FUNCTION const& f = FUNCTION())
{
    IncrementalStatistics s;
    for(int i = 0; i < n; ++i)
    {
        Vector<double> point(box.getSize());
        for(int j = 0; j < box.getSize(); ++j)
            point[j] = GlobalRNG().uniform(box[j].first, box[j].second);
        if(isInside(point)) s.addValue(f(point));
    }
    double regionVolume = boxVolume(box) * s.n/n;
    return make_pair(regionVolume * s.getMean(),
        regionVolume * s.getStandardErrorSummary().error95());
}

}//end namespace
#endif
