#ifndef IGMDK_GLOBAL_NUMERICAL_OPTIMIZATION_H
#define IGMDK_GLOBAL_NUMERICAL_OPTIMIZATION_H
#include <cmath>
#include <complex>
#include <iomanip>
#include "../Utils/Bits.h"
#include "../Utils/Vector.h"
#include "../Utils/Utils.h"
#include "../Sorting/Sort.h"
#include "../Utils/Queue.h"
#include "../Heaps/Heap.h"
#include "../RandomNumberGeneration/Random.h"
#include "../RandomTreap/Treap.h"
#include "../ComputationalGeometry/Point.h"
#include "../Optimization/Metaheuristics.h"
#include "../NumericalMethods/Matrix.h"
#include "NumericalOptimization.h"
#include "DifferentialEvolution.h"
#include "MetaheuristicWrappers.h"
namespace igmdk{

struct AgnosticStepSampler
{
    Vector<double> operator()(Vector<double> const& x)const
    {//must ensure finite samples in all components
        assert(isfinite(norm(x)));
        int maxTries = 10;//reasonable infinity protection
        while(maxTries--)
        {//correct on first attempt unless near very large numbers
            Vector<double> u = GlobalRNG().randomUnitVector(x.getSize()),
                result =  x + u *
                (findDirectionScale(x, u)/10 * GlobalRNG().Levy());
            if(isfinite(norm(result))) return result;
        }
        return x;//can't avoid infinities
    }
};

struct BoxConstrainedStepSampler
{
    Vector<pair<double, double> > box;
    AgnosticStepSampler s;
    BoxConstrainedStepSampler(Vector<pair<double, double> > const& theBox):
        box(theBox){}
    Vector<double> operator()(Vector<double> const& x)const
        {return boxTrim(s(x), box);}
};

Vector<pair<double, double> > makeAgnosticBox(int D)
{
    double inf = numeric_limits<double>::infinity();
    return Vector<pair<double, double> >(D, make_pair(-inf, inf));
}

template<typename FUNCTION, typename M_SAMPLER> struct ContinuousMMove
{
    typedef Vector<double> X;
    typedef Vector<double> MOVE;
    FUNCTION const& f;
    M_SAMPLER const& s;
    ContinuousMMove(FUNCTION const& theF, M_SAMPLER const& theS): f(theF),
        s(theS){}
    double getScore(X const& x)const{return f(x);}
    pair<MOVE, double> proposeMove(X const& x)const
    {
        Vector<double> xNext = s(x);
        assert(isfinite(norm(xNext)));
        double yNext = f(xNext), currentScore = getScore(x);
        //if started at nan allow all moves
        return make_pair(xNext, isnan(currentScore) ? 0 : yNext - getScore(x));
    }
    void applyMove(X& x, MOVE const& move)const{x = move;}
};
template<typename FUNCTION, typename M_SAMPLER> pair<Vector<double>,
    double> markovianMinimize(FUNCTION const& f, M_SAMPLER const& s,
    Vector<double> x0, int maxEvals = 1000000)
{
    assert(maxEvals > 0 && isfinite(norm(x0)));
    Vector<double> x = localSearch(ContinuousMMove<FUNCTION, M_SAMPLER>(f, s),
        x0, maxEvals);
    return make_pair(x, f(x));
}

template<typename FUNCTION, typename SAMPLER, typename M_SAMPLER>
pair<Vector<double>, double> smallBudgetHybridMinimize(FUNCTION const& f,
    SAMPLER const& s, M_SAMPLER const& ms, Vector<double> x0,
    int maxEvals = 100)
{
    assert(maxEvals > 2);//no point to have fewer
    return RCDGeneral(f, s, markovianMinimize(f, ms, x0,
        maxEvals * 0.5).first, maxEvals * 0.5);
}

template<typename FUNCTION, typename SAMPLER, typename M_SAMPLER>
pair<Vector<double>, double> hybridBeforeLocalMinimize(
    FUNCTION const& f, SAMPLER const& s, M_SAMPLER const& ms,
    Vector<pair<double, double> > const& box, int maxEvals = 1000000,
    double yPrecision = highPrecEps)
{
    assert(maxEvals > 1000);//no point to have fewer
    pair<Vector<double>, double> glsSolution = geneticLocalSearchContinuous(f,
        s, maxEvals * 0.3), deSolution = differentialEvolutionMinimize(f, s,
        box, maxEvals * 0.3);
    return hybridLocalMinimize(markovianMinimize(f, ms, (deSolution.second <
        glsSolution.second ? deSolution : glsSolution).first, maxEvals * 0.3
        ).first, f, maxEvals * 0.1);
}

//presented in numerical algs because of include cycle
template<typename FUNCTION> struct NormFunction
{
    FUNCTION f;
    NormFunction(FUNCTION const& theF): f(theF){}
    double operator()(Vector<double> const& x)const{return norm(f(x));}
};
template<typename FUNCTION> pair<Vector<double>, double> solveByOptimization(
    FUNCTION const& f, Vector<double> const& x0, int maxEvals = 1000000)
{//use scaled grad as error estimate
    NormFunction<FUNCTION> nf(f);
    pair<Vector<double>, double> xy = hybridBeforeLocalMinimize(nf,
        UnboundedSampler(x0), AgnosticStepSampler(),
        makeAgnosticBox(x0.getSize()), maxEvals - 2);
    double errorEstimate = max(normInf(estimateGradientCD(xy.first, nf))/
        max(1.0, abs(xy.second)), defaultPrecEps);
    return make_pair(xy.first, errorEstimate);
}

template<typename FUNCTION> pair<Vector<double>, double>
hybridEquationSolve(FUNCTION const& f, int D,
    double xEps = highPrecEps, int maxEvals = 1000000)
{//first opt, then Broyden local search on opt result to improve precision
    int broydenEvals = max(1000, maxEvals/4),
        optEvals = maxEvals/2 - broydenEvals;
    pair<Vector<double>, double> result = solveByOptimization(f,
        Vector<double>(D), optEvals), result2 =
        solveBroydenHybrid(f, result.first, xEps, broydenEvals);
    //keep opt error estimate if Broyden did nothing
    if(!isfinite(result2.second)) result2.second = result.second;
    //do random Broyden with the other half of evals
    result = solveBroydenLevy(f, D, xEps, maxEvals/2/1000);
    if(normInf(f(result.first)) < normInf(f(result2.first)))
        result = result2;
    return result;
}

}//end namespace igmdk
#endif
