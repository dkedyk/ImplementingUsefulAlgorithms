#ifndef IGMDK_METAHEURISTIC_WRAPPERS_H
#define IGMDK_METAHEURISTIC_WRAPPERS_H
#include <cmath>
#include "../Utils/Vector.h"
#include "../RandomNumberGeneration/Random.h"
#include "../Optimization/Metaheuristics.h"
namespace igmdk{

class UnboundedSampler
{
    Vector<double> center;
    double scaleFactor;
public:
    UnboundedSampler(Vector<double> const& theCenter, double theScaleFactor =
        10): center(theCenter), scaleFactor(theScaleFactor)
        {assert(isfinite(norm(theCenter)));}
    Vector<double> operator()()const
    {
        Vector<double> next = center;
        for(int i = 0; i < next.getSize(); ++i) next[i] = (*this)(i);
        return next;
    }
    double operator()(int i)const
    {//ensure finite samples
        double result = numeric_limits<double>::infinity();
        while(!isfinite(result)) result = center[i] + max(1.0, abs(center[i]))/
            10 * scaleFactor * GlobalRNG().Levy() * GlobalRNG().sign();
        return result;
    }
};
class BoxSampler
{
    Vector<pair<double, double> > box;
public:
    BoxSampler(Vector<pair<double, double> > const& theBox): box(theBox){}
    Vector<double> operator()()const
    {
        Vector<double> next(box.getSize());
        for(int i = 0; i < box.getSize(); ++i) next[i] = (*this)(i);
        return next;
    }
    double operator()(int i)const
        {return GlobalRNG().uniform(box[i].first, box[i].second);}
};

template<typename INCREMENTAL_FUNCTION, typename SAMPLER>struct ContinuousSMove
{
    typedef Vector<double> X;
    typedef pair<double, int> MOVE;
    INCREMENTAL_FUNCTION& f;
    SAMPLER const& s;
    ContinuousSMove(INCREMENTAL_FUNCTION& theF, SAMPLER const& theS):
        f(theF), s(theS){}
    double getScore(X const& x)const{return f(f.getXi());}
    pair<MOVE, double> proposeMove(X const& x)const
    {
        int i = GlobalRNG().mod(x.getSize());
        f.setCurrentDimension(i);
        double xi = s(i), yNext = f(xi), y = f(x[i]);
        return make_pair(MOVE(xi, i),//if started at nan allow all moves
            isnan(y) ? 0 : yNext - y);
    }
    void applyMove(X& x, MOVE const& move)const
    {
        x[move.second] = move.first;
        f.setCurrentDimension(move.second);
        f.bind(move.first);
    }
};

template<typename INCREMENTAL_FUNCTION, typename SAMPLER>
    pair<Vector<double>, double> randomCoordinateDescent(
    INCREMENTAL_FUNCTION& f, SAMPLER const& s, int maxEvals = 1000000)
{
    assert(maxEvals > 0);
    Vector<double> x = localSearch(ContinuousSMove<
        INCREMENTAL_FUNCTION, SAMPLER>(f, s), f.getX(), maxEvals);
    return make_pair(x, f(f.getXi()));
}
template<typename FUNCTION, typename SAMPLER> pair<Vector<double>, double>
    RCDGeneral(FUNCTION const &f, SAMPLER const& s,
    Vector<double> x0 = Vector<double>(), int maxEvals = 1000000)
{
    if(x0.getSize() == 0) x0 = s();
    IncrementalWrapper<FUNCTION> iw(f, x0);
    return randomCoordinateDescent(iw, s, maxEvals);
}

template<typename INCREMENTAL_FUNCTION, typename SAMPLER> pair<Vector<double>,
    double> simulatedAnnealingSMinimizeIncremental(INCREMENTAL_FUNCTION&
    f, SAMPLER const& s, Vector<double> x0, int maxEvals = 1000000)
{
    assert(maxEvals > 0 && isfinite(norm(x0)));
    Vector<double> x = selfTunedSimulatedAnnealing(ContinuousSMove<
        INCREMENTAL_FUNCTION, SAMPLER>(f, s), x0, maxEvals);
    return make_pair(x, f(f.getXi()));
}
template<typename INCREMENTAL_FUNCTION, typename SAMPLER> double
    incrementalSABeforeLocalMinimize(INCREMENTAL_FUNCTION &f, SAMPLER const& s,
    int maxEvals = 1000000, double xPrecision = highPrecEps)
{
    simulatedAnnealingSMinimizeIncremental(f, s, f.getX(), maxEvals * 0.9);
    //f remembers evals so far
    return unimodalCoordinateDescent(f, maxEvals, xPrecision);
}

template<typename FUNCTION, typename SAMPLER> class GAContinuousProblem
{
    FUNCTION const& f;
    SAMPLER const& s;
public:
    typedef Vector<double> X;
    GAContinuousProblem(FUNCTION const& theF, SAMPLER const& theS): f(theF),
        s(theS){}
    X generate()const{return s();}
    void crossover(X& x1, X& x2)const
    {//uniform crossover
        assert(x1.getSize() == x2.getSize());
        for(int k = 0; k < x1.getSize(); ++k) if(GlobalRNG().mod(2))
            swap(x1[k], x2[k]);
    }
    X localSearch(X x, int nLocalMoves)const//RBCD for the same sampler
        {return RCDGeneral(f, s, x, nLocalMoves).first;}
    double evaluate(X const& x)const{return f(x);}
};
template<typename FUNCTION, typename SAMPLER> pair<Vector<double>, double>
    geneticLocalSearchContinuous(FUNCTION const& f, SAMPLER const& s,
    int maxEvals = 1000000)
{
    int nLocalMoves = int(pow(maxEvals, 1.0/3)),populationSize = nLocalMoves;
    return geneticLocalSearch(GAContinuousProblem<FUNCTION, SAMPLER>(f, s),
        populationSize, nLocalMoves, maxEvals);
}

}//end namespace igmdk
#endif
