#ifndef IGMDK_DISCRETE_NUMERICAL_OPTIMIZATION_H
#define IGMDK_DISCRETE_NUMERICAL_OPTIMIZATION_H
#include <cmath>
#include "../Utils/Vector.h"
#include "../RandomNumberGeneration/Random.h"
namespace igmdk{

template<typename FUNCTION> Vector<double> gridMinimize(
    Vector<Vector<double> > const& sets, FUNCTION const& f = FUNCTION())
{
    assert(sets.getSize() > 0);
    long long total = 1;
    for(int i = 0; i < sets.getSize(); ++i)
    {
        assert(sets[i].getSize() > 0);
        total *= sets[i].getSize();
    }
    Vector<double> best;
    double bestScore;
    for(long long i = 0; i < total; ++i)
    {//unrank and eval
        long long rank = i;
        Vector<double> next;
        for(int j = 0; j < sets.getSize(); ++j)
        {
            next.append(sets[j][rank % sets[j].getSize()]);
            rank /= sets[j].getSize();
        }
        double score = f(next);
        if(best.getSize() == 0 || score < bestScore)
        {
            bestScore = score;
            best = next;
        }
    }
    return best;
}

template<typename FUNCTION> Vector<double> randomDiscreteMinimize(
    Vector<Vector<double> > const& sets, FUNCTION const& f = FUNCTION(),
    int evals = 20)
{
    assert(sets.getSize() > 0);
    Vector<double> best;
    double bestScore;
    for(long long i = 0; i < evals; ++i)
    {
        Vector<double> next;
        for(int j = 0; j < sets.getSize(); ++j)
            next.append(sets[j][GlobalRNG().mod(sets[j].getSize())]);
        double score = f(next);
        if(best.getSize() == 0 || score < bestScore)
        {
            bestScore = score;
            best = next;
        }
    }
    return best;
}

template<typename FUNCTION> Vector<double> compassDiscreteMinimize(
    Vector<Vector<double> > const& sets, FUNCTION const& f = FUNCTION(),
    int remainingEvals = 100)
{//use median in each set as initial solution
    Vector<int> current;
    for(int i = 0; i < sets.getSize(); ++i)
    {
        assert(sets[i].getSize() > 0);
        current.append(sets[i].getSize()/2);
    }
    return compassDiscreteMinimizeHelper(sets, current, f,
        remainingEvals).first;
}
//assumes set values are in sorted (or reverse sorted) order!
template<typename FUNCTION> pair<Vector<double>, pair<double, int> >
    compassDiscreteMinimizeHelper(Vector<Vector<double> > const& sets,
    Vector<int> current, FUNCTION const& f = FUNCTION(),
    int remainingEvals = 100)
{//start with medians
    Vector<double> best;
    for(int i = 0; i < sets.getSize(); ++i)
    {
        assert(0 <= current[i] && current[i] < sets[i].getSize());
        best.append(sets[i][current[i]]);
    }
    double bestScore = f(best);
    Vector<int> preferredSign(sets.getSize(), 1);
    for(bool isOpt = false, changedSign = false; !isOpt;)
    {
        isOpt = true;
        for(int i = 0; i < sets.getSize(); ++i)
        {
            int next = current[i] + preferredSign[i];
            if(0 <= next && next < sets[i].getSize())
            {
                if(remainingEvals-- < 1)
                    return make_pair(best, make_pair(bestScore, 0));
                best[i] = sets[i][next];
                double score = f(best);
                if(score < bestScore)
                {
                    current[i] = next;
                    bestScore = score;
                    isOpt = false;
                }
                else
                {
                    best[i] = sets[i][current[i]];
                    preferredSign[i] *= -1;
                }
            }
            else preferredSign[i] *= -1;
        }
        if(isOpt){if(!changedSign) isOpt = false;}
        else changedSign = false;
    }
    return make_pair(best, make_pair(bestScore, remainingEvals));
}

}//end namespace igmdk
#endif
