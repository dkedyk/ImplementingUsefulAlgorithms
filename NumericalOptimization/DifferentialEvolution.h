#ifndef IGMDK_DIFFERENTIAL_EVOLUTION_H
#define IGMDK_DIFFERENTIAL_EVOLUTION_H
#include <cmath>
#include "../Utils/Vector.h"
#include "../RandomNumberGeneration/Random.h"
namespace igmdk{

Vector<double> boxTrim(Vector<double> x,
    Vector<pair<double, double> > const& box)
{
    for(int i = 0; i < x.getSize(); ++i)
    {
        if(x[i] < box[i].first) x[i] = box[i].first;
        else if(isnan(x[i]) || x[i] > box[i].second) x[i] = box[i].second;
    }
    return x;
}

template<typename FUNCTION, typename SAMPLER> pair<Vector<double>, double>
    differentialEvolutionMinimize(FUNCTION const& f, SAMPLER const& s,
    Vector<pair<double, double> > const& box, int maxEvals = 1000000)
{
    assert(maxEvals > 0);
    int n = pow(maxEvals, 1.0/3);
    Vector<pair<Vector<double>, double> > population(n);
    for(int i = 0; i < n; ++i)
    {
        population[i].first = s();
        population[i].second = f(population[i].first);
    }
    maxEvals -= n;
    while(maxEvals > 0)
    {
        for(int i = 0; i < n && maxEvals-- > 0; ++i)
        {//mutate new point
            Vector<int> jkl = GlobalRNG().randomCombination(3, n);
            Vector<double> xiNew = population[i].first, xiMutated =
                boxTrim(population[jkl[0]].first + (population[jkl[1]].first -
                population[jkl[2]].first) * 0.9, box);
            //crossover with mutated point
            int D = xiNew.getSize(), randK = GlobalRNG().mod(D);
            for(int k = 0; k < D; ++k) if(GlobalRNG().mod(2) || k == randK)
                xiNew[k] = xiMutated[k];
            if(!isfinite(norm(xiNew)))
            {//enforce finite samples
                ++maxEvals;
                continue;
            }
            //select best of original and mutated
            double yiNew = f(xiNew);
            if(yiNew < population[i].second)
            {
                population[i].first = xiNew;
                population[i].second = yiNew;
            }
        }
    }
    pair<Vector<double>, double>& best = population[0];
    for(int i = 1; i < n; ++i)
        if(best.second < population[i].second) best = population[i];
    return best;
}

}//end namespace igmdk
#endif
