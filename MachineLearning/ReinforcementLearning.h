#ifndef IGMDK_REINFORCEMENT_LEARNING_H
#define IGMDK_REINFORCEMENT_LEARNING_H

#include "../Utils/Vector.h"
#include <cmath>
namespace igmdk{

double UCB1(double averageValue, int nTries, int totalTries)
    {return averageValue + sqrt(2 * log(totalTries)/nTries);}

template<typename PROBLEM> void TDLearning(PROBLEM& p)
{
    while(p.hasMoreEpisodes())
    {
        double valueCurrent = p.startEpisode();
        while(!p.isInFinalState())
        {
            double valueNext = p.pickNextState();
            p.updateCurrentStateValue(p.learningRate() * (p.reward() +
                p.discountRate() * valueNext - valueCurrent));
            p.goToNextState();
            valueCurrent = valueNext;
        }
        p.updateCurrentStateValue(p.learningRate() *
            (p.reward() - valueCurrent));
    }
}

struct DiscreteValueFunction
{
    Vector<pair<double, int> > values;
    double learningRate(int state){return 1.0/values[state].second;}
    void updateValue(int state, double delta)
    {
        ++values[state].second;
        values[state].first += delta;
    }
    DiscreteValueFunction(int n): values(n, make_pair(0.0, 1)){}
};

struct LinearCombinationValueFunction
{
    Vector<double> weights;
    int n;
    double learningRate(){return 1.0/n;}
    void updateWeights(Vector<double> const& stateFeatures, double delta)
    {//set one of the state features to 1 to have a bias weight
        assert(stateFeatures.getSize() == weights.getSize());
        for(int i = 0; i < weights.getSize(); ++i)
            weights[i] += delta * stateFeatures[i];
        ++n;
    }
    LinearCombinationValueFunction(int theN): weights(theN, 0), n(1) {}
};

}//end namespace
#endif

