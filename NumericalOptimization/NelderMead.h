#ifndef NELDER_MEAD_H
#define NELDER_MEAD_H
#include <cmath>
#include "../Utils/Vector.h"
#include "../RandomNumberGeneration/Random.h"
#include "../NumericalMethods/NumericalCommon.h"
namespace igmdk{

template<typename FUNCTION> struct NelderMead
{
    FUNCTION f;
    int D;
    Vector<double> vertexSum;//incremental centroid
    typedef pair<Vector<double>, double> P;
    Vector<P> simplex;
    double scale(P& high, double factor, int& maxEvals)
    {
        P result = high;
        //affine combination of the high point and the
        //centroid of the remaining vertices
        //centroid = (vertexSum - high)/D and
        //result = centroid * (1 - factor) + high * factor
        double centroidFactor = (1 - factor)/D;
        result.first = vertexSum * centroidFactor +
            high.first * (factor - centroidFactor);
        result.second = f(result.first);
        --maxEvals;
        if(result.second < high.second)
        {//accept scaling if improving
            vertexSum += result.first - high.first;
            high = result;
        }
        return result.second;
    }
public:
    NelderMead(int theD, FUNCTION const& theFunction = FUNCTION()):
        D(theD), f(theFunction), simplex(D + 1){assert(D > 1);}

    P minimize(Vector<double> const& initialGuess, int maxEvals = 1000000,
        double yPrecision = highPrecEps, double step = 1)
    {//initialize the simplex
        vertexSum = initialGuess;
        for(int i = 0; i < D; ++i) vertexSum[i] = 0;
        for(int i = 0; i <= D; ++i)
        {
            simplex[i].first = initialGuess;
            if(i > 0) simplex[i].first[i - 1] += GlobalRNG().uniform01() *
                step * max(1.0, abs(initialGuess[i - 1]))/10;
            simplex[i].second = f(simplex[i].first);
            --maxEvals;
            vertexSum += simplex[i].first;
        }
        for(;;)
        {//calculate high, low, and nextHigh, which must be all different
            int high = 0, nextHigh = 1, low = 2;
            if(simplex[high].second < simplex[nextHigh].second)
                swap(high, nextHigh);
            if(simplex[nextHigh].second < simplex[low].second)
            {
                swap(low, nextHigh);
                if(simplex[high].second < simplex[nextHigh].second)
                    swap(high, nextHigh);
            }
            for(int i = 3; i <= D; ++i)
            {
                if(simplex[i].second < simplex[low].second) low = i;
                else if(simplex[i].second > simplex[high].second)
                {
                    nextHigh = high;
                    high = i;
                }
                else if(simplex[i].second > simplex[nextHigh].second)
                    nextHigh = i;
            }//check if already converged
            if(maxEvals <= 0 || !isELess(simplex[low].second,
                simplex[high].second, yPrecision)) return simplex[low];
            //try to reflect
            double value = scale(simplex[high], -1, maxEvals);
            //try to double if better than low
            if(value <= simplex[low].second) scale(simplex[high], 2, maxEvals);
            else if(value >= simplex[nextHigh].second)
            {//try reflected/unrefrected halving if accepted/rejected value
                double yHi = simplex[high].second;
                if(scale(simplex[high], 0.5, maxEvals) >= yHi)
                {//contract all to get rid of the high point
                    vertexSum = simplex[low].first;
                    for(int i = 0; i <= D; ++i) if(i != low)
                    {
                        vertexSum += simplex[i].first = (simplex[i].first +
                            simplex[low].first) * 0.5;
                        simplex[i].second = f(simplex[i].first);
                        --maxEvals;
                    }
                }
            }
        }
    }

    P restartedMinimize(Vector<double> const& initialGuess,
        int maxEvals = 100000, double yPrecision = highPrecEps,
        int maxRepeats = 10, double step = 1)
    {
        P result(initialGuess, numeric_limits<double>::infinity());
        while(maxRepeats--)
        {
            double yOld = result.second;
            result = minimize(result.first, maxEvals, yPrecision, step);
            if(!isELess(result.second, yOld, yPrecision)) break;
        }
        return result;
    }
};

}//end namespace igmdk
#endif
