#ifndef IGMDK_MCMC_H
#define IGMDK_MCMC_H
#include "Random.h"
#include "../Utils/Vector.h"
#include <cmath>
namespace igmdk{

template<typename PDF> class GridRWM
{
    PDF f;
    double x, fx, aFrom, aTo;
    int from, to;
    double sampleHelper(double a)
    {
        double xNew = x + GlobalRNG().uniform(-a, a), fxNew = f(xNew);
        if(fx * GlobalRNG().uniform01() <= fxNew)
        {
            x = xNew;
            fx = fxNew;
        }
        return x;
    }
public:
    GridRWM(double x0 = 0, PDF const& theF = PDF(), int from = -10,
        int to = 20): x(x0), f(theF), fx(f(x)), aFrom(pow(2, from)),
        aTo(pow(2, to)) {}
    double sample()
    {
        for(double a = aFrom; a < aTo; a *= 2) sampleHelper(a);
        return x;
    }
};

template<typename PDF> class MultidimGridRWM
{
    PDF f;
    Vector<double> x;
    double fx, aFrom, aTo, factor;
    Vector<double> sampleHelper(double a)
    {
        Vector<double> xNew = x;
        for(int i = 0; i < xNew.getSize(); ++i)
            xNew[i] += GlobalRNG().uniform(-a, a);
        double fxNew = f(xNew);
        if(fx * GlobalRNG().uniform01() <= fxNew)
        {
            x = xNew;
            fx = fxNew;
        }
        return x;
    }
public:
    MultidimGridRWM(Vector<double> const& x0, PDF const& theF = PDF(),
        int from = -10, int to = 20): x(x0), f(theF), fx(f(x)), aFrom(pow(2,
        from)), aTo(pow(2, to)), factor(pow(2, 1.0/x.getSize())) {}
    Vector<double> sample()
    {
        for(double a = aFrom; a < aTo; a *= factor) sampleHelper(a);
        return x;
    }
};

}//end namespace
#endif
