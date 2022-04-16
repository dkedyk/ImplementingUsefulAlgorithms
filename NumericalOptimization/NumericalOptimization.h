#ifndef IGMDK_NUMERICAL_OPTIMIZATION_H
#define IGMDK_NUMERICAL_OPTIMIZATION_H
#include <cmath>
#include "../Utils/Vector.h"
#include "GoldenSection.h"
#include "NelderMead.h"
#include "CoordinateDescent.h"
#include "LBFGS.h"
namespace igmdk{

double RMRate(int i){return 1/pow(i + 1, 0.501);}

template<typename FUNCTION> pair<Vector<double>, double> hybridLocalMinimize(
    Vector<double> const& x0, FUNCTION const& f, int maxEvals = 1000000,
    double yPrecision = highPrecEps)
{
    GradientFunctor<FUNCTION> g(f);
    DirectionalDerivativeFunctor<FUNCTION> dd(f);
    int D = x0.getSize(), LBFGSevals = D < 200 ? maxEvals/2 : maxEvals;
    pair<Vector<double>, double> result = LBFGSMinimize(x0, f, g, dd,
       LBFGSevals, yPrecision);
    if(D > 1 && D < 200)
    {
        int nRestarts = 30;
        NelderMead<FUNCTION> nm(x0.getSize(), f);
        result = nm.restartedMinimize(result.first,
            (maxEvals - LBFGSevals)/nRestarts, highPrecEps, nRestarts);
    }
    return result;
}

}//end namespace igmdk
#endif
