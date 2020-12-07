#ifndef IGMDK_ALL_ROOTS_FINDER_H
#define IGMDK_ALL_ROOTS_FINDER_H
#include <cmath>
#include "Matrix.h"
#include "Interpolation.h"
#include "EquationSolving.h"
namespace igmdk{

Vector<complex<double> > findAllRoots(Vector<double> const& lowerCoefs)
{
    int n = lowerCoefs.getSize();
    Matrix<double> companion(n, n);
    for(int r = 0; r < n; ++r)
    {
        if(r > 0) companion(r, r - 1) = 1;
        companion(r, n - 1) = -lowerCoefs[r];
    }
    return QREigenHessenberg(companion);
}

template<typename FUNCTION> Vector<double> findAllRealRootsCheb(
    FUNCTION const& f, double a, double b, int maxDegree = 32,
    double duplicateXEps = highPrecEps)
{
    Vector<pair<pair<double, double>, ScaledChebAB> > pieces =
        interpolateAdaptiveHeap<IntervalCheb>(f, a, b, maxDegree).first.
        getPieces();
    PiecewiseData<EMPTY> resultFilter(duplicateXEps);
    for(int i = 0; i < pieces.getSize(); ++i)
    {
        Vector<double> rootsI = pieces[i].second.findAllRealRoots();
        for(int j = 0; j < rootsI.getSize(); ++j)
        {//range and finiteness check
            double polishedRoot =
                solveSecant(f, rootsI[j], duplicateXEps).first;
            if(isfinite(polishedRoot) && a <= polishedRoot &&
                polishedRoot <= b) resultFilter.insert(polishedRoot, EMPTY());
        }
    }
    Vector<pair<double, EMPTY> > tempResult = resultFilter.getPieces();
    Vector<double> result(tempResult.getSize());
    for(int i = 0; i < tempResult.getSize(); ++i)
        result[i] = tempResult[i].first;
    return result;
}

}//end namespace
#endif
