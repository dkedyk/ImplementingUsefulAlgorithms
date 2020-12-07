#ifndef IGMDK_COORDINATE_DESCENT_H
#define IGMDK_COORDINATE_DESCENT_H
#include <cmath>
#include "../Utils/Vector.h"
#include "../RandomNumberGeneration/Random.h"
#include "GoldenSection.h"
namespace igmdk{

template<typename FUNCTION> struct IncrementalWrapper
{
    FUNCTION f;
    mutable Vector<double> xBound;
    int i;
    mutable int evalCount;
public:
    IncrementalWrapper(FUNCTION const& theF, Vector<double> const& x0):
        f(theF), xBound(x0), i(0), evalCount(0) {}
    void setCurrentDimension(double theI)
    {
        assert(theI >= 0 && theI < xBound.getSize());
        i = theI;
    }
    int getEvalCount()const{return evalCount;}
    int getSize()const{return xBound.getSize();}
    Vector<double> const& getX()const{return xBound;}
    double getXi()const{return xBound[i];}
    double operator()(double xi)const
    {
        double oldXi = xBound[i];
        xBound[i] = xi;
        double result = f(xBound);
        ++evalCount;
        xBound[i] = oldXi;
        return result;
    }
    void bind(double xi)const{xBound[i] = xi;}
};

template<typename INCREMENTAL_FUNCTION> double unimodalCoordinateDescent(
    INCREMENTAL_FUNCTION &f, int maxEvals = 1000000,
    double xPrecision = highPrecEps)
{
    int D = f.getSize();
    Vector<int> order(D);
    for(int i = 0; i < D; ++i) order[i] = i;
    double y = f(f.getXi()), relStep = 0.1;
    while(f.getEvalCount() < maxEvals)//may be exceeded but ok
    {//use random order full cycles
        GlobalRNG().randomPermutation(order.getArray(), D);
        double yPrev = y, maxRelXStep = 0;
        for(int i = 0; i < D && f.getEvalCount() < maxEvals; ++i)
        {
            int j = order[i];
            f.setCurrentDimension(j);
            pair<double, double> resultJ = minimizeGSBracket(f, f.getXi(),
                y, true, relStep, xPrecision);
            maxRelXStep = max(maxRelXStep, abs(resultJ.first - f.getXi())/
                max(1.0, abs(f.getXi())));
            f.bind(resultJ.first);
            y = resultJ.second;
        }//done if no improvement in x or y
        if(maxRelXStep < xPrecision || !isELess(y, yPrev)) break;
        relStep = min(0.1, maxRelXStep);//take smaller first steps as converge
    }
    return y;
}
template<typename FUNCTION> pair<Vector<double>, double>
    unimodalCoordinateDescentGeneral(FUNCTION const& f, Vector<double> const&
    x0, int maxEvals = 1000000, double xPrecision = highPrecEps)
{
    IncrementalWrapper<FUNCTION> iw(f, x0);
    double y = unimodalCoordinateDescent(iw, maxEvals, xPrecision);
    return make_pair(iw.xBound, y);
}

}//end namespace igmdk
#endif
