#ifndef IGMDK_ODE_SOLVING_H
#define IGMDK_ODE_SOLVING_H
#include <cmath>
#include <iomanip>
#include "../Utils/Bits.h"
#include "../Utils/Vector.h"
#include "../Utils/Utils.h"
#include "../Sorting/Sort.h"
#include "../Utils/Queue.h"
#include "../Heaps/Heap.h"
//#include "../Optimization/Metaheuristics.h"
#include "NumericalCommon.h"
#include "EquationSolving.h"
namespace igmdk{

template<typename TWO_VAR_FUNCTION>
double RungeKutta4Step(TWO_VAR_FUNCTION const& f, double x, double y,
    double h, double f0 = numeric_limits<double>::quiet_NaN())
{
    if(isnan(f0)) f0 = f(x, y);
    double k1 = h * f0, k2 = h * f(x + h/2, y + k1/2),
        k3 = h * f(x + h/2, y + k2/2), k4 = h * f(x + h, y + k3);
    return y + (k1 + 2 * k2 + 2 * k3 + k4)/6;
}

template<typename TWO_VAR_FUNCTION> pair<pair<double, double>, double>
    RungeKuttaDormandPrinceStep(TWO_VAR_FUNCTION const& f, double x, double y,
    double h, double f0)
{
    double k1 = h * f0,
        k2 = h * f(x + h/5, y + k1/5),
        k3 = h * f(x + h * 3/10, y + k1 * 3/40 + k2 * 9/40),
        k4 = h * f(x + h * 4/5, y + k1 * 44/45 + k2 * -56/15 + k3 * 32/9),
        k5 = h * f(x + h * 8/9, y + k1 * 19372/6561 + k2 * -25360/2187 +
            k3 * 64448/6561 + k4 * -212/729),
        k6 = h * f(x + h, y + k1 * 9017/3168 + k2 * -355/33 + k3 * 46732/5247 +
            k4 * 49/176 + k5 * -5103/18656),
        yNew = y + k1 * 35/384 + k3 * 500/1113 + k4 * 125/192 +
            k5 * -2187/6784 + k6 * 11/84, f1 = f(x + h, yNew),
        k7 = h * f1;
    return make_pair(make_pair(yNew, f1), y + k1 * 5179/57600 + k3 * 7571/16695
        + k4 * 393/640 + k5 * -92097/339200 + k6 * 187/2100 + k7/40);
}
template<typename TWO_VAR_FUNCTION> pair<double, double>
    adaptiveRungeKuttaDormandPrice(TWO_VAR_FUNCTION const& f, double x0,
    double xGoal, double y0, double localERelAbs = defaultPrecEps,
    int maxIntervals = 100000, int minIntervals = -1, int upSkip = 5)
{
    if(minIntervals == -1) minIntervals = sqrt(maxIntervals);
    assert(xGoal > x0 && minIntervals > 0 && upSkip > 0);
    double hMax = (xGoal - x0)/minIntervals, hMin = (xGoal - x0)/maxIntervals,
        linearError = 0, h1 = hMax, y = y0, f0 = f(x0, y);
    bool last = false;
    int stepCounter = 0;
    for(double x = x0; !last;)
    {
        if(x + h1 > xGoal)
        {//make last step accurate
            h1 = xGoal - x;
            last = true;
        }
        pair<pair<double, double>, double> yfye =
            RungeKuttaDormandPrinceStep(f, x, y, h1, f0);
        double h2 = h1/2, xFraction = h1/(xGoal - x0);
        if(h2 < hMin || isEEqual(yfye.first.first, yfye.second,
            max(highPrecEps, localERelAbs * sqrt(xFraction))))
        {//accept step
            x += h1;
            y = yfye.first.first;
            f0 = yfye.first.second;//reuse last eval
            linearError += abs(y - yfye.second);
            if(++stepCounter == upSkip && h2 >= hMin)
            {//use larger step after few consecutive accepted steps
                h1 = min(hMax, h1 * 2);
                stepCounter = 0;
            }
        }
        else
        {//use half step
            h1 = h2;
            last = false;
            stepCounter = 0;
        }
    }
    return make_pair(y, linearError);
}

template<typename YX_FUNCTION>
Vector<double> evalYX(double x, Vector<double> const& y, YX_FUNCTION const& f)
{//assume last arg is x
    Vector<double> yAugmented = y;
    yAugmented.append(x);
    Vector<double> fAugmented = f(yAugmented);
    fAugmented.removeLast();
    return fAugmented;
}
template<typename YX_FUNCTION> struct RadauIIA5Function
{
    Vector<double> y;
    double x, h;
    YX_FUNCTION f;
    Vector<double> operator()(Vector<double> fSumDiffs)const
    {
        assert(fSumDiffs.getSize() % 3 == 0);
        int D = fSumDiffs.getSize()/3;
        double s6 = sqrt(6), ci[3] = {(4 - s6)/10, (4 + s6)/10, 1}, A[3][3] =
        {
            {(88 - 7 * s6)/360, (296 - 169 * s6)/1800, (-2 + 3 * s6)/225},
            {(296 + 169 * s6)/1800, (88 + 7 * s6)/360, (-2 - 3 * s6)/225},
            {(16 - s6)/36, (16 + s6)/36, 1.0/9}
        };
        Vector<Vector<double> > ki(3);
        for(int i = 0; i < 3; ++i)
        {
            Vector<double> yi(y);
            for(int j = 0; j < D; ++j) yi[j] += h * fSumDiffs[i * D + j];
            ki[i] = evalYX(x + h * ci[i], yi, f);
        }
        for(int i = 0; i < 3; ++i)
        {
            Vector<double> fSumi(D, 0);
            for(int j = 0; j < 3; ++j) fSumi += ki[j] * A[i][j];
            for(int j = 0; j < D; ++j)//convert fixed point to 0
                fSumDiffs[i * D + j] = fSumi[j] - fSumDiffs[i * D + j];
        }
        return fSumDiffs;
    }
};
struct RadauIIA5StepF
{
    template<typename YX_FUNCTION> Vector<double> operator()(
        YX_FUNCTION const& f, double x, Vector<double> y, double h,
        double solveERelAbs)const
    {
        RadauIIA5Function<YX_FUNCTION> r5f = {y, x, h, f};
        int D = y.getSize();
        Vector<double> fSumDiffs(3 * D, 0);
        fSumDiffs = solveBroyden(r5f, fSumDiffs, max(solveERelAbs,
            numeric_limits<double>::epsilon() * normInf(y)/h)).first;
        for(int j = 0; j < D; ++j) y[j] += h * fSumDiffs[2 * D + j];
        return y;
    }
};

template<typename YX_FUNCTION, typename STEPF> pair<Vector<double>, double>
    adaptiveStepper(YX_FUNCTION const& f, STEPF const& s, double x0,
    double xGoal, Vector<double> y0, double localERelAbs = defaultPrecEps,
    int maxIntervals = 100000, int minIntervals = -1,
    double upFactor = pow(2, 0.2))
{//assume no reuse of f0
    if(minIntervals == -1) minIntervals = sqrt(maxIntervals);
    assert(xGoal > x0 && minIntervals > 0 && upFactor > 1);
    int D = y0.getSize();
    double hMax = (xGoal - x0)/minIntervals, hMin = (xGoal - x0)/maxIntervals,
        linearError = 0, h1 = hMax;
    Vector<double> y = y0,
        y1 = Vector<double>(D, numeric_limits<double>::quiet_NaN()), f0;
    bool last = false;
    for(double x = x0; !last;)
    {
        if(x + h1 > xGoal)
        {//make last step accurate
            h1 = xGoal - x;
            last = true;
        }
        double h2 = h1/2, xFraction = h1/(xGoal - x0),
            tolERelAbs = max(highPrecEps, localERelAbs * sqrt(xFraction)),
            solveERelAbs = tolERelAbs/10;
        if(isnan(normInf(y1))) y1 = s(f, x, y, h1, solveERelAbs);
        Vector<double> y2 = s(f, x, y, h2, solveERelAbs),
            firstY2 = y2;
        y2 = s(f, x + h2, y2, h2, solveERelAbs);
        double normError = normInf(y2 - y1), normY2 = normInf(y2);
        if(h2 < hMin || isEEqual(normY2 + normError, normY2, tolERelAbs))
        {//accept step
            x += h1;
            y = y2;
            linearError += normError;
            y1 = Vector<double>(D, numeric_limits<double>::quiet_NaN());
            if(h2 >= hMin) h1 = min(hMax, h1 * upFactor);//use larger step
        }
        else
        {//use half step
            y1 = firstY2;
            h1 = h2;
            last = false;
        }
    }
    return make_pair(y, linearError);
}

template<typename YX_FUNCTION, typename BOUNDARY_FUNCTION>
struct BoundaryFunctor
{
    YX_FUNCTION const& f;
    BOUNDARY_FUNCTION const& bf;
    double x0, xGoal;
    double operator()(double b)const
    {
        return bf.evaluateGoal(adaptiveStepper(f, RadauIIA5StepF(),
            x0, xGoal, bf.getInitial(b)).first);
    }
};
template<typename YX_FUNCTION, typename BOUNDARY_FUNCTION>
    Vector<Vector<double> >solveBoundaryValue(YX_FUNCTION const& f, double x0,
    double xGoal, Vector<double> const& xPoints, BOUNDARY_FUNCTION const& bf,
    double b0 = 0)
{
    BoundaryFunctor<YX_FUNCTION, BOUNDARY_FUNCTION> fu = {f, bf, x0, xGoal};
    double bFound = solveSecant(fu, b0).first;
    Vector<Vector<double> > result;
    if(isfinite(bFound))
    {
        Vector<double> y0 = bf.getInitial(bFound);
        for(int i = 0; i < xPoints.getSize(); ++i)
        {
            y0 = adaptiveStepper(f, RadauIIA5StepF(), x0, xPoints[i],
                y0).first;
            x0 = xPoints[i];
            result.append(y0);
        }
    }
    return result;
}

}//end namespace
#endif
