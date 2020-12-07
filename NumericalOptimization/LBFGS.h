#ifndef IGMDK_LBFGS_H
#define IGMDK_LBFGS_H
#include <cmath>
#include "../Utils/Vector.h"
#include "../Utils/Utils.h"
#include "GoldenSection.h"
#include "../NumericalMethods/Differentiation.h"
namespace igmdk{

template<typename FUNCTION, typename DIRECTIONAL_DERIVATIVE> bool
    strongWolfeLineSearchMoreThuente(FUNCTION const& f,
    Vector<double> const& gradient, DIRECTIONAL_DERIVATIVE const& dd,
    Vector<double>& x, double& fx, int& maxEvals, Vector<double> const& dx,
    double yEps)
{
    double dd0 = dotProduct(dx, gradient), sLo = 0, fLo = fx, s = 1,
        sHi = numeric_limits<double>::infinity(), temp = 0.0001 * dd0;
    if(!isfinite(dd0) || dd0 >= 0) return true;
    while(maxEvals > 0 && isfinite(s) &&
        (!isfinite(sHi) || isELess(fx + dd0 * abs(sHi - sLo), fx, yEps)))
    {
        double fNew = f(x + dx * s);
        if(!isfinite(fNew)) break;
        --maxEvals;
        if(fNew - s * temp > fLo - sLo * temp) sHi = s;//case 1
        else
        {
            double ddS = dd(x + dx * s, dx);
            maxEvals -= dd.fEvals();
            if(abs(ddS) <= -0.1 * dd0)//check for early termination
            {
                sLo = s;
                fLo = fNew;
                break;
            }
            if((ddS - temp) * (sLo - s) <= 0) sHi = sLo;//case 3
            //case 2 and 3
            sLo = s;
            fLo = fNew;
        }
        if(isfinite(sHi)) s = (sLo + sHi)/2;//zooming
        else s *= 2;//case 2 init doubling
    }
    if(sLo > 0)
    {//any non-0 guaranteed to have sufficient descent
        x += dx * sLo;
        double fxFirst = fx;
        fx = fLo;
        return !isELess(fx, fxFirst, yEps);//must make good progress
    }
    return true;
}

template<typename FUNCTION1D> struct EvalCountWrapper
{//to keep track of used evaluations, which golden section doesn't
    FUNCTION1D f;
    mutable int evalCount;
    EvalCountWrapper(FUNCTION1D const& theF): f(theF), evalCount(0){}
    double operator()(double x)const
    {
        ++evalCount;
        return f(x);
    }
};
template<typename FUNCTION, typename DIRECTIONAL_DERIVATIVE> bool
    goldenSectionLineSearch(FUNCTION const& f, Vector<double> const& gradient,
    DIRECTIONAL_DERIVATIVE const& dd, Vector<double>& x, double& fx,
    int& maxEvals, Vector<double> const& dx, double yEps, bool useExact = true)
{
    if(!(norm(dx) > 0)) return false;//this way handle NaN also
    double step = 1;//well-scaled for most algorithms
    if(useExact)
    {
        EvalCountWrapper<ScaledDirectionFunction<FUNCTION> > f2(
            ScaledDirectionFunction<FUNCTION>(f, x, dx));
        pair<double, double> result = minimizeGSBracket(f2, f2.f.getS0(), fx,
            false);
        maxEvals -= f2.evalCount;//may be exceeded but ok
        if(isELess(result.second, fx, yEps))
            step = result.first - f2.f.getS0();
    }//ensure strong Wolfe conditions
    return strongWolfeLineSearchMoreThuente(f, gradient, dd, x, fx, maxEvals,
        dx * step, yEps);
}

template<typename FUNCTION, typename GRADIENT, typename DIRECTIONAL_DERIVATIVE>
    pair<Vector<double>, double> LBFGSMinimize(Vector<double> const& x0,
    FUNCTION const& f, GRADIENT const& g, DIRECTIONAL_DERIVATIVE const& dd,
    int maxEvals = 1000000, double yPrecision = highPrecEps,
    int historySize = 8, bool useExact = true)
{
    typedef Vector<double> V;
    Queue<pair<V, V> > history;
    pair<V, double> xy(x0, f(x0));
    V grad = g(xy.first), d;
    int D = xy.first.getSize(), gEvals = g.fEvals(D);
    maxEvals -= 1 + gEvals;
    double lastGoodStepSize = max(1.0, norm(x0))/10;
    while(maxEvals > 0)
    {//backtrack using d to get sufficient descent
        if(history.getSize() == 0) d = grad * (-lastGoodStepSize/norm(grad));
        pair<V, double> xyOld = xy;
        if(goldenSectionLineSearch(f, grad, dd, xy.first, xy.second, maxEvals,
            d, yPrecision, useExact))
        {//failure
            if(history.getSize() > 0)
            {//try to "restart" by purging history one at a time
                history.pop();
                continue;
            }
            else break;//gradient descent step failed
        }
        else lastGoodStepSize = norm(xy.first - xyOld.first);//success
        if((maxEvals -= gEvals) < 1) break;//out of evals
        V newGrad = g(xy.first);
        if(history.getSize() >= historySize) history.pop();
        history.push(make_pair(xy.first - xyOld.first, newGrad - grad));
        //"double recursion" algorithm to update d
        d = grad = newGrad;
        Vector<double> a, p;
        int last = history.getSize() - 1;
        for(int i = last; i >= 0; --i)
        {
            double pi = 1/dotProduct(history[i].first, history[i].second),
                ai = dotProduct(history[i].first, d) * pi;
            d -= history[i].second * ai;
            a.append(ai);
            p.append(pi);
        }//initial Hessian is scaled diagonal
        d *= 1/(dotProduct(history[last].second, history[last].second) *
            p[last]);
        for(int i = 0; i < history.getSize(); ++i)
        {
            double bi = dotProduct(history[i].second, d) * p[last - i];
            d += history[i].first * (a[last - i] - bi);
        }
        d *= -1;
    }
    return xy;
}

}//end namespace igmdk
#endif
