#ifndef IGMDK_CONJUGATE_GRADIENT_OPTIMIZATION_H
#define IGMDK_CONJUGATE_GRADIENT_OPTIMIZATION_H
#include <cmath>
#include "../Utils/Vector.h"
#include "../Utils/Utils.h"
#include "../NumericalOptimization/LBFGS.h"
namespace igmdk{

double steepestStepScale(Vector<double> const& x, double y,
Vector<double> const& grad)
{
    return max(1/max(1.0, abs(y)), max(1.0, norm(x)) *
        sqrt(numeric_limits<double>::epsilon())/norm(grad));
}
template<typename FUNCTION, typename GRADIENT, typename DIRECTIONAL_DERIVATIVE>
    pair<Vector<double>, double> conjugateGradient(Vector<double> const& x0,
    FUNCTION const& f, GRADIENT const& g, DIRECTIONAL_DERIVATIVE const& dd,
    int maxEvals = 1000000, double yPrecision = highPrecEps,
    string const& formula = "PRP+", bool useExact = true)
{
    pair<Vector<double>, double> xy(x0, f(x0));
    Vector<double> grad = g(xy.first), d;
    int D = xy.first.getSize(), gEvals = g.fEvals(D),
        restartDelta = 100, restartCountdown = 0;
    maxEvals -= gEvals;
    double a;
    while(maxEvals > 0)
    {
        if(restartCountdown <= 0)
        {
            restartCountdown = restartDelta;
            d = -grad;
            a = steepestStepScale(xy.first, xy.second, grad);
        }
        Vector<double> xOld = xy.first;
        if(goldenSectionLineSearch(f, grad, dd, xy.first, xy.second, maxEvals,
            d * a, yPrecision, useExact))
        {//failed case
            if(restartCountdown == restartDelta) break;//failed after restart
            else
            {//force restart
                restartCountdown = 0;
                continue;
            }
        }
        else --restartCountdown;
        if((maxEvals -= gEvals) < 1) break;
        Vector<double> gradNew = g(xy.first);
        double b = 0;
        if(formula == "PRP+") b = max(0.0, dotProduct(gradNew,
            (gradNew - grad))/dotProduct(grad, grad));
        else if(formula == "HZ")
        {
            Vector<double> y = gradNew - grad;
            double temp = dotProduct(y, d), b = dotProduct(y - d *
                (dotProduct(y, y) * 2/temp), gradNew)/temp;
        }
        else
        {
            DEBUG("bad formula specification");
            assert(false);
        }
        double temp = dotProduct(grad, (xy.first - xOld));
        d = -gradNew + d * b;
        grad = gradNew;
        a = dotProduct(grad, d)/temp;
    }
    return xy;
}

}//end namespace igmdk
#endif
