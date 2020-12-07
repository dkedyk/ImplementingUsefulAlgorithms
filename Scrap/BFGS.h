#ifndef IGMDK_BFGS_H
#define IGMDK_BFGS_H
#include <cmath>
#include "../Utils/Vector.h"
#include "../Utils/Utils.h"
#include "../NumericalOptimization/LBFGS.h"
#include "ConjugateGradientOptimization.h"
namespace igmdk{

Matrix<double> invertAndCorrectHessian(Matrix<double> H, bool isInverted,
    double kLimit = 100000)
{
    assert(kLimit > 1);
    int D = H.getRows();
    if(!isESymmetric(H))
    {//may have lost symmetry numerically
        for(int r = 0; r < D; ++r) for(int c = r + 1; c < D; ++c)
            H(r, c) = H(c, r) = (H(r, c) + H(c, r))/2;
        if(!isESymmetric(H))
        {//failed to correct due to inf or NaN
            if(isInverted) return H;//can't do anything else
            else assert(false);//bad input;
        }
    }
    pair<Vector<double>, Matrix<double> > eigs = QREigenSymmetric(H);
    double minEigva = max(1.0, valMax(eigs.first.getArray(), D))/kLimit;
    Matrix<double> diag(D, D);//form diagonal of the inverse
    for(int i = 0; i < D; ++i) diag(i, i) = 1/max(minEigva,
        isInverted ? 1/eigs.first[i] : eigs.first[i]);
    return eigs.second.transpose() * diag * eigs.second;
}
template<typename FUNCTION, typename GRADIENT, typename DIRECTIONAL_DERIVATIVE>
    pair<Vector<double>, double> BFGSMinimizeProper(Vector<double> const&
    initialGuess, FUNCTION const& f, GRADIENT const& g,
    DIRECTIONAL_DERIVATIVE const& dd,
    double gEps = pow(numeric_limits<double>::epsilon(), 2.0/3),
    int maxEvals = 1000000, double yPrecision = highPrecEps,
    bool useExact = true)
{
    typedef Vector<double> V;
    int D = initialGuess.getSize(), gEvals = g.fEvals(D),
        restartDelta = max(D, 100), restartCountdown = restartDelta;
    pair<V, double> xy(initialGuess, f(initialGuess));
    Matrix<double> I = Matrix<double>::identity(D),
        H = invertAndCorrectHessian(
        estimateHessianFromGradientCD(xy.first, g, gEps), false);
    V grad = g(xy.first);
    maxEvals -= 1 + gEvals + 2 * D * gEvals;
    while(maxEvals > 0)
    {//backtrack using d to get sufficient descent
        if(restartCountdown <= 0 && (maxEvals -= 2 * D * gEvals) > 0)
        {
            restartCountdown = restartDelta;
            H = invertAndCorrectHessian(
                estimateHessianFromGradientCD(xy.first, g, gEps), false);
        }
        //steepest descent if not enough evals to restart
        V d = restartCountdown <= 0 ? grad * (-1/norm(grad)) : H * -grad;
        pair<V, double> xyOld = xy;
        if(goldenSectionLineSearch(f, grad, dd, xy.first, xy.second, maxEvals, d,
            yPrecision, useExact))
        {//failed case
            if(restartCountdown <= 0) break;//failed in steepest descent
            else if(restartCountdown == restartDelta);//failed after restart
            {
                d = grad * (-1/norm(grad));//try steepest descent
                if(goldenSectionLineSearch(f, grad, dd, xy.first, xy.second,
                    maxEvals, d, yPrecision, useExact)) break;
            }//force restart
            restartCountdown = 0;
            continue;
        }
        else --restartCountdown;
        if((maxEvals -= gEvals) < 1) break;
        V newGrad = g(xy.first), s = xy.first - xyOld.first,
            y = newGrad - grad;
        double p = 1/dotProduct(y, s);
        H = (I - p * outerProduct(s, y)) * H * (I - p * outerProduct(y, s)) +
            p * outerProduct(s, s);
        grad = newGrad;
    }
    return xy;
}

}//end namespace igmdk
#endif
