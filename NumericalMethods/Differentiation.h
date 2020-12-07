#ifndef IGMDK_DIFFERENTIATION_H
#define IGMDK_DIFFERENTIATION_H
#include <cmath>
#include <complex>
#include <iomanip>
#include "../Utils/Bits.h"
#include "../Utils/Vector.h"
#include "../Utils/Utils.h"
#include "../Sorting/Sort.h"
#include "../Utils/Queue.h"
#include "../Heaps/Heap.h"
#include "../RandomNumberGeneration/Random.h"
#include "../RandomTreap/Treap.h"
#include "../ComputationalGeometry/Point.h"
//#include "../Optimization/Metaheuristics.h"
#include "Matrix.h"

namespace igmdk{

template<typename FUNCTION> double estimateDerivativeFD(FUNCTION const& f,
    double x, double fx, double fEFactor = numeric_limits<double>::epsilon())
{
    double h = sqrt(fEFactor) * max(1.0, abs(x));
    return (f(x + h) - fx)/h;
}
template<typename FUNCTION> double estimateDerivativeCD(FUNCTION const& f,
    double x, double fEFactor = numeric_limits<double>::epsilon())
{
    double h = pow(fEFactor, 1.0/3) * max(1.0, abs(x));
    return (f(x + h) - f(x - h))/(2 * h);
}
//PRESENT?
template<typename FUNCTION> struct DerivFunctor
{
    FUNCTION f;
    DerivFunctor(FUNCTION const& theF): f(theF) {}
    double operator()(double p)const{return estimateDerivativeCD(f, p);}
    int fEvals()const{return 2;}
};

double findDirectionScale(Vector<double> const& x, Vector<double> u)
{
    for(int i = 0; i < x.getSize(); ++i) u[i] *= max(1.0, abs(x[i]));
    return norm(u);
}
template<typename FUNCTION> class ScaledDirectionFunction
{
    FUNCTION f;
    Vector<double> x, d;
    double scale;
public:
    ScaledDirectionFunction(FUNCTION const& theF, Vector<double> const& theX,
        Vector<double> const& theD): f(theF), x(theX), d(theD),
        scale(findDirectionScale(x, d * (1/norm(d)))){assert(norm(d) > 0);}
    double getS0()const{return scale;}//returns x[i] if x is axis vector
    double operator()(double s)const{return f(x + d * (s - getS0()));}
};
template<typename FUNCTION> Vector<double> estimateGradientCD(
    Vector<double> const& x, FUNCTION const& f,
    double fEFactor = numeric_limits<double>::epsilon())
{
    int D = x.getSize();
    Vector<double> result(D), d(D);
    for(int i = 0; i < D; ++i)
    {
        d[i] = 1;
        ScaledDirectionFunction<FUNCTION> df(f, x, d);
        result[i] = estimateDerivativeCD(df, df.getS0(), fEFactor);
        d[i] = 0;
    }
    return result;
}
template<typename FUNCTION> struct GradientFunctor
{
    FUNCTION f;
    GradientFunctor(FUNCTION const& theF): f(theF) {}
    Vector<double> operator()(Vector<double> const& p)const
        {return estimateGradientCD(p, f);}
    int fEvals(int D)const{return 2 * D;}
};

template<typename FUNCTION> double estimateDirectionalDerivativeCD(
    Vector<double> const& x, FUNCTION const& f, Vector<double> const& d,
    double fEFactor = numeric_limits<double>::epsilon())
{//estimates grad * d if d not unit
    ScaledDirectionFunction<FUNCTION> df(f, x, d);
    return estimateDerivativeCD(df, df.getS0(), fEFactor);
}

template<typename FUNCTION> struct DirectionalDerivativeFunctor
{
    FUNCTION f;
    DirectionalDerivativeFunctor(FUNCTION const& theF): f(theF) {}
    double operator()(Vector<double> const& x, Vector<double> const& d)
        const{return estimateDirectionalDerivativeCD(x, f, d);}
    int fEvals()const{return 2;}
};

template<typename FUNCTION>  Matrix<double> estimateJacobianCD(FUNCTION const&
    f, Vector<double> x, double fEFactor = numeric_limits<double>::epsilon())
{
    int n = x.getSize();
    Matrix<double> J(n, n);
    Vector<double> dx(n);
    double temp = pow(fEFactor, 1.0/3);
    for(int c = 0; c < n; ++c)
    {
        double xc = x[c], h = max(1.0, abs(xc)) * temp;
        x[c] += h;
        Vector<double> df = f(x);
        x[c] = xc - h;
        df -= f(x);
        x[c] = xc;
        for(int r = 0; r < n; ++r) J(r, c) = df[r]/(2 * h);
    }
    return J;
}

template<typename FUNCTION> double estimate2ndDerivativeCD(FUNCTION const& f,
    double x, double fx = numeric_limits<double>::quiet_NaN(),
    double fEFactor = numeric_limits<double>::epsilon())
{
    if(!isfinite(fx)) fx = f(x);
    double h = pow(fEFactor, 1.0/4) * max(1.0, abs(x));
    return (f(x + h) - 2 * fx + f(x - h))/(h * h);
}

//if gradient is estimated by finite diff, use
//fEFactor = pow(numeric_limits<double>::epsilon(), 2.0/3)
template<typename GRADIENT> Matrix<double> estimateHessianFromGradientCD(
    Vector<double> x, GRADIENT const& g,
    double fEFactor = numeric_limits<double>::epsilon())
{
    Matrix<double> HT = estimateJacobianCD(g, x, fEFactor);
    return 0.5 * (HT + HT.transpose());//ensure symmetry
}

}//end namespace
#endif
