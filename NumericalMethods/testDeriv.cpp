#include <cassert>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <limits>
#include <memory> //for shared ptr
#include "Differentiation.h"
#include "../NumericalOptimization/GlobalNumericalOptimization.h"
#include "../RandomNumberGeneration/Statistics.h"
#include "../Utils/DEBUG.h"
#include "../ExternalMemoryAlgorithms/CSV.h"
#include "NumericalMethodsTestAuto.h"
#include "TestFunctions1D.h"
using namespace std;
using namespace igmdk;

struct DerivatorFD
{
    template<typename FUNCTION>
    double operator()(FUNCTION const& f, double x)const
        {return estimateDerivativeFD(f, x, f(x));}
};
struct DerivatorCD
{
    template<typename FUNCTION>
    double operator()(FUNCTION const& f, double x)const
        {return estimateDerivativeCD(f, x);}
};
template<typename FUNCTION> double estimateDerivativeCD4(FUNCTION const& f,
    double x, double fEFactor = numeric_limits<double>::epsilon())
{
    double h = pow(fEFactor, 1.0/5) * max(1.0, abs(x));
    return (-f(x + 2 * h) + 8 * f(x + h) - 8 * f(x - h) + f(x - 2 * h))/(12 * h);
}
struct DerivatorCD4
{
    template<typename FUNCTION>
    double operator()(FUNCTION const& f, double x)const
        {return estimateDerivativeCD4(f, x);}
};
struct DerivatorCheb
{
    ScaledChebAB pCheb;
    DerivatorCheb(ScaledChebAB const& thePCheb): pCheb(thePCheb) {}
    template<typename FUNCTION> double operator()(FUNCTION const& dummy, double x)const
        {return pCheb.evalDeriv(x);}
};
struct DerivatorChebPiecewise
{
    typedef GenericPiecewiseInterpolation<IntervalCheb::INTERPOLANT> T;
    T pCheb;
    DerivatorChebPiecewise(T const& thePCheb): pCheb(thePCheb) {}
    template<typename FUNCTION> double operator()(FUNCTION const& dummy, double x)const
        {return pCheb.evalDeriv(x);}
};

template<typename FUNCTION> void testDerivHelper(FUNCTION const& f, double a, double b,
    Vector<Vector<string> > & matrix)
{
    double length = b - a, safety = length * 0.05;
    a += safety;
    b -= safety;

    DEBUG("FD");
    matrix.lastItem().append("FD");
    testDerivator(f, DerivatorFD(), a, b, matrix);
    DEBUG("CD");
    matrix.lastItem().append("CD");
    testDerivator(f, DerivatorCD(), a, b, matrix);
    DEBUG("FPS");
    matrix.lastItem().append("FPS");
    testDerivator(f, DerivatorCD4(), a, b, matrix);
    DEBUG("Cheb Doubling");
    matrix.lastItem().append("Cheb Doubling");
    DerivatorCheb dc(ScaledChebAB(adaptiveChebEstimate(ScaledFunctionM11<FUNCTION>(a, b, f)), a, b));
    testDerivator(f, dc, a, b, matrix);
    DEBUG("Cheb64 Adaptive Piecewise");
    matrix.lastItem().append("Cheb64 Adaptive Piecewise");
    DerivatorChebPiecewise dcp(interpolateAdaptiveHeap<IntervalCheb>(f, a, b, 64).first);
    testDerivator(f, dcp, a, b, matrix);
}

template<typename FUNCTION, typename DERIVATOR> pair<double, double> testDerivator(
    FUNCTION const& f, DERIVATOR const& de, double a, double b,
    Vector<Vector<string> > & matrix, int n = 1000000)
{
    DEBUG("eval start");
    double maxRandRelErr = numeric_limits<double>::epsilon();
    for(int i = 0; i < n; ++i)
    {
        double x = GlobalRNG().uniform(a, b), answer = f.deriv(x),
            diff = de(f, x) - answer;
        maxRandRelErr = max(maxRandRelErr, abs(diff/max(1.0, answer)));
        if(i == 0)
        {
            DEBUG(TestFunctions1D::evalCount);
            matrix.lastItem().append(to_string(TestFunctions1D::evalCount));
        }
    }
    double relAbsErrorDigits = log10(maxRandRelErr);
    DEBUG(relAbsErrorDigits);
    matrix.lastItem().append(toStringDouble(relAbsErrorDigits));
    TestFunctions1D::evalCount = 0;
}

void testDeriv()
{
    Vector<Vector<string> > matrix;
    Vector<TestFunctions1D::MetaF> fs = TestFunctions1D::getFunctions();
    for(int i = 0; i < fs.getSize(); ++i)
    {
        if(!isfinite(fs[i].deriv((fs[i].getA() + fs[i].getB())/2))) continue;
        string name = fs[i].getName();
        DEBUG(name);
        matrix.append(Vector<string>());
        matrix.lastItem().append(name);
        testDerivHelper(fs[i], fs[i].getA(), fs[i].getB(), matrix);
    }
    int reportNumber = time(0);
    string filename = "reportDeriv" + to_string(reportNumber) + ".csv";
    createCSV(matrix, filename.c_str());
    Vector<string> names;
    names.append("Evals");
    names.append("Error");
    createAugmentedCSVFiles(matrix, names, filename);
}

struct GradTest
{
    double operator()(Vector<double> const& x)const
    {
        return x[0] * x[0] + x[1];
    }
    Vector<double> grad(Vector<double> const& x)const
    {
        Vector<double> result(2);
        result[0] = 2 * x[0];
        result[1] = 1;
        return result;
    }
    double dd(Vector<double> const& x, Vector<double> const& d)
    {
        return dotProduct(grad(x), d);
    }
};
void testGradDD()
{
    GradTest g;
    Vector<double> x(2, 5), d(2, 1);
    DEBUG(normInf(g.grad(x) - estimateGradientCD(x, g)));
    DEBUG(abs(g.dd(x, d) - estimateDirectionalDerivativeCD(x, g, d)));
}

int main()
{
    testDeriv();
    return 0;
    testGradDD();
    return 0;
    testELessAuto();
    return 0;

    return 0;
}
