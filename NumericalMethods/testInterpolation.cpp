#include <cassert>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <limits>
#include <memory> //for shared ptr
#include "Interpolation.h"
#include "CubicSpline.h"
#include "../NumericalOptimization/GlobalNumericalOptimization.h"
#include "../RandomNumberGeneration/Statistics.h"
#include "../Utils/DEBUG.h"
#include "../ExternalMemoryAlgorithms/CSV.h"
#include "NumericalMethodsTestAuto.h"
#include "TestFunctions1D.h"
using namespace std;
using namespace igmdk;

template<typename FUNCTION, typename INTERPOLANT> pair<double, double> testInterpolant(
    FUNCTION const& f, INTERPOLANT const& in, double a, double b,
    Vector<Vector<string> > & matrix, int n = 1000000)
{
    DEBUG("eval start");
    DEBUG(TestFunctions1D::evalCount);
    matrix.lastItem().append(to_string(TestFunctions1D::evalCount));
    double maxRandRelErr = numeric_limits<double>::epsilon();
    for(int i = 0; i < n; ++i)
    {
        double x = GlobalRNG().uniform(a, b), answer = f(x),
            diff = in(x) - answer;
        maxRandRelErr = max(maxRandRelErr, abs(diff/max(1.0, answer)));
    }
    double relAbsErrorDigits = log10(maxRandRelErr);
    DEBUG(relAbsErrorDigits);
    matrix.lastItem().append(toStringDouble(relAbsErrorDigits));
    TestFunctions1D::evalCount = 0;
}


template<typename FUNCTION>
void testInterpolGivenPointsHelper(FUNCTION const& f, Vector<pair<double, double> > const& xy,
    Vector<Vector<string> > & matrix)
{
    assert(xy.getSize() >= 2);
    double a = xy[0].first, b = xy.lastItem().first;
    DEBUG("DynLin");
    DynamicLinearInterpolation dl(xy);
    matrix.lastItem().append("DynLin");
    testInterpolant(f, dl, a, b, matrix);
    DEBUG("NAK Cub");
    NotAKnotCubicSplineInterpolation no(xy);
    matrix.lastItem().append("NAK Cub");
    testInterpolant(f, no, a, b, matrix);
}

template<typename FUNCTION> void testInterpolRandomPointsHelper(FUNCTION const& f, double a, double b,
    Vector<Vector<string> > & matrix, int n = 1000)
{
    DEBUG("Random");
    Vector<pair<double, double> > xy;
    xy.append(make_pair(a, f(a)));
    xy.append(make_pair(b, f(b)));
    n -= 2;
    for(int i = 0; i < n; ++i)
    {
        double x = GlobalRNG().uniform(a, b);
        xy.append(make_pair(x, f(x)));
    }
    quickSort(xy.getArray(), 0, xy.getSize() - 1, PairFirstComparator<double, double>());
    testInterpolGivenPointsHelper(f, xy, matrix);
}

template<typename FUNCTION> void testInterpolChosenPointsHelper(FUNCTION const& f, double a, double b,
    Vector<Vector<string> > & matrix)
{
    DEBUG("Dynamic");
    DEBUG("Cheb Adapt");
    ScaledChebAB cfa(adaptiveChebEstimate(ScaledFunctionM11<FUNCTION>(a, b, f)), a, b);
    matrix.lastItem().append("Cheb Adapt");
    testInterpolant(f, cfa, a, b, matrix);
    DEBUG("Cheb64 Adaptive");
    GenericPiecewiseInterpolation<IntervalCheb::INTERPOLANT> pCheb = interpolateAdaptiveHeap<IntervalCheb>(f, a, b, 64).first;
    matrix.lastItem().append("Cheb64 Adaptive");
    testInterpolant(f, pCheb, a, b, matrix);
}

void testInterpolChosenPoints()
{
    Vector<Vector<string> > matrix;
    bool testDynamic = true, testRandom = true;
    Vector<TestFunctions1D::MetaF> fs = TestFunctions1D::getFunctions();
    for(int i = 0; i < fs.getSize(); ++i)
    {
        string name = fs[i].getName();
        DEBUG(name);
        matrix.append(Vector<string>());
        matrix.lastItem().append(name);
        if(testDynamic) testInterpolChosenPointsHelper(fs[i], fs[i].getA(), fs[i].getB(), matrix);
        if(testRandom) testInterpolRandomPointsHelper(fs[i], fs[i].getA(), fs[i].getB(), matrix);
    }

    int reportNumber = time(0);
    string filename = "reportInterp" + to_string(reportNumber) + ".csv";
    DEBUG(matrix.getSize());
    createCSV(matrix, filename.c_str());
    Vector<string> names;
    names.append("Evals");
    names.append("Error");
    createAugmentedCSVFiles(matrix, names, filename);
}

int main()
{
    testInterpolChosenPoints();
    return 0;
    testAllAutoNumericalMethods();
    return 0;

    return 0;
}
