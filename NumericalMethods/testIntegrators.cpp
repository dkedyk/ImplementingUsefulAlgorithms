#include <cassert>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <limits>
#include <memory> //for shared ptr
#include "Integration.h"
#include "../NumericalOptimization/GlobalNumericalOptimization.h"
#include "../RandomNumberGeneration/Sobol.h"
#include "../Utils/DEBUG.h"
#include "../ExternalMemoryAlgorithms/CSV.h"
#include "NumericalMethodsTestAuto.h"
#include "TestFunctions1D.h"
using namespace std;
using namespace igmdk;

template<typename FUNCTION> struct MCIFunctionHelper
{
    FUNCTION const& f;
    MCIFunctionHelper(FUNCTION const& theF): f(theF){}
    double operator()(Vector<double>const& x)const
    {
        return f(x[0]);
    }
};

template<typename FUNCTION> pair<double, double> integrateCCDoubling(
    FUNCTION const& f, double a, double b, double maxRelAbsError = highPrecEps,
    int maxEvals = 5000, int minEvals = 17)
{//use doubling for error estimate and convergence criteria
    int n = minEvals - 1;
    assert(minEvals <= maxEvals && isPowerOfTwo(n));
    Vector<double> fx(n + 1);
    ScaledFunctionM11<FUNCTION> fM11(a, b, f);
    for(int i = 0; i <= n; ++i) fx[i] = fM11(cos(PI() * i/n));
    ScaledChebAB che(fx, a, b);
    double result = che.integrate().first,//want to double once
        oldResult = numeric_limits<double>::quiet_NaN();
    while(maxEvals >= fx.getSize() + n && (isnan(oldResult) ||
        !isEEqual(result, oldResult, maxRelAbsError)))
    {
        fx = reuseChebEvalPoints(fM11, fx);
        ScaledChebAB che2(fx, a, b);
        n *= 2;
        oldResult = result;
        result = che2.integrate().first;
    }
    return make_pair(result, abs(result - oldResult));
}

class IntervalSimpson
{
    double left, dx, fLeft, fLM, fMiddle, fMR, fRight;
    static double integrateHelper(double f0, double f1, double f2, double dx)
        {return dx * (f0 + 4 * f1 + f2)/3;}
    template<typename FUNCTION> void initHelper(FUNCTION const& f, double a,
        double b, double fa, double fm, double fb)
    {
        left = a;
        dx = (b - a)/4;
        fLeft = fa;
        fLM = f(a + dx);
        fMiddle = fm;
        fMR = f(a + 3 * dx);
        fRight = fb;
    }
    template<typename FUNCTION> IntervalSimpson(FUNCTION const& f, double a,
        double b, double fa, double fm, double fb)
            {initHelper(f, a, b, fa, fm, fb);}
public:
    template<typename FUNCTION> IntervalSimpson(FUNCTION const& f, double a,
        double b){initHelper(f, a, b, f(a), f(a + (b - a)/4), f(b));}
    double integrate()const
    {
        return integrateHelper(fLeft, fLM, fMiddle, dx) +
            integrateHelper(fMiddle, fMR, fRight, dx);
    }
    double length()const{return 4 * dx;}
    double error()const
    {
        return abs(integrate() -
            integrateHelper(fLeft, fMiddle, fRight, 2 * dx));
    }
    template<typename FUNCTION>
    Vector<IntervalSimpson> split(FUNCTION const& f)const
    {
        Vector<IntervalSimpson> result;
        double middle = left + 2 * dx, right = middle + 2 * dx;
        result.append(IntervalSimpson(f, left, middle, fLeft, fLM, fMiddle));
        result.append(IntervalSimpson(f, middle, right, fMiddle, fMR, fRight));
        return result;
    }
    static int initEvals(){return 5;}
    static int splitEvals(){return 4;}
};


template<typename FUNCTION> pair<double, double> integrateHybridSimpson(
    FUNCTION const& f, double a, double b, double eRelAbs = highPrecEps,
    int maxEvals = 1000000, int minEvals = -1)
{
    int CCEvals = min(maxEvals/2, 1000);
    pair<double, double> resultCC = integrateCC(f, a, b, CCEvals);
    if(isEEqual(resultCC.first, resultCC.first + resultCC.second, eRelAbs))
        return resultCC;
    pair<double, double> resultAS = integrateAdaptiveHeap<IntervalSimpson>(f,
        a, b, eRelAbs, maxEvals - CCEvals, minEvals);
    return resultCC.second < resultAS.second ? resultCC : resultAS;
}

template<typename FUNCTION> void testIntegrateHelper(FUNCTION const& f, double left, double right,
    Vector<Vector<string> > & matrix)
{
    double answer = f.integralAB();
    int n = 1000000;//n = 10000000//this larger n caused out of memory barycentric
    DEBUG("CC");
    matrix.lastItem().append("CC");
    printResult2(integrateCC(f, left, right), answer, matrix);
    DEBUG("CCDoubling");
    matrix.lastItem().append("CCDoubling");
    printResult2(integrateCCDoubling(f, left, right), answer, matrix);
    DEBUG("AdaptSimp");
    matrix.lastItem().append("AS");
    printResult2(integrateAdaptiveHeap<IntervalSimpson>(f, left, right), answer, matrix);
    DEBUG("AdaptGLK");
    matrix.lastItem().append("AGLK");
    printResult2(integrateAdaptiveHeap<IntervalGaussLobattoKronrod>(f, left, right), answer, matrix);
    DEBUG("Hybr");
    matrix.lastItem().append("Hybr");
    printResult2(integrateHybrid(f, left, right), answer, matrix);
    DEBUG("MC");
    matrix.lastItem().append("MC");
    printResult2(MonteCarloIntegrate(Vector<pair<double, double> >(1, make_pair(left, right)),
        n, InsideTrue(), MCIFunctionHelper<FUNCTION>(f)), answer, matrix);
    //put Sobol here?
    DEBUG("TrapData");
    matrix.lastItem().append("TrapData");
    Vector<pair<double, double> > evals;
    evals.append(pair<double, double>(left, f(left)));
    evals.append(pair<double, double>(right, f(right)));
    for(int i = 0; i < n - 2; ++i)
    {
        double x = GlobalRNG().uniform(left, right);
        evals.append(pair<double, double>(x, f(x)));
    }
    printResult2(integrateFromData(evals), answer, matrix);
}

void testIntegrators()
{
    Vector<Vector<string> > matrix;
    Vector<TestFunctions1D::MetaF> fs = TestFunctions1D::getFunctions();
    for(int i = 0; i < fs.getSize(); ++i)
    {
        if(!isfinite(fs[i].integralAB())) continue;
        string name = fs[i].getName();
        DEBUG(name);
        matrix.append(Vector<string>());
        matrix.lastItem().append(name);
        testIntegrateHelper(fs[i], fs[i].getA(), fs[i].getB(), matrix);
    }
    int reportNumber = time(0);
    string filename = "reportIntegrate" + to_string(reportNumber) + ".csv";
    createCSV(matrix, filename.c_str());
    Vector<string> names;
    names.append("Error");
    names.append("Estimate");
    names.append("Evals");
    createAugmentedCSVFiles(matrix, names, filename);
}

struct MultivarITest1
{//from https://en.wikipedia.org/wiki/Multiple_integral
    //for more see https://www.wolframalpha.com/examples/Integrals.html
    int getD()const{return 2;}
    double operator()(Vector<double> const& x)const
    {
        assert(x.getSize() == getD());
        return x[0] * x[0] + 4 * x[1];
    }
    Vector<pair<double, double> > getBox()const
    {
        Vector<pair<double, double> > box(2);
        box[0] = make_pair(11.0, 14.0);
        box[1] = make_pair(7.0, 10.0);
        return box;
    }
    double solution()const{return 1719;}
};
template<typename INTEGRATOR1D, typename FUNCTION>
testMultidimIntegrationHelper(int evalsPerDim = 33)
{
    FUNCTION f;
    DEBUG(f.solution());
    RecursiveIntegralFunction<FUNCTION, INTEGRATOR1D> rif(f.getBox(),
        evalsPerDim);
    pair<double, double> result = rif.integrate();
    DEBUG(result.first);
    DEBUG(result.second);
}
void testMultidimIntegration()
{//also try smooth functions like multidim exponential on 01 or multidim normal pdf because know answer! good choice
    int nEvals = 1000000;
    testMultidimIntegrationHelper<Cheb1DIntegrator, MultivarITest1>(
        nextPowerOfTwo(pow(nEvals, 1.0/MultivarITest1().getD())/2) + 1);//for Cheb specifically
    MultivarITest1 f;
    double resultMonte = MonteCarloIntegrate(f.getBox(), nEvals, InsideTrue(), f).first;
    DEBUG(resultMonte);
    double resultSobol = SobolIntegrate(f.getBox(), nEvals, InsideTrue(), f).first;
    DEBUG(resultSobol);
}

int main()
{
    testIntegrators();
    return 0;
    testMultidimIntegration();
    return 0;

    return 0;
}
