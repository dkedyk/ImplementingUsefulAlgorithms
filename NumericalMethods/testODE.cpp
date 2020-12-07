#include <cassert>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <limits>
#include <memory> //for shared ptr
#include "ODESolving.h"
#include "TestFunctions1D.h"

#include "../NumericalOptimization/GlobalNumericalOptimization.h"
#include "../RandomNumberGeneration/Statistics.h"
#include "../Utils/DEBUG.h"
#include "../ExternalMemoryAlgorithms/CSV.h"
#include "NumericalMethodsTestAuto.h"
using namespace std;
using namespace igmdk;

int evalCount2 = 0;

void printResult(pair<double, double> const result, double answer)
{
    double relResult = answer == 0 ? result.first : abs(answer - result.first)/answer;
    DEBUG(relResult);
    DEBUG(result.second);
    DEBUG(evalCount2);
    evalCount2 = 0;
}

template<typename TWO_VAR_FUNCTION>
double EulerSolve(TWO_VAR_FUNCTION const& f, double x0, double xGoal,
    double y0, int nIntervals = 10000)
{
    assert(nIntervals > 0);
    double h = (xGoal - x0)/nIntervals, y = y0;
    for(double x = x0; nIntervals--; x += h) y += h * f(x, y);
    return y;
}
template<typename TWO_VAR_FUNCTION>
pair<double, double> adaptiveRungeKutta4(TWO_VAR_FUNCTION const& f, double x0,
    double xGoal, double y0, double localERelAbs = defaultPrecEps,
    int maxIntervals = 100000, int minIntervals = -1,
    double upFactor = pow(2, 0.2))
{
    if(minIntervals == -1) minIntervals = sqrt(maxIntervals);
    assert(xGoal > x0 && minIntervals > 0 && upFactor > 1);
    double hMax = (xGoal - x0)/minIntervals, hMin = (xGoal - x0)/maxIntervals,
        linearError = 0, h1 = hMax, y = y0,
        y1 = numeric_limits<double>::quiet_NaN(), f0;
    bool last = false;
    for(double x = x0; !last;)
    {
        if(x + h1 > xGoal)
        {//make last step accurate
            h1 = xGoal - x;
            last = true;
        }
        if(isnan(y1))
        {
            f0 = f(x, y);
            y1 = RungeKutta4Step(f, x, y, h1, f0);
        }
        double h2 = h1/2, y2 = RungeKutta4Step(f, x, y, h2, f0), firstY2 = y2,
            xFraction = h1/(xGoal - x0);
        y2 = RungeKutta4Step(f, x + h2, y2, h2);
        if(h2 < hMin || isEEqual(y2, y1,
            max(highPrecEps, localERelAbs * sqrt(xFraction))))
        {//accept step
            x += h1;
            y = y2;
            linearError += abs(y2 - y1);
            y1 = numeric_limits<double>::quiet_NaN();
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

template<typename FUNCTION> struct DummyODE1D
{
    FUNCTION f;
    double operator()(double x, double y)const{return f(x);}
};
template<typename FUNCTION> void testODEHelper(FUNCTION const& f, double left, double right,
    Vector<Vector<string> > & matrix)
{
    DummyODE1D<FUNCTION> fxy = {f};
    double exact = f.integralAB();
    DEBUG("Euler10k");
    matrix.lastItem().append("Euler10k");
    printResult2(make_pair(EulerSolve(fxy, left, right, 0), numeric_limits<double>::infinity()), exact, matrix);
    DEBUG("DoublingRK4");
    matrix.lastItem().append("DoublingRK4");
    printResult2(adaptiveRungeKutta4(fxy, left, right, 0), exact, matrix);
    DEBUG("DP");
    matrix.lastItem().append("DP");
    printResult2(adaptiveRungeKuttaDormandPrice(fxy, left, right, 0), exact, matrix);
}

void testODE()
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
        testODEHelper(fs[i], fs[i].getA(), fs[i].getB(), matrix);
    }
    int reportNumber = time(0);
    string filename = "reportODE" + to_string(reportNumber) + ".csv";
    createCSV(matrix, filename.c_str());
    Vector<string> names;
    names.append("Error");
    names.append("Estimate");
    names.append("Evals");
    createAugmentedCSVFiles(matrix, names, filename);
}

struct D2F
{
    double operator()(double x, double y)const{++evalCount2; return x + y;}
};
//bad trap example
struct D2F2
{
    double operator()(double x, double y)const{++evalCount2; return x == 0 ? 1 : y * (-2 * x + 1/x);}
};
void testRungeKutta()
{
    DEBUG("Adaptive2");
    double x = 1;
    double exact = 3 * exp(x) - x - 1;
    DEBUG(exact);
    printResult(adaptiveRungeKutta4(D2F(), 0, 1, 2), exact);
    x = 2;
    exact = x * exp(-x * x);
    DEBUG(exact);
    printResult(adaptiveRungeKutta4(D2F2(), 0, 2, 0), exact);
    DEBUG("AdaptiveDP");
    x = 1;
    exact = 3 * exp(x) - x - 1;
    DEBUG(exact);
    printResult(adaptiveRungeKuttaDormandPrice(D2F(), 0, 1, 2), exact);
    x = 2;
    exact = x * exp(-x * x);
    DEBUG(exact);
    printResult(adaptiveRungeKuttaDormandPrice(D2F2(), 0, 2, 0), exact);
    DEBUG("E,MM");
    x = 1;
    exact = 3 * exp(x) - x - 1;
    DEBUG(exact);
    printResult(make_pair(EulerSolve(D2F(), 0, 1, 2), 0), exact);
    x = 2;
    exact = x * exp(-x * x);
    DEBUG(exact);
    printResult(make_pair(EulerSolve(D2F2(), 0, 2, 0), 0), exact);
}

//from Fausett; ignore x
struct FPStiff1 : public MultivarFuncHelper::F1DBase
{
    double operator()(Vector<double> const& yx)const{++evalCount2; return 98 * yx[0] + 198 * yx[1];}
};
struct FPStiff2 : public MultivarFuncHelper::F1DBase
{
    double operator()(Vector<double> const& yx)const{return -99 * yx[0] - 199 * yx[1];}
};
struct FPStiffDummy : public MultivarFuncHelper::F1DBase
{
    double operator()(Vector<double> const& yx)const{return 0;}
};


template<typename YX_FUNCTION> struct BackwardEulerFunction
{//for Euler don't need f basis and y is fine
    Vector<double> y;
    double x, h;
    YX_FUNCTION f;
    Vector<double> operator()(Vector<double> const& yNext)const
        {return yNext - (y + evalYX(x + h, yNext, f) * h);}
};
struct BackwardEulerStepF
{
    template<typename YX_FUNCTION> Vector<double> operator()(
        YX_FUNCTION const& f, double x, Vector<double> const& y, double h,
        double solveERelAbs)const
    {
        BackwardEulerFunction<YX_FUNCTION> bef = {y, x, h, f};
        return solveBroydenHybrid(bef, y, solveERelAbs).first;
    }
};

template<typename YX_FUNCTION> struct ImplicitTrapezoidFunction
{
    Vector<double> y, fx;
    double x, h;
    YX_FUNCTION f;
    Vector<double> operator()(Vector<double> const& fSumNext)const
        {return fSumNext - (fx + evalYX(x + h, y + fSumNext * h, f)) * 0.5;}
};
struct ImplicitTrapezoidStepF
{
    template<typename YX_FUNCTION> Vector<double> operator()(
        YX_FUNCTION const& f, double x, Vector<double> const& y, double h,
        double solveERelAbs)const
    {//ignore f0 reuse for simplicity
        ImplicitTrapezoidFunction<YX_FUNCTION> itf =
            {y, evalYX(x, y ,f), x, h, f};
        return y + solveBroydenHybrid(itf, Vector<double>(y.getSize(), 0),
            max(solveERelAbs, numeric_limits<double>::epsilon() * normInf(y)/h)
            ).first * h;
    }
};
template<typename YX_FUNCTION, typename STEPF> Vector<double>
    ImplicitFixedStepper(YX_FUNCTION const& f, STEPF const& s, double x0,
    double xGoal, Vector<double> y, int nIntervals = 100000)
{
    assert(nIntervals > 0);
    double h = (xGoal - x0)/nIntervals;
    for(double x = x0; nIntervals--; x += h) y = s(f, x, y, h, defaultPrecEps);
    return y;
}

void testStiffSolvers()
{
    MultivarFuncHelper f;
    FPStiff1 f1;
    FPStiff2 f2;
    FPStiffDummy fd;
    f.fs.append(&f1);
    f.fs.append(&f2);
    f.fs.append(&fd);
    Vector<double> y0(2);
    y0[0] = 1;
    y0[1] = 0;
    DEBUG("Imp Euler");
    Vector<double> result = ImplicitFixedStepper(f, BackwardEulerStepF(), 0, 1, y0);
    DEBUG("result");
    result.debug();
    Vector<double> correctResult(2);
    correctResult[0] = 2 * exp(-1) - exp(-100);
    correctResult[1] = -exp(-1) + exp(-100);
    DEBUG("correctResult");
    correctResult.debug();
    DEBUG(normInf(result - correctResult));
    DEBUG(evalCount2);
    evalCount2 = 0;
    DEBUG("Imp Trap");
    result = ImplicitFixedStepper(f, ImplicitTrapezoidStepF(), 0, 1, y0);
    DEBUG("result");
    result.debug();
    DEBUG(normInf(result - correctResult));
    DEBUG(evalCount2);
    evalCount2 = 0;
    pair<Vector<double>, double> result2;
    DEBUG("Adap Imp Trap");
    result2 = adaptiveStepper(f, ImplicitTrapezoidStepF(), 0, 1, y0);
    result = result2.first;
    DEBUG("result");
    result.debug();
    DEBUG(result2.second);
    DEBUG(normInf(result - correctResult));
    DEBUG(evalCount2);
    evalCount2 = 0;
    DEBUG("RadauIIA5");
    result = ImplicitFixedStepper(f, RadauIIA5StepF(), 0, 1, y0, 100);
    DEBUG("result");
    result.debug();
    DEBUG(normInf(result - correctResult));
    DEBUG(evalCount2);
    evalCount2 = 0;
    DEBUG("Adap RadauIIA5");
    result2 = adaptiveStepper(f, RadauIIA5StepF(), 0, 1, y0);
    result = result2.first;
    DEBUG("result");
    result.debug();
    DEBUG(result2.second);
    DEBUG(normInf(result - correctResult));
    DEBUG(evalCount2);
    evalCount2 = 0;
}

//from Fausett; ignore x
struct FPB1 : public MultivarFuncHelper::F1DBase
{
    double operator()(Vector<double> const& yx)const{++evalCount2; return yx[1];}
};
struct FPB2 : public MultivarFuncHelper::F1DBase
{
    double operator()(Vector<double> const& yx)const{return 2 * yx[0] * yx[1];}
};
struct FPBBV
{
    double evaluateGoal(Vector<double> const& yGoal)const{return yGoal[0] + yGoal[1] - 0.25;}
    Vector<double> getInitial(double b)const
    {
        Vector<double> y0(2);
        y0[0] = 1;
        y0[1] = b;
        return y0;
    };
};
void testBoundaryValue()
{
    MultivarFuncHelper f;
    FPB1 f1;
    FPB2 f2;
    FPStiffDummy fd;
    f.fs.append(&f1);
    f.fs.append(&f2);
    f.fs.append(&fd);
    Vector<double> xPoints(5);
    for(int i = 1; i <= 5; ++i) xPoints[i - 1] = i * 0.2;
    Vector<Vector<double> > result = solveBoundaryValue(f, 0, 1, xPoints, FPBBV());
    DEBUG("result");
    for(int i = 0; i < result.getSize(); ++i)
    {
        DEBUG(i);
        DEBUG(xPoints[i]);
        DEBUG("result[i]");
        result[i].debug();
    }
    DEBUG(evalCount2);
    evalCount2 = 0;
}

int main()
{
    testODE();
    return 0;
    testStiffSolvers();
    return 0;
    testBoundaryValue();
    return 0;
    testRungeKutta();

    return 0;
}
