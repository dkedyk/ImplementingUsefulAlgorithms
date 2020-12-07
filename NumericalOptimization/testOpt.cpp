#include <cassert>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <limits>
#include <memory> //for shared ptr
#include "NumericalOptimization.h"
#include "testOptCommon.h"
#include "SPSA.h"
using namespace std;
using namespace igmdk;

struct TestFunctionsMin
{
    struct BaseF
    {
        virtual Vector<double> operator()(Vector<double> const& x)const = 0;
        virtual string name()const = 0;
        virtual Vector<double> getX0()const = 0;
        virtual Vector<double> getAnswer()const = 0;
    };
    struct ExtendedRosenbrock: public BaseF
    {//From Dennis & Schnabel
        int n;
        ExtendedRosenbrock(int theN = 2): n(theN) {assert(theN % 2 == 0);}
        Vector<double> operator()(Vector<double> const& x)const
        {
            Vector<double> fx(n);
            for(int i = 0; i < n/2; ++i)
            {
                double i1 = 2 * i, i2 = i1 + 1;
                fx[i1] = 10 * (x[i2] - x[i1] * x[i1]);
                fx[i2] = 1 - x[i1];
            }
            return fx;
        }
        string name()const{return "ExtendedRosenbrock" + to_string(n);}
        Vector<double> getX0()const
        {
            Vector<double> x0 = Vector<double>(n, 1);
            for(int i = 0; i < n/2; ++i) x0[2 * i] = 1.2;
            return x0;
        }
        Vector<double> getAnswer()const{return Vector<double>(n, 1);}
    };
    struct ExtendedPowellSingular: public BaseF
    {//From Dennis & Schnabel, J singular at solution
        int n;
        ExtendedPowellSingular(int theN = 4): n(theN) {assert(theN % 4 == 0);}
        Vector<double> operator()(Vector<double> const& x)const
        {
            Vector<double> fx(n);
            for(int i = 0; i < n/4; ++i)
            {
                double i1 = 4 * i, i2 = i1 + 1, i3 = i2 + 1, i4 = i3 + 1;
                fx[i1] = x[i1] + 10 * x[i2];
                fx[i2] = sqrt(5) * (x[i3] - x[i4]);
                fx[i3] = (x[i2] - 2 * x[i3]) * (x[i2] - 2 * x[i3]);
                fx[i4] = sqrt(10) * (x[i1] - x[i4]) * (x[i1] - x[i4]);
            }
            return fx;
        }
        string name()const{return "ExtendedPowellSingular" + to_string(n);}
        Vector<double> getX0()const
        {
            Vector<double> x0 = Vector<double>(n, 1);
            for(int i = 0; i < n/4; ++i)
            {
                x0[4 * i] = 3;
                x0[4 * i + 1] = -1;
                x0[4 * i + 2] = 0;
                x0[4 * i + 2] = 1;
            }
            return x0;
        }
        Vector<double> getAnswer()const{return Vector<double>(n, 0);}
    };
    struct HelicalValley: public BaseF
    {//From Dennis & Schnabel
        Vector<double> operator()(Vector<double> const& x)const
        {
            Vector<double> fx(3);
            double q = 0.5/PI() * atan(x[1]/x[0]);
            if(x[0] < 0) q += 0.5;
            fx[0] = 10 * (x[2] - 10 * q);
            fx[1] = 10 * (sqrt(x[0] * x[0] + x[1] * x[1]) - 1);
            fx[2] = x[2];
            return fx;
        }
        string name()const{return "HelicalValley";}
        Vector<double> getX0()const
        {
            Vector<double> x0(3, 0);
            x0[0] = -1;
            return x0;
        }
        Vector<double> getAnswer()const{return -getX0();}
    };
    struct VariableDimensionF: public BaseF
    {//From More et al
        int n;
        VariableDimensionF(int theN = 2): n(theN) {}
        Vector<double> operator()(Vector<double> const& x)const
        {
            Vector<double> fx(n + 2);
            for(int i = 0; i < n; ++i)
            {
                fx[i] = x[i] - 1;
                fx[n] += (i + 1) * fx[i];
            }
            fx[n + 1] = fx[n] * fx[n];
            return fx;
        }
        string name()const{return "VariableDimensionF" + to_string(n);}
        Vector<double> getX0()const
        {
            Vector<double> x0(n, 1);
            for(int i = 0; i < n; ++i) x0[i] -= (i + 1.0)/n;
            return x0;
        }
        Vector<double> getAnswer()const{return Vector<double>(n, 1);}
    };
    struct LinearFFullRank: public BaseF
    {//From More et al
        int n;
        LinearFFullRank(int theN = 2): n(theN) {}
        Vector<double> operator()(Vector<double> const& x)const
        {
            double sum = 0;
            for(int i = 0; i < n; ++i) sum += x[i];
            Vector<double> fx(n);
            for(int i = 0; i < n; ++i) fx[i] = x[i] - 2 * sum/n - 1;
            return fx;
        }
        string name()const{return "LinearFFullRank" + to_string(n);}
        Vector<double> getX0()const{return Vector<double>(n, 1);}
        Vector<double> getAnswer()const{return -Vector<double>(n, 1);}
    };
    struct BrownBadScaled: public BaseF
    {//From More et al
        Vector<double> operator()(Vector<double> const& x)const
        {
            Vector<double> fx(3);
            fx[0] = x[0] - 1000000;
            fx[1] = x[1] - 2.0/1000000;
            fx[2] = x[0] * x[1] - 2;
            return fx;
        }
        string name()const{return "BrownBadScaled";}
        Vector<double> getX0()const{return Vector<double>(2, 1);}
        Vector<double> getAnswer()const
        {
            Vector<double> x0(2, 0);
            x0[0] = 1000000;
            x0[1] = 2.0/1000000;
            return x0;
        }
    };
    struct Beale: public BaseF
    {//From More et al
        Vector<double> operator()(Vector<double> const& x)const
        {
            Vector<double> fx(3);
            fx[0] = 1.5;
            fx[1] = 2.25;
            fx[2] = 2.625;
            for(int i = 0; i < fx.getSize(); ++i)
                fx[i] -= x[0] * (1 - pow(x[1], i + 1));
            return fx;
        }
        string name()const{return "Beale";}
        Vector<double> getX0()const{return Vector<double>(2, 1);}
        Vector<double> getAnswer()const
        {
            Vector<double> x0(2, 0);
            x0[0] = 3;
            x0[1] = 0.5;
            return x0;
        }
    };
    struct BiggsExp6: public BaseF
    {//From More et al;
        Vector<double> operator()(Vector<double> const& x)const
        {
            Vector<double> fx(6);
            for(int i = 0; i < fx.getSize(); ++i)
            {
                double ti = (i + 1.0)/10,
                    yi = exp(-ti) - 5 * exp(-10 * ti) + 3 * exp(-4 * ti);
                fx[i] = x[2] * exp(-ti * x[0]) - x[3] * exp(-ti * x[1]) +
                    x[5] * exp(-ti * x[4]) - yi;
            }
            return fx;
        }
        string name()const{return "BiggsExp6";}
        Vector<double> getX0()const
        {
            Vector<double> x0(6, 1);
            x0[1] = 2;
            return x0;
        }
        Vector<double> getAnswer()const
        {
            Vector<double> a(6);
            a[0] = 1;
            a[1] = 10;
            a[2] = 1;
            a[3] = 5;
            a[4] = 4;
            a[5] = 3;
            return a;
        }
    };
    static int evalCount;
    struct MetaF
    {
        shared_ptr<BaseF> f;
        template<typename F> MetaF(shared_ptr<F> const& theF): f(theF){}
        double operator()(Vector<double> const& x)const
        {
            ++evalCount;
            return norm((*f)(x));
        }
        string getName()const{return f->name();}
        Vector<double> getX0()const{return f->getX0();}
        pair<Vector<double>, double> getAnswer()const
            {return make_pair(f->getAnswer(), norm((*f)(f->getAnswer())));}
    };
    static Vector<MetaF> getFunctions()
    {
        Vector<MetaF> result;
        result.append(MetaF(make_shared<ExtendedRosenbrock>()));
        result.append(MetaF(make_shared<ExtendedPowellSingular>()));
        result.append(MetaF(make_shared<HelicalValley>()));
        result.append(MetaF(make_shared<VariableDimensionF>()));
        result.append(MetaF(make_shared<LinearFFullRank>()));
        result.append(MetaF(make_shared<BrownBadScaled>()));
        result.append(MetaF(make_shared<Beale>()));
        result.append(MetaF(make_shared<BiggsExp6>()));
        result.append(MetaF(make_shared<ExtendedRosenbrock>(10)));
        result.append(MetaF(make_shared<ExtendedPowellSingular>(12)));
        result.append(MetaF(make_shared<VariableDimensionF>(10)));
        result.append(MetaF(make_shared<LinearFFullRank>(10)));
        result.append(MetaF(make_shared<ExtendedRosenbrock>(30)));
        result.append(MetaF(make_shared<ExtendedPowellSingular>(32)));
        result.append(MetaF(make_shared<VariableDimensionF>(30)));
        result.append(MetaF(make_shared<LinearFFullRank>(30)));
        result.append(MetaF(make_shared<ExtendedRosenbrock>(100)));
        result.append(MetaF(make_shared<ExtendedPowellSingular>(100)));
        result.append(MetaF(make_shared<VariableDimensionF>(100)));
        result.append(MetaF(make_shared<LinearFFullRank>(100)));
        //large D functions
        result.append(MetaF(make_shared<ExtendedRosenbrock>(1000)));
        result.append(MetaF(make_shared<ExtendedPowellSingular>(1000)));
        result.append(MetaF(make_shared<VariableDimensionF>(1000)));
        result.append(MetaF(make_shared<LinearFFullRank>(1000)));
        result.append(MetaF(make_shared<ExtendedRosenbrock>(10000)));
        result.append(MetaF(make_shared<ExtendedPowellSingular>(10000)));
        result.append(MetaF(make_shared<VariableDimensionF>(10000)));
        result.append(MetaF(make_shared<LinearFFullRank>(10000)));
        return result;
    }
};
int TestFunctionsMin::evalCount = 0;

//USE SAME BACKTRACK FOR EQ AND OPT USING FUNCTION X *X?
template<typename FUNCTION> bool backtrackLineSearch(FUNCTION const& f,
    Vector<double> const& gradient, Vector<double>& x, double& fx,
    int& maxEvals, Vector<double> const& dx, double yEps)
{
    double fxFirst = fx, dd = -(dx * gradient), minDescent = dd * 0.0001;
    if(!isfinite(minDescent) || minDescent <= 0) return true;
    for(double s = 1; maxEvals > 0; s /= 2)
    {
        double fNew = f(x + dx * s), fGoal = fx - minDescent * s;
        --maxEvals;
        if(fGoal > fNew)
        {
            x += dx * s;
            fx = fNew;
            break;
        }//step to small if goal was with c = 1
        else if(!isELess(fxFirst - dd * s, fx, yEps)) break;
    }
    return !isELess(fx, fxFirst, yEps);//must make good progress
}

template<typename FUNCTION, typename GRADIENT, typename DIRECTIONAL_DERIVATIVE>
    pair<Vector<double>, double> LBFGSMinimizeNW(Vector<double> const& x0,
    FUNCTION const& f, GRADIENT const& g, DIRECTIONAL_DERIVATIVE const& dd,
    int maxEvals = 1000000, double yPrecision = highPrecEps,
    int historySize = 8, bool useExact = false)
{
    typedef Vector<double> V;
    Queue<pair<V, V> > history;
    pair<V, double> xy(x0, f(x0));
    V grad = g(xy.first), d = -grad;
    int D = xy.first.getSize(), gEvals = g.fEvals(D);
    maxEvals -= 1 + gEvals;
    while(maxEvals > 0)
    {//backtrack using d to get sufficient descent
        pair<V, double> xyOld = xy;
        bool failed = goldenSectionLineSearch(f, grad, dd, xy.first, xy.second,
            maxEvals, d, yPrecision, useExact);
        if(failed || (maxEvals -= gEvals) <= 0) break;
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
        d *= 1/(dotProduct(history[last].second, history[last].second) * p[last]);
        for(int i = 0; i < history.getSize(); ++i)
        {
            double bi = dotProduct(history[i].second, d) * p[last - i];
            d += history[i].first * (a[last - i] - bi);
        }
        d *= -1;
    }
    return xy;
}

template<typename FUNCTION, typename SUBGRADIENT> pair<Vector<double>, double>
    subgradientDescent(Vector<double> const& x0, FUNCTION const& f,
    SUBGRADIENT const& g, int maxEvals = 1000000)
{//TODO: scale relative to first subgrad
    pair<Vector<double>, double> xy(x0, f(x0));
    int stepCount = 1;
    while((maxEvals -= g.fEvals(xy.first.getSize()) + 1) > 0)
    {
        Vector<double> subgrad = g(xy.first), step = subgrad * (-max(1.0,
            norm(xy.first))/10/norm(subgrad)/stepCount),
            xNew = xy.first + step;
        double yNew = f(xNew);
        if(isfinite(yNew) && isfinite(normInf(xNew)))
        {
            xy.first = xNew;
            xy.second = yNew;
        }
    }
    return xy;
}

template<typename INCREMENTAL_FUNCTION> double compassMinimizeIncremental(
    INCREMENTAL_FUNCTION &f, int maxEvals = 1000000,
    double xPrecision = highPrecEps)
{//use scaled step
    f.setCurrentDimension(0);
    double yBest = f(f.getXi()), step = 0.1 * max(1.0, norm(f.getX()));
    int D = f.getSize(), nD = 2 * D, nCycleEvals = 0;
    Vector<int> order(nD);
    for(int i = 0; i < nD; ++i) order[i] = i;
    GlobalRNG().randomPermutation(order.getArray(), nD);
    for(; step > xPrecision && f.getEvalCount() < maxEvals; step /= 2)
    {//poll in all directions in random order
        for(int i = 0; i < nD && f.getEvalCount() < maxEvals; ++i)
        {
            int d = order[nCycleEvals++ % nD], j = d/2, sign = d % 2 ? 1 : -1;
            f.setCurrentDimension(j);
            double xjNew = f.getXi() + sign * step, yNew = f(xjNew);
            if((!isfinite(yBest) && isfinite(yNew)) || yNew < yBest)
            {//found good enough step
                f.bind(xjNew);
                yBest = yNew;
                step *= 4;
                break;
            }
        }
        if(nCycleEvals >= nD)//permute only when had enough evals and not
        {//during a cycle
            GlobalRNG().randomPermutation(order.getArray(), nD);
            nCycleEvals = 0;
        }
    }
    return yBest;
}
template<typename FUNCTION> pair<Vector<double>, double> compassMinimize(
    Vector<double> const& x0, FUNCTION const& f, int maxEvals = 1000000,
    double xPrecision = highPrecEps)
{
    IncrementalWrapper<FUNCTION> iw(f, x0);
    double y = compassMinimizeIncremental(iw, maxEvals, xPrecision);
    return make_pair(iw.xBound, y);
}

template<typename TEST_SET, typename FUNCTION> void debugResultHelper(
    pair<Vector<double>, double> const& result,
    FUNCTION const& f, Vector<Vector<string> > & matrix, int start)
{
    debugResultHelperBatch<TEST_SET>(Vector<pair<Vector<double>, double> >(1, result),
        f, matrix, start);
}

template<typename POINT, typename FUNCTION> void debugResultNew(
    pair<POINT, double> const& result,
    FUNCTION const& f, Vector<Vector<string> > & matrix, int start)
{
    debugResultHelper<TestFunctionsMin>(result, f, matrix, start);
}

template<typename FUNCTION> void testAllSolversLargeD(FUNCTION const& f,
    Vector<Vector<string> >& matrix)
{
    int start = 0;
    GradientFunctor<FUNCTION> g(f);
    DirectionalDerivativeFunctor<FUNCTION> dd(f);
    DEBUG("Compass");
    matrix.lastItem().append("Compass");
    debugResultNew(compassMinimize(f.getX0(), f), f, matrix, start);
    DEBUG("UnimodalCD");
    matrix.lastItem().append("UnimodalCD");
    debugResultNew(unimodalCoordinateDescentGeneral(f, f.getX0()), f, matrix, start);
    DEBUG("LBFGSMinimize");
    matrix.lastItem().append("LBFGSMinimize");
    debugResultNew(LBFGSMinimize(f.getX0(), f, g, dd), f, matrix, start);
}
//for unimodal 2 works, for powell probably just lucky step choice gives it good result!
template<typename FUNCTION> void testAllSolvers(FUNCTION const& f, Vector<Vector<string> >& matrix)
{
    GradientFunctor<FUNCTION> g(f);
    DirectionalDerivativeFunctor<FUNCTION> dd(f);
    int D = f.getX0().getSize(), start = 0;
    DEBUG("metaSPSA");
    matrix.lastItem().append("metaSPSA");
    start = clock();
    debugResultNew(metaSPSA(f.getX0(), f), f, matrix, start);
    DEBUG("Compass");
    matrix.lastItem().append("Compass");
    debugResultNew(compassMinimize(f.getX0(), f), f, matrix, start);
    DEBUG("UnimodalCD");
    matrix.lastItem().append("UnimodalCD");
    start = clock();
    debugResultNew(unimodalCoordinateDescentGeneral(f, f.getX0()), f, matrix, start);
    DEBUG("NelderMead");
    matrix.lastItem().append("NelderMead");
    NelderMead<FUNCTION> nm(D, f);
    start = clock();
    debugResultNew(nm.minimize(f.getX0()), f, matrix, start);
    DEBUG("RestartedNelderMead");
    matrix.lastItem().append("RestartedNelderMead");
    NelderMead<FUNCTION> nmr(D, f);
    start = clock();
    debugResultNew(nmr.restartedMinimize(f.getX0()), f, matrix, start);
    DEBUG("RestartedNelderMead100");
    matrix.lastItem().append("RestartedNelderMead100");
    NelderMead<FUNCTION> nmr100(D, f);
    start = clock();
    debugResultNew(nmr100.restartedMinimize(f.getX0(), 10000, highPrecEps, 100), f, matrix, start);

    DEBUG("SubgradientDescent");
    matrix.lastItem().append("SubgradientDescent");
    start = clock();
    debugResultNew(subgradientDescent(f.getX0(), f, g), f, matrix, start);

    DEBUG("LBFGSMinimizeNW");
    matrix.lastItem().append("LBFGSMinimizeNW");
    start = clock();
    debugResultNew(LBFGSMinimizeNW(f.getX0(), f, g, dd), f, matrix, start);
    DEBUG("LBFGSMinimize");
    matrix.lastItem().append("LBFGSMinimize");
    start = clock();
    debugResultNew(LBFGSMinimize(f.getX0(), f, g, dd), f, matrix, start);
    DEBUG("LBFGSMinimizeMT");
    matrix.lastItem().append("LBFGSMinimizeMT");
    start = clock();
    debugResultNew(LBFGSMinimize(f.getX0(), f, g, dd, 1000000, highPrecEps, 8, false), f, matrix, start);
    DEBUG("HybridLocalMinimize");
    matrix.lastItem().append("HybridLocalMinimize");
    start = clock();
    debugResultNew(hybridLocalMinimize(f.getX0(), f), f, matrix, start);
}

void testAllFunctions()
{
    Vector<Vector<string> > matrix;
    string name;
    Vector<TestFunctionsMin::MetaF> fs = TestFunctionsMin::getFunctions();
    for(int i = 0; i < fs.getSize(); ++i)
    {
        string name = fs[i].getName();
        DEBUG(name);
        int D = fs[i].getX0().getSize();
        if(D >= 1000)
        {
            DEBUG("large scale case");
            continue;
        }
        matrix.append(Vector<string>());
        matrix.lastItem().append(name);
        testAllSolvers(fs[i], matrix);
    }
    createMinReport("reportMin", matrix);
}

void testAllFunctionsLargeD()
{
    Vector<Vector<string> > matrix;
    string name;
    Vector<TestFunctionsMin::MetaF> fs = TestFunctionsMin::getFunctions();
    for(int i = 0; i < fs.getSize(); ++i)
    {

        string name = fs[i].getName();
        DEBUG(name);
        int D = fs[i].getX0().getSize();
        if(D < 1000)
        {
            DEBUG("small scale case");
            continue;
        }
        matrix.append(Vector<string>());
        matrix.lastItem().append(name);
        testAllSolversLargeD(fs[i], matrix);
    }
    createMinReport("reportMinLargeD", matrix);
}

int main()
{
    testAllFunctions();
    return 0;
    testAllFunctionsLargeD();
    return 0;
}
