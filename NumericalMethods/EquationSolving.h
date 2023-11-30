#ifndef IGMDK_EQUATION_SOLVING_H
#define IGMDK_EQUATION_SOLVING_H
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
#include "Matrix.h"
#include "Differentiation.h"
namespace igmdk{

bool haveDifferentSign(double a, double b){return (a < 0) != (b < 0);}
template<typename FUNCTION> pair<double, double> solveFor0(FUNCTION const& f,
    double xLeft, double xRight, double relAbsXPrecision = highPrecEps)
{
    double yLeft = f(xLeft), xMiddle = xLeft;
    assert(xRight > xLeft && haveDifferentSign(yLeft, f(xRight)));
    while(isELess(xLeft, xRight, relAbsXPrecision))
    {//below formula more robust than simple average
        xMiddle = xLeft + (xRight - xLeft)/2;
        double yMiddle = f(xMiddle);
        if(haveDifferentSign(yLeft, yMiddle)) xRight = xMiddle;
        else
        {
            xLeft = xMiddle;
            yLeft = yMiddle;
        }
    }//best guess and worst-case error
    return make_pair(xLeft + (xRight - xLeft)/2, xRight - xLeft);
}
template<typename FUNCTION> pair<double, double> find1SidedInterval0(
    FUNCTION const& f, double x0 = 0, double d = 0.001, int maxEvals = 30)
{
    assert(maxEvals >= 2 && isfinite(x0 + d) && d != 0 &&
        !isEEqual(x0, x0 + d));
    for(double xLast = x0, f0 = f(x0); --maxEvals > 0; d *= 2)
    {
        double xNext = x0 + d, fNext = f(xNext);
        if(!isfinite(xNext) || isnan(fNext)) break;
        if(haveDifferentSign(f0, fNext))
            return d < 0 ? make_pair(xNext, xLast) : make_pair(xLast, xNext);
        xLast = xNext;
    }
    return make_pair(numeric_limits<double>::quiet_NaN(),
        numeric_limits<double>::quiet_NaN());
}
template<typename FUNCTION> pair<double, double> findInterval0(
    FUNCTION const& f, double x0, double d, int maxEvals)
{
    assert(maxEvals >= 2 && isfinite(x0 + d) && isELess(x0, x0 + d));
    for(double xLast = x0, f0 = f(x0); --maxEvals > 0; d = -d)
    {
        double xNext = x0 + d, fNext = f(xNext);
        if(!isfinite(xNext) || isnan(fNext)) break;
        if(haveDifferentSign(f0, fNext)) return d < 0 ?
            make_pair(xNext, x0 - (xLast - x0)) : make_pair(xLast, xNext);
        if(d < 0)
        {
            xLast = x0 - d;
            d *= 2;
        }
    }
    return make_pair(numeric_limits<double>::quiet_NaN(),
        numeric_limits<double>::quiet_NaN());
}
template<typename FUNCTION> pair<double, double> exponentialSearch(
    FUNCTION const& f, double x0 = 0, double step = 0.001,
    double xERelAbs = highPrecEps, int maxExpEvals = 60)
{
    pair<double, double> i0 = findInterval0(f, x0, step * max(1.0, abs(x0)),
        maxExpEvals);
    return !isnan(i0.first) ? solveFor0(f, i0.first, i0.second, xERelAbs) : i0;
}
template<typename FUNCTION> pair<double, double> exponentialSearch1Sided(
    FUNCTION const& f, double x0 = 0, double step = 0.001,//negate for left
    double xERelAbs = highPrecEps, int maxExpEvals = 60)
{
    pair<double, double> i0 = find1SidedInterval0(f, x0,
        step * max(1.0, abs(x0)), maxExpEvals);
    return !isnan(i0.first) ? solveFor0(f, i0.first, i0.second, xERelAbs) : i0;
}

struct MultivarFuncHelper
{
    struct F1DBase
        {virtual double operator()(Vector<double> const& x)const = 0;};
    Vector<F1DBase*> fs;//beware storage is elsewhere
    Vector<double> operator()(Vector<double> const& x)const
    {
        Vector<double> y(fs.getSize());
        for(int i = 0; i < fs.getSize(); ++i) y[i] = (*fs[i])(x);
        return y;
    }
};

double normInf(double x){return abs(x);}
template<typename FUNCTION, typename X> bool equationBacktrack(
    FUNCTION const& f, X& x, X& fx, int& maxEvals, X const& dx, double xEps)
{
    bool failed = true;
    for(double s = 1;
        maxEvals > 0 && normInf(dx) * s > xEps * (1 + normInf(x)); s /= 2)
    {
        X fNew = f(x + dx * s);
        --maxEvals;
        if(normInf(fNew) <= (1 - 0.0001 * s) * normInf(fx))
        {
            failed = false;
            x += dx * s;
            fx = fNew;
            break;
        }
    }
    return failed;
}

template<typename X> struct BroydenSecant
{//for secant
    static int getD(double dummy){return 1;}
    static double generateUnitStep(double dummy, double infNormSize)
        {return GlobalRNG().normal(0, infNormSize);}
    class InverseOperator
    {
        double b;
    public:
        template<typename FUNCTION> InverseOperator(FUNCTION const& f,
            double x): b(estimateDerivativeCD(f, x)){}
        void addUpdate(double df, double dx){b = df/dx;}
        double operator*(double fx)const{return fx/b;}
    };
};
template<> struct BroydenSecant<Vector<double> >
{//for Broyden
    typedef Vector<double> X;
    static int getD(X const& x){return x.getSize();}
    static X generateUnitStep(X const& x, double infNormSize)
    {
        return GlobalRNG().randomUnitVector(getD(x)) *
            GlobalRNG().normal(0, infNormSize);
    }
    class InverseOperator
    {
        QRDecomposition qr;
    public:
        template<typename FUNCTION> InverseOperator(FUNCTION const& f,
            X const& x): qr(estimateJacobianCD(f, x)){}
        void addUpdate(X const& df, X dx)
        {
            double ndx2 = norm(dx);
            dx *= (1/ndx2);
            qr.rank1Update(df * (1/ndx2) - qr * dx, dx);
        }
        X operator*(X const& fx)const{return qr.solve(fx);}
    };
};
template<typename FUNCTION, typename X> bool equationTryRandomStep(
    FUNCTION const& f, X& x, X& fx, double stepNorm)
{
    X dx = BroydenSecant<X>::generateUnitStep(x, stepNorm), fNew = f(x + dx);
    bool improved = normInf(fNew) < normInf(fx);
    if(improved)
    {
        x += dx;
        fx = fNew;
    }
    return !improved;
}

template<typename FUNCTION, typename X> pair<X, double> solveBroyden(
    FUNCTION const& f, X const& x0, double xEps = highPrecEps,
    int maxEvals = 1000)
{
    //DEBUG("start");
    int D = BroydenSecant<X>::getD(x0), failCount = 0;
    X x = x0, fx = f(x);
    assert(D == BroydenSecant<X>::getD(fx) && maxEvals >= 2 * D + 1);
    typedef typename BroydenSecant<X>::InverseOperator BIO;
    BIO B(f, x);
    maxEvals -= 2 * D + 1;
    double lastGoodNorm = 1, xError = numeric_limits<double>::infinity();
    if(!isfinite(normInf(fx))) return make_pair(x, xError);
    while(maxEvals > 0)
    {
        /*DEBUG(normInf(f(x)));
        DEBUG(normInf(f(x0)));
        DEBUG(normInf(f(x)) - normInf(f(x0)));
        assert(normInf(f(x)) <= normInf(f(x0)));*/
        if(failCount > 1)
        {//after 2nd fail try random step
            --maxEvals;
            //assert(normInf(f(x)) <= normInf(f(x0)));
            if(!equationTryRandomStep(f, x, fx, lastGoodNorm))
            {
                assert(normInf(f(x)) <= normInf(f(x0)));
                if(maxEvals >= 2 * D + 1)//need enough evals for next step
                {
                    B = BIO(f, x);
                    maxEvals -= 2 * D;
                    failCount = 0;//back to normal
                }//else keep making random steps
                continue;
            }
            //assert(normInf(f(x)) <= normInf(f(x0)));
        }
        X dx = B * -fx, oldFx = fx, oldX = x;
        double ndx = normInf(dx);
        if(!isfinite(ndx))
        {//probably singular, after bad update or reestimation
            ++failCount;
            continue;
        }
        if(ndx < xEps * (1 + normInf(x))) break;//full step too small
        if(!equationBacktrack(f, x, fx, maxEvals, dx, xEps))
        {
            //assert(normInf(f(x)) <= normInf(f(x0)));
            xError = lastGoodNorm = ndx;//last successful full step
            failCount = 0;
        }
        else ++failCount;
            //assert(normInf(f(x)) <= normInf(f(x0)));
        if(failCount == 1)//after first fail reestimate J
            if(maxEvals >= 2 * D + 1)//need enough evals for next step
            {
                B = BIO(f, x);
                maxEvals -= 2 * D;
            }
            else ++failCount;//if cant do steps
        else B.addUpdate(fx - oldFx, x - oldX);
        //assert(normInf(f(x)) <= normInf(f(x0)));
    }
    return make_pair(x, xError);
}//for type safety  such such int x0 use wrapper
template<typename FUNCTION> pair<double, double> solveSecant(FUNCTION const&
    f, double const& x0, double xEps = highPrecEps, int maxEvals = 1000)
    {return solveBroyden(f, x0, xEps, maxEvals);}

class BroydenLMInverseOperator
{
    typedef Vector<double> X;
    int m;
    Queue<pair<X, X> > updates;
    struct Result
    {
        X Bfx;
        Vector<X> Bdfs, dxBs;
    };
public:
    BroydenLMInverseOperator(int theM): m(theM){}
    void addUpdate(X const& df, X const& dx)
    {
        if(updates.getSize() == m) updates.pop();
        updates.push(make_pair(df, dx));
    }
    X operator*(X const& fx)const
    {
        int n = updates.getSize();
        X Bfx = fx;//base case identity
        Vector<X> Bdfs(n), dxBs(n);
        for(int i = 0; i < n; ++i)
        {
            Bdfs[i] = updates[i].first;
            dxBs[i] = updates[i].second;
        }
        for(int i = 0; i < n; ++i)
        {
            X u = (updates[i].second - Bdfs[i]) *
                (1/dotProduct(updates[i].second, Bdfs[i]));
            if(!isfinite(normInf(u))) continue;//guard against div by 0
            Bfx += outerProductMultLeft(u, dxBs[i], fx);
            for(int j = i + 1; j < n; ++j)
            {
                Bdfs[j] += outerProductMultLeft(u, dxBs[i],updates[j].first);
                dxBs[j] += outerProductMultRight(u, dxBs[i],
                    updates[j].second);
            }
        }
        return Bfx;
    }
};
template<typename FUNCTION> pair<Vector<double>, double> solveLMBroyden(
    FUNCTION const& f, Vector<double> const& x0, double xEps = highPrecEps,
    int maxEvals = 1000, int m = 30)
{
    int D = x0.getSize();
    BroydenLMInverseOperator B(m);
    Vector<double> x = x0, fx = f(x), xBest = x, fBest = fx;
    double s = 1, xError = numeric_limits<double>::infinity(),eBest = xError;
    while(maxEvals-- > 0)
    {
        Vector<double> dx = B * -fx * s;
        double ndx = normInf(dx)/s;//norm of full step
        //something wrong or step too small
        if(!isfinite(ndx) || ndx < xEps * (1 + normInf(x))) break;
        Vector<double> fNew = f(x + dx);
        if(normInf(fNew) > 10 * normInf(fx)) s /= 10;
        else
        {
            B.addUpdate(fNew - fx, dx);
            fx = fNew;
            x += dx;
            xError = ndx;//use last full step
            s = 1;//back to normal step
            if(normInf(fx) < normInf(fBest))
            {
                xBest = x;
                fBest = fx;
                eBest = xError;
            }
        }
    }
    return make_pair(xBest, eBest);
}

template<typename FUNCTION> pair<Vector<double>, double> solveBroydenHybrid(
    FUNCTION const& f, Vector<double> const& x0, double xEps = highPrecEps,
    int maxEvals = 1000, int changeD = 200)
{
    return x0.getSize() > changeD ? solveLMBroyden(f, x0, xEps, maxEvals) :
        solveBroyden(f, x0, xEps, maxEvals);
}

template<typename FUNCTION> pair<double, double> solveSecantGlobal(
    FUNCTION const& f, double xEps = highPrecEps, int nSamples = 100)
{
    assert(nSamples > 0 && xEps >= numeric_limits<double>::epsilon());
    pair<double, double> best;
    double nBestfx;
    for(int i = 0; i < nSamples; ++i)
    {
        double x = GlobalRNG().Levy() * GlobalRNG().sign();
        pair<double, double> next = solveBroyden(f, x, xEps);
        double nNextfx = normInf(f(next.first));
        if(i == 0 || nNextfx < nBestfx)
        {
            best = next;
            nBestfx = nNextfx;
        }
    }
    return best;
}

template<typename FUNCTION> pair<Vector<double>, double> solveBroydenLevy(
    FUNCTION const& f, int D, double xEps = highPrecEps, int nSamples = 1000)
{
    assert(nSamples > 0 && xEps >= numeric_limits<double>::epsilon());
    typedef Vector<double> X;
    pair<X, double> best;
    double nBestfx;
    for(int i = 0; i < nSamples; ++i)
    {
        X x = GlobalRNG().randomUnitVector(D) * GlobalRNG().Levy();
        pair<X, double> next = solveBroydenHybrid(f, x, xEps);
        double nNextfx = normInf(f(next.first));
        if(i == 0 || nNextfx < nBestfx)
        {
            best = next;
            nBestfx = nNextfx;
        }
    }
    return best;
}

}//end namespace
#endif
