#ifndef IGMDK_INTERPOLATION_H
#define IGMDK_INTERPOLATION_H
#include <cmath>
#include <iomanip>
#include "../Utils/Bits.h"
#include "../Utils/Vector.h"
#include "../Utils/Utils.h"
#include "../Sorting/Sort.h"
#include "../Utils/Queue.h"
#include "../Heaps/Heap.h"
#include "../RandomTreap/Treap.h"
#include "Matrix.h"
#include "FFT.h"
namespace igmdk{

template<typename DATA> class PiecewiseData
{
    typedef Treap<double, DATA> POINTS;
    mutable POINTS values;
    double eRelAbs;
public:
    int getSize()const{return values.getSize();}
    PiecewiseData(double theERelAbs = numeric_limits<double>::epsilon()):
        eRelAbs(theERelAbs){}
    typedef typename POINTS::NodeType NODE;
    pair<NODE*, NODE*> findPiece(double x)const
    {
        assert(!isnan(x));
        NODE* left = values.inclusivePredecessor(x), *right = 0;
        if(!left) right = values.findMin();
        else
        {
            typename POINTS::Iterator i(left);
            ++i;
            right = i != values.end() ? &*i : 0;
        }
        return make_pair(left, right);
    }
    NODE* eFind(double x)const
    {
        assert(!isnan(x));
        pair<NODE*, NODE*> piece = findPiece(x);
        if(piece.first && isEEqual(x, piece.first->key, eRelAbs))
            return piece.first;
        if(piece.second && isEEqual(x, piece.second->key, eRelAbs))
            return piece.second;
        return 0;
    }
    NODE* findMin()const
    {
        assert(values.getSize() > 0);
        return values.findMin();
    }
    NODE* findMax()const
    {
        assert(values.getSize() > 0);
        return values.findMax();
    }
    bool isInERange(double x)const
        {return getSize() > 1 && findMin()->key <= x && x <= findMax()->key;}
    void insert(double x, DATA const& y)
    {
        NODE* node = eFind(x);
        if(node) node->value = y;
        else values.insert(x, y);
    }
    void eRemove(double x)
    {
        NODE* node = eFind(x);
        if(node) values.removeFound(node);
    }
    Vector<pair<double, DATA> > getPieces()const
    {
        Vector<pair<double, DATA> > result;
        for(typename POINTS::Iterator iter = values.begin();
            iter != values.end(); ++iter)
            result.append(make_pair(iter->key, iter->value));
        return result;
    }
};

class DynamicLinearInterpolation
{
    PiecewiseData<double> pd;
public:
    DynamicLinearInterpolation(){}
    DynamicLinearInterpolation(Vector<pair<double, double> > const& xy)
    {//filter points to ensure eps x distance, use machine eps
        assert(isSorted(xy.getArray(), 0, xy.getSize() - 1,
            PairFirstComparator<double, double>()));
        for(int i = 0, lastGood = 0; i < xy.getSize(); ++i)
            if(i == 0 || isELess(xy[lastGood].first, xy[i].first))
            {
                insert(xy[i].first, xy[i].second);
                lastGood = i;
            }
    }
    double operator()(double x)const
    {
        if(!pd.isInERange(x)) return numeric_limits<double>::quiet_NaN();
        typedef typename PiecewiseData<double>::NODE NODE;
        pair<NODE*, NODE*> segment = pd.findPiece(x);
        assert(segment.first && segment.second);//sanity check
        double ly = segment.first->value, lx = segment.first->key;
        return ly + (segment.second->value - ly) * (x - lx)/
            (segment.second->key - lx);
    }
    void eRemove(double x){pd.eRemove(x);}
    bool eContains(double x){return pd.eFind(x) != 0;}
    void insert(double x, double y){pd.insert(x, y);}
};

template<typename INTERPOLANT> class GenericPiecewiseInterpolation
{
    PiecewiseData<INTERPOLANT> pd;
    const INTERPOLANT* findInterpolant(double x)const
    {
        if(!pd.isInERange(x)) return 0;
        typename PiecewiseData<INTERPOLANT>::NODE* segment =
            pd.findPiece(x).first;
        assert(segment);//sanity check
        return &segment->value;
    }
public:
    GenericPiecewiseInterpolation(
        PiecewiseData<INTERPOLANT> const& thePd): pd(thePd){}
    double operator()(double x)const
    {
        const INTERPOLANT* i = findInterpolant(x);
        return i ? (*i)(x) : numeric_limits<double>::quiet_NaN();
    }
    double evalDeriv(double x)const
    {//unstable don't present!
        const INTERPOLANT* i = findInterpolant(x);
        return i ? (*i).evalDeriv(x) : numeric_limits<double>::quiet_NaN();
    }
    GenericPiecewiseInterpolation<INTERPOLANT> deriver()const
    {
        Vector<pair<double, INTERPOLANT> > pieces = pd.getPieces();
        PiecewiseData<INTERPOLANT> result;
        for(int i = 0; i < pieces.getSize(); ++i)//last one dummy
            result.insert(pieces[i].first, i < pieces.getSize() - 1 ?
                pieces[i].second.deriver() : pieces[i].second);
        return GenericPiecewiseInterpolation<INTERPOLANT>(result);
    }
    Vector<pair<pair<double, double>, INTERPOLANT> > getPieces()const
    {
        Vector<pair<double, INTERPOLANT> > pieces = pd.getPieces();
        assert(pieces.getSize() > 1);//1 real, 1 dummy
        Vector<pair<pair<double, double>, INTERPOLANT> > result;
        for(int i = 0; i < pieces.getSize(); ++i)
        {
            if(result.getSize() > 0)//set right boundary of prev piece
                result.lastItem().first.second = pieces[i].first;
            result.append(make_pair(make_pair(pieces[i].first, 0),
                pieces[i].second));
        }
        result.removeLast();//last piece is dummy
        return result;
    }
    double integrate()const
    {
        double sum = 0;
        Vector<pair<pair<double, double>, INTERPOLANT> > pieces = getPieces();
        for(int i = 0; i < pieces.getSize(); ++i)
            sum += pieces[i].second.integrate();
        return sum;
    }
};

class ChebFunction
{
    Vector<double> ci;
    bool converged;
    double ciAbsE()const{return numeric_limits<double>::epsilon() *
        lgCeiling(ci.getSize()) * (1 + normInf(ci));}
    void trim()
    {//remove if ci too small
        int oldN = ci.getSize();
        double cutoff = ciAbsE();
        while(ci.getSize() > 1 && abs(ci.lastItem()) < cutoff)
            ci.removeLast();
        if(oldN - ci.getSize() > 1) converged = true;
    }
    //f values must be sorted values evaled at cos(jPi/n) for 0 <= j <= n
    void ChebFunctionHelper(Vector<double> const& fValues)
    {//DCTI does most work
        ci = DCTI(fValues) * (2.0/(fValues.getSize() - 1));
        ci[0] /= 2;//half first and last
        ci.lastItem() /= 2;
        converged = false;
        trim();
    }
public:
    ChebFunction(Vector<double> const& fValues, int n)
        {ChebFunctionHelper(fValues);}
    template<typename FUNCTION> ChebFunction(FUNCTION const& f, int n)
    {
        assert(n > 0);
        Vector<double> fValues(n + 1);
        for(int i = 0; i <= n; ++i)
            fValues[i] = f(cos(PI() * i/n));
        ChebFunctionHelper(fValues);
    }
    bool hasConverged()const{return converged;}
    double operator()(double x)const
    {
        assert(x >= -1 && x <= 1);
        double d = 0, dd = 0;
        for(int i = ci.getSize() - 1; i >= 0; --i)
        {
            double temp = d;
            d = (i == 0 ? 1 : 2) * x * d - dd + ci[i];
            dd = temp;
        }
        return d;
    }
    ChebFunction integral(double FM1 = 0)
    {//special case for 0 polynomial
        if(ci.getSize() == 1 && ci.lastItem() == 0) return *this;
        Vector<double> result;
        result.append(FM1);
        for(int i = 1; i - 1 < ci.getSize(); ++i)
            result.append(((i - 1 > 0 ? 1 : 2) * ci[i - 1] -
                (i + 1 > ci.getSize() - 1 ? 0 : ci[i + 1]))/2/i);
        ChebFunction cf(*this);
        cf.ci = result;
        cf.ci[0] -= cf(-1);
        return cf;
    }
    ChebFunction derivative()const
    {
        int n = ci.getSize() - 1;
        ChebFunction result(*this);
        if(n == 0) result.ci[0] = 0;
        else
        {
            result.ci = Vector<double>(n, 0);
            for(int i = n; i > 0; --i) result.ci[i - 1] =
                (i + 1 > n - 1 ? 0 : result.ci[i + 1]) + 2 * i * ci[i];
            result.ci[0] /= 2;
        }
        return result;
    }
    double error()const{return converged ? ciAbsE() : abs(ci.lastItem());}
    pair<double, double> integrate()const
    {
        double result = 0;
        for(int i = 0; i < ci.getSize(); i += 2)
            result += 2 * ci[i]/(1 - i * i);
        return make_pair(result, 2 * error());
    }
    Vector<double> findAllRealRoots(double complexAbsE = highPrecEps)const
    {
        int n = ci.getSize() - 1;
        if(n == 0) return Vector<double>(1, 0);//all 0 case
        else if(n == 1) return Vector<double>();//no roots constant poly
        //setup colleague matrix
        Matrix<double> colleague(n, n);
        colleague(0, 1) = 1;
        for(int r = 1; r < n; ++r)
        {
            colleague(r, r - 1) = 0.5;
            if(r + 1 < n) colleague(r, r + 1) = 0.5;
            if(r == n - 1) for(int c = 0; c < n; ++c)
                colleague(r, c) -= ci[c]/(2 * ci[n]);
        }//solve + only keep real roots
        Vector<complex<double> > croots = QREigenHessenberg(colleague);
        Vector<double> result;//remove complex and extrapolated roots
        for(int i = 0; i < croots.getSize(); ++i) if(abs(croots[i].imag()) <
            complexAbsE && -1 <= croots[i].real() && croots[i].real() <= 1)
            result.append(croots[i].real());
        return result;
    }
    static double xToU(double x, double a, double b)
    {
        assert(a <= x && x <= b);
        double u = (2 * x - a - b)/(b - a);
        return u;
    }
    static double uToX(double u, double a, double b)
    {
        assert(-1 <= u && u <= 1);
        return ((b - a) * u + a + b)/2;
    }
};

template<typename FUNCTION>
Vector<double> reuseChebEvalPoints(FUNCTION const& f, Vector<double> const& fx)
{
    int n = 2 * (fx.getSize() - 1);
    assert(isPowerOfTwo(n));
    Vector<double> result;
    for(int i = 0; i <= n; ++i)
        result.append(i % 2 ? f(cos(PI() * i/n)) : fx[i/2]);
    return result;
}
template<typename FUNCTION> ChebFunction adaptiveChebEstimate(
    FUNCTION const& f, int maxEvals = 5000, int minEvals = 17)
{
    int n = minEvals - 1;
    assert(minEvals <= maxEvals && isPowerOfTwo(n));
    Vector<double> fx(n + 1);
    for(int i = 0; i <= n; ++i) fx[i] = f(cos(PI() * i/n));
    ChebFunction che(fx, n);
    while(maxEvals >= fx.getSize() + n && !che.hasConverged())
    {
        fx = reuseChebEvalPoints(f, fx);
        che = ChebFunction(fx, n *= 2);
    }
    return che;
}

template<typename FUNCTION> class ScaledFunctionM11
{//to allow [-1, 1] functions from any range
    FUNCTION f;
    double a, b;
public:
    ScaledFunctionM11(double theA, double theB, FUNCTION const& theF =
        FUNCTION()): f(theF), a(theA), b(theB) {assert(a < b);}
    double operator()(double u)const{return f(ChebFunction::uToX(u, a, b));}
};
struct ScaledChebAB
{//to eval Cheb at any range
    ChebFunction f;
    double a, b;
public:
    template<typename FUNCTION> ScaledChebAB(FUNCTION const& theF, int n,
        double theA, double theB): a(theA), b(theB),
        f(ScaledFunctionM11<FUNCTION>(theA, theB, theF), n) {assert(a < b);}
    ScaledChebAB(Vector<double> const& fValues, double theA, double theB):
        f(fValues, 0), a(theA), b(theB) {assert(a < b);}
    ScaledChebAB(ChebFunction const& theF, double theA, double theB):
        f(theF), a(theA), b(theB) {assert(a < b);}
    double operator()(double x)const{return f(ChebFunction::xToU(x, a, b));}
    pair<double, double> integrate()const
    {
        pair<double, double> result = f.integrate();
        result.first *= (b - a)/2;
        result.second *= (b - a)/2;
        return result;
    }
    double evalDeriv(double x)const
        {return 2/(b - a) * f.derivative()(ChebFunction::xToU(x, a, b));}
    Vector<double> findAllRealRoots()const
    {//default params good enough
        Vector<double> roots =  f.findAllRealRoots();
        for(int i = 0; i < roots.getSize(); ++i)
            roots[i] = ChebFunction::uToX(roots[i], a, b);
        return roots;
    }
};

class IntervalCheb
{
    ScaledChebAB cf;
    pair<double, double> ab;
    int maxEvals;
public:
    typedef ScaledChebAB INTERPOLANT;
    template<typename FUNCTION> IntervalCheb(FUNCTION const& f, double a,
        double b, int theMaxEvals = 64): ab(a, b), cf(f, theMaxEvals, a, b),
        maxEvals(theMaxEvals){}
    Vector<pair<double, ScaledChebAB> > getInterpolants()const
    {
        return Vector<pair<double, ScaledChebAB> >(1, make_pair(ab.first, cf));
    }
    double scaleEstimate()const{return 1;}
    double length()const{return 0;}
    double error()const//larger first
        {return cf.f.hasConverged() ? 0 : ab.second - ab.first;}
    template<typename FUNCTION>
    Vector<IntervalCheb> split(FUNCTION const& f)const
    {
        Vector<IntervalCheb> result;
        double middle = (ab.first + ab.second)/2;
        result.append(IntervalCheb(f, ab.first, middle, maxEvals));
        result.append(IntervalCheb(f, middle, ab.second, maxEvals));
        return result;
    }
    static int initEvals(int maxEvals){return maxEvals;}
    static int splitEvals(int maxEvals){return 2 * maxEvals;}
};

template<typename ADAPTIVE_INTERVAL> struct AdaptiveIntevalComparator
{
    double deltaLength;
    bool operator()(ADAPTIVE_INTERVAL const& lhs,
        ADAPTIVE_INTERVAL const& rhs)const
    {
        return (lhs.length() > deltaLength || rhs.length() > deltaLength) ?
            lhs.length() > rhs.length() : lhs.error() > rhs.error();
    }
};
template<typename INTERPOLATION_INTERVAL, typename FUNCTION>
    pair<GenericPiecewiseInterpolation<
    typename INTERPOLATION_INTERVAL::INTERPOLANT>,
    double> interpolateAdaptiveHeap(FUNCTION const& f, double a, double b,
    double param, double eRelAbs = highPrecEps, int maxEvals = 1000000,
    int minEvals = -1)
{
    typedef INTERPOLATION_INTERVAL II;
    typedef typename II::INTERPOLANT I;
    if(minEvals == -1) minEvals = sqrt(maxEvals);
    assert(a < b && maxEvals >= minEvals && minEvals >= II::initEvals(param));
    INTERPOLATION_INTERVAL i0(f, a, b, param);
    double scale = i0.scaleEstimate();
    AdaptiveIntevalComparator<II> ic = {(b - a)/minEvals};
    Heap<II, AdaptiveIntevalComparator<II> > h(ic);
    h.insert(i0);
    for(int usedEvals = II::initEvals(param);
        usedEvals < minEvals || (
        usedEvals + II::splitEvals(param) <= maxEvals &&
        !isEEqual(scale, scale + h.getMin().error(), eRelAbs));
        usedEvals += II::splitEvals(param))
    {
        II next = h.deleteMin();
        Vector<II> division = next.split(f);
        for(int i = 0; i < division.getSize(); ++i) h.insert(division[i]);
    }//process heap intervals
    PiecewiseData<I> pd;
    double error = h.getMin().error(),
        right = -numeric_limits<double>::infinity();
    while(!h.isEmpty())
    {
        II next = h.deleteMin();
        Vector<pair<double, I> > interpolants = next.getInterpolants();
        for(int i = 0; i < interpolants.getSize(); ++i)
            pd.insert(interpolants[i].first, interpolants[i].second);
    }//need dummy right endpoint interpolator
    pd.insert(b, i0.getInterpolants().lastItem().second);
    return make_pair(GenericPiecewiseInterpolation<I>(pd), error);
}

class BarycentricInterpolation
{
    Vector<pair<double, double> > xy;
    Vector<double> w;
public:
    BarycentricInterpolation(Vector<pair<double, double> >const& thexy)
    {//O(n^2)
        for(int i = 0; i < thexy.getSize(); ++i)
            addPoint(thexy[i].first, thexy[i].second);
    }
    void addPoint(double x, double y)
    {
        double wProduct = 1;
        for(int i = 0; i < xy.getSize(); ++i)
        {
            wProduct *= (x - xy[i].first);
            w[i] /= xy[i].first - x;//must update previous wi
            assert(isfinite(w[i]));//signal repeated point or overflow
        }
        w.append(1/wProduct);
        assert(isfinite(w.lastItem()));
        xy.append(make_pair(x, y));
    }
    void removePoint(int i)
    {
        int n = xy.getSize();
        assert(i >= 0 && i < n);
        for(int j = 0; j < n; ++j)
            if(j != i) w[j] *= (xy[j].first - xy[i].first);
        for(int j = i + 1; j < n; ++j)
        {//make generic vector remove func as nonmember?
            xy[j - 1] = xy[j];
            xy.removeLast();
            w[j - 1] = w[j];
            w.removeLast();
        }
    }
    double operator()(double x)const
    {
        assert(isfinite(x));
        double numSum = 0, denomSum = 0;
        for(int i = 0; i < xy.getSize(); ++i)
        {
            double factorI = w[i]/(x - xy[i].first);
            if(!isfinite(factorI)) return xy[i].second;//inf if x is in xy
            numSum += factorI * xy[i].second;
            denomSum += factorI;
        }
        return numSum/denomSum;
    }
    double evalDeriv(double x)const
    {//unstable don't present!
        assert(isfinite(x));
        double numSum = 0, denomSum = 0, px = (*this)(x);
        for(int i = 0; i < xy.getSize(); ++i)
        {
            double factorI = w[i]/(x - xy[i].first);
            if(!isfinite(factorI))
            {//duplicates no problem
                numSum = 0;
                denomSum = w[i];
                for(int j = 0; j < xy.getSize(); ++j) if(j != i)
                {
                    double factorJ = w[j]/(xy[i].first - xy[j].first);
                    numSum += factorJ * (xy[j].second - xy[i].second);
                }
                break;
            }
            numSum += factorI * (px - xy[i].second)/(x - xy[i].first);
            denomSum += factorI;
        }
        return numSum/denomSum;
    }
    Matrix<double> diffMatrix()const
    {//duplicate points impossible here due to weight filtering
        int n = xy.getSize();
        assert(n > 1);//need at least two points for this to make sense
        Matrix<double> diff(n, n);
        for(int r = 0; r < n; ++r) for(int c = 0; c < n; ++c) if(r != c)
            {
                diff(r, c) = w[c]/w[r]/(xy[r].first - xy[c].first);
                diff(r, r) -= diff(r, c);
            }
        return diff;
    }
    Vector<double> getY()const
    {
        int n = xy.getSize();
        Vector<double> y(n);
        for(int i = 0; i < n; ++i) y[i] = xy[i].second;
        return y;
    }
    BarycentricInterpolation overwriteY(Vector<double> const& y)const
    {
        int n = y.getSize();
        assert(xy.getSize() == n);
        BarycentricInterpolation result = *this;
        for(int i = 0; i < n; ++i) {result.xy[i].second = y[i];}
        return result;
    }
    BarycentricInterpolation deriver()const
        {return overwriteY(diffMatrix() * getY());}
};

}//end namespace
#endif
