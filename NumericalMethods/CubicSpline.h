#ifndef IGMDK_CUBIC_SPLINE_H
#define IGMDK_CUBIC_SPLINE_H

#include "Matrix.h"
#include "Interpolation.h"
namespace igmdk{

class NotAKnotCubicSplineInterpolation
{//book M is a, book d is R
    struct Data{double a, b, c;};
    PiecewiseData<Data> pd;
public:
    NotAKnotCubicSplineInterpolation(Vector<pair<double, double> > xy,
        double eRelAbs = numeric_limits<double>::epsilon()): pd(eRelAbs)
    {//filter points to ensure eps x distance
        int n = xy.getSize(), skip = 0;
        assert(isSorted(xy.getArray(), 0, n - 1,
            PairFirstComparator<double, double>()));
        for(int i = 1; i + skip < n; ++i)
        {
            if(!isELess(xy[i - 1].first, xy[i + skip].first, eRelAbs)) ++skip;
            if(i + skip < n) xy[i] = xy[i + skip];
        }
        while(skip--) xy.removeLast();
        n = xy.getSize();
        assert(n > 3);//need 4 points to fit a cubic
        //special logic for endpoints
        double h1 = xy[1].first - xy[0].first, h2 = xy[2].first - xy[1].first,
            t1 = (xy[1].second - xy[0].second)/h1,
            t2 = (xy[2].second - xy[1].second)/h2,
            hnm2 = xy[n - 2].first - xy[n - 3].first,
            hnm1 = xy[n - 1].first - xy[n - 2].first,
            tnm2 = (xy[n - 2].second - xy[n - 3].second)/hnm2,
            tnm1 = (xy[n - 1].second - xy[n - 2].second)/hnm1;
        //take out points 1 and n - 2
        for(int i = 1; i < n - 3; ++i) xy[i] = xy[i + 1];
        xy[n - 3] = xy[n - 1];
        xy.removeLast();
        xy.removeLast();
        n = xy.getSize();
        //setup and solve tridiagonal system
        TridiagonalMatrix<double> T(
            2 * TridiagonalMatrix<double>::identity(n));
        Vector<double> R(n);
        //boundary conditions
        double D0Factor = 2/(2 * h2 + h1), Dnm1FActor = 2/(2 * hnm2 + hnm1);
        T(0, 1) = (2 * h1 + h2) * D0Factor;
        R[0] = 6 * (t2 - t1) * D0Factor;
        T(n - 1, n - 2) = (2 * hnm1 + hnm2) * Dnm1FActor;
        R[n - 1] = 6 * (tnm1 - tnm2) * Dnm1FActor;
        for(int i = 1; i < n - 1; ++i)
        {
            double hk = xy[i].first - xy[i - 1].first,
                tk = (xy[i].second - xy[i - 1].second)/hk,
                hkp1 = xy[i + 1].first - xy[i].first, hSum = hk + hkp1,
                tkp1 = (xy[i + 1].second - xy[i].second)/hkp1;
            R[i] = 6 * (tkp1 - tk)/hSum;
            T(i, i + 1) = hkp1/hSum;
            T(i, i - 1) = hk/hSum;
        }
        Vector<double> a = solveTridiag(T, R);
        //compute b and c
        for(int i = 1; i < n + 1; ++i)
        {
            double bi = 0, ci = 0;
            if(i < n)
            {
                double hi = xy[i].first - xy[i - 1].first;
                bi = (xy[i].second - xy[i - 1].second)/hi -
                    (a[i] - a[i - 1]) * hi/6;
                ci = xy[i - 1].second - a[i - 1] * hi * hi/6;
            }
            Data datai = {a[i - 1], bi, ci};
            pd.insert(xy[i - 1].first, datai);
        }
    }
    double operator()(double x, int deriv = 0)const
    {
        assert(deriv >= 0 && deriv <= 2);//support 2 continuous derivatives
        if(!pd.isInERange(x)) return numeric_limits<double>::quiet_NaN();
        typedef typename PiecewiseData<Data>::NODE NODE;
        pair<NODE*, NODE*> segment = pd.findPiece(x);
        assert(segment.first && segment.second);//sanity check
        double dxl = x - segment.first->key, dxr = segment.second->key - x,
            aim1 = segment.first->value.a, ai = segment.second->value.a,
            bi = segment.first->value.b, ci = segment.first->value.c,
            hi = dxr + dxl;
        if(deriv == 2) return (aim1 * dxr + ai * dxl)/hi;
        if(deriv == 1) return (-aim1 * dxr * dxr + ai * dxl * dxl)/(hi * 2) +
            bi;
        return (aim1 * dxr * dxr * dxr + ai * dxl * dxl * dxl)/(hi * 6) +
            bi * dxl + ci;
    }
};

}
#endif
