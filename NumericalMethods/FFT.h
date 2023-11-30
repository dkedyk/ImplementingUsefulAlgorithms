#ifndef IGMDK_FFT_H
#define IGMDK_FFT_H

#include <cmath>
#include <complex>
#include "../Utils/Bits.h"
#include "../Utils/Utils.h"
#include "../Utils/Vector.h"
#include "NumericalCommon.h"

namespace igmdk{

complex<double> unityRootHelper(int j, int n)
    {return exp((j * PI()/n) * complex<double>(0, 1));}
Vector<complex<double> > FFTPower2(Vector<complex<double> > const&x)
{
    int n = x.getSize(), b = lgFloor(n);
    assert(isPowerOfTwo(n));
    typedef complex<double> C;
    Vector<C> result(n);
    for(unsigned int i = 0; i < n; ++i) result[reverseBits(i, b)] = x[i];
    for(int s = 1; s <= b; ++s)
    {
        int m = twoPower(s);
        C wm = unityRootHelper(-2, m);
        for(int k = 0; k < n; k += m)
        {
            C w(1, 0);
            for(int j = 0; j < m/2; ++j, w *= wm)
            {
                C t = w * result[k + j + m/2], u = result[k + j];
                result[k + j] = u + t;
                result[k + j + m/2] = u - t;
            }
        }
    }
    return result;
}

Vector<complex<double> > IFFTHelper(Vector<complex<double> > fftx)
{
    int n = fftx.getSize();
    fftx.reverse(1, n - 1);
    return fftx * (1.0/n);
}
Vector<complex<double> > inverseFFTPower2(Vector<complex<double> > const& x)
{
    assert(isPowerOfTwo(x.getSize()));
    return IFFTHelper(FFTPower2(x));
}

Vector<complex<double> > convolutionPower2(Vector<complex<double> > const& a,
    Vector<complex<double> > const& b)
{
    int n = a.getSize();
    assert(n == b.getSize() && isPowerOfTwo(n));
    Vector<complex<double> > fa = FFTPower2(a), fb = FFTPower2(b);
    for(int i = 0; i < n; ++i) fa[i] *= fb[i];
    return inverseFFTPower2(fa);
}
Vector<complex<double> > FFTGeneral(Vector<complex<double> > const& x)
{//Bluestein's algorithm
    int n = x.getSize(), m = nextPowerOfTwo(2 * n - 1);
    if(isPowerOfTwo(n)) return FFTPower2(x);
    Vector<complex<double> > a(m), b(m);//0-padded by default constructor
    for(int j = 0; j < n; ++j)
    {
        a[j] = x[j] * unityRootHelper(-j * j, n);
        b[j] = unityRootHelper(j * j, n);//could precompute b its fft
        if(j > 0) b[m - j] = b[j];
    }
    Vector<complex<double> > ab = convolutionPower2(a, b);
    while(ab.getSize() > n) ab.removeLast();
    for(int k = 0; k < n; ++k) ab[k] *= unityRootHelper(-k * k, n);
    return ab;
}
Vector<complex<double> > IFFTGeneral(Vector<complex<double> > const& x)
    {return IFFTHelper(FFTGeneral(x));}

pair<Vector<complex<double> >, Vector<complex<double> > > FFTReal2Seq(
    Vector<double> const& x, Vector<double> const& y)
{
    int n = x.getSize();
    assert(n == y.getSize());
    typedef complex<double> C;
    typedef Vector<C> VC;
    VC z(n);
    for(int i = 0; i < n; ++i) z[i] = C(x[i], y[i]);
    VC zf = FFTGeneral(z);
    pair<VC, VC> result;
    for(int i = 0; i < n; ++i)
    {
        C temp = conj(zf[(n - i) % n]);
        result.first.append(0.5 * (zf[i] + temp));
        result.second.append(0.5 * (zf[i] - temp) * C(0, -1));
    }
    return result;
}
Vector<complex<double> > FFTRealEven(Vector<double> const& x)
{
    int n = x.getSize(), n2 = n/2;
    assert(n % 2 == 0);
    typedef complex<double> C;
    typedef Vector<C> VC;
    Vector<double> xOdd(n2), xEven(n2);
    for(int i = 0; i < n; ++i) (i % 2 ? xOdd[i/2] : xEven[i/2]) = x[i];
    pair<VC, VC> xSplitF = FFTReal2Seq(xEven, xOdd);
    VC xF(n);
    C wn = unityRootHelper(-2, n), wi(1, 0);
    for(int i = 0; i < n2; ++i, wi *= wn)
    {
        xF[i] = xSplitF.first[i] + wi * xSplitF.second[i];
        xF[n2 + i] = xSplitF.first[i] - wi * xSplitF.second[i];
    }
    return xF;
}

Vector<double> DCTI(Vector<double> const& x)
{
    int n = x.getSize() - 1;
    assert(n > 0);
    Vector<double> y(2 * n), result(n + 1);
    for(int i = 0; i <= n; ++i) y[i] = x[i];
    for(int i = 1; i < n; ++i) y[2 * n - i] = x[i];
    Vector<complex<double> > yf = FFTRealEven(y);
    for(int i = 0; i <= n; ++i) result[i] = yf[i].real()/2;
    return result;
}
Vector<double> IDCTI(Vector<double> const& x)
    {return DCTI(x) * (2.0/(x.getSize() - 1));}

}//end namespace
#endif
