#include "FFT.h"
#include "../RandomNumberGeneration/Random.h"
#include "NumericalMethodsTestAuto.h"
using namespace std;
using namespace igmdk;

void FFTTestReal()
{
    int n = 4;
    Vector<double> x(n);
    for(int i = 0; i < n; ++i) x[i] = i;
    Vector<complex<double> > z(n);
    for(int i = 0; i < n; ++i) z[i] = complex<double>(x[i], 0);
    double normYDiff = normInf(FFTRealEven(x) - FFTGeneral(z));
    DEBUG(normYDiff);
    DEBUG("FFTRealEven(x)");
    FFTRealEven(x).debug();
    DEBUG("FFTGeneral(z))");
    FFTGeneral(z).debug();
}

Vector<double> slowDCTI(Vector<double> const& x)
{
    int n = x.getSize() - 1;
    Vector<double> result(n + 1);
    for(int i = 0; i <= n; ++i)
    {
        double ci = 0;
        for(int j = 0; j <= n; ++j) ci += cos(i * j * PI()/n)
            * x[j] * (j == 0 || j == n ? 0.5 : 1.0);
        result[i] = ci;
    }
    return result;
}
void DCTTestHelper(Vector<double> const& x, double eps = defaultPrecEps)
{
    double normXDiff = normInf(x - IDCTI(DCTI(x))),
        normYDiff = normInf(slowDCTI(x) - DCTI(x));
    if(normXDiff >= eps || normYDiff >= eps)
    {
        DEBUG("failed for x=");
        DEBUG(x.getSize());
        x.debug();
        DEBUG(normXDiff);
        DEBUG(normYDiff);
        DEBUG("IDCTI(DCTI(x))");
        IDCTI(DCTI(x)).debug();
        DEBUG("DCTI(x)");
        DCTI(x).debug();
        DEBUG("slowDCTI(x)");
        slowDCTI(x).debug();
        assert(false);
    }
}
void DCTTestAuto()
{
    int nMax = 100, nn = 1000;
    for(int n = 3; n <= nMax; ++n)//fails for 2 in bits
    {
        for(int j = 0; j < nn; ++j)
        {
            Vector<double> x(n);
            for(int i = 0; i < n; ++i) x[i] = GlobalRNG().uniform(-1, 1);
            DCTTestHelper(x);
        }
    }
    DEBUG("DCTTestAuto passed");
}

int main()
{
    FFTTestReal();
    FFTTestAuto();
    DCTTestAuto();
    return 0;
}
