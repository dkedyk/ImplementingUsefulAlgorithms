#include "MCMC.h"
#include "Statistics.h"
#include "../NumericalMethods/Matrix.h"
using namespace igmdk;

struct Normal01SemiPDF
{
    double mean, variance;
    Normal01SemiPDF(double theMean = 100, double theVariance = 10000):
        mean(theMean), variance(theVariance){}
    double operator()(double x)const{x -= mean; return exp(-x * x/2/variance);}
};

struct MultivarNormalSemiPDF
{
    double operator()(Vector<double> x)const
    {
        Matrix<double> b2(3, 3);
        b2(0, 0) = 1;
        b2(0, 1) = 4;
        b2(0, 2) = 5;
        b2(1, 0) = 4;
        b2(1, 1) = 20;
        b2(1, 2) = 32;
        b2(2, 0) = 5;
        b2(2, 1) = 32;
        b2(2, 2) = 64;
        LUP<double> lup(b2);
        Matrix<double> inv = inverse(lup, 3);
        //inv.debug();
        //for(int i = 0; i < x.getSize(); ++i) x[i] -= 100;
        //for(int i = 0; i < x.getSize(); ++i) DEBUG(x[i]);
        //DEBUG(inv * x * x/(-2));
        //DEBUG(exp(inv * x * x/(-2)));
        //system("PAUSE");
        return exp(dotProduct(inv * x, x)/(-2));
    }
};

void testVectorMagicMCMC2()
{
    MultidimGridRWM<MultivarNormalSemiPDF> g(Vector<double>(3, 0.5));
    int n = 10000;
    for(int i = 0; i < n; ++i) g.sample();

    Vector<double> sum(3, 0);
    Matrix<double> outerSum(3, 3);
    for(int i = 0; i < n; ++i)
    {
        Vector<double> x = g.sample();
        sum += x;
        outerSum += outerProduct(x, x);
    }
    Vector<double> mean = sum * (1.0/n);
    for(int i = 0; i < 3; ++i) DEBUG(mean[i]);
    Matrix<double> cov = (outerSum - outerProduct(mean, sum)) * (1.0/(n - 1));
    cov.debug();
}

void testMCMCMagic()
{
    GridRWM<Normal01SemiPDF> s;
    int n = 10;
    for(int i = 0; i < n; ++i) s.sample();
    IncrementalStatistics z;

    for(int i = 0; i < n; ++i) z.addValue(s.sample());
    DEBUG(z.getMean());
    DEBUG(z.getVariance());

}

int main(int argc, char *argv[])
{
    testVectorMagicMCMC2();
    return 0;
    testMCMCMagic();
    return 0;
}
