#ifndef IGMDK_DISTRIBUTION_TESTS_H
#define IGMDK_DISTRIBUTION_TESTS_H
#include "Statistics.h"
#include "../Sorting/Sort.h"
#include <cmath>
namespace igmdk{

double evaluateChiSquaredCdf(double chi, int n)
{
    assert(chi >= 0 && n > 0);
    double m = 5.0/6 - 1.0/9/n - 7.0/648/n/n + 25.0/2187/n/n/n,
        q2 = 1.0/18/n + 1.0/162/n/n - 37.0/11664/n/n/n, temp = chi/n,
        x = pow(temp, 1.0/6) - pow(temp, 1.0/3)/2 + pow(temp, 1.0/2)/3;
    return approxNormalCDF((x - m)/sqrt(q2));
}
double chiSquaredP(Vector<int> const& counts,
    Vector<double> const& means, int degreesOfFreedomRemoved = 0)
{
    double chiStat = 0;
    for(int i = 0; i < counts.getSize(); ++i)
    {//enforce 5 in each bin for good approximation
        assert(means[i] >= 5);
        chiStat += (counts[i] - means[i]) * (counts[i] - means[i])/means[i];
    }
    return 1 - evaluateChiSquaredCdf(chiStat,
        counts.getSize() - degreesOfFreedomRemoved);
}

template<typename CDF> double findMaxKDiff(Vector<double> x, CDF const& cdf)
{//helper to calculate max diff
    quickSort(x.getArray(), x.getSize());
    double level = 0, maxDiff = 0, del = 1.0/x.getSize();
    for(int i = 0; i < x.getSize(); ++i)
    {
        double cdfValue = cdf(x[i]);
        maxDiff = max(maxDiff, abs(cdfValue - level));
        level += del;
        while(i + 1 < x.getSize() && x[i] == x[i + 1])
        {
            level += del;
            ++i;
        }
        maxDiff = max(maxDiff, abs(cdfValue - level));
    }
    return maxDiff;
}
template<typename CDF> double DKWPValue(Vector<double> const& x,
    CDF const& cdf)
{//DKW invalid for p-value < 0.5
    double delta = findMaxKDiff(x, cdf);
    //DEBUG(delta);
    return min(0.5, 2 * exp(-2 * x.getSize() * delta * delta));
}

double findMaxKSDiff(Vector<double> a, Vector<double> b)
{//helper to calculate max diff
    quickSort(a.getArray(), a.getSize());
    quickSort(b.getArray(), b.getSize());
    double aLevel = 0, bLevel = 0, maxDiff = 0, delA = 1.0/a.getSize(),
        delB = 1.0/b.getSize();
    for(int i = 0, j = 0; i < a.getSize() || j < b.getSize();)
    {
        double x, nextX = numeric_limits<double>::infinity();
        bool useB = i >= a.getSize() || (j < b.getSize() && b[j] < a[i]);
        if(useB)
        {
            x = b[j++];
            bLevel += delB;
        }
        else
        {
            aLevel += delA;
            x = a[i++];
        }//handle equal values--process all before diff update
        if(i < a.getSize() || j < b.getSize())
        {
            useB = i >= a.getSize() || (j < b.getSize() && b[j] < a[i]);
            nextX = useB ? b[j] : a[i];
        }
        if(x != nextX) maxDiff = max(maxDiff, abs(aLevel - bLevel));
    }
    return maxDiff;
}
double KS2SamplePValue(Vector<double> const& a, Vector<double> const& b)
{//calculate the adjustment first, then find p-value of d
    double stddev = sqrt(1.0 * (a.getSize() + b.getSize())/
        (a.getSize() * b.getSize())),
        delta = findMaxKSDiff(a, b)/stddev;
    return 2 * exp(-2 * delta * delta);
}

}//end namespace
#endif
