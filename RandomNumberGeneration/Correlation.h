#ifndef IGMDK_CORRELATION_H
#define IGMDK_CORRELATION_H
#include "Statistics.h"
#include <cmath>
namespace igmdk{

double PearsonCorrelation(Vector<pair<double, double> > const& a)
{
    int n = a.getSize();
    assert(n > 1);
    IncrementalStatistics x, y;
    for(int i = 0; i < n; ++i)
    {
        x.addValue(a[i].first);
        y.addValue(a[i].second);
    }
    double covSum = 0;
    for(int i = 0; i < a.getSize(); ++i)
        covSum += (a[i].first - x.getMean()) * (a[i].second - y.getMean());
    covSum /= n - 1;
    double result = covSum/sqrt(x.getVariance() * y.getVariance());
    return isfinite(result) ? result : 0;//check for div by 0
}
pair<double, double> PearsonCorrelationConf(double corr, int n, double z = 2)
{
    assert(0 <= corr && corr <= 1 && n > 3);
    double stat = atanh(corr), std = 1/sqrt(n - 3);
    return make_pair(tanh(stat - z * std), tanh(stat + z * std));
}

double SpearmanCorrelation(Vector<pair<double, double> > a)
{
    Vector<double> x, y;
    for(int i = 0; i < a.getSize(); ++i)
    {
        x.append(a[i].first);
        y.append(a[i].second);
    }
    x = convertToRanks(x), y = convertToRanks(y);
    for(int i = 0; i < a.getSize(); ++i)
    {
        a[i].first = x[i];
        a[i].second = y[i];
    }
    return PearsonCorrelation(a);
}
pair<double, double> SpearmanCorrelationConf(double corr, int n, double z = 2)
{
    assert(0 <= corr && corr <= 1 && n > 1);
    double std = 1/sqrt(n - 1);
    return make_pair(corr - z * std, corr + z * std);
}

}//end namespace
#endif
