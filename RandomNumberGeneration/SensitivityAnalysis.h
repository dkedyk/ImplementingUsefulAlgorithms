#ifndef IGMDK_SENSITIVITY_ANALYSIS_H
#define IGMDK_SENSITIVITY_ANALYSIS_H
#include "Statistics.h"
#include "Bootstrap.h"
#include "../Utils/Vector.h"
#include <cmath>
namespace igmdk{

Vector<double> findSobolIndicesHelper(Vector<double> const& ya,
    Vector<Vector<double> > const& yc)
{//calculate S from YA and YC
    int n = ya.getSize(), D = yc.getSize();
    Vector<double> result(D);
    IncrementalStatistics s;
    for(int i = 0; i < n; ++i) s.addValue(ya[i]);
    double f02 = s.getMean() * s.getMean(), tempa = dotProduct(ya, ya)/n;
    for(int j = 0; j < D; ++j)
        result[j] = max(0.0, (dotProduct(ya, yc[j])/n - f02)/(tempa - f02));
    return result;
}
template<typename FUNCTOR> pair<Vector<double>, Vector<pair<double, double>>>
    findSobolIndicesSaltelli(Vector<pair<Vector<double>, double> > const&
    data, FUNCTOR const& f, int nBoots = 200, double a = 0.05)
{//calculate ya and yb
    int D = data[0].first.getSize(), n = data.getSize()/2;
    Vector<double> ya(n), yb(n), yaR(n);
    for(int i = 0; i < 2 * n; ++i)
        if(i < n) ya[i] = data[i].second;
        else yb[i - n] = data[i].second;
    //calculate yc
    Vector<Vector<double> > yc(D, Vector<double>(n)), ycR = yc;
    for(int j = 0; j < D; ++j)
        for(int i = 0; i < n; ++i)
        {
            Vector<double> x = data[n + i].first;
            x[j] = data[i].first[j];
            yc[j][i] = f(x);
        }
    //bootstrap to find standard deviations
    Vector<IncrementalStatistics> s(D);
    for(int k = 0; k < nBoots; ++k)
    {//resample data rows
        for(int i = 0; i < n; ++i)
        {
            int index = GlobalRNG().mod(n);
            yaR[i] = ya[index];
            for(int j = 0; j < D; ++j) ycR[j][i] = yc[j][index];
        }
        //evaluate
        Vector<double> indicesR = findSobolIndicesHelper(yaR, ycR);
        for(int j = 0; j < D; ++j) s[j].addValue(indicesR[j]);
    }
    Vector<double> indices = findSobolIndicesHelper(ya, yc);
    Vector<pair<double, double> > confs;
    double z = getMixedZ(a/D);
    for(int j = 0; j < D; ++j)
    {
        double delta = s[j].stdev() * z;
        confs.append(make_pair(indices[j] - delta, indices[j] + delta));
    }
    return make_pair(indices, confs);
}

}//end namespace
#endif
