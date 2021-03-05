#ifndef IGMDK_SOBOL_H
#define IGMDK_SOBOL_H

#include "../Utils/Bits.h"
#include "../Utils/Vector.h"
#include "Statistics.h"
#include <cstdlib>
#include <cmath>

namespace igmdk{

class Sobol
{//SobolPolys do not represent highest and lowest 1s
    enum{B = numeric_limits<double>::digits};
    unsigned long long k;//current simulation number
    Vector<unsigned long long> x, v;//current sample and precomputed data
    double factor;
    int vIndex(int d, int b){return d * B + b;}//to map to array index
public:
    static int maxD(){return 52;}
    Sobol(int d): factor(1.0/twoPower(B)), v(d * B, 0), x(d, 0), k(1)
    {
        assert(d <= maxD());
        unsigned char const SobolPolys[] = {0,1,1,2,1,4,2,4,7,11,13,14,1,13,16,
            19,22,25,1,4,7,8,14,19,21,28,31,32,37,41,42,50,55,56,59,62,14,21,
            22,38,47,49,50,52,56,67,70,84,97,103,115,122},
            SobolDegs[] = {1,2,3,3,4,4,5,5,5,5,5,5,6,6,6,6,6,6,7,7,7,7,7,7,7,7,
            7,7,7,7,7,7,7,7,7,7,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8};
        for(int i = 0; i < d; ++i)
            for(int j = 0; j < B; ++j)
            {
                unsigned long long value;
                int l = j - SobolDegs[i];
                //j < D case
                if(l < 0) value = (2 * j + 1) * twoPower(B - j - 1);
                else
                {//j >= D case
                    value = v[vIndex(i, l)];
                    value ^= value/twoPower(SobolDegs[i]);
                    for(int k = 1; k < SobolDegs[i]; ++k)
                        if(Bits::get(SobolPolys[i], k - 1))
                            value ^= v[vIndex(i, l + k)];
                }
                v[vIndex(i, j)] = value;
            }
        next();//generate first value
    }
    void next()
    {
        for(int i = 0, c = rightmost0Count(k++); i < x.getSize(); ++i)
            x[i] ^= v[vIndex(i, c)];
    }
    double getU01Value(int i)const{return x[i] * factor;}
    double getUValue(int i, double a, double b)const
        {return a + (b - a) * getU01Value(i);}
};

class ScrambledSobolHybrid
{
    int D;
    Vector<double> shifts;
    mutable Sobol s;
    Vector<pair<double, double> > box;
public:
    ScrambledSobolHybrid(Vector<pair<double, double> > const& theBox):
        D(theBox.getSize()), shifts(D), s(min(D, Sobol::maxD())), box(theBox)
    {
        for(int i = 0; i < D; ++i)
            shifts[i] = GlobalRNG().uniform(box[i].first, box[i].second);
    }
    Vector<double> operator()()const
    {//first get Sobol variates for the supported dimensions
        Vector<double> next(D);
        for(int i = 0; i < min(D, Sobol::maxD()); ++i)
            next[i] = s.getUValue(i, box[i].first, box[i].second);
        s.next();
        //random for remaining dimensions;
        for(int i = min(D, Sobol::maxD()); i < D; ++i)
            next[i] = GlobalRNG().uniform(box[i].first, box[i].second);
        //scramble
        for(int i = 0; i < D; ++i)
        {
            next[i] += shifts[i] - box[i].first;
            if(next[i] > box[i].second)
                next[i] -= box[i].second - box[i].first;
        }
        return next;
    }
};
template<typename SAMPLER = ScrambledSobolHybrid> class GeometricBoxWrapper
{
    static Vector<pair<double, double> > transformGeometricBox
        (Vector<pair<double, double> > box)
    {
        for(int i = 0; i < box.getSize(); ++i)
        {
            assert(box[i].first > 0 && box[i].first < box[i].second);
            box[i].first = log(box[i].first);
            box[i].second = log(box[i].second);
        }
        return box;
    }
    SAMPLER s;
public:
    GeometricBoxWrapper(Vector<pair<double, double> > const& box):
        s(transformGeometricBox(box)){}
    Vector<double> operator()()const
    {
        Vector<double> result = s();
        for(int i = 0; i < result.getSize(); ++i) result[i] = exp(result[i]);
        return result;
    }
};

template<typename TEST, typename FUNCTION> pair<double, double> SobolIntegrate(
    Vector<pair<double, double> > const& box, int n,
    TEST const& isInside = TEST(), FUNCTION const& f = FUNCTION())
{
    IncrementalStatistics s;
    ScrambledSobolHybrid so(box);
    for(int i = 0; i < n; ++i)
    {
        Vector<double> point = so();
        if(isInside(point)) s.addValue(f(point));
    }
    double regionVolume = boxVolume(box) * s.n/n;
    return make_pair(regionVolume * s.getMean(),
        regionVolume * s.getStandardErrorSummary().error95());
}

}//end namespace
#endif
