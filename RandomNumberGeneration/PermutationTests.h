#ifndef PERMUTATION_TESTS_H
#define PERMUTATION_TESTS_H

#include "Random.h"
#include "../Utils/Vector.h"
#include "../NumericalMethods/EquationSolving.h"

namespace igmdk{

template<typename PERM_TESTER> double permutationTest(PERM_TESTER& p,
    int b = 10000, int side = 0, Random<>& rng = GlobalRNG())
{//"more extreme" means larger for 1-sided, side -1 for smaller, 1 for larger
    assert(b > 1 && abs(side) <= 1);
    double f = p();
    int nLeft = 1, nRight = 1;//start with 1 not 0
    for(int i = 0; i < b; ++i)
    {
        p.exchange(rng);
        double fr = p();
        nLeft += (f <= fr);
        nRight += (f >= fr);
    }
    double leftP = nLeft * 1.0/(b + 1), rightP = nRight * 1.0/(b + 1);
    return side == 0 ? min(1.0, 2 * min(leftP, rightP)) :
        side == -1 ? leftP : rightP;
}

template<typename LOCATION_F> struct PairedTestLocationPermuter
{
    Vector<double> diffs;
    LOCATION_F f;
    PairedTestLocationPermuter(Vector<double> const& data): diffs(data){}
    PairedTestLocationPermuter(Vector<pair<double, double> > const& data)
    {
        for(int i = 0; i < data.getSize(); ++i)
            diffs.append(data[i].first - data[i].second);
    }
    double operator()()const{return f(diffs);}
    void exchange(Random<>& r)
        {for(int i = 0; i < diffs.getSize(); ++i)if(r.mod(2))diffs[i] *= -1;}
    void setShift(double shift)//must not be called after exchange
        {for(int i = 0; i < diffs.getSize(); ++i) diffs[i] += shift;}
};

template<typename PERM_TESTER> struct permConfHelper
{
    PERM_TESTER const& p;
    double a;
    int seed, b;
    permConfHelper(PERM_TESTER const& theP, double theA, int theSeed,
        int theB): p(theP), a(theA), seed(theSeed), b(theB) {}
    double operator()(double shift)const
    {
        Random<> rng(seed);
        PERM_TESTER p2 = p;
        p2.setShift(shift);
        return permutationTest(p2, b, 0, rng) - a;
    }
};
template<typename PERM_TESTER> pair<double, double> permutationConf(
    PERM_TESTER& p, double a = 0.05, int b = 10000)//assume two-sided
{
    double stat = p();
    permConfHelper<PERM_TESTER> f(p, a, time(0), b);
    double left = exponentialSearch1Sided(f, -stat, -0.001).first,
        right = exponentialSearch1Sided(f, -stat).first;
    left = isfinite(left) ? left : -stat;
    right = isfinite(right) ? right : -stat;
    return make_pair(2 * stat + left, 2 * stat + right);
}

template<typename LOCATION_F> pair<double, double> permutationLocationConf(
    Vector<double> const& data, double a = 0.05, int b = 10000)
{
    PairedTestLocationPermuter<LOCATION_F> p = {data};
    return permutationConf(p, a, b);
}

}//end namespace
#endif
