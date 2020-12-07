#ifndef IGMDK_BOOTSTRAP_H
#define IGMDK_BOOTSTRAP_H
#include "Statistics.h"
#include "../Utils/Vector.h"
#include "../Heaps/Heap.h"
#include <cstdlib>
#include <cmath>
namespace igmdk{

template<typename BOOTER> double bootstrapStde(BOOTER& booter, int b = 200)
{
    assert(b > 2);
    IncrementalStatistics s;
    for(int i = 0; i < b; ++i)
    {
        booter.boot();//resample
        s.addValue(booter.eval());
    }//beware - get standard deviation here not standard error
    return s.stdev();
}
template<typename BOOTER> pair<double, double> bootstrapNormalInterval(
    BOOTER& booter, int b = 200, double z = 2)
{
    assert(b > 2 && z > 0);
    double q = booter.eval(), stde = bootstrapStde(booter, b);
    return make_pair(q - z * stde, q + z * stde);
}
double getMixedZ(double a = 0.05, double c = 0.9)
{
    assert(a > 0 && a < 1 && c > 0 && c < 1);
    double z = find2SidedConfZ(1 - a), chebZ = sqrt(1/a);
    return z * c + chebZ * (1 - c);
}
template<typename BOOTER> pair<double, double> bootstrapMixedInterval(
    BOOTER& booter, int b = 200, double a = 0.05)
{//use 10% mix
    assert(b > 2 && a > 0 && a < 1);
    return bootstrapNormalInterval(booter, b, getMixedZ(a));
}

template<typename BOOTER> pair<double, double> bootstrapTIntervalCapped(
    BOOTER& booter, int b = 1000, int b2 = 50, double confidence = 0.95)
{
    assert(b > 2);
    double a = (1 - confidence)/2;
    int tailSize = b * a;
    if(tailSize < 1) tailSize = 1;
    if(tailSize > b/2 - 1) tailSize = b/2 - 1;
    //max heap to work with the largest value
    Heap<double, ReverseComparator<double> > left;
    Heap<double> right;//min heap for the smallest value
    double q = booter.eval();
    IncrementalStatistics s;
    for(int i = 0; i < b; ++i)
    {
        booter.boot();//resample
        double valStat = booter.eval();
        if(isfinite(valStat)) s.addValue(valStat);
        else continue;
        BOOTER booter2 = booter.cloneForNested();
        double stdeInner = bootstrapStde(booter2, b2),
            value = (valStat - q)/stdeInner;
        if(isnan(value)) continue; //possible division by 0 OK due to cap
        if(left.getSize() < tailSize)
        {//heaps not full
            left.insert(value);
            right.insert(value);
        }
        else
        {//heaps full - replace if more extreme
            if(value < left.getMin()) left.changeKey(0, value);
            else if(value > right.getMin()) right.changeKey(0, value);
        }
    }//beware - get standard deviation here not standard error
    //protect against bad data
    double normalZ = 1, chebZ = sqrt(1/(2 * a)), stde = s.n > 2 ? s.stdev() :
        numeric_limits<double>::infinity(), leftZ = right.getSize() > 0 ?
        min(max(right.getMin(), normalZ), chebZ) : chebZ, rightZ =
        left.getSize() > 0 ? min(max(-left.getMin(), normalZ), chebZ) : chebZ;
    return make_pair(q - stde * leftZ, q + stde * rightZ);
}

class BTPCFunctor
{
    Vector<Vector<double> > sets;
    bool findLeft;
    double q, aTarget;
    pair<double, double> makeInterval(int i, double a)const
    {//percentile
        assert(0 <= i && i < sets.getSize());
        double b = sets[i].getSize();
        int left = max<int>(0, a/2 * b),
        right = min<int>(b - 1, (1 - a/2) * b);
        return make_pair(sets[i][left], sets[i][right]);
    }
public:
    BTPCFunctor(double theQ, double a): q(theQ), aTarget(a), findLeft(true){}
    void addSet(){sets.append(Vector<double>());}
    void addValue(double value){sets.lastItem().append(value);}
    void closeSet()
        {quickSort(sets.lastItem().getArray(), sets.lastItem().getSize());}
    void flipSide(){findLeft = !findLeft;}
    double operator()(double a)const
    {//returns < 0 if under target
        int missCount = 0, b = sets.getSize();
        for(int i = 0; i < b; ++i)
        {//make sure that don't confuse miss left and miss right
            pair<double, double> conf = makeInterval(i, a);
            missCount += (findLeft ? q < conf.first : conf.second < q);
        }
        return aTarget/2 - missCount * 1.0/b;
    }
    double findA()const
    {
        double left = 0, right = 1, yLeft = this->operator()(left);
        return haveDifferentSign(yLeft, this->operator()(right)) ?
            solveFor0(*this, left, right).first :
            yLeft < 0 ? left : right;
    }
};
template<typename BOOTER> pair<double, double>
    bootstrapDoublePercentileInterval(BOOTER& booter, int b = 1020,
    int b2 = 98, double confidence = 0.95)
{
    assert(b > 2 && confidence > 0 && confidence < 1);
    double q = booter.eval();
    BTPCFunctor f(q, 1 - confidence);
    Vector<double> values(b);
    for(int i = 0; i < b; ++i)
    {
        booter.boot();//resample
        values[i] = booter.eval();
        BOOTER booter2 = booter.cloneForNested();
        f.addSet();
        for(int j = 0; j < b; ++j)
        {
            booter2.boot();//resample nested
            f.addValue(booter2.eval());
        }
        f.closeSet();
    }
    quickSort(values.getArray(), b);
    int left = max<int>(0, f.findA()/2 * b);
    f.flipSide();//default is left, now find right
    int right = min<int>(b - 1, (1 - f.findA()/2) * b);
    //cap at q for containment
    return make_pair(min(q, values[left]), max(q, values[right]));
}

template<typename FUNCTION, typename DATA = double> struct BasicBooter
{
    Vector<DATA> const& data;
    Vector<DATA> resample;
    FUNCTION f;
    BasicBooter(Vector<DATA> const& theData, FUNCTION const& theF = FUNCTION())
        :data(theData), f(theF), resample(theData){assert(data.getSize() > 0);}
    void boot()
    {
        for(int i = 0; i < data.getSize(); ++i)
            resample[i] = data[GlobalRNG().mod(data.getSize())];
    }
    double eval()const{return f(resample);}
    //for the bootstrap-t interval
    BasicBooter cloneForNested()const{return BasicBooter(resample, f);}
};

template<typename FUNCTION, typename DATA = double> struct MultisampleBooter
{
    Vector<Vector<DATA> const* > const& data;
    Vector<Vector<DATA> > resample;
    FUNCTION f;
    MultisampleBooter(FUNCTION const& theF = FUNCTION()): f(theF) {}
    void addSample(Vector<DATA> const& theData)
    {
        data.append(&theData);
        resample.append(theData);
    }
    void boot()
    {
        for(int i = 0; i < data.getSize(); ++i)
            for(int j = 0; j < data[i].getSize(); ++j)
                resample[i][j] = data[i][GlobalRNG().mod(data[i].getSize())];
    }
    double eval()const
    {
        assert(data.getSize() > 0);
        return f(resample);
    }
    MultisampleBooter cloneForNested()const
    {//for the bootstrap-t interval
        MultisampleBooter clone(f);
        for(int i = 0; i < data.getSize(); ++i)
            clone.addSample(*data[i]);
        return clone;
    }
};

}//end namespace
#endif
