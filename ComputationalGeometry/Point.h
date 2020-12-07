#ifndef IGMDK_POINT_H
#define IGMDK_POINT_H
#include "../Utils/Utils.h"
#include <cmath>
#include <cstdlib>//int version of abs
namespace igmdk{

template<typename KEY, int D = 2>
class Point: public ArithmeticType<Point<KEY, D> >
{
    KEY x[D];
public:
    static int const d = D;
    KEY& operator[](int i){assert(i >= 0 && i < D); return x[i];}
    KEY const& operator[](int i)const{assert(i >= 0 && i < D); return x[i];}
    int getSize()const{return D;}
    Point(){for(int i = 0; i < D; ++i) x[i] = 0;}
    Point(KEY const& x0, KEY const& x1)
    {
        assert(D == 2);//to prevent accidents for D > 2
        x[0] = x0;
        x[1] = x1;
    }
    bool operator==(Point const& rhs)const
    {
        for(int i = 0; i < D; ++i) if(x[i] != rhs.x[i]) return false;
        return true;
    }
    Point& operator+=(Point const& rhs)
    {
        for(int i = 0; i < D; ++i) x[i] += rhs.x[i];
        return *this;
    }
    Point& operator*=(double scalar)
    {
        for(int i = 0; i < D; ++i) x[i] *= scalar;
        return *this;
    }
    friend Point operator*(Point const& point, double scalar)
    {
        Point result = point;
        return result *= scalar;
    }
    Point& operator-=(Point const& rhs){return *this += rhs * -1;}
    Point operator-(){return *this * -1;}
    double friend dotProduct(Point const& a, Point const& b)
    {
        double dp = 0;
        for(int i = 0; i < D; ++i) dp += a[i] * b[i];
        return dp;
    }
};
typedef Point<double> Point2;

template<typename VECTOR> class EuclideanDistance
{
    static double iDistanceIncremental(VECTOR const& lhs, VECTOR const& rhs,
        int i)//add on a component
    {
        double x = lhs[i] - rhs[i];
        return x * x;
    }
    static double distanceIncremental(VECTOR const& lhs, VECTOR const& rhs,
        double bound = numeric_limits<double>::infinity())
    {//compute distance up to a bound
        assert(lhs.getSize() == rhs.getSize());
        double sum = 0;
        for(int i = 0; i < lhs.getSize() && sum < bound; ++i)
            sum += iDistanceIncremental(lhs, rhs, i);
        return sum;
    }
public:
    struct Distance
    {//metric functor
        double operator()(VECTOR const& lhs, VECTOR const& rhs)const
            {return sqrt(distanceIncremental(lhs, rhs));}
    };
    struct DistanceIncremental
    {//incremental functor that returns distance squared
        double operator()(VECTOR const& lhs, VECTOR const& rhs)const
            {return distanceIncremental(lhs, rhs);}
        double operator()(VECTOR const& lhs, VECTOR const& rhs, int i)const
            {return iDistanceIncremental(lhs, rhs, i);}
        double operator()(double bound, VECTOR const& lhs, VECTOR const& rhs)
            const{return distanceIncremental(lhs, rhs, bound);}
    };
};

}//end namespace
#endif
