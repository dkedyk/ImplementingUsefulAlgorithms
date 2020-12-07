#ifndef IGMDK_CONVEXHULL_H
#define IGMDK_CONVEXHULL_H
#include "Point.h"
#include "../Utils/Utils.h"
#include "../Utils/Vector.h"
#include "../Sorting/Sort.h"
#include "../LargeNumbers/LargeRational.h"
namespace igmdk{

double triangleArea(Point2 const& a, Point2 const& b, Point2 const& c)
    {return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);}
bool ccw(Point2 const& a, Point2 const& b, Point2 const& c)
    {return triangleArea(a, b, c) >= 0;}//true if the points turn left

Rational robustTriangleArea(Point2 const&a, Point2 const&b, Point2 const&c)
{
    return (Rational(b[0]) - Rational(a[0])) * (Rational(c[1]) -
        Rational(a[1])) - (Rational(b[1]) - Rational(a[1])) * (Rational(c[0])
        - Rational(a[0]));
}
bool robustCcw(Point2 const& a, Point2 const& b, Point2 const& c)
    {return !robustTriangleArea(a, b, c).isMinus();}

void processPoint(Vector<Point2>& hull, Point2 const& point)
{
    hull.append(point);
    while(hull.getSize() > 2 && ccw(hull[hull.getSize() - 3],
        hull[hull.getSize() - 2], hull[hull.getSize() - 1]))
    {
        hull[hull.getSize() - 2] = hull[hull.getSize() - 1];
        hull.removeLast();
    }
}
Vector<Point2> convexHull(Vector<Point2>& points)
{
    assert(points.getSize() > 2);
    quickSort(points.getArray(), 0, points.getSize() - 1,
        LexicographicComparator<Point2>());
    //upper hull
    Vector<Point2> result;
    result.append(points[0]);//initialize with the first two points
    result.append(points[1]);
    for(int i = 2; i < points.getSize(); ++i) processPoint(result, points[i]);
    //lower hull, remove leftmost point which is added twice
    for(int i = points.getSize() - 2; i >= 0; --i)
        processPoint(result, points[i]);
    result.removeLast();
    return result;
}

}//end namespace
#endif
