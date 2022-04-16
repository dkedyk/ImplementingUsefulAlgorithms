#include "ConvexHull.h"
#include "../Utils/Debug.h"
#include <cmath>
using namespace igmdk;

void testConvexHull()
{
    Vector<Point2> points, result;
    int N = 5;
    for(int i  = 0; i < N; ++i)
    {
        points.append(Point2(GlobalRNG().uniform01(), GlobalRNG().uniform01()));
        DEBUG(points[i][0]);
        DEBUG(points[i][1]);
    }
    result = convexHull(points);
    DEBUG(result.getSize());
    for(int i = 0; i < result.getSize(); ++i)
    {
        DEBUG(result[i][0]);
        DEBUG(result[i][1]);
    }
}

int main()
{
    testConvexHull();
	return 0;
}
