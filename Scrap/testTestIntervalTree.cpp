#include <iostream>
#include <cmath>
#include "IntervalTree.h"
using namespace igmdk;

void testIntervalTree()
{
    int N = 15000;
    Vector<Point2> points;
    for(int i = 0; i < N; ++i)
    {
        double p1 = (GlobalRNG().next() % 1000), p2 = i;
        points.append(Point2(min(p1, p2), max(p1, p2)));
        //DEBUG(min(p1, p2));
        //DEBUG(max(p1, p2));
    }
    DEBUG(clock());



    IntervalTree it(points);
    for(int k = 0; k < 10000; ++k)
    {
        Vector<int> result;
        int point = (GlobalRNG().next() % 1000);//GlobalRNG.next() % 10;
        //DEBUG(point);
        it.containingIntervals(point, result);
        for(int i = 0; i < result.getSize(); ++i)
        {
            //DEBUG(result[i].key.x[0]);
            //DEBUG(result[i].key.x[1]);
            //DEBUG(result[i]);
        }
        //DEBUG(result.getSize());
    }
    DEBUG(clock());
}

int main()
{
    testIntervalTree();
	return 0;
}
