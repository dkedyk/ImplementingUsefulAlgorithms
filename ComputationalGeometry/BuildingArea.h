#ifndef IGMDK_BUILDING_AREA_H
#define IGMDK_BUILDING_AREA_H

#include <map>
#include <vector>
#include <algorithm>//sort
#include <cassert>
namespace igmdk{

using namespace std;

struct RoofCorner
{
    double x, y;
    bool isLeft;
    RoofCorner(double theX, double theY, bool theIsLeft): x(theX), y(theY),
        isLeft(theIsLeft){}//comparison on x + give priority to left over right
    bool operator<(RoofCorner const& rhs)const
        {return x == rhs.x ? isLeft > rhs.isLeft : x < rhs.x;}
};
double buildingArea(vector<RoofCorner> corners)
{
    sort(corners.begin(), corners.end());
    double result = 0, currentHeight = 0, lastX = corners[0].x;
    map<double, int> openBuildings;//indexed and sorted by height, note that
    //don't have numerical issues with double key as no calculations are done
    for(unsigned int i = 0; i < corners.size(); ++i)
    {
        double x = corners[i].x, y = corners[i].y;
        result += currentHeight * (x - lastX);
        lastX = x;
        //manage count of open buildings
        openBuildings[y] += corners[i].isLeft ? 1 : -1;
        if(openBuildings[y] == 0) openBuildings.erase(y);
        //current height is that of tallest open building
        currentHeight = openBuildings.size() > 0 ?
            openBuildings.rbegin()->first : 0;
    }
    return result;
}
void testBuildingArea()
{
    vector<RoofCorner> points;
    points.push_back(RoofCorner(0, 1, true));
    points.push_back(RoofCorner(2, 1, false));
    points.push_back(RoofCorner(1, 2, true));
    points.push_back(RoofCorner(3, 2, false));
    assert(buildingArea(points) == 5);
}

}//end namespace
#endif
