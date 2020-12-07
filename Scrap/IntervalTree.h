#ifndef IGMDK_INTERVAL_TREE_H
#define IGMDK_INTERVAL_TREE_H

#include "../Sorting/Sort.h"
#include "../ComputationalGeometry/Point.h"

namespace igmdk{

class IntervalTree
{
    Vector<Point2 >& points;
    struct Node{int median, point, left, right;};
    int root;
    Vector<Node> nodes;
    int findRemoveMin(Vector<int>& indices)
    {
        int result = 0;
        for(int i = 1; i < indices.getSize(); ++i)
            if(points[indices[i]][0] < points[indices[result]][0]) result = i;
        int point = indices[result];
        //remove at result and keep order
        for(int j = result; j < indices.getSize() - 1; ++j)
            indices[j] = indices[j + 1];
        indices.removeLast();
        return point;
    }
    void construct(int node, Vector<int>& indices)
    {
        if(indices.getSize() == 1)
        {
            nodes[node].left = nodes[node].right = -1;
            nodes[node].point = indices[0];
            return;
        }
        nodes[node].point = findRemoveMin(indices);
        nodes[node].median = indices[indices.getSize()/2];
        Vector<int> left, right;
        for(int i = 0; i < indices.getSize(); ++i)
        {
            if(i <= indices.getSize()/2)
                left.append(indices[i]);
            else right.append(indices[i]);
        }
        nodes[node].left = nodes.getSize();
        nodes.append(Node());
        construct(nodes[node].left, left);
        if(right.getSize() > 0)
        {
            nodes[node].right = nodes.getSize();
            nodes.append(Node());
            construct(nodes[node].right, right);
        }
        else nodes[node].right = -1;
    }
    struct DComparator
    {
        Vector<Point2>& points;
        DComparator(Vector<Point2>& thePoints): points(thePoints){}
        bool operator()(int lhs, int rhs)const
            {return points[lhs][1] < points[rhs][1];}
        bool isEqual(int lhs, int rhs)const
            {return points[lhs][1] == points[rhs][1];}
    };
public:
    IntervalTree(Vector<Point2>& thePoints):points(thePoints), root(0)
    {
        assert(points.getSize() > 0);
        nodes.append(Node());
        Vector<int> indices;
        for(int i = 0; i < points.getSize(); ++i) indices.append(i);
        quickSort(indices.getArray(), 0, indices.getSize()-1,
            DComparator(points));
        construct(root, indices);
    }
    void containingIntervals(int x, Vector<int>& result, int node = 0)
    {
        if(points[nodes[node].point][0] <= x)
        {
            if(points[nodes[node].point][1] >= x)
                result.append(nodes[node].point);
            if(nodes[node].left != -1)
            {
                if(x <= points[nodes[node].median][1])
                    containingIntervals(x, result, nodes[node].left);
                if(nodes[node].right != -1)
                    containingIntervals(x, result, nodes[node].right);
            }
        }
    }
};

}
#endif
