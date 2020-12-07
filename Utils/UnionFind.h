#ifndef IGMDK_UNION_FIND_H
#define IGMDK_UNION_FIND_H

#include "../Utils/Utils.h"
#include "../Utils/Vector.h"

namespace igmdk{

class UnionFind
{
    mutable Vector<int> parent;//parent or negated size of the tree
public:
    UnionFind(int size): parent(size, -1){}
    bool isRoot(int n)const{return parent[n] < 0;};
    int find(int n)const
        {return isRoot(n) ? n : (parent[n] = find(parent[n]));}
    void join(int i, int j)
    {
        int parentI = find(i), parentJ = find(j);
        if(parentI != parentJ)
        {//parent[parentI] and parent[parentJ] are negative sizes
            if(parent[parentI] > parent[parentJ]) swap(parentI, parentJ);
            parent[parentI] += parent[parentJ];
            parent[parentJ] = parentI;
        }
    }
    bool areEquivalent(int i, int j)const{return find(i) == find(j);}
    int subsetSize(int i)const{return -parent[find(i)];}
    void addSubset(){parent.append(-1);}
};

}//end namespace
#endif
