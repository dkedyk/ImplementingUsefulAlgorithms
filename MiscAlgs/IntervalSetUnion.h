#ifndef IGMDK_SCRAP_UNION_FIND_H
#define IGMDK_SCRAP_UNION_FIND_H

#include "../Utils/Utils.h"
#include "../Utils/Vector.h"
#include "../RandomTreap/Treap.h"

namespace igmdk{

class IntervalSetUnion
{
    Treap<int, bool> treap;
public:
    int find(int i){return treap.getSize() == 0 ? -1: treap.successor(i)->key;}
    void merge(int i){return treap.remove(i);}
    void split(int i){treap.insert(i, 0);}
};

}
#endif
