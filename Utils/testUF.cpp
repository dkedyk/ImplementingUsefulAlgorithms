#include "UnionFind.h"
#include <iostream>
#include "DEBUG.h"
#include "../RandomTreap/Treap.h"//for set union
using namespace igmdk;

//NOT MENTIONED IN TEXT
class IntervalSetUnion
{
    Treap<int, bool> treap;
public:
    int find(int i)
        {return treap.getSize() == 0 ? -1 : treap.successor(i)->key;}
    void merge(int i){return treap.remove(i);}
    void split(int i){treap.insert(i, 0);}
};
int main()
{
	IntervalSetUnion iu;

	iu.split(5);
	DEBUG(iu.find(2));
	iu.split(3);
	DEBUG(iu.find(2));
	iu.merge(3);
	DEBUG(iu.find(2));
	return 0;
}
