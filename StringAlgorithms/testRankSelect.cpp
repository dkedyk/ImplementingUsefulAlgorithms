#include <cassert>
#include <iostream>
#include "../Utils/Debug.h"
#include "RankSelect.h"
using namespace std;
using namespace igmdk;

int main()
{
    unsigned long long x = 1063;
    DEBUG(popCountWord(x));
    for(int i = 0; i < 64; ++i)
    {
        DEBUG(rank64(x, i));
        int s = select64(x, i);
        DEBUG(s);
        if(s != -1) assert(rank64(x, s) == i);

    }
    RankSelect rs;
    rs.getBitset().appendValue(x, 64);
    rs.getBitset().appendValue(x, 64);
    rs.getBitset().debug();
    rs.finalize();
    DEBUG(rs.rank(0));
    DEBUG(rs.rank(1));
    DEBUG(rs.rank(2));
    DEBUG(rs.rank(3));
    DEBUG(rs.rank(4));
    DEBUG(rs.rank0(0));
    DEBUG(rs.rank0(1));
    DEBUG(rs.rank0(2));
    DEBUG(rs.rank0(3));
    DEBUG(rs.rank0(4));
    DEBUG(rs.select(0));
    DEBUG(rs.select(1));
    DEBUG(rs.select(2));
    DEBUG(rs.select(3));
    DEBUG(rs.select(4));
    DEBUG(rs.select(5));
    DEBUG(rs.select(6));
    DEBUG(rs.select0(0));
    DEBUG(rs.select0(1));
    DEBUG(rs.select0(2));
    DEBUG(rs.select0(3));
    DEBUG(rs.select0(4));
    DEBUG(rs.select0(5));
    DEBUG(rs.select0(64));


    BinaryTree bt;
    bt.addNodeInLevelOrder(true);
    bt.addNodeInLevelOrder(true);
    bt.addNodeInLevelOrder(true);
    bt.addNodeInLevelOrder(true);
    bt.addNodeInLevelOrder(false);
    bt.addNodeInLevelOrder(false);
    bt.addNodeInLevelOrder(true);
    bt.addNodeInLevelOrder(true);
    bt.addNodeInLevelOrder(true);
    bt.addNodeInLevelOrder(true);
    bt.addNodeInLevelOrder(false);
    bt.addNodeInLevelOrder(false);
    bt.addNodeInLevelOrder(false);
    bt.addNodeInLevelOrder(true);
    bt.addNodeInLevelOrder(false);
    bt.addNodeInLevelOrder(false);
    bt.addNodeInLevelOrder(false);
    bt.addNodeInLevelOrder(false);
    bt.addNodeInLevelOrder(false);
    bt.finalize();
    DEBUG(bt.parent(0));
    DEBUG(bt.parent(1));
    DEBUG(bt.parent(2));
    DEBUG(bt.parent(3));
    DEBUG(bt.parent(4));
    DEBUG(bt.parent(5));
    DEBUG(bt.parent(6));
    DEBUG(bt.parent(7));
    DEBUG(bt.parent(8));
    //DEBUG(bt.parent(9));
    DEBUG(bt.leftChild(0));
    DEBUG(bt.rightChild(0));
    DEBUG(bt.leftChild(1));
    DEBUG(bt.rightChild(1));
    DEBUG(bt.leftChild(2));
    DEBUG(bt.rightChild(2));

    OrdinalTree ot;
    ot.addNodeInLevelOrder(4);
    ot.addNodeInLevelOrder(0);
    ot.addNodeInLevelOrder(2);
    ot.addNodeInLevelOrder(0);
    ot.addNodeInLevelOrder(1);
    ot.addNodeInLevelOrder(0);
    ot.addNodeInLevelOrder(1);
    ot.addNodeInLevelOrder(0);
    ot.addNodeInLevelOrder(0);
    ot.finalize();
    DEBUG(ot.parent(0));
    DEBUG(ot.parent(1));
    DEBUG(ot.parent(2));
    DEBUG(ot.parent(3));
    DEBUG(ot.parent(4));
    DEBUG(ot.parent(5));
    DEBUG(ot.parent(6));
    DEBUG(ot.parent(7));
    DEBUG(ot.parent(8));

    DEBUG(ot.firstChild(0));
    DEBUG(ot.firstChild(1));
    DEBUG(ot.firstChild(2));
    DEBUG(ot.firstChild(3));
    DEBUG(ot.firstChild(4));
    DEBUG(ot.firstChild(5));
    DEBUG(ot.firstChild(6));
    DEBUG(ot.firstChild(7));
    DEBUG(ot.firstChild(8));

    DEBUG(ot.nextChild(0));
    DEBUG(ot.nextChild(1));
    DEBUG(ot.nextChild(2));
    DEBUG(ot.nextChild(3));
    DEBUG(ot.nextChild(4));
    DEBUG(ot.nextChild(5));
    DEBUG(ot.nextChild(6));
    DEBUG(ot.nextChild(7));
    DEBUG(ot.nextChild(8));

	return 0;
}
