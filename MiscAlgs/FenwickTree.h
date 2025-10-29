#ifndef IGMDK_FENWICK_TREE_H
#define IGMDK_FENWICK_TREE_H
#include "../Utils/Vector.h"
namespace igmdk{

template<typename ITEM> class FenwickTree
{
    Vector<ITEM> nodes;
public:
    explicit FenwickTree(int n): nodes(n, 0){}
    explicit FenwickTree(std::initializer_list<ITEM> args):
        nodes(args.size(), 0)
        {for(int i = 0; i < args.size(); ++i) addValue(i, args[i]);}
    int getSize()const{return nodes.getSize();}
    ITEM getCumulativeValue(int index)const
    {
        assert(index >= 0 && index <= getSize());
        ITEM sum = nodes[index];
        for(;index > 0;)
        {
            ++index;//convert to 1-based from 0-based
            index = index & (index - 1);
            --index;//convert back
            if(index < 0) break;
            sum += nodes[index];
        }
        return sum;
    }
    ITEM getValue(int index)const
    {
        assert(index >= 0 && index <= getSize());
        return getCumulativeValue(index) -
            (index > 0 ? getCumulativeValue(index - 1) : 0);
    }
    void addValue(int index, ITEM value)
    {
        assert(index >= 0 && index <= getSize());
        nodes[index] += value;
        for(;;)
        {
            ++index;//convert to 1-based from 0-based
            index += index & -index;
            --index;//convert back
            if(index >= getSize()) break;
            nodes[index] += value;
        }
    }
};

}//end namespace
#endif
