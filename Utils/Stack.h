#ifndef IGMDK_STACK_H
#define IGMDK_STACK_H
#include "Utils.h"
#include "Vector.h"
namespace igmdk{

template<typename ITEM, typename VECTOR = Vector<ITEM> > struct Stack
{
    VECTOR storage;
    void push(ITEM const& item){storage.append(item);}
    ITEM pop()
    {
        assert(!isEmpty());
        ITEM result = storage.lastItem();
        storage.removeLast();
        return result;
    }
    ITEM& getTop()
    {
        assert(!isEmpty());
        return storage.lastItem();
    }
    bool isEmpty(){return !storage.getSize();}
};

}//end namespace
#endif
