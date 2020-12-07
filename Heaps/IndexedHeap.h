#ifndef IGMDK_INDEXED_HEAP_H
#define IGMDK_INDEXED_HEAP_H
#include "Heap.h"
#include "../HashTable/ChainingHashTable.h"
namespace igmdk{

template<typename ITEM, typename COMPARATOR = DefaultComparator<ITEM>,
    typename HANDLE = int, typename HASHER = EHash<BUHash> > class IndexedHeap
{
    typedef ChainingHashTable<HANDLE, int, HASHER> MAP;
    MAP map;
    typedef typename MAP::NodeType* POINTER;
    typedef pair<ITEM, POINTER> Item;
    typedef PairFirstComparator<ITEM, POINTER, COMPARATOR> Comparator;
    struct Reporter
        {void operator()(Item& item, int i){item.second->value = i;}};
    Heap<Item, Comparator, Reporter> h;
public:
    IndexedHeap(COMPARATOR const& theC = COMPARATOR()): h(Comparator(theC)){}
    int getSize()const{return h.getSize();}
    ITEM const* find(HANDLE handle)
    {
        int* index = map.find(handle);
        return index ? &h[*index].first : 0;
    }
    bool isEmpty()const{return h.isEmpty();}
    void insert(ITEM const& item, HANDLE handle)
    {
        assert(!find(handle));//else map will fail with duplicate
        h.insert(Item(item, map.insert(handle, h.getSize())));
    }
    pair<ITEM, HANDLE> getMin()const
    {
        Item temp = h.getMin();
        return make_pair(temp.first, temp.second->key);
    }
    pair<ITEM, HANDLE> deleteMin()
    {
        Item temp = h.deleteMin();
        pair<ITEM, HANDLE> result = make_pair(temp.first, temp.second->key);
        map.remove(temp.second->key);
        return result;
    }
    void changeKey(ITEM const& item, HANDLE handle)
    {
        POINTER p = map.findNode(handle);
        if(p) h.changeKey(p->value, Item(item, p));
        else insert(item, handle);
    }
    void deleteKey(HANDLE handle)
    {
        int* index = map.find(handle);
        assert(index);
        h.remove(*index);
        map.remove(handle);
    }
};

template<typename ITEM, typename COMPARATOR = DefaultComparator<ITEM> >
class IndexedArrayHeap
{
    Vector<int> map;
    typedef pair<ITEM, int> Item;
    typedef PairFirstComparator<ITEM, int, COMPARATOR> Comparator;
    struct Reporter
    {
        Vector<int>& pmap;
        Reporter(Vector<int>& theMap): pmap(theMap) {}
        void operator()(Item& item, int i){pmap[item.second] = i;}
    };
    Heap<Item, Comparator, Reporter> h;
public:
    typedef Item ITEM_TYPE;
    IndexedArrayHeap(COMPARATOR const& theC = COMPARATOR()):
        h(Comparator(theC), Reporter(map)) {}
    int getSize()const{return h.getSize();}
    ITEM const* find(int handle)
    {
        assert(handle >= 0);
        return handle >= map.getSize() || map[handle] == -1 ? 0 :
            &h[map[handle]].first;
    }
    bool isEmpty()const{return h.isEmpty();}
    void insert(ITEM const& item, int handle)
    {
        assert(handle >= 0);
        if(handle >= map.getSize())
            for(int i = map.getSize(); i <= handle; ++i) map.append(-1);
        h.insert(Item(item, handle));
    }
    pair<ITEM, int> const& getMin()const{return h.getMin();}
    pair<ITEM, int> deleteMin()
    {
        Item result = h.deleteMin();
        map[result.second] = -1;
        return result;
    }
    void changeKey(ITEM const& item, int handle)
    {
        assert(handle >= 0);
        if(handle >= map.getSize() || map[handle] == -1) insert(item, handle);
        else h.changeKey(map[handle], Item(item, handle));
    }
    void deleteKey(int handle)
    {
        assert(handle >= 0 && handle < map.getSize());
        int pointer = map[handle];
        assert(pointer != -1);
        h.remove(pointer);
        map[handle] = -1;
    }
};

}//end namespace
#endif
