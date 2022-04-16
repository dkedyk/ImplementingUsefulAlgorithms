#ifndef IGMDK_LRU_CACHE_H
#define IGMDK_LRU_CACHE_H
#include "../Utils/GCFreeList.h"
#include "../HashTable/LinearProbingHashTable.h"

namespace igmdk{

template<typename KEY, typename VALUE,
    typename HASHER = EHash<BUHash> > class LRUCache
{
    typedef KVPair<KEY, VALUE> ITEM;
    typedef SimpleDoublyLinkedList<ITEM> LIST;
    typedef typename LIST::Iterator I;
    LIST l;
    int size, capacity;
    LinearProbingHashTable<KEY, I, HASHER> h;
public:
    LRUCache(int theCapacity): size(0), capacity(theCapacity)
        {assert(capacity > 0);}
    VALUE* read(KEY const& k)
    {
        I* np = h.find(k);
        if(np)
        {//put in front on access
            l.moveBefore(*np, l.begin());
            return &(*np)->value;
        }
        return 0;
    }
    typedef I Iterator;
    Iterator begin(){return l.begin();}
    Iterator end(){return l.end();}
    Iterator evicteeOnWrite(KEY const& k)//none if not full or item in cache
        {return size < capacity || h.find(k) ? end() : l.rBegin();}
    void write(KEY const& k, VALUE const& v)
    {
        VALUE* oldV = read(k);//first check if already inserted
        if(oldV) *oldV = v;//found, update
        else
        {
            Iterator evictee = evicteeOnWrite(k);
            if(evictee != end())
            {
                h.remove(evictee->key);
                l.moveBefore(evictee, l.begin());//recycle evictee
                evictee->key = k;
                evictee->value = v;
            }
            else
            {
                ++size;
                l.prepend(ITEM(k, v));
            }
            h.insert(k, l.begin());
        }
    }
};

template<typename KEY, typename VALUE, typename RESOURCE,
    typename HASHER = EHash<BUHash> > class DelayedCommitLRUCache
{
    RESOURCE& r;
    typedef LRUCache<KEY, pair<VALUE, bool>, HASHER> LRU;
    typedef typename LRU::Iterator I;
    LRU c;
    void commit(I i)
    {
        if(i->value.second) r.write(i->key, i->value.first);
        i->value.second = false;
    }
    void writeHelper(KEY const& k, VALUE const& v, bool fromWrite)
    {//first commit evictee if any
        I i = c.evicteeOnWrite(k);
        if(i != c.end()) commit(i);
        c.write(k, pair<VALUE, bool>(v, fromWrite));
    }
    DelayedCommitLRUCache(DelayedCommitLRUCache const&);//no copying allowed
    DelayedCommitLRUCache& operator=(DelayedCommitLRUCache const&);
public:
    DelayedCommitLRUCache(RESOURCE& theR, int capacity): r(theR), c(capacity)
        {assert(capacity > 0);}
    VALUE const& read(KEY const& k)
    {//first check if in cache
        pair<VALUE, bool>* mv = c.read(k);
        if(!mv)
        {//if not then read from resource and put in cache
            writeHelper(k, r.read(k), false);
            mv = c.read(k);
        }
        return mv->first;
    }
    void write(KEY const& k, VALUE const& v){writeHelper(k, v, true);}
    void flush(){for(I i = c.begin(); i != c.end(); ++i) commit(i);}
    ~DelayedCommitLRUCache(){flush();}
};

}//end namespace
#endif
