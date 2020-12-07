#ifndef IGMDK_EMFREELIST_H
#define IGMDK_EMFREELIST_H

#include "EMVector.h"

namespace igmdk{

template<typename POD, typename SERIALIZER = CastSerializer<POD> >
class EMFreelist
{
    EMVector<POD, SERIALIZER> nodes;
    EMVector<long long, IntegralSerializer<long long> > returned;
    //disallow copying
	EMFreelist(EMFreelist const&);
	EMFreelist& operator=(EMFreelist const&);
public:
    EMFreelist(string const& filenameSuffix, int cacheSize = 2):
        nodes("Nodes" + filenameSuffix, cacheSize),
        returned("Returned" + filenameSuffix){}
    long long allocate(POD const& item = POD())
    {
        if(returned.getSize() > 0)
        {//reuse the last deallocated node
            long long result = returned[returned.getSize() - 1];
            returned.removeLast();
            nodes.set(item, result);
            return result;
        }
        else
        {
            nodes.append(item);
            return nodes.getSize() - 1;
        }
    }
    void deallocate(long long i)
    {//no efficient way to check if already deallocated
        assert(i >= 0 && i < nodes.getSize());
        returned.append(i);
    }
    POD operator[](long long i)
    {//no efficient way to check if already deallocated
        assert(i >= 0 && i < nodes.getSize());
        return nodes[i];
    }
    void set(POD const& item, long long i)
    {//no efficient way to check if already deallocated
        assert(i >= 0 && i < nodes.getSize());
        nodes.set(item, i);
    }
};

}//end namespace
#endif
