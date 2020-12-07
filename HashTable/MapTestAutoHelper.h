#ifndef IGMDK_MAP_TEST_AUTO_HELPER_H
#define IGMDK_MAP_TEST_AUTO_HELPER_H

#include "../Utils/Bitset.h"

namespace igmdk{

struct Struct10_2
{
	enum{SIZE = 10};
	int array[SIZE];
	Struct10_2(int last)
	{
		for(int i = 1; i < SIZE; ++i)
		{
			array[i] = i;
		}
		array[0] = last;
	}
	bool operator==(Struct10_2 const& rhs)const
	{
		for(int i = 0; i < SIZE; ++i)
		{
			if(array[i] != rhs.array[i]) return false;
		}
		return true;
	}
	bool operator<(Struct10_2 const& rhs)const
	{
		for(int i = 0; i < SIZE; ++i)
		{
			if(array[i] < rhs.array[i]) return true;
			if(array[i] > rhs.array[i]) return false;
		}
		return false;
	}
	int getSize()const{return SIZE;}
	int const operator[](int i)const{return array[i];}
};

template<typename MAP_II> void testMapAutoHelper(int n = 100000)
{
    MAP_II m;
	for(int i = 0; i < n; ++i) m.insert(i, -i);
    Bitset<> seen(n), allSet(n);
    seen.setAll(false);
    allSet.setAll(true);
    for(typename MAP_II::Iterator e = m.end(), i = m.begin(); i != e; ++i)
    {
        assert(!seen[i->key]);
        seen.set(i->key, true);
        assert(i->value == -i->key);
    }
    assert(seen == allSet);
    for(int i = 0; i < n; ++i)
    {
        assert(m.find(i));
        assert(*m.find(i) == -i);
        m.remove(i);
        assert(!m.find(i));
    }
}

}//end namespace
#endif
