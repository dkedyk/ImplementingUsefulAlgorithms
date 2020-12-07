#include "CuckooHashTable.h"
#include "../RandomNumberGeneration/Statistics.h"
using namespace igmdk;
#include <iostream>
#include <cmath>
using namespace std;

class Xorshift64HashSeeded
{
    unsigned long long seed;
public:
    Xorshift64HashSeeded():seed(GlobalRNG().next()){}
    unsigned int hash(unsigned long long x)
        {return QualityXorshift64::transform(seed+x);}
	template<typename NUMBER> unsigned int hash(NUMBER* array, int size)
	{
	    unsigned long long sum = seed;
		for(int i = 0; i < size; ++i) sum =
            QualityXorshift64::transform(sum + array[i]);
		return sum;
	}
};

struct Fat2
{
	enum{SIZE = 10};
	int array[SIZE];
	Fat2(int last)
	{
		for(int i = 1; i < SIZE; ++i)
		{
			array[i] = i;
		}
		array[0] = last;
	}
	bool operator==(Fat2 const& rhs)const
	{
		for(int i = 0; i < SIZE; ++i)
		{
			if(array[i] != rhs.array[i]) return false;
		}
		return true;
	}
	bool operator<(Fat2 const& rhs)const
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

struct Action
{
    void operator()(int const& key, int const& item)
    {
        cout << key << endl;
    }
};

void timeRT()
{
	CuckooHashTable<int, int> t;
	int N = 1500000;
	for(int i = 0; i < N; ++i)
	{
		t.insert(i,i);
	}
    /*CuckooHashTable<int, int>::Iterator iter(t);
    while(iter.hasNext())
    {
        cout << iter.next()->key << endl;
    }

    Action action;
    t.forEach(action);*/
	for(int j = 0; j < 1; ++j)
	{
		for(int i = 0; i < N; ++i)
		{
			assert(t.find(i));
			//cout << *t.find(i) << endl;
			assert(*t.find(i) == i);
			t.remove(i);
		}
	}
    DEBUG("done");
	/*for(int i = 0; i < 1500000; ++i)
	{
		t.remove(i);
	}*/
}

int main()
{
	timeRT();
	return 0;
}
