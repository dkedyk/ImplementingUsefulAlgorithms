#include "Treap.h"
#include "LCPTreap.h"
#include "SkipList.h"
#include "Trie.h"
#include "../Utils/Debug.h"
#include "DynamicSortedSequenceTestAuto.h"
using namespace igmdk;

struct Struct10
{
	enum{SIZE = 10};
	int array[SIZE];
	Struct10(int last)
	{
		for(int i = 0; i < SIZE-1; ++i)
		{
			array[i] = i;
		}
		array[SIZE-1] = last;
	}
	bool operator==(Struct10 const& rhs)const
	{
		for(int i = 0; i < SIZE; ++i)
		{
			if(array[i] != rhs.array[i]) return false;
		}
		return true;
	}
	bool operator<(Struct10 const& rhs)const
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

struct Struct10_4
{
	enum{SIZE = 40};
	unsigned char array[SIZE];
	Struct10_4(unsigned int last)
	{
		for(int i = 0; i < SIZE; ++i)
		{
			array[i] = i;
		}
		for(int i = 0; i < 4; ++i)
		{
            array[SIZE-1 - i] = last % 256;
            last /= 256;
		}
	}
	bool operator==(Struct10_4 const& rhs)const
	{
		for(int i = 0; i < SIZE; ++i)
		{
			if(array[i] != rhs.array[i]) return false;
		}
		return true;
	}
	bool operator<(Struct10_4 const& rhs)const
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

struct Struct10_5
{
	enum{SIZE = 40};
	unsigned char array[SIZE];
	Struct10_5(int last)
	{
		for(int i = 0; i < SIZE; ++i)
		{
			array[i] = i;
		}
		for(int i = 0; i < 4; ++i)
		{
            array[i] = last % 256;
            last /= 256;
		}
	}
	bool operator==(Struct10_5 const& rhs)const
	{
		for(int i = 0; i < SIZE; ++i)
		{
			if(array[i] != rhs.array[i]) return false;
		}
		return true;
	}
	bool operator<(Struct10_5 const& rhs)const
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



void timeRT()
{
	//LcpTreap<Struct10_5, int> trie;
	//Treap<int, int> trie;
	SkipList<int, int> trie;
    //TrieWrapper<int, int> trie;
	int N = 1000000;
	for(int i = 0; i < N; ++i)
	{
		trie.insert(i, i);
	}
	for(int j = 0; j < 10; ++j)
	{
		for(int i = 0; i < N; ++i)
		{
			assert(trie.find(i));
			assert(*trie.find(i) == i);
			//trie.remove(i);
		}
	}
}
void DDDSkiplist()
{
    SkipList<int, int> Skiplist0to5;
    for(int i = 0; i < 5; ++i)
	{
		Skiplist0to5.insert(i, i);
	}
	cout << "breakpoint" << endl;
}
void DDDTreap()
{
    Treap<int, int> Treap0to9;
    int n = 10;
    for(int i = 0; i < n; ++i)
	{
		Treap0to9.insert(i, i);
	}
	cout << "breakpoint" << endl;
}

void DDDLCPTreap()
{
    typedef LCPTreap<Vector<char>, int > T;
    T t;
    Vector<char> a, b, c, d, e, f;
    a.append('h');
    a.append('i');

    b.append('h');
    b.append('e');
    b.append('y');

    c.append('h');
    c.append('e');
    c.append('l');
    c.append('l');

    d.append('b');
    d.append('y');
    d.append('e');

    f.append('b');
    f.append('y');

    e.append('h');
    e.append('i');
    e.append('h');
    e.append('o');

    t.insert(a,0);
    t.insert(b,1);
    t.insert(c,2);
    t.insert(d,3);
    t.insert(e,4);
    t.insert(f,5);
    cout << "breakpoint" << endl;
}

void DDDTTT()
{
    typedef TernaryTreapTrie<int> T;
    T TTT;

    TTT.insert((unsigned char*)"hi", 2, 0);
    TTT.insert((unsigned char*)"hey", 3, 1);
    TTT.insert((unsigned char*)"hell", 4, 2);
    TTT.insert((unsigned char*)"bye", 3, 3);
    TTT.insert((unsigned char*)"by", 2, 4);
    TTT.insert((unsigned char*)"hiho", 4, 5);
    cout << "breakpoint" << endl;
}

int main()
{
    testAllAutoDynamicSortedSequence();
    //return 0;
    //NEED PRED/SUCC TEST CASES MAYBE USING SUCC ITERATION AND PRED REVERSE ITERARION
    DDDTreap();
    //return 0;
    DDDSkiplist();
    DDDLCPTreap();
    DDDTTT();

	clock_t start = clock();
	timeRT();
	int tFL = (clock() - start);
    cout << "FL: "<<tFL << endl;
	return 0;
}
