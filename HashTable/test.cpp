#include "ChainingHashTable.h"
#include "LinearProbingHashTable.h"
#include "BloomFilter.h"
#include "HashTableTestAuto.h"
#include "../RandomNumberGeneration/Statistics.h"
#include "../ExternalMemoryAlgorithms/CSV.h"
#include <iostream>
#include <cmath>
#include <functional>
using namespace igmdk;

template<typename T> void timeRT(int N)
{
    T t;
	for(int i = 0; i < N; ++i)
	{
		t.insert(i,i);
	}
	for(int j = 0; j < 5; ++j)
	{
		for(int i = 0; i < N; ++i)
		{
			assert(t.find(i));
			assert(*t.find(i) == i);
			//t.remove(i);
		}
	}
}

template<typename H> struct FunctionTesterCI
{
    void operator()()const
    {
        timeRT<ChainingHashTable<int, int, H> >(1000000);
    }
};
template<typename H> struct FunctionTesterLI
{
    void operator()()const
    {
        timeRT<LinearProbingHashTable<int, int, H> >(1000000);
    }
};
template<typename H> struct FunctionTesterCF
{
    void operator()()const
    {
        timeRT<ChainingHashTable<Struct10_2, int, H> >(100000);
    }
};
template<typename H> struct FunctionTesterLF
{
    void operator()()const
    {
        timeRT<LinearProbingHashTable<Struct10_2, int, H> >(100000);
    }
};

template<typename H> double testInt(H const& h)
{
    int now = clock();
    unsigned int sum = 0;
    for(int i = 0; i < 1000000000; ++i)//
	{
		sum += h(i);
	}
	DEBUG(sum);
	return (clock() - now) * 1.0/CLOCKS_PER_SEC;
}

template<typename H> double testStruct10(H const& h)
{
    int now = clock();
    unsigned int sum = 0;
    for(int i = 0; i < 100000000; ++i)//
	{
		sum += h(Struct10_2(i));
	}
	DEBUG(sum);
	return (clock() - now) * 1.0/CLOCKS_PER_SEC;
}

template<typename HASHER = EHash<BUHash> > class PODHash
{//test use only, not safe in general as explained in the book
    HASHER h;
public:
    PODHash(HASHER const& theH): h(theH){}
    PODHash(unsigned long long m): h(m){}
    typedef typename HASHER::WORD_TYPE WORD_TYPE;
    unsigned long long max()const{return h.max();}
    template<typename POD> unsigned long long operator()(POD const& x)const
        {return h((unsigned char*)&x, sizeof(x));}
    typedef EMPTY Builder;
};

template<typename H> void testSpeedHHelper(string const& name, H const& h,
    Vector<Vector<string> >& matrix)
{
    Vector<string> titles, row;
    DEBUG(name);
    titles.append("Hasher");
    row.append(name);
    double intSpeed = testInt(h);
    DEBUG(intSpeed);
    titles.append("Int");
    row.append(toStringDouble(intSpeed));
    Vector<std::function<void(void)> > functors;
    Vector<string> names;
    functors.append(FunctionTesterCI<H>());
    names.append("Ch");
    functors.append(FunctionTesterLI<H>());
    names.append("LP");
    double Struct10Speed = testStruct10(PODHash<H>(h));
    DEBUG(Struct10Speed);
    titles.append("Struct10_10");
    row.append(toStringDouble(Struct10Speed));
    functors.append(FunctionTesterCF<PODHash<H> >());
    names.append("Ch");
    functors.append(FunctionTesterLF<PODHash<H> >());
    names.append("LP");
    for(int i = 0; i < functors.getSize(); ++i)
    {
        DEBUG(names[i]);
        IncrementalStatistics si = MonteCarloSimulate(
            SpeedTester<std::function<void(void)> >(functors[i]), 100);
        DEBUG(si.getMean());
        titles.append(names[i]);
        row.append(toStringDouble(si.getMean()));
        DEBUG(si.getStandardErrorSummary().error95());
        DEBUG(si.minimum);
        DEBUG(si.maximum);
        titles.append("+-");
        titles.append("Min");
        titles.append("Max");
        row.append(toStringDouble(si.getStandardErrorSummary().error95()));
        row.append(toStringDouble(si.minimum));
        row.append(toStringDouble(si.maximum));
    }
    if(matrix.getSize() == 0) matrix.append(titles);
    matrix.append(row);
}

void sillyTest()
{
    EHash<BHash<PrimeHash> > h(64);
    for(int i = 0; i < 100; ++i)
    {
        DEBUG(h(i));
        DEBUG(i & 63);
    }
    int m = twoPower(20);
    DEBUG(testStruct10(PODHash<EHash<BUHash> >(m)));
    DEBUG(testInt(EHash<BHash<PrimeHash> >(m)));
    DEBUG(testStruct10(PODHash<EHash<BHash<PrimeHash> > >(m)));
    DEBUG(testStruct10(PODHash<EHash<BHash<PrimeHash2> > >(m)));

}

void testSpeedH()
{
    Vector<Vector<string> > matrix;
    int m = twoPower(20);
    testSpeedHHelper("E-BU", EHash<BUHash>(m), matrix);
    testSpeedHHelper("E-B-Prime", EHash<BHash<PrimeHash> >(m), matrix);
    testSpeedHHelper("E-B-Prime2", EHash<BHash<PrimeHash2> >(m), matrix);
    testSpeedHHelper("E-B-FNV", EHash<BHash<FNVHash> >(m), matrix);
    testSpeedHHelper("E-M-FNV", EHash<MHash<FNVHash> >(m), matrix);
    testSpeedHHelper("E-B-FNV64", EHash<BHash<FNVHash> >(m), matrix);
    testSpeedHHelper("E-B-X64", EHash<BHash<Xorshift64Hash> >(m), matrix);
    testSpeedHHelper("E-B-Table", EHash<BHash<TableHash> >(m), matrix);
    createCSV(matrix, "TestResults.csv");
}

void DDDChaining()
{
    ChainingHashTable<int, int> chainingH0to9;
    for(int i = 0; i < 10; ++i)
	{
		chainingH0to9.insert(i, i);
	}
    cout << "breakpoint" << endl;
}

void DDDLinearProbing()
{
    LinearProbingHashTable<int, int> linearProbingH0to9;
    for(int i = 0; i < 10; ++i)
	{
		linearProbingH0to9.insert(i, i);
	}
    cout << "breakpoint" << endl;
}

void DDDBloomFilter()
{
    BloomFilter<int> bF16_3_0to9(16, 3);
    for(int i = 0; i < 10; ++i)
	{
		bF16_3_0to9.insert(i);
	}
}

int main()
{
    testAllAutoHashTable();
    //return 0;
    testSpeedH();
    return 0;
    sillyTest();
    return 0;
    DDDChaining();
    DDDLinearProbing();
    DDDBloomFilter();
	return 0;
}
