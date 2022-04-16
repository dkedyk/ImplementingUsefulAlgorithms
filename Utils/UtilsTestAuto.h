#ifndef IGMDK_UTILS_TEST_AUTO_H
#define IGMDK_UTILS_TEST_AUTO_H

using namespace std;
#include "Utils.h"
#include "Stack.h"
#include "Queue.h"
#include "Bitset.h"
#include "UnionFind.h"
#include "../RandomNumberGeneration/Random.h"

namespace igmdk{

void testComparatorsAuto()
{
    DEBUG("testComparatorsAuto");
    double a = 0, b = 1;
    DefaultComparator<double> dc;
    assert(dc(a, b));
    assert(!dc.isEqual(a, b));
    ReverseComparator<double> rc;
    assert(rc(b, a));
    assert(!rc.isEqual(b, a));
    PointerComparator<double> pc;
    assert(pc(&a, &b));
    assert(!pc.isEqual(&a, &b));
    double ab[] = {a, b};
    IndexComparator<double> ic(ab);
    assert(ic(0, 1));
    assert(!ic.isEqual(0, 1));
    PairFirstComparator<double, int> pfc;
    assert(pfc(make_pair(a, 0), make_pair(b, 0)));
    assert(!pfc.isEqual(make_pair(a, 0), make_pair(b, 0)));
    DEBUG("testComparatorsAuto passed");
}

void testStackAuto()
{
    DEBUG("testStackAuto");
    int n = 1000000;
    Stack<int> s;
    for(int i = 0; i < n/2; ++i) s.push(i);
    for(int i = 0; i < n/2; ++i)
    {
        s.push(n/2 + i);
        assert(s.pop() == n/2 + i);
        assert(s.pop() == n/2 - 1 - i);
    }
    assert(s.isEmpty());
    DEBUG("testStackAuto passed");
}

void testQueueAuto()
{
    DEBUG("testQueueAuto");
    int n = 1000000;
    Queue<int> q;
    for(int i = 0; i < n/2; ++i) q.push(i);
    for(int i = 0; i < n/2; ++i) assert(q[i] == i);
    for(int i = 0; i < n/2; ++i)
    {
        q.push(n/2 + i);
        assert(q.pop() == i * 2);
        assert(q.pop() == i * 2 + 1);
    }
    assert(q.isEmpty());
    DEBUG("testQueueAuto passed");
}

void testLgAuto()
{
    assert(lgFloor(2) == 1);
    assert(lgFloor(3) == 1);
    assert(lgCeiling(3) == 2);
    assert(lgCeiling(4) == 2);
    assert(nextPowerOfTwo(7) == 8);
    assert(nextPowerOfTwo(8) == 8);
}

template<typename WORD> void testBitsAutoHelper()
{
    assert(!is_signed<WORD>::value);
    int b = numeric_limits<WORD>::digits;
    WORD w = 0;
    for(int i = 0; i < b; ++i)
    {
        assert(!Bits::get(w, i));
        Bits::set(w, i, true);
        assert(Bits::get(w, i));
        Bits::set(w, i, Bits::flip(w, i));
        assert(!Bits::get(w, i));
    }
    assert(w == 0);
    w = Bits::upperMask(4);//expect 1*11110000
    for(int i = 0; i < 4; ++i) assert(!Bits::get(w, i));
    for(int i = 4; i < b; ++i) assert(Bits::get(w, i));
    w = Bits::lowerMask(4);//expect 0*00001111
    for(int i = 0; i < 4; ++i) assert(Bits::get(w, i));
    for(int i = 4; i < b; ++i) assert(!Bits::get(w, i));
    w = Bits::middleMask(3, 4);//expect *01111000
    for(int i = 0; i < 3; ++i) assert(!Bits::get(w, i));
    for(int i = 3; i < 7; ++i) assert(Bits::get(w, i));
    for(int i = 7; i < b; ++i) assert(!Bits::get(w, i));
    w = 0;
    Bits::setValue(w, 7, 5, 3);//expect 0*11100000
    for(int i = 0; i < 5; ++i) assert(!Bits::get(w, i));
    for(int i = 5; i < 8; ++i) assert(Bits::get(w, i));
    for(int i = 8; i < b; ++i) assert(!Bits::get(w, i));
    w = Bits::getValue(w, 5, 3);//expect 0*00000111
    for(int i = 0; i < 3; ++i) assert(Bits::get(w, i));
    for(int i = 3; i < b; ++i) assert(!Bits::get(w, i));
}
void testBitsAuto()
{
    DEBUG("testBitsAuto");
    testBitsAutoHelper<unsigned char>();
    testBitsAutoHelper<unsigned long long>();
    DEBUG("testBitsAuto passed");
}

template<typename WORD> void testPopCountAutoHelper()
{
    assert(!is_signed<WORD>::value);
    int b = numeric_limits<WORD>::digits;
    for(int i = 1; i < b; ++i)
    {
        WORD w = 0;
        Vector<int> setBits = GlobalRNG().randomCombination(i, b);
        for(int j = 0; j < i; ++j) Bits::set(w, setBits[j], true);
        assert(popCountWord(w) == i);
    }
}
void testPopCountAuto()
{
    DEBUG("testPopCountAuto");
    testPopCountAutoHelper<unsigned char>();
    testPopCountAutoHelper<unsigned long long>();
    DEBUG("testPopCountAuto passed");
}

template<typename WORD> void testReverseBitsAutoHelper()
{
    assert(!is_signed<WORD>::value);
    int b = numeric_limits<WORD>::digits;
    for(int i = 1; i < b; ++i)
    {
        WORD w = 0;
        Vector<int> setBits = GlobalRNG().randomCombination(i, b);
        for(int j = 0; j < i; ++j) Bits::set(w, setBits[j], true);
        assert(reverseBits(reverseBits(w)) == w);
    }
    assert(reverseBits<WORD>(7, 3) == 7);
}
void testReverseBitsAuto()
{
    DEBUG("testReverseBitsAuto");
    testReverseBitsAutoHelper<unsigned char>();
    testReverseBitsAutoHelper<unsigned long long>();
    DEBUG("testReverseBitsAuto passed");
}

void testBitsetAuto()
{
    DEBUG("testBitsetAuto");
    Bitset<> b(60);
    b.setValue(15, 14, 4);
    assert(b.getValue(14, 4) == 15);
    b.appendValue(14, 4);
    assert(b.getValue(60, 4) == 14);
    //test 2
    Bitset<unsigned char> bs(21);
	for(int i = 0; i < 21; i+=3)
	{
		bs.set(i, true);
	}
	bs.debug();
	for(int i = 0; i < 21; ++i) assert(bs[i] == !(i % 3));
	bs <<= 9;
	bs.debug();
	for(int i = 0; i < 9; ++i) assert(!bs[i]);
	for(int i = 9; i < 21; ++i) assert(bs[i] == !(i % 3));
	bs >>= 9;//fails 4th bit set
	bs.debug();
	for(int i = 0; i < 12; ++i) assert(bs[i] == !(i % 3));
	for(int i = 12; i < 21; ++i) assert(!bs[i]);
	bs.flip();
	for(int i = 0; i < 12; ++i){assert(bs[i] == bool(i % 3));}
	for(int i = 12; i < 21; ++i) assert(bs[i]);
	bs.debug();
	Bitset<unsigned char> bsOld = bs;
	bs.reverse();
	bs.debug();
	bs.reverse();
	assert(bs == bsOld);
    DEBUG("testBitsetAuto passed");
}

void testUFAuto()
{
    DEBUG("testUFAuto");
    UnionFind uf(10);
	assert(!uf.areEquivalent(4, 8));
	uf.join(4,8);
	assert(uf.areEquivalent(4, 8));
	assert(!uf.areEquivalent(4, 9));
	for(int i = 0; i < 20; ++i) uf.addSubset();
	assert(uf.areEquivalent(4, 8));
	assert(!uf.areEquivalent(4, 9));
	DEBUG("testUFAuto passed");
}


void testAllAutoUtils()
{
    DEBUG("testAllAutoUtils");
    testComparatorsAuto();
    testStackAuto();
    testQueueAuto();
    testLgAuto();
    testBitsAuto();
    testPopCountAuto();
    testReverseBitsAuto();
    testBitsetAuto();
    testUFAuto();
}

}//end namespace
#endif
