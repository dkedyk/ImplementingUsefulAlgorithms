#include "Bitset.h"
using namespace igmdk;
#include <bitset>

template<int n, typename WORD> void printBits(WORD x)
{
    bitset<n> b(x);
    cout << b << endl;
}

void testBits()
{
    printBits<8>(Bits::upperMask(4));//expect 11110000
}

int popCountWord2(unsigned long long x)
{
    int n = 0;
    while(x)
    {
        ++n;
        x &= x - 1;
    }
    return n;
}
void timePopCount()
{
    unsigned long long x = 3435345355555ull;
    clock_t start = clock();
    int p = popCountWord(x);
    for(int i = 0; i < 1000000000; ++i) p += popCountWord(x) + 1;
	//timeRT();
	int tFL = (clock() - start);
    cout << p << "PC: "<<tFL << endl;
    start = clock();
    p = popCountWord2(x);
    for(int i = 0; i < 1000000000; ++i) p += popCountWord2(x) + 1;
	//timeRT();
	tFL = (clock() - start);
    cout << p << "PC2: "<<tFL << endl;
}

void DDDBitset()
{
    Bitset<unsigned char> BitsetChar19Every4(19);
	for(int i = 0; i < 19; i += 4)
	{
		BitsetChar19Every4.set(i, true);
	}
	cout << "breakpoint" << endl;
}

int main()
{
    DDDBitset();
    //timePopCount();
    testBits();
	return 0;
}
