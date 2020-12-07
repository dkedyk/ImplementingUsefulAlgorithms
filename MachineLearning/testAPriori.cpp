#include "APriori.h"
#include "../Utils/Debug.h"
using namespace igmdk;

void testAPriori()
{
    Vector<Vector<int> > baskets;
    Vector<int> b1, b2, b3, b4;
    b1.append(0);
    b1.append(1);
    b1.append(2);
    b1.append(3);
    baskets.append(b1);
    b2.append(5);
    b2.append(1);
    b2.append(2);
    b2.append(4);
    baskets.append(b2);
    b3.append(7);
    b3.append(1);
    b3.append(2);
    b3.append(6);
    baskets.append(b3);
    b4.append(1);
    b4.append(0);
    b4.append(4);
    b4.append(6);
    baskets.append(b4);
    APriori ap;
    ap.noCutProcess(baskets, 3);
    for(LCPTreap<Vector<int>, int>::Iterator i(ap.counts.begin()); i != ap.counts.end(); ++i)
    {
        for(int j = 0; j < i->key.getSize(); ++j)
        {
            DEBUG(i->key[j]);
        }
        DEBUG(i->value);
    }
}

int main(int argc, char *argv[])
{
    testAPriori();
}


