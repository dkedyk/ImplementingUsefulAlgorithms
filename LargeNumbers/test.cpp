#include "LargeNumberTestAuto.h"
#include "../Utils/Debug.h"
using namespace igmdk;

void DDDNumber()
{
    Number n(2);
    Number TwoPow100 = power(n, Number(100));

    cout << "breakpoint" << endl;
}

void testRationalizeHelper(double x)
{
    DEBUG(numeric_limits<double>::digits);
    pair<long long, int> me = rationalize(x);
    DEBUG(me.first);
    DEBUG(me.second);
    DEBUG(me.first * pow(2, me.second));
}

void testRationalize()
{
    testRationalizeHelper(10);
    testRationalizeHelper(0);
    testRationalizeHelper(1.0/3);
}

int main()
{
    testAllAutoLargeNumber();
    testRationalize();
    //return 0;
    DDDNumber();

	return 0;
}
