#include "LinearProgramming.h"
#include "../Utils/DEBUG.h"
using namespace std;
using namespace igmdk;

void testSimplex()
{
    Matrix<double> B = Matrix<double>::identity(3), N(3, 2);
    N(0, 0) = -2;
    N(0, 1) = 1;
    N(1, 0) = -1;
    N(1, 1) = 2;
    N(2, 0) = 1;
    N(2, 1) = 0;
    Vector<double> b, cB(3, 0), cN;
    b.append(2);
    b.append(7);
    b.append(3);
    cN.append(-1);
    cN.append(-2);
    LinearProgrammingSimplex s(B, N, cB, cN, b);
    Vector<pair<int, double> > result = s.solve(), expectedResult;
    expectedResult.append(make_pair(0, 5.0));
    expectedResult.append(make_pair(2, 3.0));
    expectedResult.append(make_pair(1, 3.0));
    for(int i = 0; i < result.getSize(); ++i)
    {
        DEBUG(result[i].first);
        DEBUG(result[i].second);
        assert(result[i] == expectedResult[i]);
    }
}

int main()
{
    testSimplex();
    return 0;
}
