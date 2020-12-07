#include <cassert>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <limits>
#include <memory> //for shared ptr
#include "AllRootsFinder.h"
#include "../NumericalOptimization/GlobalNumericalOptimization.h"
#include "../RandomNumberGeneration/Statistics.h"
#include "../Utils/DEBUG.h"
#include "../ExternalMemoryAlgorithms/CSV.h"
#include "NumericalMethodsTestAuto.h"
#include "TestFunctions1D.h"
using namespace std;
using namespace igmdk;

void testAllRootSolver()
{
    Vector<double> coefs(2);//x^2 - 1 = 0, roots: 1, -1
    coefs[0] = -1;
    coefs[1] = 0;
    Vector<complex<double> > roots = findAllRoots(coefs);
    DEBUG("roots");
    roots.debug();
    Vector<double> coefs2(3);//MATLAB EXAMPLE x^3 -7x + 6 = 0, roots = 1, 2, 3
    coefs2[0] = 6;
    coefs2[1] = -7;
    coefs2[2] = 0;
    Vector<complex<double> > roots2 = findAllRoots(coefs2);
    DEBUG("roots2");
    roots2.debug();
    Vector<double> coefs3(3);//0 EXAMPLE x^3 = 0
    coefs3[0] = 0;
    coefs3[1] = 0;
    coefs3[2] = 0;
    Vector<complex<double> > roots3 = findAllRoots(coefs3);
    DEBUG("roots3");
    roots3.debug();
}

struct SQUARE2
{
    double operator()(double x)const{return x * x - 0.25;}
};
void testChebRoots()
{
    ChebFunction cf(SQUARE2(), 16);
    Vector<double> rroots = cf.findAllRealRoots();
    DEBUG("roots");
    rroots.debug();//expect 0.5 and -0.5
    rroots = findAllRealRootsCheb(SQUARE2(), -1, 1);
    DEBUG("roots adaptive");
    rroots.debug();//expect 0.5 and -0.5
    rroots = findAllRealRootsCheb(TestFunctions1D::Sin(), -10, 10);
    DEBUG("sin roots adaptive");
    rroots.debug();//expect 0 and +- iPi for i 1 to 3
}

int main()
{
    testChebRoots();
    return 0;
    testAllRootSolver();
    return 0;

    return 0;
}
