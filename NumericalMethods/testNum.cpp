#include <cassert>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <limits>
#include <memory> //for shared ptr
#include "NumericalMethods.h"
#include "../NumericalOptimization/GlobalNumericalOptimization.h"
#include "../RandomNumberGeneration/Statistics.h"
#include "../Utils/DEBUG.h"
#include "../ExternalMemoryAlgorithms/CSV.h"
#include "NumericalMethodsTestAuto.h"
#include "TestFunctions1D.h"
using namespace std;
using namespace igmdk;

int main()
{
    testELessAuto();
    return 0;

    DEBUG(numeric_limits<double>::min());
    DEBUG(numeric_limits<double>::max());
    DEBUG(numeric_limits<double>::epsilon());

    return 0;
}
