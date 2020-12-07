#ifndef IGMDK_DEBUG_H
#define IGMDK_DEBUG_H
#include <iostream>
#include <iomanip>
using namespace std;
namespace igmdk{
//print the expression, a space, and a new line
#define DEBUG(var) cout << #var " "<< setprecision(17) << (var) << endl;

}//end namespace
#endif
