#include "ContraintProcessing.h"
#include "../Utils/Debug.h"
using namespace igmdk;
using namespace std;

void testAC3()
{
    //Example: trace execution on (0)<->(0,1)<->(0,1,2)<->(0,1,2,3) with <->
    //denoting alldiff constraint
    AllDifferent ad;

    ConstraintGraph<AllDifferent::Handle> cg;
    for(int i = 0; i < 4; ++i)
    {
        cg.addVariable(i+1);
        cg.variables[i].setAll(true);
        ad.addVariable(i);
    }
    for(int i = 0; i < 4; ++i) cg.variables[i].debug();
    for(int i = 0; i < 3; ++i) cg.addConstraint(i, i+1, ad.handle);
        //for(int j = i + 1; j < 4; ++j)
        //    cg.addConstraint(i, j, AllDifferent());
    DEBUG(cg.SAC3());
    for(int i = 0; i < 4; ++i) cg.variables[i].debug();
    //Vector<int> result;
    //backtrackFind(cg.variables, result, AllDifferent());
    //for(int i = 0; i < result.getSize(); ++i) DEBUG(result[i]);
}


void testSudoku()
{
    int easy[] =
    {
        2,0,1,7,9,5,0,0,0,
        0,0,9,0,0,8,0,1,0,
        0,0,0,3,0,1,0,0,7,
        0,2,0,5,0,0,1,7,8,
        0,8,0,0,0,0,0,9,0,
        1,5,7,0,0,4,0,3,0,
        6,0,0,8,0,2,0,0,0,
        0,9,0,6,0,0,5,0,0,
        0,0,0,1,7,9,3,0,6
    };
    int medium[] =
    {
        0,0,5,0,0,3,2,9,0,
        9,0,0,2,0,0,0,3,4,
        0,0,0,0,1,0,8,0,0,
        0,0,0,0,9,0,0,7,1,
        0,0,0,6,0,5,0,0,0,
        7,3,0,0,2,0,0,0,0,
        0,0,7,0,6,0,0,0,0,
        6,8,0,0,0,9,0,0,2,
        0,5,2,8,0,0,6,0,0
    };
    int hard[] =
    {
        0,1,2,0,6,0,8,0,0,
        0,0,0,0,3,0,0,5,0,
        6,0,0,4,0,0,0,0,7,
        0,0,6,0,0,0,0,1,0,
        0,9,7,0,0,0,6,4,0,
        0,8,0,0,0,0,7,0,0,
        8,0,0,0,0,1,0,0,3,
        0,4,0,0,5,0,0,0,0,
        0,0,1,0,2,0,9,7,0
    };
    {
        Sudoku sd(easy);
        DEBUG("Easy Solution");
        sd.printSolution();
    }
    {
        Sudoku sd(medium);
        DEBUG("Medium Solution");
        sd.printSolution();
    }
    {
        Sudoku sd(hard);
        DEBUG("Hard Solution");
        sd.printSolution();
    }

}


int main()
{
    //testAC3();
    testSudoku();
	return 0;
}
