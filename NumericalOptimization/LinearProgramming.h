#ifndef IGMDK_LINEAR_PROGRAMMING_H
#define IGMDK_LINEAR_PROGRAMMING_H
#include "../Utils/Utils.h"
#include "../NumericalMethods/Matrix.h"
namespace igmdk{

struct LinearProgrammingSimplex
{
    Matrix<double> B, N;
    Vector<double> b, cB, cN, x;
    Vector<int> p;
    bool isUnbounded;
    bool performIteration()
    {
        LUP<double> lup(B), lupT(B.transpose());
        x = lup.solve(b);
        //check if x is optimal or find entering variable
        Vector<double> y = cN - lupT.solve(cB) * N;
        int entering = 0;
        double bestValue = y[0];
        for(int i = 1; i < y.getSize(); ++i) if(y[i] < bestValue)
            {
                bestValue = y[i];
                entering = i;
            }
        if(bestValue >= 0) return false;
        //find leaving variable
        Vector<double> a;
        for(int i = 0; i < N.rows; ++i) a.append(N(i, entering));
        a = lup.solve(a);
        int leaving = -1;
        double minRatio, maxA = -1;
        for(int i = 0; i < x.getSize(); ++i) if(a[i] > 0)
            {
                double newRatio = x[i]/a[i];
                if(leaving == -1 || minRatio > newRatio)
                {
                    leaving = i;
                    maxA = max(maxA, a[i]);
                    minRatio = newRatio;
                }
            }
        if(maxA <= 0){isUnbounded = true; return false;}
        //swap variables
        for(int i = 0; i < N.rows; ++i){swap(B(i, leaving), N(i, entering));}
        swap(p[leaving], p[entering]);
        swap(cB[leaving], cN[entering]);
        return true;
    }
    LinearProgrammingSimplex(Matrix<double>const&B0, Matrix<double>
        const& N0, Vector<double>const& cB0, Vector<double>const& cN0,
        Vector<double> const& b0): isUnbounded(false), B(B0), N(N0), cB(cB0),
        cN(cN0), b(b0), x(b)
        {for(int i = 0; i < cB.getSize() + cN.getSize(); ++i) p.append(i);}
    Vector<pair<int, double> > solve()
    {
        while(performIteration());
        Vector<pair<int, double> > result;
        if(!isUnbounded)
            for(int i = 0; i < x.getSize(); ++i)
                result.append(make_pair(p[i], x[i]));
        return result;
    }
};

}//end namespace igmdk
#endif
