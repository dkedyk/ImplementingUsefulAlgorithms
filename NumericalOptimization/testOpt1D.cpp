#include <cassert>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <limits>
#include <memory> //for shared ptr
#include "NumericalOptimization.h"
#include "testOptCommon.h"
#include "SPSA.h"
using namespace std;
using namespace igmdk;

struct OptTestFunctions1D
{
    struct BaseF
    {
        virtual double operator()(double const& x)const = 0;
        virtual string name()const = 0;
        virtual double getX0()const{return 0;}
        virtual double getAnswer(){return 3;}
    };
    struct Linear: public BaseF
    {//easy
        double operator()(double const& x)const{return x - 3;}
        string name()const{return "Linear";}
    };
    struct QuadE: public BaseF
    {//hard, differentiable
        double operator()(double const& x)const
        {
            double temp = x - 3;
            return temp * temp - 0.001;
        }
        string name()const{return "QuadE";}
    };
    struct Poly6: public BaseF
    {//differentiable, ill-conditioned
        double operator()(double const& x)const{return pow(x - 3, 6);}
        string name()const{return "Poly6";}
    };
    struct Sqrt2F: public BaseF
    {//differentiable, ill-conditioned
        double a, x0;
        Sqrt2F(double theA = 2, double theX0 = 0): a(theA), x0(theX0){}
        double operator()(double const& x)const{return x * x - a;}
        string name()const{return "Sqrt2F_" + toStringDouble(a) + "_" + toStringDouble(x0);}
        double getX0()const{return 0;}
        double getAnswer()const{return sqrt(x0);}
    };
    static int evalCount;
    struct MetaF
    {
        shared_ptr<BaseF> f;
        template<typename F> MetaF(shared_ptr<F> const& theF): f(theF){}
        double operator()(double const& x)const
        {
            ++evalCount;
            double temp = (*f)(x);
            return temp * temp;
        }
        string getName()const{return f->name();}
        double getX0()const{return f->getX0();}
        double getAnswer()const{return f->getAnswer();}
    };
    static Vector<MetaF> getFunctions()
    {
        Vector<MetaF> result;
        result.append(MetaF(make_shared<Linear>()));
        result.append(MetaF(make_shared<QuadE>()));
        result.append(MetaF(make_shared<Poly6>()));
        result.append(MetaF(make_shared<Sqrt2F>()));
        result.append(MetaF(make_shared<Sqrt2F>(2, 2)));
        return result;
    }
};
int OptTestFunctions1D::evalCount = 0;

void debugOpt1DResult(pair<double, double> const& result, double answer,
    Vector<Vector<string> >& matrix)
{
    DEBUG(normInf(result.first - answer));
    DEBUG(result.second);
    DEBUG(OptTestFunctions1D::evalCount);
    matrix.lastItem().append(toStringDouble(normInf(result.first - answer)));
    matrix.lastItem().append(toStringDouble(result.second));
    matrix.lastItem().append(toStringDouble(OptTestFunctions1D::evalCount));
    OptTestFunctions1D::evalCount = 0;
}
template<typename FUNCTION> void testOptHelper1D(FUNCTION const& f, double x0,
    Vector<Vector<string> >& matrix)
{
    DEBUG("BracketGSUnimodal");
    matrix.lastItem().append("BracketGSUnimodal");
    debugOpt1DResult(minimizeGSBracket(f, x0), f.getAnswer(), matrix);
}

void testOpt1D()
{
    Vector<Vector<string> > matrix;
    Vector<OptTestFunctions1D::MetaF> fs = OptTestFunctions1D::getFunctions();
    for(int i = 0; i < fs.getSize(); ++i)
    {
        string name = fs[i].getName();
        DEBUG(name);
        matrix.append(Vector<string>());
        matrix.lastItem().append(name);
        testOptHelper1D(fs[i], fs[i].getX0(), matrix);
    }
    int reportNumber = time(0);
    string filename = "reportOpt1D" + to_string(reportNumber) + ".csv";
    createCSV(matrix, filename.c_str());
}

int main()
{
    testOpt1D();
    return 0;
}
