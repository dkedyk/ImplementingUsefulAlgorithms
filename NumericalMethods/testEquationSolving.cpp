#include <cassert>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <limits>
#include <memory> //for shared ptr
#include "EquationSolving.h"
#include "../NumericalOptimization/GlobalNumericalOptimization.h"
#include "../RandomNumberGeneration/Statistics.h"
#include "../Utils/DEBUG.h"
#include "../ExternalMemoryAlgorithms/CSV.h"
#include "NumericalMethodsTestAuto.h"
using namespace std;
using namespace igmdk;

struct TestFunctionsMultiD
{//for integral use http://www.emathhelp.net/calculators/calculus-2/definite-integral-calculator
    struct BaseF
    {
        virtual Vector<double> operator()(Vector<double> const& x)const = 0;
        virtual string name()const = 0;
        virtual Vector<double> getX0()const = 0;
        virtual Vector<double> getAnswer()const = 0;
    };
    struct TestFunc0
    {//from Fausett
        struct Func1 : public MultivarFuncHelper::F1DBase
        {
            double operator()(Vector<double> const& x)const
                {return x[0] + pow(x[0], 3) + 10 * x[0] - x[1] - 5;}
        };
        struct Func2 : public MultivarFuncHelper::F1DBase
        {
            double operator()(Vector<double> const& x)const
                {return x[1] + x[0] + pow(x[1], 3) - 10 * x[1] + 1;}
        };
        struct Func: public BaseF
        {
            MultivarFuncHelper f;
            Func1 f1;
            Func2 f2;
            Func()
            {
                f.fs.append(&f1);
                f.fs.append(&f2);
            }
            Vector<double> operator()(Vector<double> const& x)const{return f(x);}
            string name()const{return "TestFunc0";}
            Vector<double> getX0()const{return Vector<double>(2, 0.6);}
            Vector<double> getAnswer()const{return Vector<double>(2, 0);}//???
        };
    };
    struct TestFunc1
    {//from Fausett
        struct Func1 : public MultivarFuncHelper::F1DBase
        {
            double operator()(Vector<double> const& x)const
                {return x[0] * x[0] + x[1] * x[1] + x[2] * x[2] - 1;}
        };
        struct Func2 : public MultivarFuncHelper::F1DBase
        {
            double operator()(Vector<double> const& x)const
                {return x[0] * x[0] + x[2] * x[2] - 0.25;}
        };
        struct Func3 : public MultivarFuncHelper::F1DBase
        {
            double operator()(Vector<double> const& x)const
                {return x[0] * x[0] + x[1] * x[1] - 4 * x[2];}
        };
        struct Func: public BaseF
        {
            MultivarFuncHelper f;
            Func1 f1;
            Func2 f2;
            Func3 f3;
            Func()
            {
                f.fs.append(&f1);
                f.fs.append(&f2);
                f.fs.append(&f3);
            }
            Vector<double> operator()(Vector<double> const& x)const{return f(x);}
            string name()const{return "TestFunc1";}
            Vector<double> getX0()const{return Vector<double>(3, 1);}
            Vector<double> getAnswer()const{return Vector<double>(3, 0);}//???
        };
    };
    struct ExtendedRosenbrock: public BaseF
    {//From Dennis & Schnabel
        int n;
        ExtendedRosenbrock(int theN = 2): n(theN) {assert(theN % 2 == 0);}
        Vector<double> operator()(Vector<double> const& x)const
        {
            Vector<double> fx(n);
            for(int i = 0; i < n/2; ++i)
            {
                double i1 = 2 * i, i2 = i1 + 1;
                fx[i1] = 10 * (x[i2] - x[i1] * x[i1]);
                fx[i2] = 1 - x[i1];
            }
            return fx;
        }
        string name()const{return "ExtendedRosenbrock" + to_string(n);}
        Vector<double> getX0()const
        {
            Vector<double> x0 = Vector<double>(n, 1);
            for(int i = 0; i < n/2; ++i) x0[2 * i] = 1.2;
            return x0;
        }
        Vector<double> getAnswer()const{return Vector<double>(n, 1);}
    };
    struct ExtendedPowellSingular: public BaseF
    {//From Dennis & Schnabel, J singular at solution
        int n;
        ExtendedPowellSingular(int theN = 4): n(theN) {assert(theN % 4 == 0);}
        Vector<double> operator()(Vector<double> const& x)const
        {
            Vector<double> fx(n);
            for(int i = 0; i < n/4; ++i)
            {
                double i1 = 4 * i, i2 = i1 + 1, i3 = i2 + 1, i4 = i3 + 1;
                fx[i1] = x[i1] + 10 * x[i2];
                fx[i2] = sqrt(5) * (x[i3] - x[i4]);
                fx[i3] = (x[i2] - 2 * x[i3]) * (x[i2] - 2 * x[i3]);
                fx[i4] = sqrt(10) * (x[i1] - x[i4]) * (x[i1] - x[i4]);
            }
            return fx;
        }
        string name()const{return "ExtendedPowellSingular" + to_string(n);}
        Vector<double> getX0()const
        {
            Vector<double> x0 = Vector<double>(n, 1);
            for(int i = 0; i < n/4; ++i)
            {
                x0[4 * i] = 3;
                x0[4 * i + 1] = -1;
                x0[4 * i + 2] = 0;
            }
            return x0;
        }
        Vector<double> getAnswer()const{return Vector<double>(n, 0);}
    };
    struct Trig: public BaseF
    {//From Dennis & Schnabel, they give no solution, but x = 0 works
        int n;
        Trig(int theN = 2): n(theN) {}
        Vector<double> operator()(Vector<double> const& x)const
        {

            double cSum = 0;
            for(int i = 0; i < n; ++i) cSum += cos(x[i]);
            Vector<double> fx(n, n - cSum);
            for(int i = 0; i < n; ++i)
                fx[i] += (i + 1) * (1 - cos(x[i])) - sin(x[i]);
            return fx;
        }
        string name()const{return "Trig" + to_string(n);}
        Vector<double> getX0()const{return Vector<double>(n, 1.0/n);}
        Vector<double> getAnswer()const{return Vector<double>(n, 0);}
    };
    struct HelicalValley: public BaseF
    {//From Dennis & Schnabel
        Vector<double> operator()(Vector<double> const& x)const
        {
            Vector<double> fx(3);
            double q = 0.5/PI() * atan(x[1]/x[0]);
            if(x[0] < 0) q += 0.5;
            fx[0] = 10 * (x[2] - 10 * q);
            fx[1] = 10 * (sqrt(x[0] * x[0] + x[1] * x[1]) - 1);
            fx[2] = x[2];
            return fx;
        }
        string name()const{return "HelicalValley";}
        Vector<double> getX0()const
        {
            Vector<double> x0(3, 0);
            x0[0] = -1;
            return x0;
        }
        Vector<double> getAnswer()const{return -getX0();}
    };
    struct Wood: public BaseF
    {//From Dennis & Schnabel; More et al have different definition
        Vector<double> operator()(Vector<double> const& x)const
        {
            Vector<double> fx(4);
            double t0 = x[0] * x[0] - x[1], t01 = 1 - x[0],
                t1 = x[2] * x[2] - x[3], t11 = 1 - x[2],
                t2 = 1 - x[1], t3 = 1 - x[3];
            fx[0] = 100 * t0 * t0 + t01 * t01;
            fx[1] = 90 * t1 * t1 + t11 * t11;
            fx[2] = 10.1 * (t2 * t2 + t3 * t3);
            fx[3] = 19.8 * t2 * t3;
            return fx;
        }
        string name()const{return "Wood";}
        Vector<double> getX0()const
        {
            Vector<double> x0(4, 0);
            x0[0] = -3;
            x0[1] = -1;
            x0[2] = -3;
            x0[3] = -1;
            return x0;
        }
        Vector<double> getAnswer()const{return Vector<double>(4, 1);}
    };
    struct MultiArctan: public BaseF
    {//swings Newton without backtrack into inf per Kelley
        int n;
        MultiArctan(int theN = 2): n(theN) {}
        Vector<double> operator()(Vector<double> const& x)const
        {
            Vector<double> fx(n);
            for(int i = 0; i < n; ++i) fx[i] = atan(x[i]);
            return fx;
        }
        string name()const{return "MultiArctan" + to_string(n);}
        Vector<double> getX0()const{return Vector<double>(n, 10);}
        Vector<double> getAnswer()const{return Vector<double>(n, 0);}
    };
    struct LinearFFullRank: public BaseF
    {//From More et al
        int n;
        LinearFFullRank(int theN = 2): n(theN) {}
        Vector<double> operator()(Vector<double> const& x)const
        {
            double sum = 0;
            for(int i = 0; i < n; ++i) sum += x[i];
            Vector<double> fx(n);
            for(int i = 0; i < n; ++i) fx[i] = x[i] - 2 * sum/n - 1;
            return fx;
        }
        string name()const{return "LinearFFullRank" + to_string(n);}
        Vector<double> getX0()const{return Vector<double>(n, 1);}
        Vector<double> getAnswer()const{return -Vector<double>(n, 1);}
    };
    struct GulfRND: public BaseF
    {//From More et al; paper has mistake "mi" should be "-" not m * (i + 1)
        Vector<double> operator()(Vector<double> const& x)const
        {
            Vector<double> fx(3);
            for(int i = 0; i < fx.getSize(); ++i)
            {
                double ti = (i + 1.0)/100, yi = 25 + pow(-50 * log(ti), 2.0/3);
                fx[i] = exp(-pow(abs(yi - x[1]), x[2])/x[0]) - ti;
            }
            return fx;
        }
        string name()const{return "GulfRND";}
        Vector<double> getX0()const
        {
            Vector<double> x0(3);
            x0[0] = 5;
            x0[1] = 2.5;
            x0[2] = 0.15;
            return x0;
        }
        Vector<double> getAnswer()const{return getX0() * 10;}
    };
    struct BiggsExp6: public BaseF
    {//From More et al;
        Vector<double> operator()(Vector<double> const& x)const
        {
            Vector<double> fx(6);
            for(int i = 0; i < fx.getSize(); ++i)
            {
                double ti = (i + 1.0)/10,
                    yi = exp(-ti) - 5 * exp(-10 * ti) + 3 * exp(-4 * ti);
                fx[i] = x[2] * exp(-ti * x[0]) - x[3] * exp(-ti * x[1]) +
                    x[5] * exp(-ti * x[4]) - yi;
            }
            return fx;
        }
        string name()const{return "BiggsExp6";}
        Vector<double> getX0()const
        {
            Vector<double> x0(6, 1);
            x0[1] = 2;
            return x0;
        }
        Vector<double> getAnswer()const
        {
            Vector<double> a(6);
            a[0] = 1;
            a[1] = 10;
            a[2] = 1;
            a[3] = 5;
            a[4] = 4;
            a[5] = 3;
            return a;
        }
    };
    static int evalCount;
    struct MetaF
    {
        shared_ptr<BaseF> f;
        template<typename F> MetaF(shared_ptr<F> const& theF): f(theF){}
        Vector<double> operator()(Vector<double> const& x)const
        {
            ++evalCount;
            return (*f)(x);
        }
        string getName()const{return f->name();}
        virtual Vector<double> getX0()const{return f->getX0();}
        virtual Vector<double> getAnswer()const{return f->getAnswer();}
    };
    static Vector<MetaF> getFunctions()
    {
        Vector<MetaF> result;
        //result.append(MetaF(make_shared<TestFunc0::Func>()));
        //result.append(MetaF(make_shared<TestFunc1::Func>()));
        result.append(MetaF(make_shared<ExtendedRosenbrock>()));
        result.append(MetaF(make_shared<ExtendedPowellSingular>()));
        result.append(MetaF(make_shared<Trig>()));
        result.append(MetaF(make_shared<MultiArctan>()));
        result.append(MetaF(make_shared<HelicalValley>()));
        result.append(MetaF(make_shared<LinearFFullRank>()));
        result.append(MetaF(make_shared<Wood>()));
        result.append(MetaF(make_shared<GulfRND>()));
        result.append(MetaF(make_shared<BiggsExp6>()));
        result.append(MetaF(make_shared<ExtendedRosenbrock>(10)));
        result.append(MetaF(make_shared<ExtendedPowellSingular>(12)));
        result.append(MetaF(make_shared<Trig>(10)));
        result.append(MetaF(make_shared<MultiArctan>(10)));
        result.append(MetaF(make_shared<LinearFFullRank>(10)));
        result.append(MetaF(make_shared<ExtendedRosenbrock>(100)));
        result.append(MetaF(make_shared<ExtendedPowellSingular>(100)));
        result.append(MetaF(make_shared<Trig>(100)));
        result.append(MetaF(make_shared<MultiArctan>(100)));
        result.append(MetaF(make_shared<LinearFFullRank>(100)));
        return result;
    }
};
int TestFunctionsMultiD::evalCount = 0;

template<typename FUNCTION>
void printResultNonliearEq(pair<Vector<double>, double> const& result,
    FUNCTION const& f, Vector<Vector<string> >& matrix)
{
    int evalCount = TestFunctionsMultiD::evalCount;
    double relAbsXErrorDigits = log10(max(numeric_limits<double>::epsilon(),
        normInf(result.first - f.getAnswer())));
    DEBUG(relAbsXErrorDigits);
    matrix.lastItem().append(toStringDouble(relAbsXErrorDigits));
    double absEstErrorDigits = log10(max(numeric_limits<double>::epsilon(),
        abs(result.second)));
    DEBUG(absEstErrorDigits);
    matrix.lastItem().append(toStringDouble(absEstErrorDigits));
    double absYErrorDigits = log10(max(numeric_limits<double>::epsilon(),
        normInf(f(result.first))));
    DEBUG(absYErrorDigits);
    matrix.lastItem().append(toStringDouble(absYErrorDigits));
    DEBUG(evalCount);
    matrix.lastItem().append(to_string(evalCount));
    TestFunctionsMultiD::evalCount = 0;
}
template<typename FUNCTION> void testNonlinearEqHelper(FUNCTION const& f,
    Vector<double> const& x0, Vector<Vector<string> >& matrix)
{
    DEBUG("Broyden QR");
    matrix.lastItem().append("Broyden QR");
    printResultNonliearEq(solveBroyden(f, x0), f, matrix);
    DEBUG("LMBroyden");
    matrix.lastItem().append("LMBroyden");
    printResultNonliearEq(solveLMBroyden(f, x0), f, matrix);
    DEBUG("BroydenLevy");
    matrix.lastItem().append("BroydenLevy");
    printResultNonliearEq(solveBroydenLevy(f, x0.getSize()), f, matrix);
    DEBUG("Opt");
    matrix.lastItem().append("Opt");
    printResultNonliearEq(solveByOptimization(f, x0), f, matrix);
    DEBUG("Hybrid");
    matrix.lastItem().append("Hybrid");
    printResultNonliearEq(hybridEquationSolve(f, x0.getSize()), f, matrix);
}

void testNonlinearEq()
{
    Vector<Vector<string> > matrix;
    Vector<TestFunctionsMultiD::MetaF> fs = TestFunctionsMultiD::getFunctions();
    for(int i = 0; i < fs.getSize(); ++i)
    {
        string name = fs[i].getName();
        DEBUG(name);
        matrix.append(Vector<string>());
        matrix.lastItem().append(name);
        testNonlinearEqHelper(fs[i], fs[i].getX0(), matrix);
    }
    int reportNumber = time(0);
    string filename = "reportNonlinearSolve" + to_string(reportNumber) + ".csv";
    createCSV(matrix, filename.c_str());
    Vector<string> names;
    names.append("XError");
    names.append("XEstimate");
    names.append("YError");
    names.append("Evals");
    createAugmentedCSVFiles(matrix, names, filename);
}

struct SolveTestFunctions1D
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
            return (*f)(x);
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
int SolveTestFunctions1D::evalCount = 0;


void debugEq1DResult(pair<double, double> const& result, double answer,
    Vector<Vector<string> >& matrix)
{
    DEBUG(normInf(result.first - answer));
    double relAbsXErrorDigits = log10(max(numeric_limits<double>::epsilon(),
        abs(answer - result.first)/max(1.0, answer)));
    DEBUG(relAbsXErrorDigits);
    double absYErrorDigits = log10(max(numeric_limits<double>::epsilon(),
        abs(result.second)));
    DEBUG(relAbsXErrorDigits);
    DEBUG(absYErrorDigits);
    DEBUG(SolveTestFunctions1D::evalCount);
    matrix.lastItem().append(toStringDouble(relAbsXErrorDigits));
    matrix.lastItem().append(toStringDouble(absYErrorDigits));
    matrix.lastItem().append(toStringDouble(SolveTestFunctions1D::evalCount));
    SolveTestFunctions1D::evalCount = 0;
}
template<typename FUNCTION> void testNonlinearEqHelper1D(FUNCTION const& f,
    double x0, Vector<Vector<string> >& matrix)
{
    DEBUG("ExpSearch");
    matrix.lastItem().append("ExpSearch");
    debugEq1DResult(exponentialSearch(f, x0), f.getAnswer(), matrix);
    DEBUG("Secant");
    matrix.lastItem().append("Secant");
    debugEq1DResult(solveSecant(f, x0), f.getAnswer(), matrix);
    DEBUG("SecantGlobal");
    matrix.lastItem().append("SecantGlobal");
    debugEq1DResult(solveSecantGlobal(f), f.getAnswer(), matrix);
}

void testNonlinearEq1D()
{
    Vector<Vector<string> > matrix;
    Vector<SolveTestFunctions1D::MetaF> fs = SolveTestFunctions1D::getFunctions();
    for(int i = 0; i < fs.getSize(); ++i)
    {
        string name = fs[i].getName();
        DEBUG(name);
        matrix.append(Vector<string>());
        matrix.lastItem().append(name);
        testNonlinearEqHelper1D(fs[i], fs[i].getX0(), matrix);
    }
    int reportNumber = time(0);
    string filename = "reportNonlinearEq1D" + to_string(reportNumber) + ".csv";
    createCSV(matrix, filename.c_str());
    Vector<string> names;
    names.append("XError");
    names.append("YError");
    names.append("Evals");
    createAugmentedCSVFiles(matrix, names, filename);
}

int main()
{
    testNonlinearEq();
    return 0;
    testNonlinearEq1D();//fails here check again
    return 0;

    return 0;
}
