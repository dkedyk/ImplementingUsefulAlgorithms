#include <cassert>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <limits>
#include <memory> //for shared ptr
#include "ConjugateGradientOptimization.h"
#include "DSRS.h"
#include "BFGS.h"
#include "../RandomNumberGeneration/Statistics.h"
#include "../Utils/DEBUG.h"
#include "../ExternalMemoryAlgorithms/CSV.h"
using namespace std;
using namespace igmdk;

struct TestFunctionsMin
{
    struct BaseF
    {
        virtual Vector<double> operator()(Vector<double> const& x)const = 0;
        virtual string name()const = 0;
        virtual Vector<double> getX0()const = 0;
        virtual Vector<double> getAnswer()const = 0;
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
                x0[4 * i + 2] = 1;
            }
            return x0;
        }
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
    struct VariableDimensionF: public BaseF
    {//From More et al
        int n;
        VariableDimensionF(int theN = 2): n(theN) {}
        Vector<double> operator()(Vector<double> const& x)const
        {
            Vector<double> fx(n + 2);
            for(int i = 0; i < n; ++i)
            {
                fx[i] = x[i] - 1;
                fx[n] += (i + 1) * fx[i];
            }
            fx[n + 1] = fx[n] * fx[n];
            return fx;
        }
        string name()const{return "VariableDimensionF" + to_string(n);}
        Vector<double> getX0()const
        {
            Vector<double> x0(n, 1);
            for(int i = 0; i < n; ++i) x0[i] -= (i + 1.0)/n;
            return x0;
        }
        Vector<double> getAnswer()const{return Vector<double>(n, 1);}
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
    struct BrownBadScaled: public BaseF
    {//From More et al
        Vector<double> operator()(Vector<double> const& x)const
        {
            Vector<double> fx(3);
            fx[0] = x[0] - 1000000;
            fx[1] = x[1] - 2.0/1000000;
            fx[2] = x[0] * x[1] - 2;
            return fx;
        }
        string name()const{return "BrownBadScaled";}
        Vector<double> getX0()const{return Vector<double>(2, 1);}
        Vector<double> getAnswer()const
        {
            Vector<double> x0(2, 0);
            x0[0] = 1000000;
            x0[1] = 2.0/1000000;
            return x0;
        }
    };
    struct Beale: public BaseF
    {//From More et al
        Vector<double> operator()(Vector<double> const& x)const
        {
            Vector<double> fx(3);
            fx[0] = 1.5;
            fx[1] = 2.25;
            fx[2] = 2.625;
            for(int i = 0; i < fx.getSize(); ++i)
                fx[i] -= x[0] * (1 - pow(x[1], i + 1));
            return fx;
        }
        string name()const{return "Beale";}
        Vector<double> getX0()const{return Vector<double>(2, 1);}
        Vector<double> getAnswer()const
        {
            Vector<double> x0(2, 0);
            x0[0] = 3;
            x0[1] = 0.5;
            return x0;
        }
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
        double operator()(Vector<double> const& x)const
        {
            ++evalCount;
            return norm((*f)(x));
        }
        string getName()const{return f->name();}
        Vector<double> getX0()const{return f->getX0();}
        pair<Vector<double>, double> getAnswer()const
            {return make_pair(f->getAnswer(), norm((*f)(f->getAnswer())));}
    };
    static Vector<MetaF> getFunctions()
    {
        Vector<MetaF> result;
        result.append(MetaF(make_shared<ExtendedRosenbrock>()));
        result.append(MetaF(make_shared<ExtendedPowellSingular>()));
        result.append(MetaF(make_shared<HelicalValley>()));
        result.append(MetaF(make_shared<VariableDimensionF>()));
        result.append(MetaF(make_shared<LinearFFullRank>()));
        result.append(MetaF(make_shared<BrownBadScaled>()));
        result.append(MetaF(make_shared<Beale>()));
        result.append(MetaF(make_shared<BiggsExp6>()));
        result.append(MetaF(make_shared<ExtendedRosenbrock>(10)));
        result.append(MetaF(make_shared<ExtendedPowellSingular>(12)));
        result.append(MetaF(make_shared<VariableDimensionF>(10)));
        result.append(MetaF(make_shared<LinearFFullRank>(10)));
        result.append(MetaF(make_shared<ExtendedRosenbrock>(30)));
        result.append(MetaF(make_shared<ExtendedPowellSingular>(32)));
        result.append(MetaF(make_shared<VariableDimensionF>(30)));
        result.append(MetaF(make_shared<LinearFFullRank>(30)));
        result.append(MetaF(make_shared<ExtendedRosenbrock>(100)));
        result.append(MetaF(make_shared<ExtendedPowellSingular>(100)));
        result.append(MetaF(make_shared<VariableDimensionF>(100)));
        result.append(MetaF(make_shared<LinearFFullRank>(100)));
        //large D functions
        result.append(MetaF(make_shared<ExtendedRosenbrock>(1000)));
        result.append(MetaF(make_shared<ExtendedPowellSingular>(1000)));
        result.append(MetaF(make_shared<VariableDimensionF>(1000)));
        result.append(MetaF(make_shared<LinearFFullRank>(1000)));
        result.append(MetaF(make_shared<ExtendedRosenbrock>(10000)));
        result.append(MetaF(make_shared<ExtendedPowellSingular>(10000)));
        result.append(MetaF(make_shared<VariableDimensionF>(10000)));
        result.append(MetaF(make_shared<LinearFFullRank>(10000)));
        return result;
    }
};
int TestFunctionsMin::evalCount = 0;

template<typename TEST_SET, typename FUNCTION> void debugResultHelper(
    pair<Vector<double>, double> const& result,
    FUNCTION const& f, Vector<Vector<string> > & matrix, int start)
{
    //for(int i = 0; i < result.first.getSize(); ++i) DEBUG(result.first[i]);
    //DEBUG(result.second);
    //DEBUG(f(result.first));
    double timediff = 1.0 * (clock() - start)/CLOCKS_PER_SEC;
    pair<Vector<double>, double> answer = f.getAnswer();
    double eps = numeric_limits<double>::epsilon(),
        relAbsXNorm = max(eps, normInf(result.first - answer.first)/
        max(1.0, normInf(answer.first))),
        relAbsYNorm = max(eps, abs(result.second - answer.second)/
            max(1.0, abs(answer.second)));
    DEBUG(relAbsXNorm);
    DEBUG(relAbsYNorm);
    DEBUG(TEST_SET::evalCount);
    matrix.lastItem().append(toStringDouble(relAbsXNorm));
    matrix.lastItem().append(toStringDouble(relAbsYNorm));
    matrix.lastItem().append(to_string(TEST_SET::evalCount));
    TEST_SET::evalCount = 0;
    double gradNorm = norm(estimateGradientCD(result.first, f)),
        normalizedGradNorm = gradNorm/max(1.0, abs(result.second));
    TEST_SET::evalCount = 0;
    DEBUG(normalizedGradNorm);
    matrix.lastItem().append(toStringDouble(normalizedGradNorm));
    DEBUG(timediff);
    matrix.lastItem().append(toStringDouble(timediff));
}

template<typename POINT, typename FUNCTION> void debugResultNew(
    pair<POINT, double> const& result,
    FUNCTION const& f, Vector<Vector<string> > & matrix, int start)
{
    debugResultHelper<TestFunctionsMin>(result, f, matrix, start);
}

//for unimodal 2 works, for powell probably just lucky step choice gives it good result!
template<typename FUNCTION> void testAllSolvers(FUNCTION const& f, Vector<Vector<string> >& matrix)
{
    GradientFunctor<FUNCTION> g(f);
    DirectionalDerivativeFunctor<FUNCTION> dd(f);
    int D = f.getX0().getSize(), start = 0;
    DEBUG("BFGSMinimizeProper");
    matrix.lastItem().append("BFGSMinimizeProper");
    debugResultNew(BFGSMinimizeProper(f.getX0(), f, g, dd), f, matrix, start);
    DEBUG("ConjugateGradientPRP+");
    matrix.lastItem().append("ConjugateGradientPRP+");
    start = clock();
    debugResultNew(conjugateGradient(f.getX0(), f, g, dd), f, matrix, start);
    DEBUG("ConjugateGradientHZ");
    matrix.lastItem().append("ConjugateGradientHZ");
    start = clock();
    debugResultNew(conjugateGradient(f.getX0(), f, g, dd, 1000000, highPrecEps, "HZ"), f, matrix, start);
    DEBUG("DSRS");
    debugResultNew(DSRS(f.getX0(), f), f, matrix, start);
}

void createMinReport(string const& prefix,
    Vector<Vector<string> > const& matrix, int nRepeats)
{
    int reportNumber = time(0);
    string filename = prefix + to_string(reportNumber) + ".csv";
    createCSV(matrix, filename.c_str());
    Vector<string> names;
    names.append("XError");
    names.append("YError");
    names.append("NEvals");
    names.append("ScaledGradNorm");
    names.append("TimeSeconds");
    createAugmentedCSVFiles(matrix, names, filename, nRepeats);
}

void testAllFunctions()
{
    Vector<Vector<string> > matrix;
    string name;
    Vector<TestFunctionsMin::MetaF> fs = TestFunctionsMin::getFunctions();
    for(int i = 0; i < fs.getSize(); ++i)
    {
        string name = fs[i].getName();
        DEBUG(name);
        int D = fs[i].getX0().getSize();
        if(D >= 1000)
        {
            DEBUG("large scale case");
            continue;
        }
        matrix.append(Vector<string>());
        matrix.lastItem().append(name);
        testAllSolvers(fs[i], matrix);
    }
    createMinReport("reportMin", matrix, 1);
}

int main()
{
    testAllFunctions();
    return 0;
}
