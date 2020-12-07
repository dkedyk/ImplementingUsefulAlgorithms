#ifndef IGMDK_TEST_FUNCTIONS_1D_H
#define IGMDK_TEST_FUNCTIONS_1D_H

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <limits>
#include <memory> //for shared ptr
#include "../Utils/Vector.h"
#include "Integration.h"
#include "../ExternalMemoryAlgorithms/CSV.h"

namespace igmdk{

struct TestFunctions1D
{//for integral use http://www.emathhelp.net/calculators/calculus-2/definite-integral-calculator
    struct BaseF
    {
        virtual double operator()(double x)const = 0;
        virtual string name()const{return "";}
        virtual double getA()const{return -1;}
        virtual double getB()const{return 1;}
        virtual double deriv(double x)const//not needed in most cases
            {return numeric_limits<double>::quiet_NaN();}
        virtual double integralAB()const//not needed in most cases
            {return numeric_limits<double>::quiet_NaN();}
    };
    struct Step: public BaseF//discontinuity
    {
        double operator()(double x)const{return x > 0;}
        string name()const{return "Step";}
        double deriv(double x)const{return 0;}
        double integralAB()const{return 1;}
    };
    struct Abs: public BaseF//non-differentiable
    {
        double operator()(double x)const{return abs(x);}
        string name()const{return "Abs";}
        double deriv(double x)const{return x = 0 ? 0 : x < 0 ? -1 : 1;}
        double integralAB()const{return 1;}
    };
    struct Lin: public BaseF//exact for order 1
    {
        double operator()(double x)const{return x;}
        string name()const{return "Lin";}
        double deriv(double x)const{return 1;}
        double integralAB()const{return 0;}
    };
    struct Square: public BaseF//exact for order 2
    {
        double operator()(double x)const{return x * x;}
        string name()const{return "Square";}
        double deriv(double x)const{return 2 * x;}
        double integralAB()const{return 2.0/3;}
    };
    struct Cube: public BaseF//exact for order 3
    {
        double operator()(double x)const{return x * x * x;}
        string name()const{return "Cube";}
        double deriv(double x)const{return 3 * x * x;}
        double integralAB()const{return 0;}
    };
    struct Quad: public BaseF//exact for order 2
    {
        double operator()(double x)const{return x * x * x * x;}
        string name()const{return "Quad";}
        double deriv(double x)const{return 4 * x * x * x;}
        double integralAB()const{return 0.4;}
    };
    struct Exp: public BaseF//analytic
    {
        double operator()(double x)const{return exp(x);}
        string name()const{return "Exp";}
        double deriv(double x)const{return exp(x);}
        double integralAB()const
        {
            double e = exp(1);
            return e - 1/e;
        }
    };
    struct SqrtAbs: public BaseF
    {
        double operator()(double x)const{return sqrt(abs(x));}
        string name()const{return "SqrtAbs";}
        //NEED DERIV!
        double integralAB()const{return 4.0/3;}
    };
    struct Runge: public BaseF//differentiable nonlinear
    {//deriv not given - too complex
        double operator()(double x)const{return 1/(1 + 25 * x * x);}
        string name()const{return "Runge";}
        double getA()const{return -5;}
        double getB()const{return 5;}
        double integralAB()const{return 0.4 * atan(25);}
    };
    struct Log: public BaseF//analytic slow growth
    {
        double operator()(double x)const{return log(x);}
        string name()const{return "Log";}
        double deriv(double x)const{return 1/x;}
        double getA()const{return 1;}
        double getB()const{return 2;}
        double integralAB()const{return 2 * log(2) - 1;}
    };
    struct XSinXM1: public BaseF//continuous not Lipschitz
    {
        double operator()(double x)const
        {
            double temp = x * sin(1/x);
            return isfinite(temp) ? temp : 0;
        }
        string name()const{return "XSinXM1";}
        double deriv(double x)const
        {
            double temp = 2 * x * sin(1/x) - cos(1/x);
            return isfinite(temp) ? temp : 0;
        }
        //exact integrator struggles here - use 0.00000000000001 to 1
        double integralAB()const{return 2 * 0.378530017124161;}
    };
    struct Sin: public BaseF//periodic analytic
    {
        double operator()(double x)const{return sin(x);}
        string name()const{return "Sin";}
        double deriv(double x)const{return cos(x);}
        double integralAB()const{return 0;}
    };
    struct AbsIntegral: public BaseF//single continuous derivative
    {
        double operator()(double x)const{return (x > 0 ? 1 : -1) * x * x/2;}
        string name()const{return "AbsIntegral";}
        double deriv(double x)const{return abs(x);}
        double integralAB()const{return 0;}
    };
    struct Tanh: public BaseF//bounded range slow growth
    {
        double operator()(double x)const{return tanh(x);}
        string name()const{return "Tanh";}
        double deriv(double x)const{return 1 - tanh(x) * tanh(x);}
        double integralAB()const{return 0;}
    };
    struct NormalPDF: public BaseF//analytic thin tails
    {
        double operator()(double x)const{return exp(-x*x/2)/sqrt(2 * PI());}
        string name()const{return "NormalPDF";}
        double deriv(double x)const{return -x * (*this)(x);}
        double getA()const{return -10;}
        double getB()const{return 10;}
        double integralAB()const{return 1;}
    };
    struct DeltaPDF: public BaseF//non-diff continuous, 0 tails
    {
        double operator()(double x)const{return x > -1 && x < 1 ? 1 - abs(x) : 0;}
        string name()const{return "DeltaPDF";}
        double deriv(double x)const{return x > -1 && x < 1 ? -(x = 0 ? 0 : x < 0 ? -1 : 1) : 0;}
        double getA()const{return -20;}
        double getB()const{return 20;}
        double integralAB()const{return 1;}
    };
    struct Kahaner21: public BaseF//hard to notice bump
    {
        double operator()(double x)const{return 1/pow(cosh(10 * x - 2), 2) +
            1/pow(cosh(100 * x - 40), 4) + 1/pow(cosh(1000 * x - 600), 6);}
        string name()const{return "Kahaner21";}
        double getA()const{return 0;}
        double integralAB()const{return 0.210802735500549277;}
    };
    struct F575: public BaseF//singularity at 0
    {
        struct F575Helper
        {//singularity at 0
            double operator()(double x)const
                {return exp(-0.4 * x) * cos(2 * x)/pow(x, 0.7);}
        };
        SingularityWrapper<F575Helper> fh;
        double operator()(double x)const{return fh(x);}
        string name()const{return "F575";}
        double getA()const{return 0;}
        double getB()const{return 100;}
        double integralAB()const{return 2.213498276272980295056;}
    };
    static int evalCount;
    struct MetaF
    {
        shared_ptr<BaseF> f;
        template<typename F> MetaF(shared_ptr<F> const& theF): f(theF){}
        double operator()(double x)const
        {
            assert((*f).getA() <= x && x <= (*f).getB());//do nan not assert?
            ++evalCount;
            return (*f)(x);
        }
        double deriv(double x)const
        {
            assert((*f).getA() <= x && x <= (*f).getB());//do nan not assert?
            return f->deriv(x);
        }
        string getName()const{return f->name();}
        double integralAB()const{return f->integralAB();}
        double getA()const{return f->getA();}
        double getB()const{return f->getB();}
    };
    static Vector<MetaF> getFunctions()
    {
        Vector<MetaF> result;
        result.append(MetaF(make_shared<Step>()));
        result.append(MetaF(make_shared<Abs>()));
        result.append(MetaF(make_shared<Lin>()));
        result.append(MetaF(make_shared<Square>()));
        result.append(MetaF(make_shared<Cube>()));
        result.append(MetaF(make_shared<Quad>()));
        result.append(MetaF(make_shared<Exp>()));
        result.append(MetaF(make_shared<SqrtAbs>()));
        result.append(MetaF(make_shared<Runge>()));
        result.append(MetaF(make_shared<Log>()));
        result.append(MetaF(make_shared<XSinXM1>()));
        result.append(MetaF(make_shared<Sin>()));
        result.append(MetaF(make_shared<AbsIntegral>()));
        result.append(MetaF(make_shared<Tanh>()));
        result.append(MetaF(make_shared<NormalPDF>()));
        result.append(MetaF(make_shared<DeltaPDF>()));
        result.append(MetaF(make_shared<F575>()));
        return result;
    }
};
int TestFunctions1D::evalCount = 0;//HANDLE INTERGRATION FUNCS NEXT!

void printResult2(pair<double, double> const result, double answer, Vector<Vector<string> > & matrix)
{//assume make rounding error so at least 1 bit off even if 0 error
    double relAbsErrorDigits = log10(max(numeric_limits<double>::epsilon(),
        abs(answer - result.first)/max(1.0, answer)));
    DEBUG(relAbsErrorDigits);
    double relAbsEstimateDigits = isfinite(result.second) ? log10(max(numeric_limits<double>::epsilon(),
        result.second/max(1.0, answer))) : 0;
    DEBUG(relAbsEstimateDigits);
    DEBUG(TestFunctions1D::evalCount);
    //matrix.lastItem().append(toStringDouble(relResult));
    //matrix.lastItem().append(toStringDouble(result.first));
    matrix.lastItem().append(toStringDouble(relAbsErrorDigits));
    matrix.lastItem().append(toStringDouble(relAbsEstimateDigits));
    matrix.lastItem().append(toStringDouble(TestFunctions1D::evalCount));
    TestFunctions1D::evalCount = 0;
}

}//end namespace
#endif
