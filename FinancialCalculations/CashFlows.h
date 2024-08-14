#ifndef CASH_FLOWS_H
#define CASH_FLOWS_H
#include "../Utils/Vector.h"
#include "../NumericalMethods/EquationSolving.h"
#include "../NumericalMethods/Differentiation.h"
#include "../NumericalMethods/Matrix.h"
#include "../RandomNumberGeneration/TimeSeries.h"
#include "../MachineLearning/Lasso.h"
#include "../HashTable/ChainingHashTable.h"
#include <memory>
namespace igmdk{

double calculatePriceGeneral(Vector<double> const& cashFlow,
    Vector<double> const& datesFrom0, double r)
{//r is discount rate for the period, with negative values allowed
    assert(isfinite(r) && r > -1 && cashFlow.getSize() > 0 &&
        datesFrom0.getSize() > 0);
    double result = 0;
    for(int i = 0; i < cashFlow.getSize(); ++i)
    {
        assert(isfinite(cashFlow[i]) && isfinite(datesFrom0[i]) &&
            datesFrom0[i] >= 0);
        result += cashFlow[i]/pow(1 + r, datesFrom0[i]);
    }
    return result;
}

double calculateYieldGeneral(Vector<double> const& cashFlow,
    Vector<double> const& datesFrom0, double price)
{//r is discount rate for the period, 0 index is earliest payment
    assert(isfinite(price) && cashFlow.getSize() > 0);
    //solve for r such that price equals to calculated price
    auto priceFunctor = [price, &cashFlow, &datesFrom0](double r)
        {return calculatePriceGeneral(cashFlow, datesFrom0, r) - price;};
    pair<double, double> result = exponentialSearch(priceFunctor, 0);
    return result.first;
}

Vector<double> getDatesFrom0(int n, double periodLength = 1,
    double timeToFirstPayment = 0)
{//offset to adjust date to start between payment periods
    assert(n > 0 && periodLength > 0 && periodLength <= 1 &&
        abs(timeToFirstPayment) < 1);
    Vector<double> datesFrom0;
    for(int i = 0; i < n; ++i)
        datesFrom0.append(timeToFirstPayment + i * periodLength);
    return datesFrom0;
}

double calculatePrice(Vector<double> const& cashFlow, double r)
{
    assert(isfinite(r) && cashFlow.getSize() > 0);
    return calculatePriceGeneral(cashFlow, getDatesFrom0(cashFlow.getSize()),r);
}

double calculateYield(Vector<double> const& cashFlow, double price)
{//r is discount rate for the period, 0 index is earliest payment
    assert(isfinite(price) && cashFlow.getSize() > 0);
    Vector<double> datesFrom0;
    for(int i = 0; i < cashFlow.getSize(); ++i) datesFrom0.append(i);
    return calculateYieldGeneral(cashFlow, getDatesFrom0(cashFlow.getSize()),
        price);
}

struct GeneralBond
{//both negative and positive numbers allowed
    double coupon, principal, timeToFirstPayment;
    int nPayments, nPaymentsPerYear;//last payment is principal + coupon
    GeneralBond(double theCoupon, double thePrincipal, int theNPayments,
        double theTimeToFirstPayment = 0, int theNPaymentsPerYear = 1):
        coupon(theCoupon), principal(thePrincipal), nPayments(theNPayments),
        timeToFirstPayment(theTimeToFirstPayment),
        nPaymentsPerYear(theNPaymentsPerYear)
    {
        assert(isfinite(theCoupon) && isfinite(thePrincipal) && theNPayments > 0
            && abs(theTimeToFirstPayment) < 1 && theNPaymentsPerYear > 0);
    }
    pair<Vector<double>, Vector<double> > toCashFlow()const
    {
        Vector<double> cashFlow(nPayments, coupon/nPaymentsPerYear);
        cashFlow[nPayments - 1] += principal;
        return {cashFlow, getDatesFrom0(cashFlow.getSize(),
            1.0/nPaymentsPerYear, timeToFirstPayment)};
    }
    double calculatePrice(double r)const
    {
        assert(isfinite(r));
        auto cashFlow = toCashFlow();
        return calculatePriceGeneral(cashFlow.first, cashFlow.second, r);
    }
    double calculateYield(double price)const
    {
        assert(price > 0);
        auto cashFlow = toCashFlow();
        return calculateYieldGeneral(cashFlow.first, cashFlow.second, price);
    }
    double calculateModifiedDuration(double price)const
    {//derivative of price with respect to r, evaluated at current yield
        double derivative = estimateDerivativeCD(
            [this](double r){return calculatePrice(r);}, calculateYield(price));
        //by conversion change sign for relative price
        return -derivative/price;
    }
    double calculateConvexity(double price)const
    {//2nd derivative of price with respect to r, evaluated at current yield
        double derivative2 = estimate2ndDerivativeCD(
            [this](double r){return calculatePrice(r);}, calculateYield(price));
        //by conversion for relative price
        return derivative2/price;
    }
};

struct Mortgage
{//also works for reverse mortgage
    double yield, principal;
    int nYears;//assume full years with first payment one month after
    Mortgage(double theYield, double thePrincipal, int theNYears):
        yield(theYield), principal(thePrincipal), nYears(theNYears)
    {
        assert(isfinite(theYield) && isfinite(thePrincipal) &&
            theNYears > 0);
    }
    pair<Vector<double>, Vector<double> > toCashFlow(double monthlyPayment)const
    {//first coupon due month after grant
        Vector<double> cashFlow(nYears * 12 + 1, monthlyPayment);
        cashFlow[0] = -principal;
        return {cashFlow, getDatesFrom0(cashFlow.getSize(), 1.0/12)};
    }
    double getMonthlyPayment() const
    {
        auto priceFunctor = [this](double monthlyPayment)
        {
            auto cashFlow = toCashFlow(monthlyPayment);
            return calculatePriceGeneral(cashFlow.first, cashFlow.second,
                yield);
        };
        pair<double, double> result = exponentialSearch(priceFunctor, 0);
        return result.first;
    }
};

}//end namespace
#endif
