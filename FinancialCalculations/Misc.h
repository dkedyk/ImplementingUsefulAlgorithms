#ifndef MISC_H
#define MISC_H
#include "../Utils/Vector.h"
#include "../NumericalMethods/EquationSolving.h"
#include "../NumericalMethods/Differentiation.h"
#include "../NumericalMethods/Matrix.h"
#include "../RandomNumberGeneration/TimeSeries.h"
#include "../MachineLearning/Lasso.h"
#include "../HashTable/ChainingHashTable.h"
#include <memory>
namespace igmdk{

class ReturnSpecifier
{//returns for funds of securities predicted
public:
    //all rates are percentages, not multiples
    double stockReturn, stockStd, bondReturn, bondStd, stockBondCorrelation,
        riskFreeRate, taxRateCapitalGains, taxRateFederal, taxRateLocal,
        stockDividend, bondCoupon, bondTreasuryFraction, inflationRate;
    bool isInRange(double x, double a,
        double b = numeric_limits<double>::infinity()) const //a included b not
        {return x >= a && x < b;}
    double getTotalCapitalGainsTax()const
        {return taxRateCapitalGains + taxRateLocal;}
    double getTotalIncomeTax()const{return taxRateFederal + taxRateLocal;}
    double getStockAfterTaxDividend()const
        {return stockDividend * (1 - getTotalCapitalGainsTax());}
    double getBondAfterTaxCoupon()const
    {
        return bondCoupon * (1 - taxRateFederal -
            taxRateLocal * (1 - bondTreasuryFraction));
    }
    double getStockReturn()const
        {return stockReturn - stockDividend * getTotalCapitalGainsTax();}
    double getBondReturn()const
    {
        return (bondReturn - bondCoupon) * (1 - getTotalCapitalGainsTax()) +
            getBondAfterTaxCoupon();
    }
    double getRiskFreeRate()const
        {return riskFreeRate * (1 - getTotalIncomeTax());}
    double getStockStd()const{return stockStd;}
    double getBondStd()const{return bondStd;}
    double getStockBondCorrelation()const{return stockBondCorrelation;}
    double getInflationRate()const{return inflationRate;}
    //defaults need updating - calculated as of late 2023
    ReturnSpecifier(double theTaxRateCapitalGains = 0,
        double theTaxRateFederal = 0, double theTaxRateLocal = 0,
        double theBondTreasuryFraction = 0, double theStockReturn = 0.086,
        double theStockStd = 0.17, double theBondReturn = 0.055,
        double theBondStd = 0.09, double theStockBondCorrelation = 0,
        double theRiskFreeRate = 0.048, double theStockDividend = 0.02,
        double theBondCoupon = 0.025, double theInflationRate = 0.023):
        stockReturn(theStockReturn), stockStd(theStockStd),
        bondReturn(theBondReturn), bondStd(theBondStd),
        stockBondCorrelation(theStockBondCorrelation),
        riskFreeRate(theRiskFreeRate), taxRateCapitalGains(
        theTaxRateCapitalGains), taxRateFederal(theTaxRateFederal),
        taxRateLocal(theTaxRateLocal), stockDividend(theStockDividend),
        bondCoupon(theBondCoupon),bondTreasuryFraction(theBondTreasuryFraction),
        inflationRate(theInflationRate)
    {
        assert(isInRange(theTaxRateCapitalGains, 0, 1) &&
            isInRange(theTaxRateFederal, 0, 1) &&
            isInRange(theTaxRateLocal, 0, 1) &&
            isInRange(theStockReturn, 0, 1) &&
            isInRange(theStockStd, 0, 1) &&
            isInRange(theBondReturn, 0, 1) &&
            isInRange(theBondStd, 0, 1) &&
            isInRange(theStockBondCorrelation, -1, 1) &&
            isInRange(theRiskFreeRate, 0, 1) &&
            isInRange(theStockDividend, 0, 1) &&
                theStockDividend <= theStockReturn &&
            isInRange(theBondCoupon, 0, 1) && theBondCoupon <= theBondReturn &&
            isInRange(theBondTreasuryFraction, 0, 1) &&
            isInRange(theInflationRate, 0, 1)
        );
    }
};

pair<double, double> estimateLogNormalParameters(double mean, double stdev)
{//given observed mean and stdev of lognormal, estimate its parameters
    assert(mean > 0 && stdev > 0 && isfinite(mean) && isfinite(stdev));
//formula from https://en.wikipedia.org/wiki/Log-normal_distribution
    double m2 = mean * mean, v = stdev * stdev,
        mu = log(m2/sqrt(m2 + v)), q = sqrt(log(1 + v/m2));
    return {mu, q};
//some explanation of algebra with equivalent formula in
//https://www.johndcook.com/blog/2022/02/24/find-log-normal-parameters/
}
class LognormalDistribution
{
    double mu{0}, q{0};
public:
    LognormalDistribution(){}
    LognormalDistribution(double mean, double stdev)
    {//method of moments estimator, mean must be 1 + rate
        pair<double, double> muq = estimateLogNormalParameters(mean, stdev);
        mu = muq.first;
        q = muq.second;
    }
    LognormalDistribution(Vector<double> const& samples,
        bool skipInvalid = false)
    {//maximum likelihood estimator
        IncrementalStatistics s;
        for(int i = 0; i < samples.getSize(); ++i)
        {//0 and negative impossible
            if(samples[i] <= 0)
            {
                if(skipInvalid) continue;
                else assert(samples[i] > 0);
            }
            s.addValue(log(samples[i]));
        }
        mu = s.getMean();
        q = s.stdev();
    }
    double getMu()const{return mu;}
    double getQ()const{return q;}
    void setMu(double theMu){mu = theMu;}
    void setQ(double theQ){assert(theQ >= 0); q = theQ;}
    double getMean()const{return exp(mu + q*q/2);}
    double getStdev()const{return sqrt((exp(q*q) - 1) * exp(2 * mu + q*q));}
    double getMedian()const{return exp(mu);}
    double sample()const{return GlobalRNG().logNormal(mu, q);}
    LognormalDistribution& operator+=(LognormalDistribution const& rhs)
    {//parameters are additive in log
        mu += rhs.mu;
        q = sqrt(q * q + rhs.q * rhs.q);
        return *this;
    }
};

pair<double, double> estimateLogNormalParametersFromMedian(
    double median, double stdev)
{//given observed median and stdev of lognormal, estimate its parameters
    assert(median > 0 && stdev > 0 && isfinite(median) && isfinite(stdev));
    double mu = log(median);
    auto qFunctor = [mu, stdev](double q)
    {
        LognormalDistribution l;
        l.setMu(mu);
        l.setQ(q);
        return l.getStdev() - stdev;
    };
    //0 natural boundary for q, and stdev a clear upper bound
    double q = solveFor0(qFunctor, 0, stdev).first;
    return {mu, q};
}
double estimateArithmeticNominalStockReturn(double peRatio, double inflation,
    double stdev)//for current PE use IShares ITOT, inflation 10 year vs TIPS
{//assume e/p is geometric real return and lognormal returns
    assert(peRatio > 0 && stdev > 0);
    double geometricNominalReturn = 1 + 1/peRatio + inflation;
    //this is slightly more accurate than adding variance drag
    pair<double, double> muq =
        estimateLogNormalParametersFromMedian(geometricNominalReturn, stdev);
    LognormalDistribution l;
    l.setMu(muq.first);
    l.setQ(muq.second);
    return l.getMean();
}
double estimateArithmeticNominalBondReturn(double yield, double riskFreeRate,
    double stdev)//for current yield use IShares AGG and risk-free 10-year
{//assume 1/2 of credit risk premium realized
    assert(riskFreeRate > 0 && riskFreeRate < yield && stdev > 0);
    double geometricNominalReturn = 1 + riskFreeRate + (yield - riskFreeRate)/2;
    //this is slightly more accurate than adding variance drag
    pair<double, double> muq =
        estimateLogNormalParametersFromMedian(geometricNominalReturn, stdev);
    LognormalDistribution l;
    l.setMu(muq.first);
    l.setQ(muq.second);
    return l.getMean();
}

struct PercentileManager
{
    Vector<double> values;
    PercentileManager(Vector<double> const& theValues): values(theValues)
    {
        assert(theValues.getSize() > 0);
        quickSort(values.getArray(), values.getSize());
    }
    double getPercentile(double percentile)const
    {
        assert(percentile > 0 && percentile < 1);
        return values[int(percentile * values.getSize())];
    }
    double getTrimmedMean()const
        {return trimmedMean(values, 0.2, true);}
};

LognormalDistribution multigoalAdjustedLognormal(ReturnSpecifier const&
    returnSpecifier, Vector<pair<double, int>> const& goalWeights, int nYears)
{//action takes place at the end of year with growth and payments
    LognormalDistribution current, annual(1 +
        returnSpecifier.getStockReturn(), returnSpecifier.getStockStd());
    double remainingFraction = 1;
    for(int year = 1; year <= nYears; ++year)
    {//update with time
        current += annual;
        double goalAsW0Fraction = 0;
        for(int i = 0; i < goalWeights.getSize(); ++i)
        {//simple inefficient search
            if(year == goalWeights[i].second)
            {
                goalAsW0Fraction = goalWeights[i].first;
                break;
            }
        }
        if(goalAsW0Fraction != 0)
        {//recalculate based on mean-variance matching
            remainingFraction -= goalAsW0Fraction;
            double newMean = current.getMean() - goalAsW0Fraction *
                pow(1 + returnSpecifier.getRiskFreeRate(), year);
            assert(newMean > 0);
            current = LognormalDistribution(newMean, current.getStdev());
        }
    }
    assert(remainingFraction > 0);
    //scale for implicit leverage
    current = LognormalDistribution(current.getMean()/remainingFraction,
        current.getStdev()/remainingFraction);
    //adjust for geometric for single year
    current.setMu(current.getMu()/nYears);
    current.setQ(current.getQ()/sqrt(nYears));
    return current;
}

//Below is experimental code not presented in the book
double estimateStockBondSalesForRebalance(ReturnSpecifier const&
    returnSpecifier = ReturnSpecifier(), double bondFraction = 0.1,
    int nSimulations = 1000000)
{
    assert(0 < bondFraction && bondFraction < 1 && nSimulations > 0);
    LognormalDistribution stocks(1 + returnSpecifier.getStockReturn(),
        returnSpecifier.getStockStd()), bonds(1 +
        returnSpecifier.getBondReturn(), returnSpecifier.getBondStd());
    IncrementalStatistics percentSales;
    for(int i = 0; i < nSimulations; ++i)
    {
        double cash = (1 - bondFraction) *
            returnSpecifier.getStockAfterTaxDividend() +
            bondFraction * returnSpecifier.getBondAfterTaxCoupon(),
            stockValue = (1 - bondFraction) *
            (stocks.sample() - returnSpecifier.getStockAfterTaxDividend()),
            bondValue = bondFraction *
            (bonds.sample() - returnSpecifier.getBondAfterTaxCoupon()),
            total = stockValue + bondValue,
            saleAmount = max(0.0,fabs(bondValue - bondFraction * total) - cash);
        percentSales.addValue(saleAmount/total);
    }
    return percentSales.getMean();
}

double CRRAUtility(double wealth, double a)
{
    if(wealth <= 0) return -numeric_limits<double>::infinity();
    assert(a > 0);
    double temp = 1 - a;
    return abs(temp) < numeric_limits<double>::epsilon() ?
        log(wealth) : pow(wealth, temp)/temp;
}
double inverseCRRAUtility(double utility, double a)
{
    assert(a > 0);
    double temp = 1 - a;
    return abs(temp) < numeric_limits<double>::epsilon() ?
        exp(utility) : pow(utility * temp, 1.0/temp);
}
double expectedCRRACertaintyEquivalent(Vector<double> const& values, double a)
{
    IncrementalStatistics s;
    for(int i = 0; i < values.getSize(); ++i)
        s.addValue(CRRAUtility(values[i], a));
    return inverseCRRAUtility(s.getMean(), a);
}
double trimmedCRRACertaintyEquivalent(Vector<double> const& values, double a)
{
    IncrementalStatistics s;
    for(int i = 0.2 * values.getSize(); i < 0.8 * values.getSize(); ++i)
        s.addValue(CRRAUtility(values[i], a));
    return inverseCRRAUtility(s.getMean(), a);
}

Vector<double> getMSCIWorldReturns()
{
    double returns[] =
    {
        -2.36,4.16,22.80,16.13,-7.35,
        12.67,34.63,-42.19,11.66,20.95,
        10.84,15.23,33.99,-19.32,-16.21,
        -14.21,26.44,24.43,15.76,13.48,
        20.72,5.08,22.50,-5.23,18.28,
        -17.02,16.61,23.29,16.16,41.89,
        40.56,4.72,21.93,9.71,-4.79,
        25.67,10.95,16.52,0.68,13.40,
        -25.47,-15.24,22.48,18.36,-3.09
    };
    Vector<double> result;
    for(int i = 0; i < sizeof(returns)/sizeof(returns[0]); ++i)
        result.append(returns[i]/100.0);
    return result;
}

}//end namespace
#endif
