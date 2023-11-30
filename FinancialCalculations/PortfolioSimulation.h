#ifndef PORTFOLIO_SIMULATION_H
#define PORTFOLIO_SIMULATION_H
#include "../Utils/Vector.h"
#include "../NumericalMethods/EquationSolving.h"
#include "../NumericalMethods/Differentiation.h"
#include "../NumericalMethods/Matrix.h"
#include "../RandomNumberGeneration/TimeSeries.h"
#include "../MachineLearning/Lasso.h"
#include "../HashTable/ChainingHashTable.h"
#include "Misc.h"
#include <memory>
namespace igmdk{

struct PortfolioSimulationResult
{
    double ruinChance;
    PercentileManager percentiles;
    PortfolioSimulationResult(Vector<double> const& values): percentiles(values)
    {
        int ruinCount = 0;
        for(int i = 0; i < values.getSize(); ++i) ruinCount += (values[i] <= 0);
        ruinChance = ruinCount * 1.0/values.getSize();
    }
    double getMedian()const{return percentiles.getPercentile(0.50);}
    void debug(double inflationFactor = 1) const
    {
        DEBUG(ruinChance);
        DEBUG(percentiles.getPercentile(0.05)/inflationFactor);
        DEBUG(percentiles.getPercentile(0.25)/inflationFactor);
        DEBUG(getMedian()/inflationFactor);
        DEBUG(percentiles.getPercentile(0.75)/inflationFactor);
        DEBUG(percentiles.getPercentile(0.95)/inflationFactor);
    }
    double riskFreeRank(double riskFreeReturn) const
    {//find position of risk-free-return with modified binary search
        int left = 0, right = percentiles.values.getSize() - 1;
        while(right - left > 1)
        {
            int middle = left + (right-left)/2;
            if(percentiles.values[middle] < riskFreeReturn) left = middle;
            else right = middle;
        }
        return right * 1.0/percentiles.values.getSize();
    }
};

class StockBondAsset
{
    double value;
    LognormalDistribution returns;
    static LognormalDistribution makeReturnDistribution(ReturnSpecifier const&
        returnSpecifier, double bondFraction)
    {
        MeanVariancePortfolio mvp(makeStockBondMVP(returnSpecifier));
        Vector<double> weights;
        weights.append(1 - bondFraction);
        weights.append(bondFraction);
        pair<double, double> ms = mvp.evaluate(weights);
        return LognormalDistribution(1 + ms.first, ms.second);
    }
public:
    StockBondAsset(ReturnSpecifier const& returnSpecifier = ReturnSpecifier(),
        double initialValue = 0, double bondFraction = 0): value(initialValue),
        returns(makeReturnDistribution(returnSpecifier, bondFraction)) {}
    void simulateStep(double netSavings)
        {value = value * returns.sample() + netSavings;}
    double getValue()const{return value;}
};
class RiskFreeAsset
{
    double value, riskFreeRate;
public:
    RiskFreeAsset(ReturnSpecifier const& returnSpecifier = ReturnSpecifier(),
        double initialValue = 0): value(initialValue),
        riskFreeRate(returnSpecifier.getRiskFreeRate()) {}
    void simulateStep(double netSavings)
        {value = value * (1 + riskFreeRate) + netSavings;}
    double getValue()const{return value;}
};

template<typename FINANCIAL_ASSET>
PortfolioSimulationResult performSimulation(FINANCIAL_ASSET const&
    initialFinancialAsset, Vector<double> const& netSavings,
    int nSimulations = 1000000)
{
    assert(netSavings.getSize() > 0 && nSimulations > 0);
    Vector<double> values;
    for(int i = 0; i < nSimulations; ++i)
    {
        FINANCIAL_ASSET finacialAsset = initialFinancialAsset;
        for(int step = 0; step < netSavings.getSize(); ++step)
        {//check if already ruined
            if(finacialAsset.getValue() < 0) break;
            finacialAsset.simulateStep(netSavings[step]);
        }//0 if ruined
        values.append(max(0.0, finacialAsset.getValue()));
    }
    return PortfolioSimulationResult(values);
}

Vector<double> generateDCASavings(double initialSavings, double growthRate,
    int term)
{//first one not growing, rest are
    Vector<double> netSavings;
    for(int i = 0; i < term; ++i)
        netSavings.append(initialSavings * pow(1 + growthRate, i));
    return netSavings;
}
Vector<double> generateRetirementExpenses(double initialExpenses,
    double growthRate, int term)
    {return generateDCASavings(initialExpenses, growthRate, term) * -1;}
Vector<double> lifetimeNetSavings(double initialSavings, double initialExpenses,
    double growthRate, int term, int retirementYear)
{//save until retirementYear
    Vector<double> netSavings = generateDCASavings(initialSavings, growthRate,
        retirementYear);
    //spend after that
    netSavings.appendVector(generateRetirementExpenses(initialExpenses * pow(1 +
        growthRate, retirementYear), growthRate, term - retirementYear));
    return netSavings;
}

}//end namespace
#endif
