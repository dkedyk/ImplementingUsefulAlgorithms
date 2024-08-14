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
#include "Annuity.h"
#include <memory>
#include "../NumericalMethods/EquationSolving.h"
namespace igmdk{

struct PortfolioSimulationResult
{
    double ruinChance;
    PercentileManager percentiles, ruinPercentiles;
    PortfolioSimulationResult(Vector<double> const& values,
        Vector<double> const& ruinValues): percentiles(values),
        ruinPercentiles(ruinValues.getSize() > 0 ? ruinValues :
        Vector<double>(1, 0))
    {
        int ruinCount = 0;
        for(int i = 0; i < values.getSize(); ++i) ruinCount += (values[i] <= 0);
        ruinChance = ruinCount * 1.0/values.getSize();
    }
    double getMedian()const{return percentiles.getPercentile(0.50);}
    void debug(double inflationFactor = 1) const
    {
        if(ruinChance > 0)
        {
            DEBUG(ruinChance);
            double targetPercentile = 0.05/ruinChance;
            if(targetPercentile <= 1)
            {
                double percentile5TimeToRuin =
                    ruinPercentiles.getPercentile(targetPercentile);
                DEBUG(percentile5TimeToRuin);
            }
        }
        else
            DEBUG(expectedCRRACertaintyEquivalent(percentiles.values, 3)/
                inflationFactor);
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
    void join(PortfolioSimulationResult const& other)
    {
        percentiles.join(other.percentiles);
    }
};

class StockBondAsset
{
    double value;
    LognormalDistribution returns;
public:
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
    StockBondAsset(ReturnSpecifier const& returnSpecifier = ReturnSpecifier(),
        double initialValue = 0, double bondFraction = 0): value(initialValue),
        returns(makeReturnDistribution(returnSpecifier, bondFraction)) {}
    void simulateStep(double netSavings)
        {value = value * returns.sample() + netSavings;}
    double getValue()const{return value;}
    void setValue(double newValue){value = newValue;}
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
    void setValue(double newValue){value = newValue;}
};

template<typename FINANCIAL_ASSET>
PortfolioSimulationResult performActuarialSimulation(FINANCIAL_ASSET const&
    initialFinancialAsset, Vector<double> const& netSavings,
    Vector<double> const& nextYearSurvivalProbabilities,
    int nSimulations = 1000000)
{
    assert(netSavings.getSize() > 0 && nSimulations > 0 &&
        netSavings.getSize() == nextYearSurvivalProbabilities.getSize());
    Vector<double> values, ruinValues;
    for(int i = 0; i < nSimulations; ++i)
    {
        FINANCIAL_ASSET financialAsset = initialFinancialAsset;
        for(int step = 0; step < netSavings.getSize(); ++step)
        {//check if already ruined
            if(financialAsset.getValue() < 0)
            {
                ruinValues.append(step);
                break;
            }
            if(nextYearSurvivalProbabilities[step] < GlobalRNG().uniform01())
            {//died, assume continuation till the end with no cash flows to
             //calculate normalized inheritance
                for(; step < netSavings.getSize(); ++step)
                    financialAsset.simulateStep(0);
                break;
            }
            financialAsset.simulateStep(netSavings[step]);
        }//0 if ruined
        values.append(max(0.0, financialAsset.getValue()));
    }
    return PortfolioSimulationResult(values, ruinValues);
}
template<typename FINANCIAL_ASSET, typename SURVIVAL_ESTIMATOR>
PortfolioSimulationResult performActuarialSimulation(FINANCIAL_ASSET const&
    initialFinancialAsset, Vector<double> const& netSavings,
    SURVIVAL_ESTIMATOR const& e, int nSimulations = 1000000)
{
    Vector<double> nextYearSurvivalProbabilities;
    for(int i = 0; i < netSavings.getSize(); ++i)
        nextYearSurvivalProbabilities.append(
            getFutureSurvivalProbability(e, i));
    return performActuarialSimulation(initialFinancialAsset, netSavings,
        nextYearSurvivalProbabilities, nSimulations);
}
template<typename FINANCIAL_ASSET>
PortfolioSimulationResult performSimulation(FINANCIAL_ASSET const&
    initialFinancialAsset, Vector<double> const& netSavings,
    int nSimulations = 1000000)
{
    return performActuarialSimulation(initialFinancialAsset, netSavings,
        Vector<double>(netSavings.getSize(), 1), nSimulations);
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
    double growthRate, int term, int delay = 0)
{
    assert(delay >= 0 && delay < term);
    Vector<double> result(delay, 0);
    result.appendVector(generateDCASavings(initialExpenses *
        pow(1 + growthRate, delay), growthRate, term - delay) * -1);
    return result;
}

template<typename FINANCIAL_ASSET>
double findMaximumExpenseRatio(FINANCIAL_ASSET const& financialAsset,
    double growthRate, int term, double maximumRuinProbability = 0.05)
{
    auto f = [&](double expenseRatio)->double
    {
        Vector<double> retirementExpenses = generateRetirementExpenses(
            expenseRatio * 100, growthRate, term);
        PortfolioSimulationResult result =
            performSimulation(financialAsset, retirementExpenses);
        return result.ruinChance - maximumRuinProbability;
    };//start with 0.1% and take 5% steps at first, accuracy 10^-4
    pair<double, double> result = exponentialSearch1Sided(f, 0.001, 0.05,
        0.0001);
    return result.first;
}

template<typename FINANCIAL_ASSET>
PortfolioSimulationResult performSequentialAssetSimulation(
    Vector<FINANCIAL_ASSET> const& initialFinancialAssets,
    Vector<double> const& netSavings, int nSimulations = 1000000)
{
    assert(netSavings.getSize() > 0 && nSimulations > 0 &&
        initialFinancialAssets.getSize() == netSavings.getSize());
    Vector<double> values, ruinValues;
    for(int i = 0; i < nSimulations; ++i)
    {
        Vector<FINANCIAL_ASSET> financialAssets = initialFinancialAssets;
        int step = 0;
        for(; step < netSavings.getSize(); ++step)
        {//check if already ruined
            if(financialAssets[step].getValue() < 0)
            {
                ruinValues.append(step);
                break;
            }
            financialAssets[step].simulateStep(netSavings[step]);
            if(step + 1 < financialAssets.getSize())
                financialAssets[step + 1].setValue(
                    financialAssets[step].getValue());
        }//0 if ruined
        values.append(max(0.0, financialAssets.lastItem().getValue()));
    }
    return PortfolioSimulationResult(values, ruinValues);
}
Vector<StockBondAsset> makeSequantialAssets(ReturnSpecifier const&
    returnSpecifier, Vector<double> const& bondFractions, double initialValue)
{
    Vector<StockBondAsset> result;
    for(int i = 0; i < bondFractions.getSize(); ++i)
        result.append(StockBondAsset(returnSpecifier, i == 0 ? initialValue : 0,
            bondFractions[i]));
    return result;
}
struct CRRA3Functor
{
    double initialValue;
    ReturnSpecifier const& returnSpecifier;
    Vector<double> const& netSavings;
    double finalHomeValue;
    double operator()(Vector<double> const& bondFractions) const
    {
        PortfolioSimulationResult result = performSequentialAssetSimulation(
            makeSequantialAssets(returnSpecifier, bondFractions, initialValue),
            netSavings);
        double inflationFactor = pow(1 + returnSpecifier.getInflationRate(),
            netSavings.getSize());
        //add final home value
        result.percentiles.values +=
            Vector<double>(result.percentiles.values.getSize(), finalHomeValue);
        return expectedCRRACertaintyEquivalent(result.percentiles.values, 3)/
            inflationFactor;
    }
};
pair<Vector<double>, double> findOptimalBondFractions(int term,
    ReturnSpecifier const& returnSpecifier, double initialValue,
    double initialIncome, double homeValue = 0, double maxBondFraction = 0.4,
    double stepDelta = 0.05, int maxEvals = 300)
{
    Vector<double> dcaSavings = generateDCASavings(initialIncome,
        returnSpecifier.getInflationRate(), term);
    CRRA3Functor functor = {initialValue, returnSpecifier, dcaSavings,
        homeValue * pow(1 + returnSpecifier.getInflationRate(), term)};
    //first pick best constant allocation
    double bestScore = 0;
    Vector<double> bestbondFractions(term, 0);
    for(double bondFraction = 0; bondFraction <= maxBondFraction + 0.001;
        bondFraction += stepDelta)
    {
        Vector<double> nextBondFractions(term, bondFraction);
        double nextScore = functor(nextBondFractions);
        --maxEvals;
        if(nextScore > bestScore)
        {//2nd eval to make sure
            nextScore = functor(nextBondFractions);
            --maxEvals;
            if(nextScore > bestScore)
            {
                bestScore = nextScore;
                bestbondFractions = nextBondFractions;
            }
        }
    }
    //do random search with monotonicity for half of budget
    int nPartitions = 1 + (maxBondFraction + 0.001)/stepDelta;
    assert(nPartitions >= 2);
    for(int i = 0; i < maxEvals/2; ++i)
    {
        Vector<double> nextBondFractions(term, 0);
        //each index is last item of prev partition
        Vector<int> splitIndexes;
        for(int j = 0; j < nPartitions - 1; ++j)
            splitIndexes.append(GlobalRNG().mod(term));
        quickSort(splitIndexes.getArray(), splitIndexes.getSize());
        int choice = 0;
        for(int j = 0; j < term; ++j)
        {
            nextBondFractions[j] = choice * stepDelta;
            if(choice < splitIndexes.getSize() && splitIndexes[choice] == j)
                ++choice;//advance if hit partition bound
        }
        double nextScore = functor(nextBondFractions);
        if(nextScore > bestScore)
        {//2nd eval to make sure
            nextScore = functor(nextBondFractions);
            if(nextScore > bestScore)
            {
                bestScore = nextScore;
                bestbondFractions = nextBondFractions;
            }
        }
    }
    //finish with local search
    for(int i = 0; i < maxEvals/2; ++i)
    {//replace random index by different fraction
        Vector<double> nextBondFractions = bestbondFractions;
        int changedAt = 0;
        do
        {
            int j = GlobalRNG().mod(term);
            int choice = GlobalRNG().mod(nPartitions);
            nextBondFractions[j] = choice * stepDelta;
            changedAt = j;
            //avoid null move
        }while (nextBondFractions == bestbondFractions);
        //monotonic correction backward
        for(int j = changedAt - 1; j >= 0; --j) nextBondFractions[j] =
            min(nextBondFractions[j], nextBondFractions[changedAt]);
        //monotonic correction forward
        for(int j = changedAt + 1; j < term; ++j) nextBondFractions[j] =
            max(nextBondFractions[j], nextBondFractions[changedAt]);
        double nextScore = functor(nextBondFractions);
        if(nextScore > bestScore)
        {//2nd eval to make sure
            nextScore = functor(nextBondFractions);
            if(nextScore > bestScore)
            {
                bestScore = nextScore;
                bestbondFractions = nextBondFractions;
            }
        }
    }
    return {bestbondFractions, bestScore};
}

}//end namespace
#endif
