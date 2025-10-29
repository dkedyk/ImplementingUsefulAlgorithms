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
    PercentileManager percentiles, ruinPercentiles, maxDrawdownPercentiles;
    PortfolioSimulationResult(Vector<double> const& values,
        Vector<double> const& ruinValues,
        Vector<double> const& maxDrawdownValues): percentiles(values),
        ruinPercentiles(ruinValues.getSize() > 0 ? ruinValues :
        Vector<double>(1, 0)), maxDrawdownPercentiles(maxDrawdownValues)
    {
        int ruinCount = 0;
        for(int i = 0; i < values.getSize(); ++i) ruinCount += (values[i] <= 0);
        ruinChance = ruinCount * 1.0/values.getSize();
    }
    double getMedian()const{return percentiles.getPercentile(0.50);}
    void debug() const
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
        {
            DEBUG(expectedCRRACertaintyEquivalent(percentiles.values, 3));
            DEBUG(expectedCRRACertaintyEquivalent(percentiles.values, 6));
        }
        DEBUG(percentiles.getPercentile(0.05));
        DEBUG(percentiles.getPercentile(0.25));
        DEBUG(getMedian());
        DEBUG(percentiles.getPercentile(0.75));
        DEBUG(percentiles.getPercentile(0.95));
        DEBUG(maxDrawdownPercentiles.getPercentile(0.95));
    }
    double riskFreeRank(double realRiskFreeReturn) const
    {//find position of risk-free-return with modified binary search
        int left = 0, right = percentiles.values.getSize() - 1;
        while(right - left > 1)
        {
            int middle = left + (right-left)/2;
            if(percentiles.values[middle] < realRiskFreeReturn) left = middle;
            else right = middle;
        }
        return right * 1.0/percentiles.values.getSize();
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
        return LognormalDistribution(1 + ms.first -
            returnSpecifier.getInflationRate(), ms.second);
    }
    StockBondAsset(ReturnSpecifier const& returnSpecifier = ReturnSpecifier(),
        double initialValue = 0, double bondFraction = 0): value(initialValue),
        returns(makeReturnDistribution(returnSpecifier, bondFraction)) {}
    void simulateStep(double netSavings, int step)
        {value = value * returns.sample() + netSavings;}
    double getValue()const{return value;}
    void setValue(double newValue){value = newValue;}
};

class RiskFreeAsset
{
    double value, realRiskFreeRate;
public:
    RiskFreeAsset(ReturnSpecifier const& returnSpecifier = ReturnSpecifier(),
        double initialValue = 0): value(initialValue), realRiskFreeRate(
        returnSpecifier.getRealRiskFreeRate())
    {}
    void simulateStep(double netSavings, int step)
        {value = value * (1 + realRiskFreeRate) + netSavings;}
    double getValue()const{return value;}
    void setValue(double newValue){value = newValue;}
};
template<typename RISKY_ASSET> class MixedAsset
{
    RISKY_ASSET riskyAsset;
    RiskFreeAsset riskFreeAsset;
    double riskFreeFraction;
    bool isUpRebalanced;
public:
    MixedAsset(RISKY_ASSET const& theRiskyAsset,
        ReturnSpecifier const& returnSpecifier =
        ReturnSpecifier(), double initialRiskFreeValue = 0,
        double theRiskFreeFraction = 0, bool theIsUpRebalanced = false):
        riskyAsset(theRiskyAsset), riskFreeAsset(returnSpecifier,
        initialRiskFreeValue), riskFreeFraction(theRiskFreeFraction),
        isUpRebalanced(theIsUpRebalanced) {}
    void simulateStep(double netSavings, int step)
    {
        riskyAsset.simulateStep(netSavings * (1 - riskFreeFraction), step);
        riskFreeAsset.simulateStep(netSavings * riskFreeFraction, step);
        if(!isUpRebalanced || riskFreeAsset.getValue() <
           getValue() * riskFreeFraction) setValue(getValue());
    }
    double getValue()const
        {return riskyAsset.getValue() + riskFreeAsset.getValue();}
    void setValue(double newValue)
    {
        riskyAsset.setValue(newValue * (1 - riskFreeFraction));
        riskFreeAsset.setValue(newValue * riskFreeFraction);
    }
    double getRiskFreeFraction()const{return riskFreeFraction;}
};

MixedAsset<StockBondAsset> makeMVOptimalRiskyAsset(double value,
    double stockTarget, double tangency,
    ReturnSpecifier const& returnSpecifier = ReturnSpecifier())
{//mix with risk-free to get target stock
    assert(stockTarget >= 0 && stockTarget <= 1 && tangency >= 0 && tangency);
    double riskFreeFraction = tangency <= stockTarget ? 0 :
        1 - stockTarget/tangency;
    return MixedAsset<StockBondAsset>(StockBondAsset(returnSpecifier,
        value * (1 - riskFreeFraction), 1 - max(stockTarget, tangency)),
        returnSpecifier, value * riskFreeFraction, riskFreeFraction, false);
}
//construct from a raw allocation such as 25% stocks/25% bonds/25% risk-free/
//25% up-rebalanced risk-free, where the weight add up to 1
MixedAsset<MixedAsset<StockBondAsset>> makeGeneralAsset(double initialValue,
    double bondFraction, double riskFreeFraction,
    double riskFreeFractionUpRebalanced = 0,
    ReturnSpecifier const& returnSpecifier = ReturnSpecifier())
{
    assert(initialValue > 0 && bondFraction + riskFreeFraction +
        riskFreeFractionUpRebalanced <= 1);
    //convert weights to recursive
    riskFreeFraction /= 1 - riskFreeFractionUpRebalanced;
    bondFraction /= 1 - riskFreeFraction;
    double nondedicatedInitialValue =
        initialValue * (1 - riskFreeFractionUpRebalanced);
    return MixedAsset<MixedAsset<StockBondAsset>>(MixedAsset<StockBondAsset>(
        StockBondAsset(returnSpecifier, nondedicatedInitialValue *
        (1 - riskFreeFraction), bondFraction), returnSpecifier,
        nondedicatedInitialValue * riskFreeFraction, riskFreeFraction, false),
        returnSpecifier, initialValue * riskFreeFractionUpRebalanced,
        riskFreeFractionUpRebalanced, true);
}

struct TotalAssetPolicy
{
    double totalStockTarget, nonFinancialStockFraction, netSavings,
        realRiskFreeRate, minStock, maxStock;
    int finalStep;
    double getNonFinancialPV(int step) const
    {
        assert(step <= finalStep);
        return ordinaryAnnuityPV(netSavings, finalStep - step,
            realRiskFreeRate);
    }
public:
    double operator()(double value, int step) const
    {//formula from the Smart Investing book
        if(value <= 0) return 0;
        double financialAssetFraction = value/(value + getNonFinancialPV(step)),
            financialStockFraction = (totalStockTarget -
            nonFinancialStockFraction * (1 - financialAssetFraction))/
            financialAssetFraction;
        return max(minStock, min(maxStock, financialStockFraction));
    }
};
template<typename POLICY> class DynamicallyRebalancedAsset
{
    double value;
    double tangency;
    POLICY policy;
    ReturnSpecifier returnSpecifier;
    double riskFreeUpRelancedFraction;
public:
    DynamicallyRebalancedAsset(double theValue, double theTangency,
        POLICY const& thePolicy, double theRiskFreeUpRelancedFraction,
        ReturnSpecifier const& theReturnSpecifier = ReturnSpecifier()): value(
            theValue), policy(thePolicy), returnSpecifier(theReturnSpecifier),
            riskFreeUpRelancedFraction(theRiskFreeUpRelancedFraction) {}
    void simulateStep(double netSavings, int step)
    {//create temporary static asset and use it to simulate
        double riskyValue = value * (1 - riskFreeUpRelancedFraction);
        double stockTarget = policy(riskyValue, step);
        MixedAsset<StockBondAsset> asset = makeMVOptimalRiskyAsset(riskyValue,
            stockTarget, tangency, returnSpecifier);
        asset.simulateStep(netSavings, step);
        value = asset.getValue() + value * riskFreeUpRelancedFraction *
            (1 + returnSpecifier.getRealRiskFreeRate());
    }
    double getValue()const{return value;}
    void setValue(double newValue){value = newValue;}
};

template<typename FINANCIAL_ASSET, typename SIM_RECORDER>
PortfolioSimulationResult performActuarialSimulation(FINANCIAL_ASSET const&
    initialFinancialAsset, SIM_RECORDER& simRecorder, Vector<double> const&
    netSavings, Vector<double> const& nextYearSurvivalProbabilities,
    int nSimulations = 1000000)
{
    assert(netSavings.getSize() > 0 && nSimulations > 0 &&
        netSavings.getSize() == nextYearSurvivalProbabilities.getSize());
    Vector<double> values, ruinValues, maxDrawdownValues, breachValues;
    double inflationRate = ReturnSpecifier().getInflationRate();
    for(int i = 0; i < nSimulations; ++i)
    {
        FINANCIAL_ASSET financialAsset = initialFinancialAsset;
        double initialValue = financialAsset.getValue(), maxDrawdown = 0,
            maxValue = initialValue;
        for(int step = 0; step < netSavings.getSize(); ++step)
        {//check if ruined
            if(financialAsset.getValue() < 0)
            {
                ruinValues.append(step);
                break;
            }
            if(nextYearSurvivalProbabilities[step] < GlobalRNG().uniform01())
            {//died, assume continuation till the end with no cash flows to
             //calculate normalized inheritance
                for(; step < netSavings.getSize(); ++step)
                    financialAsset.simulateStep(0, step);
                break;
            }
            financialAsset.simulateStep(netSavings[step], step);
            double realValue = financialAsset.getValue();
            simRecorder(i, step, realValue);
            if(realValue > maxValue) maxValue = realValue;
            else maxDrawdown = max(maxDrawdown, 1 - realValue/maxValue);
        }//0 if ruined
        simRecorder(i, netSavings.getSize(), financialAsset.getValue());
        values.append(max(0.0, financialAsset.getValue()));
        maxDrawdownValues.append(maxDrawdown);
    }
    return PortfolioSimulationResult(values, ruinValues, maxDrawdownValues);
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
    auto dummySimRecorder = [](int, int, double){};
    return performActuarialSimulation(initialFinancialAsset, dummySimRecorder,
        netSavings, nextYearSurvivalProbabilities, nSimulations);
}
template<typename FINANCIAL_ASSET>
PortfolioSimulationResult performSimulation(FINANCIAL_ASSET const&
    initialFinancialAsset, Vector<double> const& netSavings,
    int nSimulations = 1000000)
{
    auto dummySimRecorder = [](int, int, double){};
    return performActuarialSimulation(initialFinancialAsset, dummySimRecorder,
        netSavings, Vector<double>(netSavings.getSize(), 1), nSimulations);
}

Vector<double> generateRetirementExpenses(double initialExpenses,
    int term, int delay = 0)
{
    assert(delay >= 0 && delay < term);
    Vector<double> result(delay, 0);
    result.appendVector(Vector<double>(term - delay, initialExpenses) * -1);
    return result;
}

template<typename FINANCIAL_ASSET>
double findMaximumExpenseRatio(FINANCIAL_ASSET const& financialAsset,
    int term, double maximumRuinProbability = 0.05)
{
    auto f = [&](double expenseRatio)->double
    {
        Vector<double> retirementExpenses = generateRetirementExpenses(
            expenseRatio * 100, term);
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
            financialAssets[step].simulateStep(netSavings[step], step);
            if(step + 1 < financialAssets.getSize())
                financialAssets[step + 1].setValue(
                    financialAssets[step].getValue());
        }//0 if ruined
        values.append(max(0.0, financialAssets.lastItem().getValue()));
    }
    return PortfolioSimulationResult(values, ruinValues, Vector(1, 0.0));
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
        //add final home value
        result.percentiles.values +=
            Vector<double>(result.percentiles.values.getSize(), finalHomeValue);
        return expectedCRRACertaintyEquivalent(result.percentiles.values, 3);
    }
};
pair<Vector<double>, double> findOptimalBondFractions(int term,
    ReturnSpecifier const& returnSpecifier, double initialValue,
    double initialIncome, double homeValue = 0, double maxBondFraction = 0.4,
    double stepDelta = 0.05, int maxEvals = 300)
{
    Vector<double> dcaSavings(term, initialIncome);
    CRRA3Functor functor = {initialValue, returnSpecifier, dcaSavings,
        homeValue};
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
