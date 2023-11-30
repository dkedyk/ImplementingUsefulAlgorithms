#ifndef FINANCIAL_CALCULATIONS_TEST_AUTO_H
#define FINANCIAL_CALCULATIONS_TEST_AUTO_H
#include <string>
#include "CashFlows.h"
#include "Annuity.h"
#include "MeanVarianceOptimization.h"
#include "PortfolioSimulation.h"
#include "Misc.h"
#include "OptionPricing.h"
#include "../RandomNumberGeneration/testCommon.h"
#include "../RandomNumberGeneration/DistributionTests.h"
#include "../ExternalMemoryAlgorithms/CSV.h"

using namespace std;

namespace igmdk{

void testMSCIDitributionMatch(bool useLog = true)
{//both not rejected by DKW test
    Vector<double> samples = getMSCIWorldReturns();
    if(useLog)
    {//for lognormal
        for(int i = 0; i < samples.getSize(); ++i)
        {
            samples[i] = log(1 + samples[i]/100);
        }
    }
    IncrementalStatistics s;
    for(int i = 0; i < samples.getSize(); ++i)
    {
        s.addValue(samples[i]);
    }
    double mean = s.getMean(), stdev = s.stdev();
    auto normalCDF = [mean, stdev](double x){return approxNormalCDF((x - mean)/stdev);};
    DEBUG(DKWPValue(samples, normalCDF));
    //overwrite with random for compare
    for(int i = 0; i < samples.getSize(); ++i)
    {
        samples[i] = GlobalRNG().normal(mean, stdev);
    }
    DEBUG(DKWPValue(samples, normalCDF));
}

void testCashFlowAuto()
{
    DEBUG("testCashFlowAuto");
    Vector<double> cashFlowSimple;
    cashFlowSimple.append(1000);
    cashFlowSimple.append(1000);
    DEBUG(calculatePrice(cashFlowSimple, 0));
    DEBUG(calculateYield(cashFlowSimple, 2000));
    DEBUG("testCashFlowAuto passed");
}

void testBondAuto()
{
    DEBUG("testBondAuto");
    GeneralBond simple(50, 1000, 10, 0, 2);
    DEBUG(simple.calculateYield(1000));
    DEBUG("testBondAuto passed");
}

void testBondDurationConvexityAuto()
{
    DEBUG("testBondDurationConvexityAuto");
    GeneralBond simple(50, 1000, 10);//5% 10 year trading at par
    double currentYield = simple.calculateYield(1000);
    DEBUG(currentYield);
    DEBUG(simple.calculatePrice(currentYield));
    double increment = 0.03;
    DEBUG(simple.calculatePrice(currentYield + increment));
    double modifiedDuration = simple.calculateModifiedDuration(1000);
    DEBUG(modifiedDuration);
    double convexity = simple.calculateConvexity(1000);
    DEBUG(convexity);
    DEBUG(-modifiedDuration * increment + convexity/2 * increment * increment);
    DEBUG("testBondDurationConvexityAuto passed");
}

void testMortgageAuto()
{
    DEBUG("testMortgageAuto");
    Mortgage mortgage(0.03945, 100000, 23);
    DEBUG(mortgage.getMonthlyPayment());
    DEBUG("testMortgageAuto passed");
}

void testMortgageReverseAuto()
{//reverse Mortgage is annuity
    DEBUG("testMortgageReverseAuto");
    Mortgage mortgage(0.035, -100000, 23);
    DEBUG(mortgage.getMonthlyPayment());
    DEBUG("testMortgageReverseAuto passed");
}

void testAnnuityAuto()
{
    DEBUG("testAnnuityAuto");
    //annuity for 62 year old with 100K giving 550 per month:
    GeneralBond annuity(550, 0, 12 * 22);//84 life expectancy
    DEBUG("expected life model");
    DEBUG(annuity.calculateYield(100000) * 12);
    DEBUG("Female 62 based on SSA tables");
    Annuity<> annuity2(62, convertToSurvivalProbabilities(getSSAFemaleDeathProbabilities()));
    DEBUG(annuity2.calculateYield(100000, 550));

    DEBUG("testAnnuityAuto passed");
}

void testAnnuityAuto2()
{
    DEBUG("testAnnuityAuto2");
    for(int age = 35; age < 65; age += 5)
    {
        DEBUG(age);
        //annuity for x year old male giving 10K per month at 1.5 real yield:
        DEBUG("Male x based on SSA tables");
        Annuity<> annuity2(age, convertToSurvivalProbabilities(getSSAMaleDeathProbabilities()));
        DEBUG(annuity2.calculatePrice(7000, 0.015));
    }
    DEBUG("testAnnuityAuto passed");
}

void testAnnuityAutoJoint()
{
    DEBUG("testAnnuityAutoJoint");
    //annuity for 62 year old with 100K giving 550 per month:

    DEBUG("Male 65 Female 60 based on SSA tables");
    JointSurvivalEstimator e;
    e.addPerson(convertToSurvivalProbabilities(getSSAFemaleDeathProbabilities()), 60);
    e.addPerson(convertToSurvivalProbabilities(getSSAMaleDeathProbabilities()), 65);
    Annuity<JointSurvivalEstimator> annuity2(0, e);
    DEBUG(annuity2.calculateYield(100000, 500));

    DEBUG("testAnnuityAutoJoint passed");
}

void testAnnuityAutoJoint2()
{
    DEBUG("testAnnuityAutoJoint");
    //annuity for 62 year old with 100K giving 550 per month:
    for(int age = 35; age < 65; age += 5)
    {
        DEBUG(age);
        DEBUG("Male 35 Female 34 based on SSA tables");
        JointSurvivalEstimator e;
        e.addPerson(convertToSurvivalProbabilities(getSSAMaleDeathProbabilities()), age - 1);
        e.addPerson(convertToSurvivalProbabilities(getSSAFemaleDeathProbabilities()), age);
        Annuity<JointSurvivalEstimator> annuity2(0, e);
        DEBUG(annuity2.calculatePrice(7000, 0.015));
    }
    DEBUG("testAnnuityAutoJoint2 passed");
}

void testAnnuityAutoJointSpending()
{
    DEBUG("testAnnuityAutoJointSpending");
    Vector<Vector<string> > matrix;
    for(int age = 35; age <= 95; age += 1)
    {
        DEBUG(age);
        JointSurvivalEstimator e;
        e.addPerson(convertToSurvivalProbabilities(getSSAMaleDeathProbabilities()), age);
        e.addPerson(convertToSurvivalProbabilities(getSSAFemaleDeathProbabilities()), age);
        Annuity<JointSurvivalEstimator> annuity2(0, e);
        DEBUG("actuarial present value of spending 1K a month at 2.5% real rate");
        DEBUG(annuity2.calculatePrice(1000, 0.025));
        Vector<string> row(1, to_string(annuity2.calculatePrice(1000, 0.025)));
        matrix.append(row);
    }
    //createCSV(matrix, "spending_.csv");
}

void testAnnuityAutoSingleFemaleSpending()
{
    DEBUG("testAnnuityAutoSingleFemaleSpending");
    for(int age = 60; age <= 100; age += 5)
    {
        DEBUG(age);
        JointSurvivalEstimator e;
        e.addPerson(convertToSurvivalProbabilities(getSSAFemaleDeathProbabilities()), age);
        Annuity<JointSurvivalEstimator> annuity2(0, e);
        //DEBUG("actuarial present value of spending 1K a month at 2.5% real rate");
        //DEBUG(annuity2.calculatePrice(1000, 0.025));
        DEBUG("actuarial present value of spending 10K a month at 6% nominal rate");
        DEBUG(annuity2.calculatePrice(10000, 0.06));
    }
}

void testOptionPricer()
{
    DEBUG("testOptionPricer");
    double price = 62, strikePrice = 60, r = 0.1, q = 0.2, strikeTime = 5.0/12;
    assert(isEEqual(priceEuropeanCall(price, r, q, strikePrice, strikeTime), 5.7977812415148975));
    double riskFreeRate = 0.03;
    assert(isEEqual(priceEuropeanPut(price, riskFreeRate, r, q, strikePrice, strikeTime), 3.0524492711477933));
    DEBUG("testOptionPricer passed");
}

void testEstimateLogNormalParameters()
{
    double mean = 1.08, stdev = 0.17;
    DEBUG(mean);
    DEBUG(stdev);
    pair<double, double> result = estimateLogNormalParameters(mean, stdev);
    double mu = result.first;
    double q = result.second, q2 = q * q;
    DEBUG(mu);
    DEBUG(q);
    DEBUG(exp(mu));
    DEBUG(mean - stdev * stdev/2);
    DEBUG(exp(mu + q2/2));
    DEBUG((exp(q2) - 1) * exp(2 * mu + q2));

    DEBUG(mean - exp(mu + q2/2));
    DEBUG(stdev * stdev - (exp(q2) - 1) * exp(2 * mu + q2));
}

void testReturnEstimator()
{
    DEBUG("Return Estimator");
    double pe = 20, bondYield = 0.055, inflation = 0.023, stdevStock = 0.17, stdevBond = 0.09, riskFreeRate = 0.048;
    double stockReturn = estimateArithmeticNominalStockReturn(pe, inflation, stdevStock);
    DEBUG(stockReturn - 1);//accurate
    DEBUG(1/pe + inflation + stdevStock*stdevStock/2);//simple
    LognormalDistribution ld(stockReturn, stdevStock);
    double bondReturn = estimateArithmeticNominalBondReturn(bondYield, riskFreeRate, stdevBond);
    DEBUG(bondReturn - 1);//accurate
    DEBUG(riskFreeRate + (bondYield - riskFreeRate)/2 + stdevBond*stdevBond/2);//simple
}

void testMultigoal()
{
    {//bond ladder
        double spendingFraction = 0.5;
        int nYears = 1;
        double riskFreeRate = 0.03;
        DEBUG("return rate on remainder");
        ReturnSpecifier returnSpecifier;


        Vector<pair<double, int>> goals;
        goals.append({spendingFraction, nYears});
        LognormalDistribution muq = multigoalAdjustedLognormal(returnSpecifier, goals, nYears);
        DEBUG("return rate on remainder with leverage");
        DEBUG(muq.getMedian());
    }
}

void testCRRA()
{
    Vector<double> values;
    values.append(100);
    values.append(200);
    DEBUG(expectedCRRACertaintyEquivalent(values, 0.001));
    DEBUG(expectedCRRACertaintyEquivalent(values, 1));
    DEBUG(expectedCRRACertaintyEquivalent(values, 1.5));
    DEBUG(expectedCRRACertaintyEquivalent(values, 2));
    DEBUG(expectedCRRACertaintyEquivalent(values, 2.5));
    DEBUG(expectedCRRACertaintyEquivalent(values, 3));
    DEBUG(expectedCRRACertaintyEquivalent(values, 10));
}

void testGoalCalculator()
{
    double amount = 100000, rate = 0.03;
    int term = 10;
    DEBUG(calculateGoalMonthlySavings(amount, term, rate));
    DEBUG(calculateGoalMonthlySavings(amount, term, rate) * 12 * term/amount);
}

void testPortfolioSimulationLifetime()
{
    //ReturnSpecifier returnSpecifier(0.2, 0.3);//tax
    ReturnSpecifier returnSpecifier;//no tax
    for(int retirementStep = 30; retirementStep <= 30; retirementStep += 5)
    {
        DEBUG(retirementStep);
        for(int term = 60; term <= 60; term += 10)
        {
            DEBUG(term);
            double inflationFactor = pow(1 + returnSpecifier.getInflationRate(), term);
            double riskFreeReturn = 0;
            Vector<double> netSavings = lifetimeNetSavings(100, 100, returnSpecifier.getInflationRate(), term, retirementStep);
            {
                PortfolioSimulationResult safe = performSimulation(RiskFreeAsset(returnSpecifier), netSavings);
                riskFreeReturn = safe.getMedian();
                DEBUG("risk-free");
                safe.debug(inflationFactor);
            }




            for(int percentStock = 50; percentStock <= 100; percentStock += 10)
            {
                DEBUG(percentStock);
                PortfolioSimulationResult result = performSimulation(StockBondAsset(returnSpecifier, 0, 1 - percentStock/100.0), netSavings);
                result.debug(inflationFactor);
                //DEBUG(result.sharpeRatio(riskFreeReturn));
                double riskFreeRank = result.riskFreeRank(riskFreeReturn);
                DEBUG(riskFreeRank);
            }
        }
    }


}

void testPortfolioSimulationRetirement()
{
    //ReturnSpecifier returnSpecifier(0.2, 0.3);//tax
    ReturnSpecifier returnSpecifier;//no tax
    for(int expensesNowPercent = 2; expensesNowPercent <= 4; ++expensesNowPercent)
    {
        DEBUG(expensesNowPercent);
        int term = 50;
        double inflationFactor = pow(1 + returnSpecifier.getInflationRate(), term);
        double riskFreeReturn = 0;
        Vector<double> retirementExpenses = generateRetirementExpenses(expensesNowPercent, returnSpecifier.getInflationRate(), term);
        {
            PortfolioSimulationResult safe = performSimulation(RiskFreeAsset(returnSpecifier, 100), retirementExpenses);
            riskFreeReturn = safe.getMedian();
            DEBUG("risk-free");
            safe.debug(inflationFactor);
        }
        for(int percentStock = 50; percentStock <= 100; percentStock += 10)
        {
            DEBUG(percentStock);
            {
                PortfolioSimulationResult result = performSimulation(StockBondAsset(returnSpecifier, 100, 1 - percentStock/100.0), retirementExpenses);
                result.debug(inflationFactor);
                //DEBUG(result.sharpeRatio(riskFreeReturn));
                double riskFreeRank = result.riskFreeRank(riskFreeReturn);
                DEBUG(riskFreeRank);
            }
        }
    }

}

void testPortfolioSimulationStatic()
{
    //ReturnSpecifier returnSpecifier(0.2, 0.3);//tax
    ReturnSpecifier returnSpecifier;//no tax
        int term = 30;
        double inflationFactor = pow(1 + returnSpecifier.getInflationRate(), term);
        double riskFreeReturn = 0;
        Vector<double> noSavings(term, 0);
        {
            PortfolioSimulationResult safe = performSimulation(RiskFreeAsset(returnSpecifier, 100), noSavings);
            riskFreeReturn = safe.getMedian();
            DEBUG("risk-free");
            safe.debug(inflationFactor);
        }
        for(int percentStock = 50; percentStock <= 100; percentStock += 10)
        {
            DEBUG(percentStock);
            {
                PortfolioSimulationResult result = performSimulation(StockBondAsset(returnSpecifier, 100, 1 - percentStock/100.0), noSavings);
                result.debug(inflationFactor);
                //DEBUG(result.sharpeRatio(riskFreeReturn));
                double riskFreeRank = result.riskFreeRank(riskFreeReturn);
                DEBUG(riskFreeRank);
            }
        }

}

void testPortfolioSimulationDCA()
{
    //ReturnSpecifier returnSpecifier(0.2, 0.3);//tax
    ReturnSpecifier returnSpecifier;//no tax
        int term = 30;
        double inflationFactor = pow(1 + returnSpecifier.getInflationRate(), term);
        double riskFreeReturn = 0;
        Vector<double> dcaSavings = generateDCASavings(100, returnSpecifier.getInflationRate(), term);
        {
            PortfolioSimulationResult safe = performSimulation(RiskFreeAsset(returnSpecifier), dcaSavings);
            riskFreeReturn = safe.getMedian();
            DEBUG("risk-free");
            safe.debug(inflationFactor);
        }

        for(int percentStock = 50; percentStock <= 100; percentStock += 10)
        {
            DEBUG(percentStock);
            {//note that first year contribution doesn't grow but amount saved in first year is reduced by inflation
                PortfolioSimulationResult result = performSimulation(StockBondAsset(returnSpecifier, 0, 1 - percentStock/100.0), dcaSavings);
                result.debug(inflationFactor);
                //DEBUG(result.sharpeRatio(riskFreeReturn));
                double riskFreeRank = result.riskFreeRank(riskFreeReturn);
                DEBUG(riskFreeRank);
            }
        }
}

void testMinVarCalculationCurrent()
{
    DEBUG("testMinVarCalculationCurrent");
    ReturnSpecifier returnSpecifier;
    MeanVariancePortfolio mvp(makeStockBondMVP(returnSpecifier));
    Vector<pair<double, Vector<double> > > frontier =
        mvp.findOptimalPortfolioWeightsRange(21);
    DEBUG(returnSpecifier.getRiskFreeRate());
    for(int i = 0; i < frontier.getSize(); ++i)
    {
        DEBUG(i);
        DEBUG(frontier[i].first);
        frontier[i].second.debug();
        pair<double, double> ms = mvp.evaluate(frontier[i].second);
        DEBUG(ms.first);
        DEBUG(ms.second);
        DEBUG(ms.first - 2 * ms.second);
        DEBUG(ms.first + 2 * ms.second);

        LognormalDistribution ld(1 + ms.first, ms.second);
        DEBUG(ld.getMedian() - 1);

        DEBUG("Sharpe");
        DEBUG((ms.first - returnSpecifier.getRiskFreeRate())/ms.second);
        DEBUG("CRRA 2");
        DEBUG(ms.first - ms.second * ms.second);
    }

    DEBUG("Optimal Sharpe");
    pair<double, Vector<double> > optimal = mvp.findOptimalSharpeWeights(returnSpecifier.getRiskFreeRate());
    DEBUG(optimal.second[0]);
    DEBUG(-optimal.first);
    pair<double, double> ms = mvp.evaluate(optimal.second);
    DEBUG(ms.first);
    DEBUG(ms.second);
    DEBUG((ms.first - returnSpecifier.getRiskFreeRate())/ms.second);
}

void testMinVarCalculationCurrentTaxes()
{
    DEBUG("testMinVarCalculationCurrentTaxes");
    ReturnSpecifier returnSpecifier(0.2, 0.32, 0.06);
    MeanVariancePortfolio mvp(makeStockBondMVP(returnSpecifier));
    Vector<pair<double, Vector<double> > > frontier =
        mvp.findOptimalPortfolioWeightsRange(11);
    DEBUG(returnSpecifier.getRiskFreeRate());
    for(int i = frontier.getSize()/2; i < frontier.getSize(); ++i)
    {
        DEBUG(i);
        DEBUG(frontier[i].first);
        frontier[i].second.debug();
        pair<double, double> ms = mvp.evaluate(frontier[i].second);
        DEBUG(ms.first);
        DEBUG(ms.second);
        DEBUG(ms.first - 2 * ms.second);
        DEBUG(ms.first + 2 * ms.second);

        LognormalDistribution ld(1 + ms.first, ms.second);
        DEBUG(ld.getMedian() - 1);

        DEBUG("Sharpe");
        DEBUG((ms.first - returnSpecifier.getRiskFreeRate())/ms.second);
        DEBUG("CRRA 2");
        DEBUG(ms.first - ms.second * ms.second);
        if(i == frontier.getSize() - 1)
        {
            DEBUG("with 10% risk-free");
            DEBUG(ms.first * 0.9 + returnSpecifier.getRiskFreeRate() * 0.1);
            DEBUG(ms.second * 0.9);
            DEBUG(0.9 * (ms.first - 2 * ms.second));
            DEBUG(0.9 * (ms.first + 2 * ms.second));
            LognormalDistribution ld2(1 + ms.first * 0.9 + returnSpecifier.getRiskFreeRate() * 0.1, ms.second * 0.9);
            DEBUG(ld2.getMedian() - 1);
        }
    }

    DEBUG("Optimal Sharpe");
    pair<double, Vector<double> > optimal = mvp.findOptimalSharpeWeights(returnSpecifier.getRiskFreeRate());
    DEBUG(optimal.second[0]);
    DEBUG(-optimal.first);
    pair<double, double> ms = mvp.evaluate(optimal.second);
    DEBUG(ms.first);
    DEBUG(ms.second);
    DEBUG((ms.first - returnSpecifier.getRiskFreeRate())/ms.second);

    LognormalDistribution ld(1 + ms.first, ms.second);
    DEBUG(ld.getMedian() - 1);

    DEBUG("with 10% risk-free");
    DEBUG(ms.first * 0.9 + returnSpecifier.getRiskFreeRate() * 0.1);
    DEBUG(ms.second * 0.9);
    LognormalDistribution ld2(1 + ms.first * 0.9 + returnSpecifier.getRiskFreeRate() * 0.1, ms.second * 0.9);
    DEBUG(ld2.getMedian() - 1);
}

void testMultiAccountReturns()
{
    DEBUG("total with all accounts");
    ReturnSpecifier taxable(0.2, 0.32, 0.06), taxEfficient;
    MultiAccountReturns mar;
    mar.addAccount(300, 0, 1.0/3, taxable);
    mar.addAccount(30, 0.75, 0, taxEfficient);
    mar.addAccount(500, 0.25, 0, taxEfficient);
    DEBUG(mar.getTotalGeometricReturnRate());
    DEBUG(mar.getTotalStockProportion());
    MultiAccountReturns allStock;
    allStock.addAccount(300, 0, 0, taxable);
    allStock.addAccount(30, 0, 0, taxEfficient);
    allStock.addAccount(500, 0, 0, taxEfficient);
    DEBUG(allStock.getTotalGeometricReturnRate());
    MultiAccountReturns lessBonds;
    lessBonds.addAccount(300, 0, 1.0/3, taxable);
    lessBonds.addAccount(30, 0.25, 0, taxEfficient);
    lessBonds.addAccount(500, 0.10, 0, taxEfficient);
    DEBUG(lessBonds.getTotalGeometricReturnRate());
    DEBUG(lessBonds.getTotalStockProportion());
    MultiAccountReturns lessBondsNoEM;
    lessBondsNoEM.addAccount(300, 0, 0, taxable);
    lessBondsNoEM.addAccount(30, 0.25, 0, taxEfficient);
    lessBondsNoEM.addAccount(500, 0.10, 0, taxEfficient);
    DEBUG(lessBondsNoEM.getTotalGeometricReturnRate());
    DEBUG(lessBondsNoEM.getTotalStockProportion());
    MultiAccountReturns targetBondsNoEM;
    targetBondsNoEM.addAccount(300, 0, 0, taxable);
    targetBondsNoEM.addAccount(30, 0.75, 0, taxEfficient);
    targetBondsNoEM.addAccount(500, 0.25, 0, taxEfficient);
    DEBUG(targetBondsNoEM.getTotalGeometricReturnRate());
    DEBUG(targetBondsNoEM.getTotalStockProportion());
}

void testMultiAccountReturnsImbalance()
{
    DEBUG("25-75 vs 50-50");
    ReturnSpecifier taxEfficient;
    MultiAccountReturns mar;
    mar.addAccount(1, 0.75, 0, taxEfficient);
    mar.addAccount(1, 0.25, 0, taxEfficient);
    DEBUG(mar.getTotalGeometricReturnRate());
    MultiAccountReturns optimal;
    optimal.addAccount(1, 0.5, 0, taxEfficient);
    DEBUG(optimal.getTotalGeometricReturnRate());
}

void testMultiAccountReturnsEmergencyFund()
{
    DEBUG("with em vs without");
    ReturnSpecifier taxable(0.2, 0.32, 0.06);
    MultiAccountReturns mar;
    mar.addAccount(1, 0, 2.0/3, taxable);
    DEBUG(mar.getTotalGeometricReturnRate());
    MultiAccountReturns optimal;
    optimal.addAccount(1, 0, 0, taxable);
    DEBUG(optimal.getTotalGeometricReturnRate());
}

void testEstimateStockBondRebalance()
{
    DEBUG("no taxes case");
    double stockBondPercents = estimateStockBondSalesForRebalance();
    DEBUG(stockBondPercents);
    DEBUG("taxes case");
    stockBondPercents = estimateStockBondSalesForRebalance(ReturnSpecifier(0.2, 0.32, 0.06));
    DEBUG(stockBondPercents);
}

void testLiabilityPVDurationCalculation()
{
    DEBUG("emergency fund 100 and 10 years of spending 100 at 5% and 2.5% inflation")
    Vector<pair<double, double>> liabilities;
    liabilities.append({100, 0});//emf
    for(int i = 0; i < 10; ++i)
    {
        liabilities.append({100 * pow(1.025, i), i});
    }
    pair<double, double> result = calculatePVAndDuration(liabilities, 0.05);
    DEBUG(result.first);
    DEBUG(result.second);
    DEBUG("same but for next year");
    Vector<pair<double, double>> liabilities2;
    liabilities2.append({100 * 1.025, 0});//emf
    for(int i = 0; i < 9; ++i)
    {
        liabilities2.append({100 * 1.025 * pow(1.025, i), i});
    }
    result = calculatePVAndDuration(liabilities2, 0.05);
    DEBUG(result.first);
    DEBUG(result.second);
}

void testAllAutoFinancialCalculations()
{
    DEBUG("testAllAutoFinancialCalculations");
    testMSCIDitributionMatch();
    //return;
    testOptionPricer();
    testCRRA();
    testMultigoal();
    //return;
    testGoalCalculator();
    //return;
    testCashFlowAuto();
    testBondAuto();
    testBondDurationConvexityAuto();
    //return;
    testMortgageAuto();
    testMortgageReverseAuto();
    testAnnuityAuto();
    testAnnuityAuto2();
    testAnnuityAutoJoint();
    testAnnuityAutoJoint2();
    testAnnuityAutoSingleFemaleSpending();
    testAnnuityAutoJointSpending();
    //return;
    testEstimateLogNormalParameters();
    testReturnEstimator();
    //return;
    testMinVarCalculationCurrent();
    testMinVarCalculationCurrentTaxes();
    testMultiAccountReturns();
    testMultiAccountReturnsImbalance();
    testMultiAccountReturnsEmergencyFund();
    testEstimateStockBondRebalance();
    testLiabilityPVDurationCalculation();
    //testPortfolioSimulationLifetime();
    //testPortfolioSimulationRetirement();
    //testPortfolioSimulationStatic();
    testPortfolioSimulationDCA();
}

}//end namespace
#endif
