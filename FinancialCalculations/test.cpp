#include "FinancialCalculationsTestAuto.h"
#include <cassert>
#include <iostream>
#include <cstdlib>
#include <ctime>
using namespace igmdk;

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

void testMakeReturnFile()
{
    DEBUG("testMakeReturnFile");
    Vector<Vector<string> > matrix;
    ReturnSpecifier returnSpecifier;

    LognormalDistribution stockReturns = StockBondAsset::makeReturnDistribution(
        returnSpecifier, 0), bondReturns =
        StockBondAsset::makeReturnDistribution(returnSpecifier, 1);

    int n = 20;
    Vector<string> labels = Vector<string>(n, "StockReturn");
    labels.appendVector(Vector<string>(n, "BondReturn"));
    matrix.append(labels);
    Vector<string> selections = Vector<string>(n, "0");
    selections.appendVector(Vector<string>(n, ""));
    matrix.append(selections);
    matrix.append(Vector<string>(2*n, ""));
    matrix.append(Vector<string>(2*n, ""));
    Vector<string> constants(2 * n, "");
    constants[0] = "Geometric:";
    constants[1] = to_string(stockReturns.getMedian()-1);
    constants[2] = "Treasury:";
    constants[3] = to_string(returnSpecifier.getRiskFreeRate());
    constants[4] = "Inflation:";
    constants[5] = to_string(returnSpecifier.getInflationRate());
    constants[6] = "TangencyStocks:";
    MeanVariancePortfolio mvp(makeStockBondMVP(returnSpecifier));
    constants[7] = to_string(mvp.findOptimalSharpeWeights(
        returnSpecifier.getRiskFreeRate()).second[0]);
    constants[8] = "Bonds:";
    constants[9] = to_string(bondReturns.getMedian()-1);
    constants[10] = "GeometricReturns100By5:";


    pair<double, Vector<double> > tangency = mvp.findOptimalSharpeWeights(returnSpecifier.getRiskFreeRate());
    pair<double, double> mqTangency = mvp.evaluate(tangency.second);
    Vector<pair<double, Vector<double> > > frontier =
        mvp.findOptimalPortfolioWeightsRange(21);
    for(int i = 0; i < frontier.getSize(); ++i)
    {
        int frontierIndex = frontier.getSize() - 1 - i;
        pair<double, double> mq;
        if(tangency.second[0] > frontier[frontierIndex].second[0])
        {//stock below tangency, will mix with risk-free to get target stock
            double tangencyFraction = max(frontier[frontierIndex].second[0]/tangency.second[0],
                numeric_limits<double>::epsilon());//to avoid numerical issues
            mq.first = mqTangency.first * tangencyFraction +
                returnSpecifier.getRiskFreeRate() * (1 - tangencyFraction);
            mq.second = mqTangency.second * tangencyFraction;
        }
        else
        {//stocks at or above tangency
            mq = mvp.evaluate(frontier[frontierIndex].second);
        }
        LognormalDistribution lognormal(1 + mq.first, mq.second);
        constants[11 + i] = to_string(lognormal.getMedian() - 1);
    }
    matrix.append(constants);
    for(int i = 0; i <= 100; ++i)
    {
        Vector<string> row;
        for(int i = 0; i < n; ++i)
            row.append(to_string(stockReturns.sample()-1));
        for(int i = 0; i < n; ++i)
            row.append(to_string(bondReturns.sample()-1));
        matrix.append(row);
    }
    //createCSV(matrix, "return_simulations.csv");
}

void testPrintJointSurvivalProbabilities()
{
    DEBUG("testPrintJointSurvivalProbabilities");
    Vector<Vector<string> > matrix;
    int age = 40;
    JointSurvivalEstimator e;
    e.addPerson(convertToSurvivalProbabilities(getSSAMaleDeathProbabilities()), age);
    e.addPerson(convertToSurvivalProbabilities(getSSAFemaleDeathProbabilities()), age);
    for(int i = 0; i <= 102; ++i)//spreadsheet needs + 2
    {
        DEBUG(i);
        DEBUG(e.getSurvivalProbability(0, i));
        Vector<string> row(1, to_string(e.getSurvivalProbability(0, i)));
        matrix.append(row);
    }
    //createCSV(matrix, "survival_probs.csv");
}

void testAnnuityAutoSingleFemaleSpending()
{
    DEBUG("testAnnuityAutoSingleFemaleSpending");
    for(int age = 70; age <= 70; age += 5)
    {
        DEBUG(age);
        JointSurvivalEstimator e;
        e.addPerson(convertToSurvivalProbabilities(getSSAFemaleDeathProbabilities()), age);
        Annuity<JointSurvivalEstimator> annuity2(0, e);
        DEBUG("actuarial present value of spending 5K a month at 0.5% real rate");
        DEBUG(annuity2.calculatePrice(5000, 0.005));
        DEBUG("actuarial present value of spending 5K a month at 2% real rate");
        DEBUG(annuity2.calculatePrice(5000, 0.02));
        Annuity<JointSurvivalEstimator> annuity3(0, e, 12, 0, 0, 0.1, 0.01, 2);
        DEBUG("commercial present value of spending 5K a month at 0.5% real rate");
        DEBUG(annuity3.calculatePrice(5000, 0.005));
        DEBUG("commercial present value of spending 5K a month at 2% real rate");
        DEBUG(annuity3.calculatePrice(5000, 0.02));
        Annuity<JointSurvivalEstimator> annuity5(0, e, 12, 0.03, 0, 0.1, 0.01, 2);
        DEBUG("commercial present value of spending 5K a month at 3.7% rate 3% COLA");
        DEBUG(annuity5.calculatePrice(5000, 0.037));
        DEBUG("commercial present value of spending 5K a month at 4.7% rate 3% COLA");
        DEBUG(annuity5.calculatePrice(5000, 0.055));
        e.addPerson(convertToSurvivalProbabilities(getSSAMaleDeathProbabilities()), age);
        Annuity<JointSurvivalEstimator> annuity4(0, e, 12, 0, 0, 0.1, 0.01, 2);
        DEBUG("couple commercial present value of spending 5K a month at 0.5% real rate");
        DEBUG(annuity4.calculatePrice(5000, 0.005));
        DEBUG("couple commercial present value of spending 5K a month at 2% real rate");
        DEBUG(annuity4.calculatePrice(5000, 0.02));
        Annuity<JointSurvivalEstimator> annuity6(0, e, 12, 0.03, 0, 0.1, 0.01, 2);
        DEBUG("couple commercial present value of spending 5K a month at 3.7% rate 3% COLA");
        DEBUG(annuity6.calculatePrice(5000, 0.037));
        DEBUG("couple commercial present value of spending 5K a month at 4.7% rate 3% COLA");
        DEBUG(annuity6.calculatePrice(5000, 0.047));
        Annuity<JointSurvivalEstimator> annuity7(0, e, 12, 0.03);
        DEBUG("couple fair present value of spending 5K a month at 3.7% rate 3% COLA");
        DEBUG(annuity7.calculatePrice(5000, 0.037));
        DEBUG("couple fair present value of spending 5K a month at 4.7% rate 3% COLA");
        DEBUG(annuity7.calculatePrice(5000, 0.047));
    }
}


/* From early 2024:
    double pe = 25.07, bondYield = 0.0485, inflation = 0.0235, stdevStock = 0.17, stdevBond = 0.09, riskFreeRate = 0.042, dp = 0.014, earningsGrowth = 0.04;
    double internationalPE = 15.38, usFraction = 0.6;
*/

void testReturnEstimator()
{
    DEBUG("Return Estimator");
    //From early 2025:
    double pe = 28.2, bondYield = 0.05, inflation = 0.024, stdevStock = 0.17, stdevBond = stdevStock/2, riskFreeRate = 0.046, dp = 0.012, earningsGrowth = 0.04;
    double internationalPE = 15.78, usFraction = 0.66;
    double stockReturnUS = estimateArithmeticNominalStockReturn(pe, inflation, stdevStock, dp, earningsGrowth);
    double stockReturnEXUS = estimateArithmeticNominalStockReturn(internationalPE, inflation, stdevStock);
    double totalStockReturn = stockReturnUS * usFraction + stockReturnEXUS * (1 - usFraction);
    DEBUG(totalStockReturn - 1);

    double bondReturn = estimateArithmeticNominalBondReturn(bondYield, riskFreeRate, stdevBond);
    DEBUG(bondReturn - 1);//accurate
}

void testPortfolioSimulationNoCashFlow()
{
    ReturnSpecifier returnSpecifier;//no tax
    double tangencyStockFraction = getTangencyStockFraction(returnSpecifier);
        int term = 30;
        DEBUG(term);
        double realRiskFreeReturn = 0;
        Vector<double> noSavings(term, 0);
        for(int percentStock = 0; percentStock <= 100; percentStock += 20)
        {
            DEBUG(percentStock);
            PortfolioSimulationResult result = performSimulation(makeMVOptimalRiskyAsset(100, percentStock/100.0, tangencyStockFraction, returnSpecifier), noSavings, 1000000);
            result.debug();
            if(percentStock == 0)
            {
                realRiskFreeReturn = result.getMedian();
            }
            else
            {
                double riskFreeRank = result.riskFreeRank(realRiskFreeReturn);
                DEBUG(riskFreeRank);
            }
        }

}

void testPortfolioSimulationRetirement()
{
    ReturnSpecifier returnSpecifier;//no tax
    double tangencyStockFraction = getTangencyStockFraction(returnSpecifier);

    Vector<Vector<string> > matrix;
    Vector<string> constants;
    constants.append("Initial Expense Ratio %");
    constants.append("Stock %");
    constants.append("Median Units");
    constants.append("Ruin Probability");
    constants.append("5th % Time to Ruin");
    matrix.append(constants);

    for(int expensesNowPercent = 2; expensesNowPercent <= 6; ++expensesNowPercent)
    {
        DEBUG(expensesNowPercent);
        int term = 30;
        //int term = 50;
        //int term = 1000;
        //DEBUG(term);
        Vector<double> retirementExpenses = generateRetirementExpenses(expensesNowPercent, term);
        double realRiskFreeReturn = 0;
        for(int percentStock = 0; percentStock <= 100; percentStock += 20)
        {
            DEBUG(percentStock);
            PortfolioSimulationResult result = performSimulation(makeMVOptimalRiskyAsset(100, percentStock/100.0, tangencyStockFraction, returnSpecifier), retirementExpenses, 1000000);
            //result.debug();
            if(percentStock == 0)
            {
                realRiskFreeReturn = result.getMedian();
            }
            else
            {
                double riskFreeRank = result.riskFreeRank(realRiskFreeReturn);
                //DEBUG(riskFreeRank);
            }

            Vector<string> row;
            row.append(to_string(expensesNowPercent));
            row.append(to_string(percentStock));
            row.append(to_string(result.getMedian()));
            row.append(to_string(result.ruinChance));
            if(result.ruinChance > 0)
            {
                double targetPercentile = 0.05/result.ruinChance;
                if(targetPercentile <= 1)
                {
                    double percentile5TimeToRuin =
                        result.ruinPercentiles.getPercentile(targetPercentile);
                    row.append(to_string(percentile5TimeToRuin));
                }
                else row.append("N/A");
            }
            else row.append("N/A");
            matrix.append(row);
        }
    }

    //createCSV(matrix, "expense_ratio_simulation.csv");
    //createCSV(matrix, "expense_ratio_simulation_long.csv");
}

void testPortfolioSimulationRetirementTaxable()
{
    ReturnSpecifier returnSpecifier(0.15, 0.15, 0.00);//tax
    int term = 100;
    double taxableFraction = 0.5;

    Vector<Vector<string> > matrix;
    int n = 2;
    Vector<string> constants(n, "");
    constants[0] = "Years";
    constants[1] = "Expense Ratio";
    matrix.append(constants);
    for(; term >= 0; term -= 5)
    {
        if(term == 0) term = 1;
        DEBUG(term);
        double fraction = findMaximumExpenseRatio(StockBondAsset(returnSpecifier, 100, 0), term, 0.05);
        double taxAdjustedFraction = fraction/(1 + taxableFraction * (returnSpecifier.taxRateCapitalGains + returnSpecifier.taxRateLocal));
        DEBUG(taxAdjustedFraction);

        Vector<string> row;
        row.append(to_string(term));
        row.append(to_string(taxAdjustedFraction));
        matrix.append(row);
    }
    createCSV(matrix, "early_retirement_amounts.csv");

}

void testPortfolioSimulationRetirementDelay()
{//for partial ladder
    //ReturnSpecifier returnSpecifier(0.2, 0.3);//tax
    ReturnSpecifier returnSpecifier;//no tax
    double tangencyStockFraction = getTangencyStockFraction(returnSpecifier);
    for(int expensesNowPercent = 2; expensesNowPercent <= 6; ++expensesNowPercent)
    {
        DEBUG(expensesNowPercent);
        int term = 50;
        DEBUG(term);
        int delay = 20;
        int initialBudget = 100;
        double TIPSRate = returnSpecifier.getRealRiskFreeRate();
        initialBudget -= ordinaryAnnuityPV(expensesNowPercent, delay, TIPSRate);
        DEBUG(ordinaryAnnuityPV(expensesNowPercent, delay, TIPSRate));
        double yearsToDepleteWithTIPS = log(expensesNowPercent/100.0/(expensesNowPercent/100.0 - TIPSRate))/TIPSRate;
        DEBUG(yearsToDepleteWithTIPS);

        JointSurvivalEstimator e;
        e.addPerson(convertToSurvivalProbabilities(getSSAMaleDeathProbabilities()), 70);
        e.addPerson(convertToSurvivalProbabilities(getSSAFemaleDeathProbabilities()), 70);
        Annuity<JointSurvivalEstimator> annuity5(0, e, 1, 0, 0, 0.1, 0.01, 2);
        DEBUG("commercial present value TIPS annuity");
        DEBUG(annuity5.calculatePrice(expensesNowPercent, returnSpecifier.getRealRiskFreeRate()));

        assert(initialBudget > 0);
        double realRiskFreeReturn = 0;
        Vector<double> retirementExpenses = generateRetirementExpenses(
            expensesNowPercent, term, delay);
        for(int percentStock = 0; percentStock <= 100; percentStock += 20)
        {
            DEBUG(percentStock);
            PortfolioSimulationResult result = performActuarialSimulation(makeMVOptimalRiskyAsset(initialBudget, percentStock/100.0, tangencyStockFraction, returnSpecifier), retirementExpenses, e);
            result.debug();
            if(percentStock == 0)
            {
                realRiskFreeReturn = result.getMedian();
            }
            else
            {
                double riskFreeRank = result.riskFreeRank(realRiskFreeReturn);
                DEBUG(riskFreeRank);
            }
        }
    }

}

void testPortfolioSimulationRetirementActurial()
{
    ReturnSpecifier returnSpecifier;//no tax
    double tangencyStockFraction = getTangencyStockFraction(returnSpecifier);
    JointSurvivalEstimator e;
    e.addPerson(convertToSurvivalProbabilities(getSSAMaleDeathProbabilities()), 70);
    e.addPerson(convertToSurvivalProbabilities(getSSAFemaleDeathProbabilities()), 70);

    Vector<Vector<string> > matrix;
    Vector<string> constants;
    constants.append("Initial Expense Ratio %");
    constants.append("Stock %");
    constants.append("Median Units");
    constants.append("Ruin Probability");
    constants.append("5th % Time to Ruin");
    matrix.append(constants);

    for(int expensesNowPercent = 2; expensesNowPercent <= 6; ++expensesNowPercent)
    {
        DEBUG(expensesNowPercent);
        int term = 50;
        //DEBUG(term);
        double riskFreeReturn = 0;
        Vector<double> retirementExpenses = generateRetirementExpenses(expensesNowPercent, term);
        for(int percentStock = 0; percentStock <= 100; percentStock += 20)
        {
            DEBUG(percentStock);
            PortfolioSimulationResult result = performActuarialSimulation(makeMVOptimalRiskyAsset(100, percentStock/100.0, tangencyStockFraction, returnSpecifier), retirementExpenses, e);
            if(percentStock == 0)
            {
                riskFreeReturn = result.getMedian();
            }
            else
            {
                double riskFreeRank = result.riskFreeRank(riskFreeReturn);
                //DEBUG(riskFreeRank);
            }

            Vector<string> row;
            row.append(to_string(expensesNowPercent));
            row.append(to_string(percentStock));
            row.append(to_string(result.getMedian()));
            row.append(to_string(result.ruinChance));
            if(result.ruinChance > 0)
            {
                double targetPercentile = 0.05/result.ruinChance;
                if(targetPercentile <= 1)
                {
                    double percentile5TimeToRuin =
                        result.ruinPercentiles.getPercentile(targetPercentile);
                    row.append(to_string(percentile5TimeToRuin));
                }
                else row.append("N/A");
            }
            else row.append("N/A");
            matrix.append(row);
        }
    }
    //createCSV(matrix, "expense_ratio_simulation_actuarial.csv");
}

void testPortfolioSimulationNoCashFlowUpRebalanced()
{
    DEBUG("testPortfolioSimulationNoCashFlowUpRebalanced");
    Vector<Vector<string> > matrix;
    Vector<string> constants;
    constants.append("Stock Fraction");
    for(int i = 0; i < 3; ++i)
    {
        constants.append("Median");
        constants.append("0.1th %");
        constants.append("5th %");
        constants.append("95th %");
        constants.append("CRRA(3)");
        constants.append("CRRA(6)");
        constants.append("Risk-free Rank");
        constants.append("MaxDrawdown 5%");
        constants.append("MaxDrawdown 50%");
        constants.append("MaxDrawdown 95%");
    }
    constants.append("Put cost");
    matrix.append(constants);


    ReturnSpecifier returnSpecifier;//no tax
        int term = 30;
        DEBUG(term);
        double realRiskFreeReturn = 0;
        Vector<double> noSavings(term, 0);
        {
            PortfolioSimulationResult safe = performSimulation(RiskFreeAsset(returnSpecifier, 100), noSavings, 10000000);
            realRiskFreeReturn = safe.getMedian();
            DEBUG("risk-free");
            safe.debug();
        }

        double tangencyStockFraction = getTangencyStockFraction(returnSpecifier);


        for(int percentStock = 5; percentStock <= 100; percentStock += 5)
        {
            double stockFraction = percentStock/100.0;
            DEBUG(percentStock);
            Vector<PortfolioSimulationResult> results;
            {//optimal
                PortfolioSimulationResult result = performSimulation(makeMVOptimalRiskyAsset(100, stockFraction, tangencyStockFraction, returnSpecifier), noSavings, 10000000);
                results.append(result);
            }
            {//risk-free rebalanced
                PortfolioSimulationResult result = performSimulation(makeGeneralAsset(100, 0, 1 - stockFraction, 0, returnSpecifier), noSavings, 10000000);
                results.append(result);
            }
            {//risk-free up-rebalanced
                PortfolioSimulationResult result = performSimulation(makeGeneralAsset(100, 0, 0, 1 - stockFraction, returnSpecifier), noSavings, 10000000);
                results.append(result);
            }

            Vector<string> row;
            row.append(to_string(stockFraction));
            for(int i = 0; i < 3; ++i)
            {
                PortfolioSimulationResult result = results[i];
                double riskFreeRank = result.riskFreeRank(realRiskFreeReturn);
                row.append(to_string(result.getMedian()));
                row.append(to_string(result.percentiles.getPercentile(0.001)));
                row.append(to_string(result.percentiles.getPercentile(0.05)));
                row.append(to_string(result.percentiles.getPercentile(0.95)));
                row.append(to_string(expectedCRRACertaintyEquivalent(result.percentiles.values, 3)));
                row.append(to_string(expectedCRRACertaintyEquivalent(result.percentiles.values, 6)));
                row.append(to_string(riskFreeRank));
                row.append(to_string(result.maxDrawdownPercentiles.getPercentile(0.95)));
            }
            matrix.append(row);

        }
    createCSV(matrix, "protection_losses.csv");
}

void testPortfolioSimulationDCA()
{
    int term = 30;
    {
        DEBUG(term);
        DEBUG("401(k) DCA");
        ReturnSpecifier returnSpecifier;//no tax
        double tangencyStockFraction = getTangencyStockFraction(returnSpecifier);
        double initialValue = 0;
        double realRiskFreeReturn = 0;
        Vector<double> dcaSavings(term, 100);
        int joinIndex = 0;
        for(int percentStock = 0; percentStock <= 100; percentStock += 20)
        {
            DEBUG(percentStock);
            PortfolioSimulationResult result = performSimulation(makeMVOptimalRiskyAsset(initialValue, percentStock/100.0, tangencyStockFraction, returnSpecifier), dcaSavings);
            result.debug();
            if(percentStock == 0)
            {
                realRiskFreeReturn = result.getMedian();
            }
            else
            {
                double riskFreeRank = result.riskFreeRank(realRiskFreeReturn);
                DEBUG(riskFreeRank);
            }
        }
    }
}

void testPortfolioSimulationDCAVarious()
{
    Vector<Vector<string> > matrix;
    Vector<string> constants;
    constants.append("Strategy");
    constants.append("Stock Fraction");
    constants.append("Median");
    constants.append("5th %");
    constants.append("95th %");
    constants.append("CRRA(3)");
    constants.append("CRRA(6)");
    constants.append("MaxDrawdown 95%");
    matrix.append(constants);

    int term = 30;
    DEBUG(term);
    DEBUG("401(k) DCA");
    ReturnSpecifier returnSpecifier;//no tax
    double tangencyStockFraction = getTangencyStockFraction(returnSpecifier);
    double initialValue = 0;
    double units = ordinaryAnnuityPV(100, term, returnSpecifier.getRealRiskFreeRate());
    DEBUG(units);
    DEBUG(returnSpecifier.getRealRiskFreeRate());
    Vector<double> dcaSavings(term, 100), noSavings(term, 0);

    for(int percentStock = 0; percentStock <= 100; percentStock += 20)
    {
        DEBUG(percentStock);
        Vector<pair<string, PortfolioSimulationResult>> results;
        {
            DEBUG("Paid PV");
            PortfolioSimulationResult result = performSimulation(makeMVOptimalRiskyAsset(units, percentStock/100.0, tangencyStockFraction, returnSpecifier), noSavings);
            results.append({"Paid PV", result});
        }
        {
            DEBUG("Dynamic");
            TotalAssetPolicy policy = {percentStock/100.0, 0, dcaSavings[0], returnSpecifier.getRealRiskFreeRate(), 0, 1, term};
            PortfolioSimulationResult result = performSimulation(DynamicallyRebalancedAsset<TotalAssetPolicy>(initialValue, tangencyStockFraction, policy, 0, returnSpecifier), dcaSavings);
            results.append({"Dynamic", result});
        }
        {
            DEBUG("Constrained0209");
            TotalAssetPolicy policy = {percentStock/100.0, 0, dcaSavings[0], returnSpecifier.getRealRiskFreeRate(), 0.20, 0.9, term};
            PortfolioSimulationResult result = performSimulation(DynamicallyRebalancedAsset<TotalAssetPolicy>(initialValue, tangencyStockFraction, policy, 0, returnSpecifier), dcaSavings);
            results.append({"Constrained0209", result});
        }
        {
            DEBUG("Static");
            PortfolioSimulationResult result = performSimulation(makeMVOptimalRiskyAsset(initialValue, percentStock/100.0, tangencyStockFraction, returnSpecifier), dcaSavings);
            results.append({"Static", result});
        }

        for(int i = 0; i < results.getSize(); ++i)
        {
            auto const& result = results[i].second;
            Vector<string> row;
            row.append(results[i].first);
            row.append(to_string(percentStock/100.0));
            row.append(to_string(result.getMedian()));
            row.append(to_string(result.percentiles.getPercentile(0.05)));
            row.append(to_string(result.percentiles.getPercentile(0.95)));
            row.append(to_string(expectedCRRACertaintyEquivalent(result.percentiles.values, 3)));
            row.append(to_string(expectedCRRACertaintyEquivalent(result.percentiles.values, 6)));
            row.append(to_string(result.maxDrawdownPercentiles.getPercentile(0.95)));
            matrix.append(row);
        }

    }

    //createCSV(matrix, "dca_comparison.csv");
}

void testPortfolioSimulationDCAPaths()
{
    DEBUG("testPortfolioSimulationDCAPaths - slow");
    ReturnSpecifier returnSpecifier;
    int term = 30;
    pair<Vector<double>, double> best30 = findOptimalBondFractions(term,
        returnSpecifier, 0, 100, 0);
    DEBUG(best30.second);
    best30.first.debug();
}

double ARVASpending(double currentAmount, int yearToUse, double rate)
{
    return OrdinaryAnnuityCalculator::PMTBeginning(currentAmount, yearToUse, rate);
}
double MertonSpending(double currentAmount, int yearToUse, double rate)
{
    ReturnSpecifier returnSpecifier;
    //to find stock weight from geometric rate solve quadratic equation
    //pg = wp - (wq)^2/2


    double pg = rate - returnSpecifier.getRealRiskFreeRate();
    //DEBUG(rate);
    //DEBUG(returnSpecifier.getRealRiskFreeRate());
    double p = returnSpecifier.getStockReturn() - returnSpecifier.getRiskFreeRate();
    double q = returnSpecifier.getStockStd();
    //DEBUG(pg);
    //DEBUG(p);
    //DEBUG(q);
    double w = (p - sqrt(p * p - 4 * (q*q/2) * pg))/(2 * (q*q/2));
    //DEBUG(w);
    //DEBUG(w*p - (w*q)*(w*q)/2);

    double s = p/q;
    double kelly = s/q;
    double a = kelly/w;
    //DEBUG(a);
    if(!isfinite(a)) a = 100000000;//large sentinel
    double rho = 0;
    double v = (rho + (a - 1)*(returnSpecifier.getRealRiskFreeRate() + s*s/a/2))/a;
    //DEBUG(v);
    double merton = v > 0 ? v/(1 - exp(-v * yearToUse)) : 1.0/yearToUse;
    //DEBUG(merton);
    //if(pg > 0) system("PAUSE");
    return merton * currentAmount;
}
template<typename FINANCIAL_ASSET>
PortfolioSimulationResult performActuarialARVASimulation(FINANCIAL_ASSET const&
    initialFinancialAsset, double realDiscountRate,
    Vector<double> const& nextYearSurvivalProbabilities,
    bool useARVA, int nSimulations = 1000000)
{
    assert(nextYearSurvivalProbabilities.getSize() > 0 && nSimulations > 0);
    double initialSpending = useARVA ? ARVASpending(initialFinancialAsset.getValue(),
        nextYearSurvivalProbabilities.getSize(), realDiscountRate) :
        MertonSpending(initialFinancialAsset.getValue(),
        nextYearSurvivalProbabilities.getSize(), realDiscountRate);
    IncrementalStatistics result;
    //DEBUG(initialSpending);
    Vector<double> values;
    for(int i = 0; i < nSimulations; ++i)
    {
        FINANCIAL_ASSET financialAsset = initialFinancialAsset;
        for(int step = 0; step < nextYearSurvivalProbabilities.getSize();++step)
        {
            if(financialAsset.getValue() < 0)
            {//artificially ruined (sim first, substract later)
                continue;
            }
            //DEBUG(financialAsset.getValue());
            //DEBUG(nextYearSurvivalProbabilities.getSize()-step);
            double spending = useARVA ? ARVASpending(financialAsset.getValue(),
                nextYearSurvivalProbabilities.getSize()-step,
                realDiscountRate) :
                MertonSpending(financialAsset.getValue(),
                nextYearSurvivalProbabilities.getSize()-step,
                realDiscountRate);
            //DEBUG(spending);
            //values.append(spending/initialSpending);
            values.append(spending);
            financialAsset.setValue(financialAsset.getValue());
            financialAsset.simulateStep(-spending, step);
            if(nextYearSurvivalProbabilities[step] < GlobalRNG().uniform01())
            {//died, done spending
                break;
            }
        }
    }
    return PortfolioSimulationResult(values, Vector<double>(), Vector(1, 0.0));
}

void testARVARetirement()
{
    bool useARVA = true;
    ReturnSpecifier returnSpecifier;//no tax
    double tangencyStockFraction = getTangencyStockFraction(returnSpecifier);
    {
        int term = 50;
        //int term = 1000;
        DEBUG(term);
        double riskFreeReturn = 0;
        for(int percentStock = 0; percentStock <= 100; percentStock += 20)
        {
            DEBUG(percentStock);
            double discountRate;
            if(percentStock/100.0 < tangencyStockFraction)
            {
                LognormalDistribution ld = StockBondAsset::makeReturnDistribution(returnSpecifier, 1 - tangencyStockFraction);
                double tangecyFraction = (percentStock/100.0)/tangencyStockFraction;
                double mean = tangecyFraction * (ld.getMean()-1) + (1 - tangecyFraction) * returnSpecifier.getRealRiskFreeRate();
                double stdev = tangecyFraction * ld.getStdev();
                discountRate = mean - stdev * stdev/2;
            }
            else discountRate = StockBondAsset::makeReturnDistribution(returnSpecifier, 1 - percentStock/100.0).getMedian() - 1;
            DEBUG(discountRate);
            double initialSpending;
            if(useARVA)
            {
                initialSpending = ARVASpending(100, term, discountRate);
            }
            else
            {
                initialSpending = MertonSpending(100, term, discountRate);
            }
            DEBUG(initialSpending);
            PortfolioSimulationResult result = performActuarialARVASimulation(
                makeMVOptimalRiskyAsset(100, percentStock/100.0, tangencyStockFraction, returnSpecifier),
                discountRate, Vector<double>(term, 1), useARVA);
            result.debug();
            if(percentStock == 0)
            {
                riskFreeReturn = result.getMedian();
            }
            else
            {
                double riskFreeRank = result.riskFreeRank(riskFreeReturn);
                DEBUG(riskFreeRank);
            }
        }
    }
}

void testMinVarCalculationTotal(ReturnSpecifier returnSpecifier =
    ReturnSpecifier(), bool adjustWithTangency = true, double minStock = 0)
{
    DEBUG("testMinVarCalculationTotal");

    //2 less
    //returnSpecifier.stockReturn -= 0.02;
    //2 more
    //returnSpecifier.stockReturn += 0.02;
    MeanVariancePortfolio mvp(makeStockBondMVP(returnSpecifier));

    DEBUG("Optimal Sharpe");
    pair<double, Vector<double> > tangency = mvp.findOptimalSharpeWeights(returnSpecifier.getRiskFreeRate());
    DEBUG(tangency.second[0]);
    pair<double, double> mqTangency = mvp.evaluate(tangency.second);
    DEBUG(mqTangency.first);
    DEBUG(mqTangency.second);
    DEBUG((mqTangency.first - returnSpecifier.getRiskFreeRate())/mqTangency.second);
    double tangencySharpe = (mqTangency.first - returnSpecifier.getRiskFreeRate())/mqTangency.second;

    Vector<pair<double, Vector<double> > > frontier =
        mvp.findOptimalPortfolioWeightsRange(21);
    DEBUG(returnSpecifier.getRiskFreeRate());
    double stockSharpe = (returnSpecifier.getStockReturn() - returnSpecifier.getRiskFreeRate())/returnSpecifier.getStockStd();
    //DEBUG(stockSharpe);
    for(int i = 0; i < frontier.getSize(); ++i)
    {
        if(frontier[i].second[0] < minStock) continue;
        frontier[i].second.debug();
        pair<double, double> mq;
        if(tangency.second[0] > frontier[i].second[0])
        {//stock below tangency, will mix with risk-free to get target stock
            double tangencyFraction = max(frontier[i].second[0]/tangency.second[0],
                numeric_limits<double>::epsilon());//to avoid numerical issues
            mq.first = mqTangency.first * tangencyFraction +
                returnSpecifier.getRiskFreeRate() * (1 - tangencyFraction);
            mq.second = mqTangency.second * tangencyFraction;
        }
        else
        {//stocks at or above tangency
            mq = mvp.evaluate(frontier[i].second);
        }
        DEBUG(mq.first);
        DEBUG(mq.second);
        DEBUG(mq.second/mqTangency.second);
        DEBUG(mq.first - 2 * mq.second);
        DEBUG(mq.first + 2 * mq.second);

        LognormalDistribution lognormal(1 + mq.first, mq.second);
        DEBUG(lognormal.getMedian() - 1);
        pair<double, double> predictionInterval10 = (lognormal*=10).getMedianPredictionInterval(10);
        DEBUG(predictionInterval10.first - 1);
        DEBUG(predictionInterval10.second - 1);
        pair<double, double> predictionInterval30 = (lognormal*=3).getMedianPredictionInterval(30);
        DEBUG(predictionInterval30.first - 1);
        DEBUG(predictionInterval30.second - 1);
        double CRRA3Score = mq.first - 3 * mq.second * mq.second/2;
        DEBUG(CRRA3Score);

        DEBUG("Sharpe");
        DEBUG((mq.first - returnSpecifier.getRiskFreeRate())/mq.second);
        double sharpe = (mq.first - returnSpecifier.getRiskFreeRate())/mq.second;
        //DEBUG(tangencySharpe - sharpe);
        DEBUG(1 - (tangencySharpe - sharpe)/(tangencySharpe - stockSharpe));
    }
}

void testMinVarCalculationMakeTable(ReturnSpecifier returnSpecifier =
    ReturnSpecifier(), bool adjustWithTangency = true, double minStock = 0.01)
{
    DEBUG("testMinVarCalculationMakeTable");

    Vector<Vector<string> > matrix;
    Vector<string> constants;
    constants.append("Stock %");
    constants.append("Geometric Real Return");
    constants.append("5th %");
    constants.append("95th %");
    constants.append("Sharpe Ratio");
    constants.append("Sharpe Ratio Fraction To Tangency Captured");
    constants.append("10-year 5th %");
    constants.append("10-year 95th %");
    constants.append("30-year 5th %");
    constants.append("30-year 95th %");
    matrix.append(constants);

    //2 less
    //returnSpecifier.stockReturn -= 0.02;
    //2 more
    //returnSpecifier.stockReturn += 0.02;
    MeanVariancePortfolio mvp(makeStockBondMVP(returnSpecifier));

    DEBUG("Optimal Sharpe");
    pair<double, Vector<double> > tangency = mvp.findOptimalSharpeWeights(returnSpecifier.getRiskFreeRate());
    DEBUG(tangency.second[0]);
    pair<double, double> mqTangency = mvp.evaluate(tangency.second);
    DEBUG(mqTangency.first);
    DEBUG(mqTangency.second);
    DEBUG((mqTangency.first - returnSpecifier.getRiskFreeRate())/mqTangency.second);
    double tangencySharpe = (mqTangency.first - returnSpecifier.getRiskFreeRate())/mqTangency.second;

    Vector<pair<double, Vector<double> > > frontier =
        mvp.findOptimalPortfolioWeightsRange(21);
    DEBUG(returnSpecifier.getRiskFreeRate());
    double stockSharpe = (returnSpecifier.getStockReturn() - returnSpecifier.getRiskFreeRate())/returnSpecifier.getStockStd();
    //DEBUG(stockSharpe);
    for(int i = 0; i < frontier.getSize(); ++i)
    {
        if(frontier[i].second[0] < minStock) continue;
        frontier[i].second.debug();
        pair<double, double> mq;
        if(tangency.second[0] > frontier[i].second[0])
        {//stock below tangency, will mix with risk-free to get target stock
            double tangencyFraction = max(frontier[i].second[0]/tangency.second[0],
                numeric_limits<double>::epsilon());//to avoid numerical issues
            mq.first = mqTangency.first * tangencyFraction +
                returnSpecifier.getRiskFreeRate() * (1 - tangencyFraction);
            mq.second = mqTangency.second * tangencyFraction;
        }
        else
        {//stocks at or above tangency
            mq = mvp.evaluate(frontier[i].second);
        }
        DEBUG(mq.first);
        DEBUG(mq.second);
        DEBUG(mq.second/mqTangency.second);
        DEBUG(mq.first - 2 * mq.second);
        DEBUG(mq.first + 2 * mq.second);

        LognormalDistribution lognormal(1 + mq.first - returnSpecifier.getInflationRate(), mq.second);

        Vector<string> row;
        row.append(to_string(frontier[i].second[0]));
        row.append(to_string(lognormal.getMedian() - 1));

        pair<double, double> predictionInterval1 = lognormal.getMedianPredictionInterval(1);
        DEBUG(lognormal.getMedian() - 1);
        pair<double, double> predictionInterval10 = (lognormal*=10).getMedianPredictionInterval(10);
        DEBUG(predictionInterval10.first - 1);
        DEBUG(predictionInterval10.second - 1);
        pair<double, double> predictionInterval30 = (lognormal*=3).getMedianPredictionInterval(30);
        DEBUG(predictionInterval30.first - 1);
        DEBUG(predictionInterval30.second - 1);
        double CRRA3Score = mq.first - 3 * mq.second * mq.second/2;
        DEBUG(CRRA3Score);

        DEBUG("Sharpe");
        DEBUG((mq.first - returnSpecifier.getRiskFreeRate())/mq.second);
        double sharpe = (mq.first - returnSpecifier.getRiskFreeRate())/mq.second;
        //DEBUG(tangencySharpe - sharpe);
        DEBUG(1 - (tangencySharpe - sharpe)/(tangencySharpe - stockSharpe));


        row.append(to_string(predictionInterval1.first-1));
        row.append(to_string(predictionInterval1.second-1));
        row.append(to_string((mq.first - returnSpecifier.getRiskFreeRate())/mq.second));
        row.append(to_string(1 - (tangencySharpe - sharpe)/(tangencySharpe - stockSharpe)));
        row.append(to_string(predictionInterval10.first-1));
        row.append(to_string(predictionInterval10.second-1));
        row.append(to_string(predictionInterval30.first-1));
        row.append(to_string(predictionInterval30.second-1));
        matrix.append(row);
    }

    //createCSV(matrix, "allocation_table.csv");
}

void testMakeReturnChartByAllocation(ReturnSpecifier returnSpecifier =
    ReturnSpecifier(), bool adjustWithTangency = true, double minStock = 0.1)
{
    DEBUG("testMakeReturnChartByAllocation");


    Vector<Vector<string> > matrix;
    Vector<string> constants;
    constants.append("Stock %");
    constants.append("Year");
    constants.append("5th %");
    constants.append("50th %");
    constants.append("95th %");
    matrix.append(constants);

    MeanVariancePortfolio mvp(makeStockBondMVP(returnSpecifier));

    DEBUG("Optimal Sharpe");
    pair<double, Vector<double> > tangency = mvp.findOptimalSharpeWeights(returnSpecifier.getRiskFreeRate());
    pair<double, double> mqTangency = mvp.evaluate(tangency.second);
    double tangencySharpe = (mqTangency.first - returnSpecifier.getRiskFreeRate())/mqTangency.second;

    Vector<pair<double, Vector<double> > > frontier =
        mvp.findOptimalPortfolioWeightsRange(5);
    DEBUG(returnSpecifier.getRiskFreeRate());
    double stockSharpe = (returnSpecifier.getStockReturn() - returnSpecifier.getRiskFreeRate())/returnSpecifier.getStockStd();
    //DEBUG(stockSharpe);
    int years = 30;
    for(int i = 0; i < frontier.getSize(); ++i)
    {
        if(frontier[i].second[0] < minStock) continue;
        frontier[i].second.debug();
        pair<double, double> mq;
        if(tangency.second[0] > frontier[i].second[0])
        {//stock below tangency, will mix with risk-free to get target stock
            double tangencyFraction = max(frontier[i].second[0]/tangency.second[0],
                numeric_limits<double>::epsilon());//to avoid numerical issues
            mq.first = mqTangency.first * tangencyFraction +
                returnSpecifier.getRiskFreeRate() * (1 - tangencyFraction);
            mq.second = mqTangency.second * tangencyFraction;
        }
        else
        {//stocks at or above tangency
            mq = mvp.evaluate(frontier[i].second);
        }

        {
            Vector<string> row;
            row.append(to_string(frontier[i].second[0]));
            row.append("0");
            row.append("1");
            row.append("1");
            row.append("1");
            matrix.append(row);
        }



        for(int j = 1; j <= years; ++j)
        {
            LognormalDistribution lognormal(1 + mq.first - returnSpecifier.getInflationRate(), mq.second);
            lognormal*=j;

            pair<double, double> predictionInterval10 = lognormal.getMedianPredictionInterval();

            Vector<string> row;
            row.append(to_string(frontier[i].second[0]));
            row.append(to_string(j));
            row.append(to_string(predictionInterval10.first));
            row.append(to_string(lognormal.getMedian()));
            row.append(to_string(predictionInterval10.second));
            matrix.append(row);
        }
    }

    //createCSV(matrix, "path_simulations.csv");
}
void testMinVarCalculationCurrentAverageReturns()
{
    DEBUG("testMinVarCalculationCurrentHistoricalStockReturns");
    //set stock return to 6% real, adjust to nominal arithmetic
    ReturnSpecifier returnSpecifier;
    returnSpecifier.stockReturn = 0.06 + returnSpecifier.getInflationRate() +
        returnSpecifier.getStockStd() * returnSpecifier.getStockStd()/2;
    testMinVarCalculationTotal(returnSpecifier, false);
}
void testMinVarCalculationCurrent()
{
    DEBUG("testMinVarCalculationCurrent");
    ReturnSpecifier returnSpecifier;
    testMinVarCalculationTotal(returnSpecifier, false);
}

void testMinVarCalculationCurrentTaxes()
{
    DEBUG("testMinVarCalculationCurrentTaxes");
    ReturnSpecifier returnSpecifier(0.15, 0.20, 0.00);
    testMinVarCalculationTotal(returnSpecifier, false, 0.99);
}

pair<double, double> estimateRealReturnRate(double taxEfficientStock, double taxEfficientBond,
    double taxEfficientRF, ReturnSpecifier const& taxable,
    ReturnSpecifier const& taxEfficient, double taxableStock, double taxableRF)
{
    double totalStock = taxEfficientStock + taxableStock, totalRF = taxEfficientRF + taxableRF, totalBond = taxEfficientBond;
    MeanVariancePortfolio stockBondMVP = makeStockBondMVP(taxEfficient);
    Vector<double> weights;
    double bondFraction = totalBond/(totalBond + totalStock);
    weights.append(1 - bondFraction);
    weights.append(bondFraction);
    pair<double, double> meanstd = stockBondMVP.evaluate(weights);
    double rfFraction = totalRF/(totalRF + totalBond + totalStock);
    DEBUG(rfFraction);
    DEBUG(bondFraction);
    DEBUG(totalStock/(totalRF + totalBond + totalStock));
    meanstd.first = meanstd.first * (1 - rfFraction) + rfFraction * taxEfficient.getRiskFreeRate();
    meanstd.second *= (1 - rfFraction);
    //DEBUG(meanstd.first);
    DEBUG(meanstd.second);
    DEBUG((meanstd.first - taxEfficient.getRiskFreeRate())/meanstd.second);
    double geometricReturn = LognormalDistribution(1 + meanstd.first, meanstd.second).getMedian() - 1;
    //DEBUG(geometricReturn);
    double taxes = (taxableStock * (taxEfficient.getStockReturn() - taxable.getStockReturn()) +
        taxableRF * (taxEfficient.getRiskFreeRate() - taxable.getRiskFreeRate()))/
        (totalRF + totalBond + totalStock);
    return {geometricReturn - taxes - taxEfficient.getInflationRate(), meanstd.second};
}

void testMultiAccountReturns()
{
    DEBUG("total with all accounts");
    ReturnSpecifier taxable(0.15, 0.2, 0), taxEfficient;
    double totalTaxable = 500, totalTaxEfficient = 500, totalEF = 200, bondFraction = 0.2;
    DEBUG(estimateRealReturnRate((totalTaxEfficient - totalEF) * (1 - bondFraction), (totalTaxEfficient - totalEF) * bondFraction, totalEF, taxable, taxEfficient, totalTaxable, 0).first);
}

void testEstimateStockBondRebalance()
{
    DEBUG("no taxes case");
    double stockBondPercents = estimateStockBondSalesForRebalance();
    DEBUG(stockBondPercents);
    DEBUG("taxes case");
    stockBondPercents = estimateStockBondSalesForRebalance(ReturnSpecifier(0.15, 0.2, 0));
    DEBUG(stockBondPercents);
}

double evaluateLognormalCDF(double relativeReturnTarget, double mean,
    double stdev, int nPeriods)
{// mean is 1-based
    assert(relativeReturnTarget > 0 && mean > 0 && stdev > 0 && nPeriods > 0);
    LognormalDistribution ld(mean, stdev);
    ld *= nPeriods;
    //DEBUG(log(relativeReturnTarget));
    //DEBUG(ld.getMu());
    //DEBUG(ld.getQ());
    double z = (log(relativeReturnTarget) - ld.getMu())/ld.getQ();
    //DEBUG(z);
    return approxNormalCDF(z);
}

pair<double, Vector<double> > stockOfMaxHitChance(double relativeReturnTarget,
    MeanVariancePortfolio const& mvp, int nPeriods)
{
    auto functor = [relativeReturnTarget, nPeriods](double mean, double stdev)
    {//min value best functor
        return -(1 - evaluateLognormalCDF(relativeReturnTarget, 1 + mean, stdev,
            nPeriods));
    };
    return mvp.findOptimalFrontierPoint(functor);
}

void testGoalSeek()
{
    DEBUG("testGoalSeek");
    //with 100% stocks
    DEBUG(1 - evaluateLognormalCDF(1.04, 1.08, 0.17, 1));//chance to at least outperform t-bill in 1 year
    DEBUG(1 - evaluateLognormalCDF(1.5, 1.08, 0.17, 5));//chance to at least 50% in 5
    DEBUG(1 - evaluateLognormalCDF(10, 1.08, 0.17, 30));//chance to at least 10x wealth in 30
    DEBUG(1 - evaluateLognormalCDF(1, 1.08, 0.17, 10));//chance to at least gain nothing in 10
    DEBUG(1 - evaluateLognormalCDF(pow(1.04, 10), 1.08, 0.17, 10));//chance to at least outperform t-bill in 10 years
    DEBUG(1 - evaluateLognormalCDF(pow(1.04, 30), 1.08, 0.17, 30));//chance to at least outperform t-bill in 30 years
    DEBUG(evaluateLognormalCDF(0.8, 1.08, 0.17, 1));//chance to lose at least 20% in 1 year
    DEBUG(1 - evaluateLognormalCDF(1.2, 1.08, 0.17, 1));//chance to gain at least 20% in 1 year
    DEBUG(evaluateLognormalCDF(0.8, 1.08, 0.17, 5));//chance to lose at least 20% in 5 years
    DEBUG(evaluateLognormalCDF(1, 1.08, 0.17, 1));//chance to lose in 1 year
    DEBUG(evaluateLognormalCDF(1, 1.08, 0.17, 5));//chance to lose in 5 years
    DEBUG(evaluateLognormalCDF(0.5, 1.08, 0.17, 1));//chance to lose at least 50% in 1 year
    DEBUG(evaluateLognormalCDF(0.5, 1.08, 0.17, 5));//chance to lose at least 50% in 5 years
    //with 60% stocks
    DEBUG(1 - evaluateLognormalCDF(1.04, 1.06, 0.11, 1));//chance to at least outperform t-bill in 1 year
    DEBUG(1 - evaluateLognormalCDF(10, 1.06, 0.11, 30));//chance to at least 10x wealth in 30
    DEBUG(1 - evaluateLognormalCDF(1, 1.06, 0.11, 10));//chance to at least gain nothing in 10
    DEBUG(1 - evaluateLognormalCDF(pow(1.04, 10), 1.06, 0.11, 10));//chance to at least outperform t-bill in 10 years
    DEBUG(1 - evaluateLognormalCDF(pow(1.04, 30), 1.06, 0.11, 30));//chance to at least outperform t-bill in 30 years
    //max chance to double in 20 years
    MeanVariancePortfolio mvp = makeStockBondMVP(ReturnSpecifier());
    pair<double, Vector<double> > result = stockOfMaxHitChance(2, mvp, 20);
    DEBUG(result.first);
    result.second.debug();
}

void testGoalSeekTargetDatePath()
{
    DEBUG("testGoalSeekTargetDatePath");
    ReturnSpecifier returnSpecifier;

    Vector<Vector<string> > matrix;
    int n = 7;
    Vector<string> constants(n, "");
    constants[0] = "Multiple";
    constants[1] = "Years";
    constants[2] = "20/80";
    constants[3] = "40/60";
    constants[4] = "60/40";
    constants[5] = "80/20";
    constants[6] = "100/0";
    matrix.append(constants);

    MeanVariancePortfolio mvp(makeStockBondMVP(returnSpecifier));
    pair<double, Vector<double> > tangency = mvp.findOptimalSharpeWeights(returnSpecifier.getRiskFreeRate());
    pair<double, double> mqTangency = mvp.evaluate(tangency.second);

    Vector<pair<double, Vector<double> > > frontier =
        mvp.findOptimalPortfolioWeightsRange(6);
    for(double multiple = 0.5; multiple < 8.1; multiple *= 2)
    {
        for(int nYears = 1; nYears <= 30; nYears += 5)
        {
            Vector<string> row;
            row.append(to_string(multiple));
            row.append(to_string(nYears));
            double target = pow(1 + returnSpecifier.getRiskFreeRate(), nYears);
            for(int i = 1; i < frontier.getSize(); ++i)
            {
                pair<double, double> mq;
                if(tangency.second[0] > frontier[i].second[0])
                {//stock below tangency, will mix with risk-free to get target stock
                    double tangencyFraction = max(frontier[i].second[0]/tangency.second[0],
                        numeric_limits<double>::epsilon());//to avoid numerical issues
                    mq.first = mqTangency.first * tangencyFraction +
                        returnSpecifier.getRiskFreeRate() * (1 - tangencyFraction);
                    mq.second = mqTangency.second * tangencyFraction;
                }
                else
                {//stocks at or above tangency
                    mq = mvp.evaluate(frontier[i].second);
                }
                double successProbability = 1 - evaluateLognormalCDF(target * multiple, 1 + mq.first, mq.second, nYears);//chance to at least outperform t-bill multiple
                row.append(to_string(successProbability));

            }
            matrix.append(row);
            if(nYears == 1) nYears = 0;
        }
    }
    //createCSV(matrix, "target_fund_success_probs.csv");
}

void testTotalBondEffect()
{
    DEBUG("testTotalBondEffect");
    ReturnSpecifier returnSpecifier;
    MeanVariancePortfolio mvp(makeStockBondMVP(returnSpecifier));


    pair<double, Vector<double> > tangency = mvp.findOptimalSharpeWeights(returnSpecifier.getRiskFreeRate());
    pair<double, double> mqTangency = mvp.evaluate(tangency.second);

    Vector<pair<double, Vector<double> > > frontier =
        mvp.findOptimalPortfolioWeightsRange(21);

    Vector<Vector<string> > matrix;
    int n = 7;
    Vector<string> constants(n, "");
    constants[0] = "Stock Fraction";
    constants[1] = "Geometric Return Optimal";
    constants[2] = "Std Optimal Mix";
    constants[3] = "Geometric Return Stock-Bond";
    constants[4] = "Std Stock-Bond";
    constants[5] = "Geometric Return Stock-Risk-free";
    constants[6] = "Std Stock-Bond";
    matrix.append(constants);


    for(int i = 0; i < frontier.getSize(); ++i)
    {
        pair<double, double> mq;
        double stockFraction = frontier[i].second[0];
        if(tangency.second[0] > stockFraction)
        {//stock below tangency, will mix with risk-free to get target stock
            double tangencyFraction = max(stockFraction/tangency.second[0],
                numeric_limits<double>::epsilon());//to avoid numerical issues
            mq.first = mqTangency.first * tangencyFraction +
                returnSpecifier.getRiskFreeRate() * (1 - tangencyFraction);
            mq.second = mqTangency.second * tangencyFraction;
        }
        else
        {//stocks at or above tangency
            mq = mvp.evaluate(frontier[i].second);
        }
        LognormalDistribution lognormal(1 + mq.first, mq.second);

        Vector<string> row;
        row.append(to_string(stockFraction));
        //optimal
        row.append(to_string(lognormal.getMedian() - 1));
        row.append(to_string(mq.second));
        //stock-bond
        mq = mvp.evaluate(frontier[i].second);
        lognormal = LognormalDistribution(1 + mq.first, mq.second);
        row.append(to_string(lognormal.getMedian() - 1));
        row.append(to_string(mq.second));
        //stock-risk-free
        mq.first = returnSpecifier.getStockReturn() * stockFraction +
                returnSpecifier.getRiskFreeRate() * (1 - stockFraction);
        mq.second = returnSpecifier.getStockStd() * stockFraction;
        lognormal = LognormalDistribution(1 + mq.first, mq.second);
        row.append(to_string(lognormal.getMedian() - 1));
        row.append(to_string(mq.second));
        matrix.append(row);
    }
    // createCSV(matrix, "mean_var_total_bond_effect.csv");
}

void testPickedPortfolios()
{
    DEBUG("testPickedPortfolios");

    Vector<Vector<string> > matrix;
    Vector<string> constants;
    constants.append("Risk Group");
    constants.append("Stocks/Bonds/Risk-free");
    constants.append("Median");
    constants.append("5th %");
    constants.append("95th %");
    constants.append("CRRA (3)");
    constants.append("CRRA (6)");
    constants.append("Risk-free Rank");
    constants.append("95th % Max Drawdown");
    matrix.append(constants);

    ReturnSpecifier returnSpecifier;//no tax
    for(int term = 5; term <= 20; term *= 2)
    {
        DEBUG(term);
        Vector<string> termRow(constants.getSize(), "");
        termRow[0] = "Years";
        termRow[1] = to_string(term);
        matrix.append(termRow);

        double realRiskFreeReturn = 0;
        Vector<double> noSavings(term, 0);
        {
            DEBUG("risk-free");
            PortfolioSimulationResult safe = performSimulation(RiskFreeAsset(returnSpecifier, 100), noSavings, 10000000);
            realRiskFreeReturn = safe.getMedian();
        }
        for(int i = 0; i < 7; ++i)
        {
            string name = "", allocation = "";
            double bondFraction = 0, riskFreeFraction = 0;
            if(i == 0)
            {
                name = "High1_15";
                allocation = "100/0/0";
            }
            else if(i == 1)
            {
                name = "High2";
                allocation = "90/10/0";
                bondFraction = 0.1;
            }
            else if(i == 2)
            {
                name = "Start70";
                allocation = "70/30/0";
                bondFraction = 0.3;
                riskFreeFraction = 0.1;
            }
            else if(i == 3)
            {
                name = "Medium3";
                allocation = "60/30/10";
                bondFraction = 0.3;
                riskFreeFraction = 0.1;
            }
            else if(i == 4)
            {
                name = "Medium4";
                allocation = "40/20/40";
                bondFraction = 0.2;
                riskFreeFraction = 0.4;
            }
            else if(i == 5)
            {
                name = "Low6";
                allocation = "30/15/55";
                bondFraction = 0.15;
                riskFreeFraction = 0.55;
            }
            else if(i == 6)
            {
                name = "Low8";
                allocation = "20/10/70";
                bondFraction = 0.1;
                riskFreeFraction = 0.7;
            }
            DEBUG(name);
            PortfolioSimulationResult result = performSimulation(makeGeneralAsset(100, bondFraction, riskFreeFraction, 0, returnSpecifier), noSavings, 10000000);
            double riskFreeRank = result.riskFreeRank(realRiskFreeReturn);

            Vector<string> row;
            row.append(name);
            row.append(allocation);
            row.append(to_string(result.getMedian()));
            row.append(to_string(result.percentiles.getPercentile(0.05)));
            row.append(to_string(result.percentiles.getPercentile(0.95)));
            row.append(to_string(expectedCRRACertaintyEquivalent(result.percentiles.values, 3)));
            row.append(to_string(expectedCRRACertaintyEquivalent(result.percentiles.values, 6)));
            row.append(to_string(riskFreeRank));
            row.append(to_string(result.maxDrawdownPercentiles.getPercentile(0.95)));
            matrix.append(row);
        }
    }
    createCSV(matrix, "picked_portfolios_compare.csv");
}


void testOptionInsurance()
{
    DEBUG("testOptionInsurance");
    ReturnSpecifier returnSpecifier;//no tax
    double price = 1, strikePrice = 1, strikeTime = 1;
    double geometricReturn = returnSpecifier.getStockReturn() - returnSpecifier.getStockStd() * returnSpecifier.getStockStd()/2;
    for(int i = 0; i < 4; ++i)
    {
        if(i == 0) strikePrice = 1;
        else if(i == 1) strikePrice = 0.75;
        else if(i == 2) strikePrice = 0.5;
        else if(i == 3) strikePrice = 0.25;
        DEBUG(strikePrice);
        double putPrice = priceEuropeanPut(price, returnSpecifier.getRiskFreeRate(), returnSpecifier.getStockStd(), strikePrice, strikeTime);
        DEBUG(putPrice);
        DEBUG(geometricReturn - putPrice);
        double riskFreeEquivalentFraction = strikePrice;
        double stockRiskFreeReturn = geometricReturn * (1 - riskFreeEquivalentFraction) + returnSpecifier.getRiskFreeRate() * riskFreeEquivalentFraction;
        DEBUG(stockRiskFreeReturn);
    }
}

//assume 1 units of consumption
double calculateRetirementNeedTo120(int age,
    int retirementEmergencyFundUnits, double realRiskFreeRate)
{
    return ordinaryAnnuityPV(1, 120 - age, realRiskFreeRate) + retirementEmergencyFundUnits;
}
struct FinancialIndependenceSimRecorder
{
    Vector<double> stepsToIndepence;
    double realRiskFreeRate;
    int startAge, retirementEmergencyFundUnits;
    Vector<bool> simulationRecorded;
    int maxSteps;
    bool use25X;
    void operator()(int simulation, int step, double realWealth)
    {
        if(!simulationRecorded[simulation])
        {
            if(step < maxSteps)
            {
                int age = startAge + step;
                //DEBUG(age);
                //DEBUG(realWealth);
                //DEBUG(calculateRetirementNeedTo120(age, retirementEmergencyFundUnits, realRiskFreeRate));
                if(use25X ? (realWealth >= 25) : (realWealth >= calculateRetirementNeedTo120(age,
                    retirementEmergencyFundUnits, realRiskFreeRate)))
                {//got enough
                    stepsToIndepence.append(age);
                    simulationRecorded[simulation] = true;
                }
            }
            else
            {//reach end of sim
                simulationRecorded[simulation] = true;
                stepsToIndepence.append(startAge + maxSteps);
            }
        }
    }
};
PercentileManager performFinancialIndepenceSimulation(int age,
    double currentMultiple, double savingsMultiple,
    int retirementEmergencyFundUnits, double stockFraction, bool use25X, int nSimulations = 1000000)
{//calculate distribution of age hitting financial independence, assuming have
    //constant real consumption and will use guaranteed funding to 120
    ReturnSpecifier returnSpecifier;//no tax
    Vector<double> dcaSavings(120 - age, savingsMultiple);//1 units saved, can work until 120
    FinancialIndependenceSimRecorder sr;
    sr.realRiskFreeRate = returnSpecifier.getRealRiskFreeRate();
    sr.maxSteps = dcaSavings.getSize();
    sr.simulationRecorded = Vector<bool>(nSimulations, false);
    sr.startAge = age;
    use25X = use25X;
    sr.retirementEmergencyFundUnits = retirementEmergencyFundUnits;
    //ignore sim result, only use recorder output
    performActuarialSimulation(makeMVOptimalRiskyAsset(currentMultiple, stockFraction, getTangencyStockFraction(returnSpecifier), returnSpecifier), sr,//replace by optimal asset given tangency
        dcaSavings, Vector<double>(dcaSavings.getSize(), 1), nSimulations);
    return sr.stepsToIndepence;
}

void testFinancialIndependenceDistribution()
{
    ReturnSpecifier returnSpecifier;//no tax
    double tangencyStockFraction = getTangencyStockFraction(returnSpecifier);
    int age = 22;//college grad
    JointSurvivalEstimator e;//assume couple of same age for actuarial
    e.addPerson(convertToSurvivalProbabilities(getSSAMaleDeathProbabilities()), age);
    e.addPerson(convertToSurvivalProbabilities(getSSAFemaleDeathProbabilities()), age);
    double currentMultiple = 0;//save nothing so far
    double savingsMultiple = 0.5;//saving half of income
    int retirementEmergencyFundUnits = 5;//need 5 years of retirement emergency fund
    bool use25X = true;
    for(int percentStock = 0; percentStock <= 100; percentStock += 20)
    {
        DEBUG(percentStock);
        PercentileManager percentiles = performFinancialIndepenceSimulation(age, currentMultiple, savingsMultiple, retirementEmergencyFundUnits, percentStock/100.0, use25X);
        DEBUG(percentiles.getPercentile(0.05));
        DEBUG(percentiles.getPercentile(0.25));
        DEBUG(percentiles.getPercentile(0.5));
        DEBUG(percentiles.getPercentile(0.75));
        DEBUG(percentiles.getPercentile(0.95));
        //calculate total-life consumption utility
        IncrementalStatistics utility, utilityActuarial;
        Vector<double> CRRAValues;
        for(int i = 0; i < percentiles.values.getSize(); ++i)
        {//double consumption in retirement years, no actuarial adjustment
            double ageFI = percentiles.values[i];
            double lifetimeValue = (ageFI - age) + 2 * (120 - ageFI);
            utility.addValue(lifetimeValue/(120 - age));
            bool survivedToRetirement = true;
            for(int ageJ = age; ageJ < 120; ++ageJ)
            {
                if(getFutureSurvivalProbability(e, ageJ - age) < GlobalRNG().uniform01())
                {
                    if(ageJ <= ageFI) survivedToRetirement = false;
                    break;
                }
                utilityActuarial.addValue(ageJ <= ageFI ? 1 : 2);
                if(ageJ <= ageFI) CRRAValues.append(1);
            }
            //CRRA utility assuming will spend using ARVA
            if(use25X && survivedToRetirement)
            {
                //calculate survivalProbabilities
                Vector<double> nextYearSurvivalProbabilities;
                for(int ageJ = ageFI + 1; ageJ < 120; ++ageJ)
                {
                    nextYearSurvivalProbabilities.append(getFutureSurvivalProbability(e, ageJ - age));
                }
                double retirementAmount = 25;
                StockBondAsset asset(returnSpecifier, retirementAmount, 1 - percentStock/100.0);
                double discountRate = StockBondAsset::makeReturnDistribution(returnSpecifier, 1 - percentStock/100.0).getMedian() - 1;
                PortfolioSimulationResult result = performActuarialARVASimulation(asset, discountRate, nextYearSurvivalProbabilities, 1);
                CRRAValues.appendVector(result.percentiles.values * 2);//double amounts for retirement
            }
        }
        if(use25X)
        {
            //calculate CRRA utility
            DEBUG(expectedCRRACertaintyEquivalent(CRRAValues, 1));
            DEBUG(expectedCRRACertaintyEquivalent(CRRAValues, 3));
            DEBUG(expectedCRRACertaintyEquivalent(CRRAValues, 6));
        }
        else
        {
            DEBUG(utility.getMean());
            DEBUG(utilityActuarial.getMean());
        }
    }
}

double invertAllocationToCrraAWithFormula(double financialStocks, double
    weightNonfinancial, ReturnSpecifier const& returnSpecifier)
{
    double stocksSharpeRatio = (returnSpecifier.getStockReturn() -
        returnSpecifier.getRiskFreeRate())/returnSpecifier.getStockStd(),
        KellyStocks = stocksSharpeRatio/returnSpecifier.getStockStd();
    double totalStocks = financialStocks * (1 - weightNonfinancial);
    return KellyStocks/totalStocks;
}

void testCRRAInvertion()
{
    Vector<Vector<string> > matrix;
    Vector<string> constants;
    constants.append("Stock %");
    constants.append("Resulting CRRA a values");
    for(int i = 0; i < 10; ++i) constants.append(" ");
    matrix.append(constants);
    Vector<string> constants2;
    constants2.append("Nonfinancial %");
    for(int i = 0; i <= 10; ++i)
    {
        double adjustedWeightNonfinancial = i == 10 ? 0.99 : i/10.0;
        constants2.append(to_string(adjustedWeightNonfinancial));
    }
    matrix.append(constants2);
    ReturnSpecifier returnSpecifier;
    for(int financialStocks = 0; financialStocks <= 10; ++ financialStocks)
    {
        Vector<string> row;
        double adjustedFinancialStocks =
            financialStocks == 0 ? 0.01 : financialStocks/10.0;
        row.append(to_string(adjustedFinancialStocks));
        for(int weightNonfinancial = 0; weightNonfinancial <= 10;
            ++weightNonfinancial)
        {
            double adjustedWeightNonfinancial =
                weightNonfinancial == 10 ? 0.99 : weightNonfinancial/10.0,
                a = invertAllocationToCrraAWithFormula(adjustedFinancialStocks,
                adjustedWeightNonfinancial, returnSpecifier);
            row.append(to_string(a));
        }
        matrix.append(row);
    }
    createCSV(matrix, "crra_inversions.csv");
}

//returns financial weights of risk-free, bonds, and stocks
Vector<double> getRiskFreeBondStockWeightsFormula(double
    weightNonfinancial, double weightTaxEfficient, double a, ReturnSpecifier
    const& returnSpecifier, bool allowBondsInTaxable = true, double
    minRiskFreeWeight = 0)
{
    assert(weightNonfinancial >= 0 && weightNonfinancial < 1 &&
        weightTaxEfficient >= 0 && weightTaxEfficient <= 1 && a > 0 &&
        minRiskFreeWeight >= 0 && minRiskFreeWeight <= 1);
    double stocksSharpeRatio = (returnSpecifier.getStockReturn() -
        returnSpecifier.getRiskFreeRate())/returnSpecifier.getStockStd(),
        KellyStocks = stocksSharpeRatio/returnSpecifier.getStockStd(),
        CRRAStocks = KellyStocks/a,
        financialStocks = min(1.0, CRRAStocks/(1 - weightNonfinancial)),//cap at 100% stocks
        tangencyStocks = getTangencyStockFraction(returnSpecifier),
        financialRiskFree = 0,
        bonds = 1 - financialStocks;
    if(financialStocks < tangencyStocks)
    {//below tangency prorate
        financialRiskFree = 1 - financialStocks/tangencyStocks;
        financialStocks = tangencyStocks * (1 - financialRiskFree);
        bonds = (1 - tangencyStocks) * (1 - financialRiskFree);
    }
    if(!allowBondsInTaxable)
    {
        bonds = min(bonds, weightTaxEfficient);
        financialRiskFree = 1 - financialStocks - bonds;
    }
    if(financialRiskFree < minRiskFreeWeight)
    {//trim first bonds then stocks
        double gap = minRiskFreeWeight - financialRiskFree;
        financialRiskFree = minRiskFreeWeight;
        bonds -= gap;
        if(bonds < 0)
        {
            financialStocks += bonds;
            bonds = 0;
        }
    }
    Vector<double> fullResult;
    fullResult.append(financialRiskFree);
    fullResult.append(bonds);
    fullResult.append(financialStocks);
    return fullResult;
}

void testOptimizePortfolio()
{
    double weightNonfinancial = 0.0;
    ReturnSpecifier returnSpecifier;
    ReturnSpecifier returnSpecifierTaxable(0.15, 0.25, 0.05);
    double weightTaxEfficient = 0.0;
    double a = 3;

    string name = "";

    Vector<Vector<string> > matrix;
    Vector<string> constants;
    constants.append("Name");
    constants.append("Weight Nonfinancial");
    constants.append("Weight Tax-efficient ");
    constants.append("CRRA a");
    constants.append("Allow Taxable Bonds");
    constants.append("Risk-free %");
    constants.append("Bond %");
    constants.append("Stock %");
    constants.append("CRRA CE");
    constants.append("Formula Risk-free %");
    constants.append("Formula Bond %");
    constants.append("Formula Stock %");
    constants.append("Formula CRRA CE");
    constants.append("Constant Financial Allocation");
    constants.append("CFA CRRA CE");
    matrix.append(constants);

    for(int i = 0; i < 5; ++i)
    {
        if(i == 0)
        {
            name = "Taxable Retired";
            weightNonfinancial = 0.0;
            weightTaxEfficient = 0.0;
        }
        else if(i == 1)
        {
            name = "Tax-Efficient Retired";
            weightNonfinancial = 0.0;
            weightTaxEfficient = 1.0;
        }
        else if(i == 2)
        {
            name = "Late Career";
            weightNonfinancial = 0.1;
            weightTaxEfficient = 0.5;
        }
        else if(i == 3)
        {
            name = "Mid Career";
            weightNonfinancial = 0.5;
            weightTaxEfficient = 0.5;
        }
        else
        {
            name = "Early Career";
            weightNonfinancial = 0.9;
            weightTaxEfficient = 0.5;
        }

        for(int j = 0; j < 3; ++j)
        {
            if(j == 0) a = 1;
            else if(j == 1) a = 3;
            else a = 6;
            {//only relevant for high risk tolerance
                //if(a != 1) continue;
                for(int allowTaxableBonds = 0; allowTaxableBonds < 2; ++allowTaxableBonds)
                {//taxable bonds only relevant for the first case
                    if(allowTaxableBonds == 0 && i != 0) continue;

                    pair<Vector<double>, double> rbs = getRiskFreeBondStockWeights(weightNonfinancial, weightTaxEfficient, a, returnSpecifier, returnSpecifierTaxable, Vector<double>(), allowTaxableBonds);

                    Vector<double> rbsFormula = getRiskFreeBondStockWeightsFormula(weightNonfinancial, weightTaxEfficient, a, returnSpecifier, allowTaxableBonds);
                    double formulaCertaintyEquivalent = getRiskFreeBondStockWeights(weightNonfinancial, weightTaxEfficient, a, returnSpecifier, returnSpecifierTaxable, rbsFormula, allowTaxableBonds).second;

                    Vector<string> row;
                    row.append(name);
                    row.append(to_string(weightNonfinancial));
                    row.append(to_string(weightTaxEfficient));
                    row.append(to_string(a));
                    row.append(to_string(allowTaxableBonds));
                    row.append(to_string(rbs.first[0]));
                    row.append(to_string(rbs.first[1]));
                    row.append(to_string(rbs.first[2]));
                    row.append(to_string(rbs.second));

                    row.append(to_string(rbsFormula[0]));
                    row.append(to_string(rbsFormula[1]));
                    row.append(to_string(rbsFormula[2]));
                    row.append(to_string(formulaCertaintyEquivalent));

                    //benchmarks - constant financial
                    {
                        Vector<double> testX(3, 0);
                        string allocation = "";
                        if(a == 1)
                        {
                            allocation = "0/10/90";
                            testX[1] = 0.1;
                            testX[2] = 0.9;
                        }
                        else if(a == 3)
                        {
                            allocation = "10/30/60";
                            testX[0] = 0.1;
                            testX[1] = 0.3;
                            testX[2] = 0.6;
                        }
                        else if(a == 6)
                        {
                            allocation = "55/15/30";
                            testX[0] = 0.55;
                            testX[1] = 0.15;
                            testX[2] = 0.3;
                        }
                        pair<Vector<double>, double> test = getRiskFreeBondStockWeights(weightNonfinancial, weightTaxEfficient, a, returnSpecifier, returnSpecifierTaxable, testX);
                        row.append(allocation);
                        row.append(to_string(test.second));
                    }
                    matrix.append(row);
                }
            }
        }
    }

    //createCSV(matrix, "optimal_allocations.csv");
}

int main()
{
    testAllAutoFinancialCalculations();
    //return 0;
    testMSCIDitributionMatch();
    testAnnuityAutoSingleFemaleSpending();
    testPrintJointSurvivalProbabilities();
    testMakeReturnFile();
    //return 0;
    testReturnEstimator();
    testMinVarCalculationTotal();
    //return 0;
    //testMinVarCalculationCurrentAverageReturns();
    testMinVarCalculationCurrent();
    testMinVarCalculationCurrentTaxes();
    testMultiAccountReturns();
    //return 0;
    testEstimateStockBondRebalance();
    testGoalSeek();
    testGoalSeekTargetDatePath();
    testTotalBondEffect();
    testOptionInsurance();
    testCRRAInvertion();
    //return 0;
    //testPortfolioSimulationRetirement();
    //testPortfolioSimulationRetirementTaxable();
    //testPortfolioSimulationRetirementDelay();
    //testPortfolioSimulationRetirementActurial();
    //testPortfolioSimulationNoCashFlow();
    //testPortfolioSimulationNoCashFlowUpRebalanced();
    //testPickedPortfolios();
    //testPortfolioSimulationDCA();
    //testPortfolioSimulationDCAVarious();
    //testPortfolioSimulationDCAPaths();
    //testARVARetirement();
    //testFinancialIndependenceDistribution();
    //testMakeReturnChartByAllocation();
    //testMinVarCalculationMakeTable();
    //testOptimizePortfolio();
	return 0;
}
