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

void testPortfolioSimulationRetirement()
{
    ReturnSpecifier returnSpecifier;//no tax
    //for(int expensesNowPercent = 2; expensesNowPercent <= 4; ++expensesNowPercent)
    for(int expensesNowPercent = 3; expensesNowPercent <= 5; ++expensesNowPercent)
    {
        DEBUG(expensesNowPercent);
        int term = 30;
        //int term = 50;
        DEBUG(term);
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

void testPortfolioSimulationRetirementTaxable5()
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
        double fraction = findMaximumExpenseRatio(StockBondAsset(returnSpecifier, 100, 0), returnSpecifier.getInflationRate(), term);
        double taxAdjustedFraction = fraction/(1 + taxableFraction*(returnSpecifier.taxRateCapitalGains + returnSpecifier.taxRateLocal));
        DEBUG(taxAdjustedFraction);

        Vector<string> row;
        row.append(to_string(term));
        row.append(to_string(taxAdjustedFraction));
        matrix.append(row);
    }
    //createCSV(matrix, "early_retirement_amounts.csv");

}

double ordinaryAnnuityPV(double payment, int n, double r)
{
    assert(payment > 0 && n > 0 && r > 0);
    return payment * (1 - pow(1 + r, -n))/r;
}
double ordinaryAnnuityPayment(double PV, int n, double r)
{
    assert(PV > 0 && n > 0 && r > 0);
    return PV/((1 - pow(1 + r, -n))/r);
}

void testPortfolioSimulationRetirementDelay()
{
    //ReturnSpecifier returnSpecifier(0.2, 0.3);//tax
    ReturnSpecifier returnSpecifier;//no tax
    for(int expensesNowPercent = 2; expensesNowPercent <= 6; ++expensesNowPercent)
    {
        DEBUG(expensesNowPercent);
        int term = 50;
        DEBUG(term);
        int delay = 20;
        int initialBudget = 100;
        double TIPSRate = returnSpecifier.getRiskFreeRate() - returnSpecifier.getInflationRate();
        initialBudget -= ordinaryAnnuityPV(expensesNowPercent, delay, TIPSRate);
        DEBUG(ordinaryAnnuityPV(expensesNowPercent, delay, TIPSRate));
        double yearsToDepleteWithTIPS = log(expensesNowPercent/100.0/(expensesNowPercent/100.0 - TIPSRate))/TIPSRate;
        DEBUG(yearsToDepleteWithTIPS);

        JointSurvivalEstimator e;
        e.addPerson(convertToSurvivalProbabilities(getSSAMaleDeathProbabilities()), 70);
        e.addPerson(convertToSurvivalProbabilities(getSSAFemaleDeathProbabilities()), 70);
        Annuity<JointSurvivalEstimator> annuity5(0, e, 1, 0, 0, 0.1, 0.01, 2);
        DEBUG("commercial present value TIPS annuity");
        DEBUG(annuity5.calculatePrice(expensesNowPercent, returnSpecifier.getRiskFreeRate() - returnSpecifier.getInflationRate()));

        assert(initialBudget > 0);
        double inflationFactor = pow(1 + returnSpecifier.getInflationRate(), term);
        double riskFreeReturn = 0;
        Vector<double> retirementExpenses = generateRetirementExpenses(
            expensesNowPercent, returnSpecifier.getInflationRate(), term, delay);
        {
            PortfolioSimulationResult safe = performActuarialSimulation(RiskFreeAsset(returnSpecifier, initialBudget), retirementExpenses, e);
            //PortfolioSimulationResult safe = performSimulation(RiskFreeAsset(returnSpecifier, initialBudget), retirementExpenses);
            riskFreeReturn = safe.getMedian();
            DEBUG("risk-free");
            safe.debug(inflationFactor);
        }
        for(int percentStock = 50; percentStock <= 100; percentStock += 10)
        {
            DEBUG(percentStock);
            {
                PortfolioSimulationResult result = performActuarialSimulation(StockBondAsset(returnSpecifier, initialBudget, 1 - percentStock/100.0), retirementExpenses, e);
                //PortfolioSimulationResult result = performSimulation(StockBondAsset(returnSpecifier, initialBudget, 1 - percentStock/100.0), retirementExpenses);
                result.debug(inflationFactor);
                //DEBUG(result.sharpeRatio(riskFreeReturn));
                double riskFreeRank = result.riskFreeRank(riskFreeReturn);
                DEBUG(riskFreeRank);
            }
        }
    }

}

void testPortfolioSimulationRetirementActurial()
{
    ReturnSpecifier returnSpecifier;//no tax
    JointSurvivalEstimator e;
    e.addPerson(convertToSurvivalProbabilities(getSSAMaleDeathProbabilities()), 70);
    e.addPerson(convertToSurvivalProbabilities(getSSAFemaleDeathProbabilities()), 70);

    for(int expensesNowPercent = 2; expensesNowPercent <= 6; ++expensesNowPercent)
    {
        DEBUG(expensesNowPercent);
        int term = 50;
        DEBUG(term);
        double inflationFactor = pow(1 + returnSpecifier.getInflationRate(), term);
        double riskFreeReturn = 0;
        Vector<double> retirementExpenses = generateRetirementExpenses(expensesNowPercent, returnSpecifier.getInflationRate(), term);

        {
            PortfolioSimulationResult safe = performActuarialSimulation(RiskFreeAsset(returnSpecifier, 100), retirementExpenses, e);
            riskFreeReturn = safe.getMedian();
            DEBUG("risk-free");
            safe.debug(inflationFactor);
        }
        for(int percentStock = 50; percentStock <= 100; percentStock += 10)
        {
            DEBUG(percentStock);
            {
                PortfolioSimulationResult result = performActuarialSimulation(StockBondAsset(returnSpecifier, 100, 1 - percentStock/100.0), retirementExpenses, e);
                result.debug(inflationFactor);
                //DEBUG(result.sharpeRatio(riskFreeReturn));
                double riskFreeRank = result.riskFreeRank(riskFreeReturn);
                DEBUG(riskFreeRank);
            }
        }
    }

}

void testPortfolioSimulationNoCashFlow()
{
    ReturnSpecifier returnSpecifier;//no tax
        int term = 30;
        DEBUG(term);
        double inflationFactor = pow(1 + returnSpecifier.getInflationRate(), term);
        double riskFreeReturn = 0;
        Vector<double> noSavings(term, 0);
        {
            PortfolioSimulationResult safe = performSimulation(RiskFreeAsset(returnSpecifier, 100), noSavings, 1000000);
            riskFreeReturn = safe.getMedian();
            DEBUG("risk-free");
            safe.debug(inflationFactor);
        }
        for(int percentStock = 50; percentStock <= 100; percentStock += 10)
        {
            DEBUG(percentStock);
            {
                PortfolioSimulationResult result = performSimulation(StockBondAsset(returnSpecifier, 100, 1 - percentStock/100.0), noSavings, 1000000);
                result.debug(inflationFactor);
                double riskFreeRank = result.riskFreeRank(riskFreeReturn);
                DEBUG(riskFreeRank);
            }
        }

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
        double inflationFactor = pow(1 + returnSpecifier.getInflationRate(), term);
        double riskFreeReturn = 0;
        Vector<double> noSavings(term, 0);
        {
            PortfolioSimulationResult safe = performSimulation(RiskFreeAsset(returnSpecifier, 100), noSavings, 10000000);
            riskFreeReturn = safe.getMedian();
            DEBUG("risk-free");
            safe.debug(inflationFactor);
        }

        //calculate tangency
        MeanVariancePortfolio mvp(makeStockBondMVP(returnSpecifier));
        pair<double, Vector<double> > tangency = mvp.findOptimalSharpeWeights(returnSpecifier.getRiskFreeRate());
        pair<double, double> mqTangency = mvp.evaluate(tangency.second);
        double tangencyStockFraction = tangency.second[0];


        for(int percentStock = 5; percentStock <= 100; percentStock += 5)
        {
            double stockFraction = percentStock/100.0;
            DEBUG(percentStock);
            Vector<PortfolioSimulationResult> results;
            {//optimal
                double bondFraction = 1 - stockFraction, riskFreeFraction = 0;
                if(stockFraction < tangencyStockFraction)
                {
                    bondFraction = 1 - tangencyStockFraction;
                    riskFreeFraction = 1 - stockFraction/tangencyStockFraction;
                }
                PortfolioSimulationResult result = performSimulation(makeGeneralAsset(100, 1 - stockFraction, 0, 0, returnSpecifier), noSavings, 10000000);
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

            double price = 100, strikePrice = 100 - percentStock, strikeTime = 1;
            double putPrice = numeric_limits<double>::quiet_NaN();
            if(stockFraction < 1) putPrice = priceEuropeanPut(price, returnSpecifier.getRiskFreeRate(), returnSpecifier.getStockReturn(), returnSpecifier.getStockStd(), strikePrice, strikeTime);

            Vector<string> row;
            row.append(to_string(stockFraction));
            for(int i = 0; i < 3; ++i)
            {
                PortfolioSimulationResult result = results[i];
                double riskFreeRank = result.riskFreeRank(riskFreeReturn);
                row.append(to_string(result.getMedian()/inflationFactor));
                row.append(to_string(result.percentiles.getPercentile(0.001)/inflationFactor));
                row.append(to_string(result.percentiles.getPercentile(0.05)/inflationFactor));
                row.append(to_string(result.percentiles.getPercentile(0.95)/inflationFactor));
                row.append(to_string(expectedCRRACertaintyEquivalent(result.percentiles.values, 3)/inflationFactor));
                row.append(to_string(expectedCRRACertaintyEquivalent(result.percentiles.values, 6)/inflationFactor));
                row.append(to_string(riskFreeRank));
                row.append(to_string(result.maxDrawdownPercentiles.getPercentile(0.05)));
                row.append(to_string(result.maxDrawdownPercentiles.getPercentile(0.5)));
                row.append(to_string(result.maxDrawdownPercentiles.getPercentile(0.95)));
            }
            row.append(to_string(putPrice/100));
            matrix.append(row);

        }
    //createCSV(matrix, "protection_losses.csv");
}

void testPortfolioSimulationDCA()
{
    int term = 30;
    PortfolioSimulationResult* temp = 0, *tempR = 0;
    {
        DEBUG("taxable DCA");
        ReturnSpecifier returnSpecifier(0.15, 0.20);//tax
        double inflationFactor = pow(1 + returnSpecifier.getInflationRate(), term);
        double initialValue = 0;
        double riskFreeReturn = 0;
        Vector<double> dcaSavings = generateDCASavings(100, returnSpecifier.getInflationRate(), term);
        {
            PortfolioSimulationResult safe = performSimulation(RiskFreeAsset(returnSpecifier, initialValue), dcaSavings);
            riskFreeReturn = safe.getMedian();
            DEBUG("risk-free");
            safe.debug(inflationFactor);
            tempR = new PortfolioSimulationResult(safe);
        }

        for(int percentStock = 100; percentStock <= 100; percentStock += 10)
        {
            DEBUG(percentStock);
            {//note that first year contribution doesn't grow but amount saved in first year is reduced by inflation
                PortfolioSimulationResult result = performSimulation(StockBondAsset(returnSpecifier, initialValue, 1 - percentStock/100.0), dcaSavings);
                result.debug(inflationFactor);
                //DEBUG(result.sharpeRatio(riskFreeReturn));
                double riskFreeRank = result.riskFreeRank(riskFreeReturn);
                DEBUG(riskFreeRank);
                temp = new PortfolioSimulationResult(result);
            }
        }
    }
    {
        DEBUG(term);
        DEBUG("401k DCA");
        ReturnSpecifier returnSpecifier;//no tax
        double inflationFactor = pow(1 + returnSpecifier.getInflationRate(), term);
        double initialValue = 0;
        double riskFreeReturn = 0;
        Vector<double> dcaSavings = generateDCASavings(100, returnSpecifier.getInflationRate(), term);
        {
            PortfolioSimulationResult safe = performSimulation(RiskFreeAsset(returnSpecifier, initialValue), dcaSavings);
            riskFreeReturn = safe.getMedian();
            DEBUG("risk-free");
            safe.debug(inflationFactor);
            DEBUG("joined");
            safe.join(*tempR);
            safe.debug(inflationFactor);
        }

        for(int percentStock = 50; percentStock <= 100; percentStock += 10)
        {
            DEBUG(percentStock);
            {//note that first year contribution doesn't grow but amount saved in first year is reduced by inflation
                PortfolioSimulationResult result = performSimulation(StockBondAsset(returnSpecifier, initialValue, 1 - percentStock/100.0), dcaSavings);
                result.debug(inflationFactor);
                //DEBUG(result.sharpeRatio(riskFreeReturn));
                double riskFreeRank = result.riskFreeRank(riskFreeReturn);
                DEBUG(riskFreeRank);
                DEBUG("joined");
                result.join(*temp);
                result.debug(inflationFactor);

            }
        }
    }
    delete temp;
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
    //Mortgage mortgage(rate, currentAmount, yearToUse);
    //return mortgage.getMonthlyPayment() * 12;
    return ordinaryAnnuityPayment(currentAmount, yearToUse, rate);
}
template<typename FINANCIAL_ASSET>
PortfolioSimulationResult performActuarialARVASimulation(FINANCIAL_ASSET const&
    initialFinancialAsset,
    double discountRate, double inflationRate,
    Vector<double> const& nextYearSurvivalProbabilities,
    int nSimulations = 1000000)
{
    assert(nextYearSurvivalProbabilities.getSize() > 0 && nSimulations > 0);
    double initialSpending = ARVASpending(initialFinancialAsset.getValue(),
        nextYearSurvivalProbabilities.getSize(), discountRate - inflationRate);
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
                values.append(0);
                continue;
            }
            if(nextYearSurvivalProbabilities[step] < GlobalRNG().uniform01())
            {//died, done spending
                break;
            }
            //DEBUG(financialAsset.getValue());
            //DEBUG(nextYearSurvivalProbabilities.getSize()-step);
            double spending = ARVASpending(financialAsset.getValue(),
                nextYearSurvivalProbabilities.getSize()-step,
                discountRate - inflationRate);
            //DEBUG(spending);
            //values.append(spending/initialSpending);
            values.append(spending);
            financialAsset.setValue(financialAsset.getValue()/
                (1 + inflationRate));
            financialAsset.simulateStep(-spending);
        }
    }
    return PortfolioSimulationResult(values, Vector<double>(), Vector(1, 0.0));
}

void testARVARetirement()
{
    ReturnSpecifier returnSpecifier;//no tax
    {
        int term = 50;
        DEBUG(term);
        double riskFreeReturn = 0;
        {
            PortfolioSimulationResult safe = performActuarialARVASimulation(RiskFreeAsset(returnSpecifier, 100), returnSpecifier.getRiskFreeRate(), returnSpecifier.getInflationRate(), Vector<double>(term, 1));
            DEBUG("risk-free");
            safe.debug();
            riskFreeReturn = safe.getMedian();
        }
        for(int percentStock = 50; percentStock <= 100; percentStock += 10)
        {
            DEBUG(percentStock);
            {
                StockBondAsset asset(returnSpecifier, 100, 1 - percentStock/100.0);
                double discountRate = StockBondAsset::makeReturnDistribution(returnSpecifier, 1 - percentStock/100.0).getMedian() - 1;
                PortfolioSimulationResult result = performActuarialARVASimulation(asset, discountRate, returnSpecifier.getInflationRate(), Vector<double>(term, 1));
                result.debug();
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

    Vector<pair<double, Vector<double> > > frontier =
        mvp.findOptimalPortfolioWeightsRange(21);
    DEBUG(returnSpecifier.getRiskFreeRate());
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
    }
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
    createCSV(matrix, "mean_var_total_bond_effect.csv");
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
    constants.append("CRRA(3)");
    constants.append("CRRA(6)");
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

        double inflationFactor = pow(1 + returnSpecifier.getInflationRate(), term);
        double riskFreeReturn = 0;
        Vector<double> noSavings(term, 0);
        {
            DEBUG("risk-free");
            PortfolioSimulationResult safe = performSimulation(RiskFreeAsset(returnSpecifier, 100), noSavings, 10000000);
            riskFreeReturn = safe.getMedian();
        }
        for(int i = 0; i < 5; ++i)
        {
            string name = "", allocation = "";
            double bondFraction = 0, riskFreeFractionUpRebalanced = 0;
            if(i == 0)
            {
                name = "High";
                allocation = "100/0/0";
            }
            else if(i == 1)
            {
                name = "Medium-high";
                allocation = "75/25/0";
                bondFraction = 0.25;
            }
            else if(i == 2)
            {
                name = "Medium";
                allocation = "60/15/25";
                riskFreeFractionUpRebalanced = 0.25;
                bondFraction = 0.15;
            }
            else if(i == 3)
            {
                name = "Medium-low";
                allocation = "50/25/25";
                riskFreeFractionUpRebalanced = 0.25;
                bondFraction = 0.25;
            }
            else if(i == 4)
            {
                name = "Low";
                allocation = "25/10/65";
                riskFreeFractionUpRebalanced = 0.65;
                bondFraction = 0.1;
            }
            DEBUG(name);
            PortfolioSimulationResult result = performSimulation(makeGeneralAsset(100, bondFraction, 0, riskFreeFractionUpRebalanced, returnSpecifier), noSavings, 10000000);
            double riskFreeRank = result.riskFreeRank(riskFreeReturn);

            Vector<string> row;
            row.append(name);
            row.append(allocation);
            row.append(to_string(result.getMedian()/inflationFactor));
            row.append(to_string(result.percentiles.getPercentile(0.05)/inflationFactor));
            row.append(to_string(result.percentiles.getPercentile(0.95)/inflationFactor));
            row.append(to_string(expectedCRRACertaintyEquivalent(result.percentiles.values, 3)/inflationFactor));
            row.append(to_string(expectedCRRACertaintyEquivalent(result.percentiles.values, 6)/inflationFactor));
            row.append(to_string(riskFreeRank));
            row.append(to_string(result.maxDrawdownPercentiles.getPercentile(0.95)));
            matrix.append(row);
        }
    }
    //createCSV(matrix, "picked_portfolios_compare.csv");
}


void testOptionInsurance()
{
    ReturnSpecifier returnSpecifier;//no tax
    double price = 100, strikePrice = 80, strikeTime = 1;
    double geometricReturn = returnSpecifier.getStockReturn() - returnSpecifier.getStockStd() * returnSpecifier.getStockStd()/2;
    for(int i = 0; i < 6; ++i)
    {
        if(i == 0) strikePrice = 80;
        else if(i == 1) strikePrice = 75;
        else if(i == 2) strikePrice = 50;
        else if(i == 3) strikePrice = 25;
        else if(i == 4) strikePrice = 10;
        else if(i == 5) strikePrice = 100;
        DEBUG(strikePrice);
        double putPrice = priceEuropeanPut(price, returnSpecifier.getRiskFreeRate(), returnSpecifier.getStockReturn(), returnSpecifier.getStockStd(), strikePrice, strikeTime);
        DEBUG(putPrice/100);
        DEBUG(geometricReturn - putPrice/100);
        double riskFreeEquivalentFraction = strikePrice/100;//almost
        DEBUG(geometricReturn * (1 - riskFreeEquivalentFraction) + returnSpecifier.getRiskFreeRate() * riskFreeEquivalentFraction);
    }
}


int main()
{
    testAllAutoFinancialCalculations();
    //return 0;
    testMSCIDitributionMatch();
    testMultigoal();
    testAnnuityAutoSingleFemaleSpending();
    testPrintJointSurvivalProbabilities();
    testMakeReturnFile();
    //return;
    testReturnEstimator();
    testMinVarCalculationTotal();
    //return 0;
    testMinVarCalculationCurrentAverageReturns();
    testMinVarCalculationCurrent();
    testMinVarCalculationCurrentTaxes();
    testMultiAccountReturns();
    //return 0;
    testEstimateStockBondRebalance();
    testGoalSeek();
    testGoalSeekTargetDatePath();
    testTotalBondEffect();
    testOptionInsurance();
    //return 0;
    //testPortfolioSimulationRetirement();
    //testPortfolioSimulationRetirementTaxable5();
    //testPortfolioSimulationRetirementDelay();
    //testPortfolioSimulationNoCashFlow();
    //testPortfolioSimulationNoCashFlowUpRebalanced();
    //testPickedPortfolios();
    //testPortfolioSimulationDCA();
    //testPortfolioSimulationDCAPaths();
    //testPortfolioSimulationRetirementActurial();
    //testARVARetirement();
	return 0;
}
