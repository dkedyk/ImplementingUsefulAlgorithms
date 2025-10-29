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

void testCashFlowAuto()
{
    DEBUG("testCashFlowAuto");
    Vector<double> cashFlowSimple;
    cashFlowSimple.append(1000);
    cashFlowSimple.append(1000);
    assert(isEEqual(calculatePrice(cashFlowSimple, 0), 2000));
    assert(isEEqual(calculateYield(cashFlowSimple, 2000), 7.2759576141834261e-15));
    DEBUG("testCashFlowAuto passed");
}

void testBondAuto()
{
    DEBUG("testBondAuto");
    GeneralBond simple(50, 1000, 10, 0, 2);
    assert(isEEqual(simple.calculateYield(1000), 0.057161739605136958));
    DEBUG("testBondAuto passed");
}

void testOrdinaryAnnuityCalculatorAuto()
{
    DEBUG("testOrdinaryAnnuityCalculatorAuto");
    double pmt = 1000, r = 0.04;
    int n = 10;
    {
        assert(isEEqual(pmt/(1 + r), OrdinaryAnnuityCalculator::PV(pmt, 1, r), numeric_limits<float>::epsilon()));
        assert(isEEqual(pmt, OrdinaryAnnuityCalculator::FV(pmt, 1, r), numeric_limits<float>::epsilon()));
        assert(isEEqual(1, OrdinaryAnnuityCalculator::NPER(pmt, pmt/(1 + r), r), numeric_limits<float>::epsilon()));
        double pvCalculated = OrdinaryAnnuityCalculator::PV(pmt, n, r);
        assert(isEEqual(pmt, OrdinaryAnnuityCalculator::PMT(OrdinaryAnnuityCalculator::PV(pmt, n, r), n, r), numeric_limits<float>::epsilon()));
        assert(isEEqual(OrdinaryAnnuityCalculator::FV(pmt, n, r), OrdinaryAnnuityCalculator::FVFromPV(pvCalculated, n, r), numeric_limits<float>::epsilon()));
        assert(isEEqual(n, OrdinaryAnnuityCalculator::NPER(pmt, pvCalculated, r), numeric_limits<float>::epsilon()));
        assert(isEEqual(r, OrdinaryAnnuityCalculator::RATE(pmt, pvCalculated, n), numeric_limits<float>::epsilon()));
    }
    {//payment at beginning
        assert(isEEqual(pmt, OrdinaryAnnuityCalculator::PVBeginning(pmt, 1, r), numeric_limits<float>::epsilon()));
        assert(isEEqual(pmt * (1 + r), OrdinaryAnnuityCalculator::FVBeginning(pmt, 1, r), numeric_limits<float>::epsilon()));
        assert(isEEqual(1, OrdinaryAnnuityCalculator::NPERBeginning(pmt, pmt, r), numeric_limits<float>::epsilon()));
        double pvCalculatedBeginning = OrdinaryAnnuityCalculator::PVBeginning(pmt, n, r);
        assert(isEEqual(pmt, OrdinaryAnnuityCalculator::PMTBeginning(OrdinaryAnnuityCalculator::PVBeginning(pmt, n, r), n, r), numeric_limits<float>::epsilon()));
        assert(isEEqual(OrdinaryAnnuityCalculator::FVBeginning(pmt, n, r), OrdinaryAnnuityCalculator::FVFromPV(pvCalculatedBeginning, n, r), numeric_limits<float>::epsilon()));
        assert(isEEqual(n, OrdinaryAnnuityCalculator::NPERBeginning(pmt, pvCalculatedBeginning, r), numeric_limits<float>::epsilon()));
        assert(isEEqual(r, OrdinaryAnnuityCalculator::RATEBeginning(pmt, pvCalculatedBeginning, n), numeric_limits<float>::epsilon()));
    }
    DEBUG("testOrdinaryAnnuityCalculatorAuto passed");
}

double priceEuropeanPut2(double price, double r, double q, double strikePrice,
    double strikeTime, double t0 = 0)
{//Black-Scholes formula
//r, q, and time have same period unit, usually annual
    assert(price > 0 && isfinite(r) && r > 0 && q > 0 && strikePrice > 0 &&
        strikeTime > 0 && t0 <= strikeTime);
    double t = strikeTime - t0, temp = q * sqrt(t),
        d1 = (log(price/strikePrice) + (r + q * q/2) * t)/temp,
        d2 = d1 - temp, discountedStrike = strikePrice * exp(-r * t);
    DEBUG(log(price/strikePrice));
    DEBUG(r + q * q/2);
    DEBUG(d1);
    DEBUG(d2);
    DEBUG("Stock fraction");
    DEBUG(-approxNormalCDF(-d1));
    DEBUG("Risk-free fraction");
    DEBUG(approxNormalCDF(-d2));
    return -price * approxNormalCDF(-d1) + discountedStrike * approxNormalCDF(-d2);
}

void testOptionPricerAuto()
{
    DEBUG("testOptionPricerAuto");
    double price = 62, strikePrice = 60, r = 0.1, q = 0.2, strikeTime = 5.0/12;
    assert(isEEqual(priceEuropeanCall(price, r, q, strikePrice, strikeTime), 5.7977812415148975));
    assert(isEEqual(priceEuropeanPut(price, r, q, strikePrice, strikeTime), 1.3491486680631866));
    DEBUG(priceEuropeanPut(1, 0.02, 0.17, 0.75, 1));
    assert(isEEqual(priceEuropeanPut(1, 0.02, 0.17, 0.75, 1), 0.0020239428840884699));
    DEBUG(priceEuropeanPut(1, 0.02, 0.17, 0.95, 1));
    assert(isEEqual(priceEuropeanPut(1, 0.02, 0.17, 0.95, 1), 0.036654680954770535));
    assert(isEEqual(priceEuropeanPut(1, 0.02, 0.17, 0.95, 1), priceEuropeanPut2(1, 0.02, 0.17, 0.95, 1)));
    DEBUG("testOptionPricerAuto passed");
}

void testBondDurationConvexityAuto()
{
    DEBUG("testBondDurationConvexityAuto");
    GeneralBond simple(50, 1000, 10);//5% 10 year trading at par
    double currentYield = simple.calculateYield(1000);
    assert(isEEqual(currentYield, 0.057263819488078287));
    assert(isEEqual(simple.calculatePrice(currentYield), 1000.0000000000043));
    double increment = 0.03;
    assert(isEEqual(simple.calculatePrice(currentYield + increment), 824.08749056605666));
    double modifiedDuration = simple.calculateModifiedDuration(1000);
    assert(isEEqual(modifiedDuration, 6.6644536793343008));
    double convexity = simple.calculateConvexity(1000);
    assert(isEEqual(convexity, 58.832023284912111));
    assert(isEEqual(-modifiedDuration * increment + convexity/2 * increment * increment, -0.17345919990181857));
    DEBUG("testBondDurationConvexityAuto passed");
}

void testMortgageAuto()
{
    DEBUG("testMortgageAuto");
    assert(isEEqual(548.0188297427261, calculateMortgageMonthlyPayment(100000, 23, 0.03945), numeric_limits<float>::epsilon()));
    DEBUG("testMortgageAuto passed");
}

void testAnnuityAuto()
{
    DEBUG("testAnnuityAuto");
    Annuity<> annuity(70, convertToSurvivalProbabilities(getSSAFemaleDeathProbabilities()));
    assert(isEEqual(annuity.calculateYield(100000, 500), -0.0020841008604693336));
    assert(isEEqual(annuity.calculatePrice(500, annuity.calculateYield(100000, 500)), 100000, numeric_limits<float>::epsilon()));
    assert(isEEqual(annuity.calculatePayment(100000, annuity.calculateYield(100000, 500)), 500, numeric_limits<float>::epsilon()));
    Annuity<> annuity2(70, convertToSurvivalProbabilities(getSSAMaleDeathProbabilities()));
    assert(isEEqual(annuity2.calculateYield(100000, 500), -0.01801336140538479));
    assert(isEEqual(annuity2.calculatePrice(500, annuity2.calculateYield(100000, 500)), 100000, numeric_limits<float>::epsilon()));
    assert(isEEqual(annuity2.calculatePayment(100000, annuity2.calculateYield(100000, 500)), 500, numeric_limits<float>::epsilon()));
    DEBUG("testAnnuityAuto passed");
}

void testAnnuityAutoJoint()
{
    DEBUG("testAnnuityAutoJoint");
    JointSurvivalEstimator e;
    e.addPerson(convertToSurvivalProbabilities(getSSAFemaleDeathProbabilities()), 70);
    e.addPerson(convertToSurvivalProbabilities(getSSAMaleDeathProbabilities()), 70);
    Annuity<JointSurvivalEstimator> annuity(0, e);
    assert(isEEqual(annuity.calculateYield(100000, 500), 0.015618553427113512));
    assert(isEEqual(annuity.calculatePrice(500, annuity.calculateYield(100000, 500)), 100000, numeric_limits<float>::epsilon()));
    assert(isEEqual(annuity.calculatePayment(100000, annuity.calculateYield(100000, 500)), 500, numeric_limits<float>::epsilon()));
    DEBUG("testAnnuityAutoJoint passed");
}

void testEstimateLogNormalParametersAuto()
{
    DEBUG("testEstimateLogNormalParametersAuto");
    double mean = 1.08, stdev = 0.17;
    pair<double, double> result = estimateLogNormalParameters(mean, stdev);
    double mu = result.first, q = result.second, q2 = q * q;
    assert(isEEqual(mu, 0.064723482321357662));
    assert(isEEqual(q, 0.15644525441681395));
    assert(isEEqual(mean - exp(mu + q2/2), 0));
    assert(isEEqual(stdev * stdev - (exp(q2) - 1) * exp(2 * mu + q2), 0));
    DEBUG("testEstimateLogNormalParametersAuto passed");
}

void testCRRAAuto()
{
    DEBUG("testCRRAAuto");
    Vector<double> values;
    values.append(100);
    values.append(200);
    assert(isEEqual(expectedCRRACertaintyEquivalent(values, 0.001), 149.99150480145843));
    assert(isEEqual(expectedCRRACertaintyEquivalent(values, 1), 141.42135623730945));
    assert(isEEqual(expectedCRRACertaintyEquivalent(values, 1.5), 137.2583002030479));
    assert(isEEqual(expectedCRRACertaintyEquivalent(values, 2), 133.33333333333334));
    assert(isEEqual(expectedCRRACertaintyEquivalent(values, 2.5), 129.72880065637221));
    assert(isEEqual(expectedCRRACertaintyEquivalent(values, 3), 126.49110640673517));
    assert(isEEqual(expectedCRRACertaintyEquivalent(values, 10), 107.9825604906332));
    DEBUG("testCRRAAuto passed");
}

void testMVOptimalRiskyAssetAuto()
{
    DEBUG("testMVOptimalRiskyAssetAuto");
    double value = 1, tangency = 0.7;
    MixedAsset<StockBondAsset> stocks = makeMVOptimalRiskyAsset(value, 1, tangency);
    assert(isEEqual(stocks.getRiskFreeFraction(), 0));
    MixedAsset<StockBondAsset> riskFree = makeMVOptimalRiskyAsset(value, 0, tangency);
    assert(isEEqual(riskFree.getRiskFreeFraction(), 1));
    MixedAsset<StockBondAsset> tangencyAsset = makeMVOptimalRiskyAsset(value, tangency, tangency);
    assert(isEEqual(tangencyAsset.getRiskFreeFraction(), 0));
    MixedAsset<StockBondAsset> halfTangency = makeMVOptimalRiskyAsset(value, tangency/2, tangency);
    assert(isEEqual(halfTangency.getRiskFreeFraction(), 0.5));
    DEBUG("testMVOptimalRiskyAssetAuto passed");
}

void testAllAutoFinancialCalculations()
{
    DEBUG("testAllAutoFinancialCalculations");
    testCashFlowAuto();
    testBondAuto();
    testOrdinaryAnnuityCalculatorAuto();
    testOptionPricerAuto();
    testBondDurationConvexityAuto();
    testMortgageAuto();
    testAnnuityAuto();
    testAnnuityAutoJoint();
    testEstimateLogNormalParametersAuto();
    testCRRAAuto();
    testMVOptimalRiskyAssetAuto();
    DEBUG("testAllAutoFinancialCalculationsPassed");
}

}//end namespace
#endif
