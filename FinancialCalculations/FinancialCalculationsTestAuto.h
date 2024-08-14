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

void testOptionPricerAuto()
{
    DEBUG("testOptionPricerAuto");
    double price = 62, strikePrice = 60, r = 0.1, q = 0.2, strikeTime = 5.0/12;
    assert(isEEqual(priceEuropeanCall(price, r, q, strikePrice, strikeTime), 5.7977812415148975));
    double riskFreeRate = 0.03;
    assert(isEEqual(priceEuropeanPut(price, riskFreeRate, r, q, strikePrice, strikeTime), 3.0524492711477933));
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
    Mortgage mortgage(0.03945, 100000, 23);
    assert(isEEqual(mortgage.getMonthlyPayment(), 548.0188297427261));
    DEBUG("testMortgageAuto passed");
}

void testMortgageReverseAuto()
{
    DEBUG("testMortgageReverseAuto");
    Mortgage mortgage(0.035, -100000, 23);
    assert(isEEqual(mortgage.getMonthlyPayment(), -525.11858027067376));
    DEBUG("testMortgageReverseAuto passed");
}

void testAnnuityAuto()
{
    DEBUG("testAnnuityAuto");
    Annuity<> annuity(70, convertToSurvivalProbabilities(getSSAFemaleDeathProbabilities()));
    assert(isEEqual(annuity.calculateYield(100000, 500), -0.0020841008604693336));
    Annuity<> annuity2(70, convertToSurvivalProbabilities(getSSAMaleDeathProbabilities()));
    assert(isEEqual(annuity2.calculateYield(100000, 500), -0.01801336140538479));
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

void testAllAutoFinancialCalculations()
{
    DEBUG("testAllAutoFinancialCalculations");
    testCashFlowAuto();
    testBondAuto();
    testOptionPricerAuto();
    testBondDurationConvexityAuto();
    testMortgageAuto();
    testMortgageReverseAuto();
    testAnnuityAuto();
    testAnnuityAutoJoint();
    testEstimateLogNormalParametersAuto();
    testCRRAAuto();
    DEBUG("testAllAutoFinancialCalculationsPassed");
}

}//end namespace
#endif
