#ifndef ANNUITY_H
#define ANNUITY_H
#include "../Utils/Vector.h"
#include "../NumericalMethods/EquationSolving.h"
#include "../NumericalMethods/Differentiation.h"
#include "../NumericalMethods/Matrix.h"
#include "../RandomNumberGeneration/TimeSeries.h"
#include "../MachineLearning/Lasso.h"
#include "../HashTable/ChainingHashTable.h"
#include "CashFlows.h"
#include <memory>
namespace igmdk{

Vector<double> getSSAMaleDeathProbabilities()
{//from from SSA 2023 TR https://www.ssa.gov/oact/STATS/table4c6.html
    double probabilies[] = {
        0.005837,0.000410,0.000254,0.000207,0.000167,0.000141,0.000123,0.000113,
        0.000108,0.000114,0.000127,0.000146,0.000174,0.000228,0.000312,0.000435,
        0.000604,0.000814,0.001051,0.001250,0.001398,0.001524,0.001612,0.001682,
        0.001747,0.001812,0.001884,0.001974,0.002070,0.002172,0.002275,0.002368,
        0.002441,0.002517,0.002590,0.002673,0.002791,0.002923,0.003054,0.003207,
        0.003333,0.003464,0.003587,0.003735,0.003911,0.004137,0.004452,0.004823,
        0.005214,0.005594,0.005998,0.006500,0.007081,0.007711,0.008394,0.009109,
        0.009881,0.010687,0.011566,0.012497,0.013485,0.014595,0.015702,0.016836,
        0.017908,0.018943,0.020103,0.021345,0.022750,0.024325,0.026137,0.028125,
        0.030438,0.033249,0.036975,0.040633,0.044710,0.049152,0.054265,0.059658,
        0.065568,0.072130,0.079691,0.088578,0.098388,0.109139,0.120765,0.133763,
        0.148370,0.164535,0.182632,0.202773,0.223707,0.245124,0.266933,0.288602,
        0.309781,0.330099,0.349177,0.366635,0.384967,0.404215,0.424426,0.445648,
        0.467930,0.491326,0.515893,0.541687,0.568772,0.597210,0.627071,0.658424,
        0.691346,0.725913,0.762209,0.800319,0.840335,0.882352,0.926469,0.972793
    };
    Vector<double> result;
    for(int i = 0; i < sizeof(probabilies)/sizeof(probabilies[0]); ++i)
        result.append(probabilies[i]);
    return result;
}

Vector<double> getSSAFemaleDeathProbabilities()
{//from SSA 2023 TR https://www.ssa.gov/oact/STATS/table4c6.html
    double probabilies[] = {
        0.004907,0.000316,0.000196,0.000160,0.000129,0.000109,0.000100,0.000096,
        0.000092,0.000089,0.000092,0.000104,0.000123,0.000145,0.000173,0.000210,
        0.000257,0.000314,0.000384,0.000440,0.000485,0.000533,0.000574,0.000617,
        0.000655,0.000700,0.000743,0.000796,0.000851,0.000914,0.000976,0.001041,
        0.001118,0.001186,0.001241,0.001306,0.001386,0.001472,0.001549,0.001637,
        0.001735,0.001850,0.001950,0.002072,0.002217,0.002383,0.002573,0.002777,
        0.002984,0.003210,0.003476,0.003793,0.004136,0.004495,0.004870,0.005261,
        0.005714,0.006227,0.006752,0.007327,0.007926,0.008544,0.009173,0.009841,
        0.010529,0.011265,0.012069,0.012988,0.014032,0.015217,0.016634,0.018294,
        0.020175,0.022321,0.025030,0.027715,0.030631,0.033900,0.037831,0.042249,
        0.047148,0.052545,0.058685,0.065807,0.074052,0.083403,0.093798,0.104958,
        0.117435,0.131540,0.146985,0.163592,0.181562,0.200724,0.219958,0.239460,
        0.258975,0.278225,0.296912,0.314727,0.333610,0.353627,0.374844,0.397335,
        0.421175,0.446446,0.473232,0.501626,0.531724,0.563627,0.597445,0.633292,
        0.671289,0.711567,0.754261,0.799516,0.840335,0.882352,0.926469,0.972793
    };
    Vector<double> result;
    for(int i = 0; i < sizeof(probabilies)/sizeof(probabilies[0]); ++i)
        result.append(probabilies[i]);
    return result;
}

Vector<double> getSSADeathProbabilities()
{//gender - neural - weighted using male/female ratio 1:1
    return (getSSAMaleDeathProbabilities() + getSSAFemaleDeathProbabilities()) *
        0.5;
}

Vector<double> convertToSurvivalProbabilities(
    Vector<double> const& deathByNextYearProbabilities)
{
    double currentSurvivalProportion = 1;
    Vector<double> percentSurvivedByAge;
    for(int i = 0; i < deathByNextYearProbabilities.getSize(); ++i)
    {
        percentSurvivedByAge.append(currentSurvivalProportion);
        currentSurvivalProportion *= 1 - deathByNextYearProbabilities[i];
    }
    return percentSurvivedByAge;
}

class SurvivalEstimator
{
    Vector<double> percentSurvivedByAge;
public:
    double getSurvivalProbability(int ageNow, int ageUntil) const
    {
        assert(ageNow >= 0 && ageNow <= ageUntil);
        //corner cases - ageUntil and possibly ageNow past table knowledge
        if(ageUntil >= percentSurvivedByAge.getSize()) return 0;
        double percentUntil = percentSurvivedByAge[ageUntil];
        double percentNow = percentSurvivedByAge[ageNow];
        return percentUntil/percentNow;//account for current age
    }
    SurvivalEstimator(Vector<double> const& percentSurvivedByAgeTable):
        percentSurvivedByAge(percentSurvivedByAgeTable)
    {//survival is percentage and nonincreasing in age
        double prevPercent = 1;
        for(int i = 0; i < percentSurvivedByAge.getSize(); ++i)
        {
            assert(prevPercent >= percentSurvivedByAge[i] &&
                percentSurvivedByAge[i] >= 0);
            prevPercent = percentSurvivedByAge[i];
        }
    }
};

class JointSurvivalEstimator
{
    Vector<pair<SurvivalEstimator, int> > people;
public:
    double getSurvivalProbability(int ageNow, int ageUntil) const
    {
        double totalDeathProbability = 1;
        for(int i = 0; i < people.getSize(); ++i)
            totalDeathProbability *= 1 - people[i].first.getSurvivalProbability(
                ageNow + people[i].second, ageUntil + people[i].second);
        return 1 - totalDeathProbability;
    }
    void addPerson(SurvivalEstimator const& survivalEstimator, int ageOffset)
        {people.append({survivalEstimator, ageOffset});}
};

template<typename SURVIVAL_ESTIMATOR> double getFutureSurvivalProbability(
    SURVIVAL_ESTIMATOR const& e, int yearsInFuture)
{//assuming setup estimator with 0 starting age
    //at future age to next year
    //in case of joint estimator incremental probabilities
    //are correct only from the starting age
    return e.getSurvivalProbability(0, yearsInFuture + 1)/
       e.getSurvivalProbability(0, yearsInFuture);
}

template<typename SURVIVAL_ESTIMATOR = SurvivalEstimator> class Annuity
{
    int age, nPaymentsPerAgeUnit;
    SURVIVAL_ESTIMATOR survivalEstimator;
    double stepUpR;
    int deferralYears;//assume no fraction years allowed for simplicity
    double initialFee, annualFee;
    int tableDelay;
public:
    Annuity(int theAge, SURVIVAL_ESTIMATOR const& theSurvivalEstimator,
        int theNPaymentsPerAgeUnit = 12, double theStepUpR = 0,
        int theDeferralYears = 0, double theInitialFee = 0,
        double theAnnualFee = 0, int theTableDelay = 0): age(theAge),
        stepUpR(theStepUpR), deferralYears(theDeferralYears),
        nPaymentsPerAgeUnit(theNPaymentsPerAgeUnit),
        survivalEstimator(theSurvivalEstimator), initialFee(theInitialFee),
        annualFee(theAnnualFee), tableDelay(theTableDelay)
    {assert(theAge >= 0 && theNPaymentsPerAgeUnit > 0 && theDeferralYears >= 0
        && 0 <= theStepUpR && 0 <= initialFee && initialFee < 1 &&
        0 <= annualFee && annualFee < 1 && 0 <= tableDelay && tableDelay < 10);}
    Vector<double> calculateEstimatedCashFlow(double payment) const
    {
        assert(payment > 0 && isfinite(payment));
        Vector<double> cashFlow;
        for(int ageNext = age;; ++ageNext)
        {
            double survivalProbability =
                survivalEstimator.getSurvivalProbability(age - tableDelay,
                ageNext - tableDelay);
            if(!isELess(0, survivalProbability)) break;
            double expectedPayment = payment * survivalProbability;
            if(ageNext < age + deferralYears) expectedPayment = 0;
            else if(stepUpR > 0)
                expectedPayment *= pow(1 + stepUpR, ageNext - age);
            //assume same payment without smoothing for survival change
            for(int i = 0; i < nPaymentsPerAgeUnit; ++i)
                cashFlow.append(expectedPayment);
        }
        return cashFlow;
    }
    double calculatePrice(double payment, double r) const
    {
        assert(payment > 0 && isfinite(payment) && isfinite(r));
        r -= annualFee;//implicit r lower for the fee
        Vector<double> cashFlow = calculateEstimatedCashFlow(payment);
        //price higher for the fee
        return (1 + initialFee) * igmdk::calculatePriceGeneral(cashFlow,
            getDatesFrom0(cashFlow.getSize(), 1.0/nPaymentsPerAgeUnit), r);
    }
    double calculatePayment(double price, double r) const
    {
        assert(price > 0 && isfinite(r));
        price /= 1 + initialFee;//implicit price lower for the fee
        r -= annualFee;//implicit r lower for the fee
        //solve for r such that price equals to calculated price
        auto priceFunctor = [price, r, this](double payment)
            {return calculatePrice(payment, r) - price;};
        return exponentialSearch1Sided(priceFunctor, 0).first;
    }
    double calculateYield(double price, double payment) const
    {//no price and r adjustments for yield calculation
        Vector<double> cashFlow = calculateEstimatedCashFlow(payment);
        return igmdk::calculateYieldGeneral(cashFlow,
            getDatesFrom0(cashFlow.getSize(), 1.0/nPaymentsPerAgeUnit), price);
    }
};

}//end namespace
#endif
