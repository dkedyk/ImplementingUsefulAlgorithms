#ifndef MEAN_VARIANCE_OPTIMIZATION_H
#define MEAN_VARIANCE_OPTIMIZATION_H
#include "../Utils/Vector.h"
#include "../NumericalMethods/EquationSolving.h"
#include "../NumericalMethods/Differentiation.h"
#include "../NumericalMethods/Matrix.h"
#include "../RandomNumberGeneration/TimeSeries.h"
#include "../MachineLearning/Lasso.h"
#include "../HashTable/ChainingHashTable.h"
#include <memory>
#include "Misc.h"
#include "../NumericalOptimization/GoldenSection.h"
namespace igmdk{

struct MeanVariancePortfolio
{
    Vector<double> means;
    Matrix<double> covariances;
    MeanVariancePortfolio(Vector<double> const& theMeans,
        Matrix<double> const& theCovariances): means(theMeans),
        covariances(theCovariances)
    {//input checks
        int n = means.getSize();//first do basic checks
        assert(n > 0 && n == covariances.getRows() &&
            n == covariances.getRows());
        for(int row = 0; row < n; ++row)
            for(int column = 0; column < n; ++column)
            {
                if(row == column) assert(covariances(row, column) > 0);
                assert(isEEqual(covariances(row, column),
                     covariances(column, row)));
            }
    }
    pair<double, double> evaluate(Vector<double> const& weights)const
    {//return mean and std
        double mean = dotProduct(means, weights),
            variance = dotProduct(covariances * weights, weights);
        return {mean, sqrt(variance)};
    }
    Vector<double> findOptimalPortfolioWeights(double targetMean)const
    {//Example - with stock, bond 0.07, 0.2; 0.035, 0.1; cov 0.05;
        //targetMean 0.06 get A =
        //0.20^2  0.05^2 -0.07   -1
        //0.05^2  0.10^2 -0.035  -1
        //0.07    0.035   0     0
        //   1     1    0     0
        //and b = [0, 0, 0.06, 1] (transposed)
        int n = means.getSize(), m = n + 2;
        Matrix<double> A(m, m);
        for(int column = 0; column < n; ++column)
        {
            for(int row = 0; row < n; ++row)
                A(row, column) = covariances(row, column);
            A(m - 2, column) = means[column];
            A(m - 1, column) = 1;
        }
        for(int row = 0; row < n; ++row)
        {
            A(row, m - 2) = -means[row];
            A(row, m - 1) = -1;
        }
        Vector<double> b(m, 0);
        b[m - 2] = targetMean;
        b[m - 1] = 1;
        //use LUP decomposition to solve
        LUP lup(A);
        Vector<double> wlm = lup.solve(b);
        wlm.removeLast();//remove l
        wlm.removeLast();//remove m
        return wlm;
    }
    pair<double, double> findMeanRange()const
    {
        double maxMean = means[0], minMean = means[0];
        for(int i = 1; i < means.getSize(); ++i)
        {
            maxMean = max(maxMean, means[i]);
            minMean = min(minMean, means[i]);
        }
        return {minMean, maxMean};
    }
    Vector<pair<double, Vector<double> > > findOptimalPortfolioWeightsRange(
        int nFrontierPoints) const
    {//find portfolio optimal frontier in min to max mean range, with
        //nFrontierPoints
        assert(nFrontierPoints >= 2);
        pair<double, double> minMax = findMeanRange();
        double minMean = minMax.first, maxMean = minMax.second,
            deltaMean = (maxMean - minMean)/(nFrontierPoints - 1);
        //for each mean
        Vector<pair<double, Vector<double> > > result;
        for(int i = 0; i < nFrontierPoints; ++i)
        {
            double meanI = min(minMean + i * deltaMean, maxMean);
            assert(minMean <= meanI && meanI <= maxMean);
            result.append({meanI, findOptimalPortfolioWeights(meanI)});
        }
        return result;
    }
    template <typename FUNCTOR>
    pair<double, Vector<double> > findOptimalFrontierPoint(FUNCTOR const& f)
        const
    {//assume f takes portfolio mean and std and produces a double result,
        //with smaller value better
        pair<double, double> minMax = findMeanRange();
        //find optimal return
        auto functor = [this, &f](double targetMean)
        {
            Vector<double> w = findOptimalPortfolioWeights(targetMean);
            pair<double, double> mq = evaluate(w);
            return f(mq.first, mq.second);
        };
        pair<double, double> best =
            minimizeGS(functor, minMax.first, minMax.second);
        return {best.second, findOptimalPortfolioWeights(best.first)};
    }
    pair<double, Vector<double> > findOptimalSharpeWeights(double riskFreeRate)
        const
    {
        return findOptimalFrontierPoint([riskFreeRate](double mean,
            double stdev){return -(mean - riskFreeRate)/stdev;});
    }
};

MeanVariancePortfolio makeStockBondMVP(ReturnSpecifier const& returnSpecifier)
{
    Vector<double> means;
    means.append(returnSpecifier.getStockReturn());
    means.append(returnSpecifier.getBondReturn());
    //remember to square standard deviations
    Matrix<double> covariances(2, 2);
    covariances(0, 0) = pow(returnSpecifier.getStockStd(), 2);
    covariances(1, 1) = pow(returnSpecifier.getBondStd(), 2);
    covariances(0, 1) = covariances(1, 0) =
        returnSpecifier.getStockBondCorrelation() * sqrt(covariances(0, 0) *
        covariances(1, 1));
    return MeanVariancePortfolio(means, covariances);
}

double getTangencyStockFraction(ReturnSpecifier const& returnSpecifier)
{
    MeanVariancePortfolio mvp(makeStockBondMVP(returnSpecifier));
    pair<double, Vector<double> > tangency =
        mvp.findOptimalSharpeWeights(returnSpecifier.getRiskFreeRate());
    return tangency.second[0];
}

//returns financial weights of risk-free, bonds, and stocks, and the CRRA
//certainty equivalent
pair<Vector<double>, double> getRiskFreeBondStockWeights(double
    weightNonfinancial, double weightTaxEfficient, double a, ReturnSpecifier
    const& returnSpecifier, ReturnSpecifier const& returnSpecifierTaxable,
    Vector<double> const& evaluateRBS, bool allowBondsInTaxable = true, double
    minRiskFreeWeight = 0)
{
    assert(weightNonfinancial >= 0 && weightNonfinancial < 1 &&
        weightTaxEfficient >= 0 && weightTaxEfficient <= 1 && a > 0 &&
        (evaluateRBS.getSize() == 0 || evaluateRBS.getSize() == 3) &&
        minRiskFreeWeight >= 0 && minRiskFreeWeight <= 1);
    double stockTaxPenalty = returnSpecifier.getStockReturn() -
        returnSpecifierTaxable.getStockReturn(),
        bondTaxPenalty = returnSpecifier.getBondReturn() -
        returnSpecifierTaxable.getBondReturn(),
        riskFreeTaxPenalty = returnSpecifier.getRiskFreeRate() -
        returnSpecifierTaxable.getRiskFreeRate(),
        bondTaxDelta = bondTaxPenalty - stockTaxPenalty,
        riskFreeTaxDelta = riskFreeTaxPenalty - stockTaxPenalty;
    auto minObjective = [weightNonfinancial, weightTaxEfficient, a,
        stockTaxPenalty, bondTaxDelta, riskFreeTaxDelta, &returnSpecifier,
        allowBondsInTaxable, minRiskFreeWeight]
        (Vector<double> const& x)
    {
        double weightRiskFree = x[0], weightBond = x[1],
            weightStock = 1 - weightNonfinancial - weightRiskFree - weightBond,
            taxExtraPenalty = riskFreeTaxDelta *
            max(0.0, weightRiskFree + weightBond - weightTaxEfficient) +
            bondTaxDelta *
            max(0.0, weightRiskFree + weightBond - weightTaxEfficient);
        if(weightRiskFree + weightBond > weightTaxEfficient)
        {
            double bondTaxable = max(0.0, weightBond - weightTaxEfficient),
            riskFreeTaxable = weightRiskFree -
                max(0.0, weightTaxEfficient - weightBond);
            taxExtraPenalty = bondTaxable * bondTaxDelta +
                riskFreeTaxable * riskFreeTaxDelta;
        }
        double taxPenalty = stockTaxPenalty * (1 - weightTaxEfficient) +
            taxExtraPenalty;
        double expectedReturn = (weightNonfinancial + weightRiskFree) *
            returnSpecifier.getRiskFreeRate() +
            weightBond * returnSpecifier.getBondReturn() +
            weightStock * returnSpecifier.getStockReturn() -
            taxPenalty * (1 - weightNonfinancial) -
            returnSpecifier.getInflationRate();
        double var = (pow(weightBond * returnSpecifier.getBondStd(), 2) +
            weightBond * returnSpecifier.getBondStd() *
            returnSpecifier.getStockBondCorrelation() * weightStock *
            returnSpecifier.getStockStd() +
            pow(weightStock * returnSpecifier.getStockStd(), 2));
        double crra = expectedReturn - a/2 * var,
            penalty = 1000 * (max(0.0, -weightRiskFree) + max(0.0, -weightBond)
            + max(0.0, -weightStock) +
            (weightRiskFree < minRiskFreeWeight * (1 - weightNonfinancial)));
        if(!allowBondsInTaxable)
            penalty += 1000 * (weightBond > weightTaxEfficient);
        return -crra + penalty;
    };
    pair<Vector<double>, double> result;
    if(evaluateRBS.getSize() == 3)
    {//only evaluated specified risk-free-stock-bond combination
        Vector<double> x0 = evaluateRBS;
        x0.removeLast();
        x0 *= 1 - weightNonfinancial;
        result = {x0, minObjective(x0)};
    }
    else
    {//optimize
        Vector<double> x0(2, 0);// start with 0 bond and correct risk-free
        x0[0] = minRiskFreeWeight * (1 - weightNonfinancial);
        result = hybridLocalMinimize(x0, minObjective);
    }
    Vector<double> fullResult = result.first;
    fullResult.append(1 - weightNonfinancial - fullResult[0] - fullResult[1]);
    fullResult *= 1/(1 - weightNonfinancial);
    return {fullResult, -result.second};
}

}//end namespace
#endif
