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
        // and b = [0, 0, 0.06, 1] (transposed)
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

class MultiAccountReturns
{
    struct AccountData
    {
        LognormalDistribution annualReturns;
        double weight, stockProportion;
        AccountData(pair<double, double> const& meanstd, double theWeight,
            double theStockProportion): stockProportion(theStockProportion),
            annualReturns(1 + meanstd.first, meanstd.second), weight(theWeight)
        {}
    };
    Vector<AccountData> accounts;
public:
    void addAccount(double weight, double bondFraction, double riskFreeFraction,
        ReturnSpecifier const& returnSpecifier)
    {
        assert(weight > 0 && bondFraction >= 0 && bondFraction <= 1 &&
            riskFreeFraction >= 0 && riskFreeFraction <= 1);
        MeanVariancePortfolio stockBondMVP = makeStockBondMVP(returnSpecifier);
        Vector<double> weights;
        weights.append(1 - bondFraction);
        weights.append(bondFraction);
        pair<double, double> meanstd = stockBondMVP.evaluate(weights);
        meanstd.first = (1 - riskFreeFraction) * meanstd.first +
            riskFreeFraction * returnSpecifier.getRiskFreeRate();
        meanstd.second *= 1 - riskFreeFraction;
        accounts.append(AccountData(meanstd, weight,
            (1 - bondFraction) * (1 - riskFreeFraction)));
    }
    double getTotalGeometricReturnRate()const
    {
        double totalWeight = 0, totalGeometricRate = 0;
        for(int i = 0; i < accounts.getSize(); ++i)
        {
            totalWeight += accounts[i].weight;
            totalGeometricRate += accounts[i].weight *
                accounts[i].annualReturns.getMedian();
        }
        return (totalGeometricRate)/totalWeight - 1;
    }
    double getTotalStockProportion()const
    {
        double totalWeight = 0, totalStockProportion = 0;
        for(int i = 0; i < accounts.getSize(); ++i)
        {
            totalWeight += accounts[i].weight;
            totalStockProportion += accounts[i].weight *
                accounts[i].stockProportion;
        }
        return totalStockProportion/totalWeight;
    }
};

}//end namespace
#endif
