#ifndef IGMDK_TIME_SERIES_H
#define IGMDK_TIME_SERIES_H
#include "Statistics.h"
#include "../Utils/Vector.h"
#include "../Utils/Queue.h"
#include "../NumericalOptimization/NumericalOptimization.h"
#include <cmath>
namespace igmdk{

struct ARMA
{
    enum{HISTORY_LIMIT = 100};//heuristic for computational savings, can
    //also use doubling but more complicated
    Vector<double> phi, theta;
    double std;
    static pair<Vector<double>, Vector<double> > unpackParams(
        Vector<double> const& params, int p)
    {//params = phi|theta|std
        assert(p >= 0 && p <= params.getSize());
        Vector<double> thePhi(p), theTheta(params.getSize() - p);
        for(int i = 0; i < p; ++i) thePhi[i] = params[i];
        for(int i = p; i < params.getSize(); ++i)
            theTheta[i - p] = params[i];
        return make_pair(thePhi, theTheta);
    }
    Vector<double> calculatePsi(int maxSize)const
    {//first weight is 1, rest is given by relation
        Vector<double> psi(max<int>(maxSize,
            max<int>(phi.getSize(), theta.getSize() + 1)), 1);
        for(int i = 0; i < psi.getSize() - 1; ++i)
        {
            double sum = 0;
            for(int j = 0; j < min(i + 1, phi.getSize()); ++j)
                sum += phi[j] * psi[i - j];
            double thetai = i < theta.getSize() ? theta[i] : 0;
            psi[i + 1] = thetai + sum;
        }
        return psi;
    }
    double calculateR0()const
    {//infinite sum of squared psi, converges or diverges
        //1000 is typically more than enough, or can use doubling with
        //sum monitoring
        Vector<double> weights = calculatePsi(1000);
        double sum = 0;//loop order below matters for numerical stability
        for(int i = weights.getSize() - 1; i >= 0; --i)
            sum += weights[i] * weights[i];
        return sum;//for divergent this is huge or not finite
    }
    int desiredHistorySize()const
    {//p for AR models, and at least HISTORY_LIMIT for MA and ARMA models
        return max(phi.getSize() + theta.getSize(),
            theta.getSize() == 0 ? 0 : HISTORY_LIMIT);
    }
    Vector<double> calculateTruncatedInnovations(Vector<double> const& x)const
    {
        Vector<double> innovations = x;
        for(int i = 1; i < x.getSize(); ++i)
        {
            int nPast = min(i, desiredHistorySize());
            Vector<double> xPast(nPast);//nPast values through i - 1
            for(int j = 0; j < nPast; ++j) xPast[j] = x[i - nPast + j];
            innovations[i] -= operator()(xPast);
        }
        return innovations;
    }
    struct LLFunctor
    {
        Vector<double> const& x;
        int phiSize;
        double operator()(Vector<double> const& params)const
        {//std doesn't matter here as long as not 0
            pair<Vector<double>, Vector<double> > result =
                unpackParams(params, phiSize);
            ARMA a(result.first, result.second, 1);
            return a.concentratedLL(x);
        }
    };
    Vector<double> gammaAll(int h)const
    {//calculate using numerical approximation from psi
        Vector<double> gamma(h, 0), psi = calculatePsi(1000);//enough
        for(int i = 0; i < h; ++i)//order below sums smaller numbers first
            for(int j = psi.getSize() - 1 - i; j >= 0; --j)
                gamma[i] += psi[j] * psi[j + i];
        return gamma * (std * std);
    }
    Vector<double> gammaAllEmpirical(int n, Vector<double> const& x)const
    {//first n empirical autocovariances assuming mean 0
        assert(n > 0 && x.getSize() > 0);
        Vector<double> result(n);
        for(int i = 0; i < min(n, x.getSize()); ++i)
        {
            double sum = 0;
            for(int j = 0; j + i < x.getSize(); ++j) sum += x[j + i] * x[j];
            result[i] = sum/x.getSize();
        }
        return result;
    }

    double rho(int h, Vector<double> const& gamma)const
    {
        assert(h >= 0 && h < gamma.getSize());
        return gamma[h]/gamma[0];
    }
    Matrix<double> calculatePACFHelper(int n, Vector<double> const& gamma)const
    {//Durbin-Levinson to solve matrix equation
        assert(n >= 0 && gamma.getSize() == n + 1);
        Matrix<double> PACF(n + 1, n + 1);
        PACF(0, 0) = 0;
        for(int i = 1; i < n + 1; ++i)
        {
            double sum1 = 0, sum2 = 0;
            for(int j = 1; j < i; ++j)
            {
                sum1 += PACF(i - 1, j) * rho(i - j, gamma);
                sum2 += PACF(i - 1, j) * rho(j, gamma);
            }
            PACF(i, i) = (rho(i, gamma) - sum1)/(1 - sum2);
            for(int j = 1; j < i; ++j) PACF(i, j) =
                PACF(i - 1, j) - PACF(i, i) * PACF(i - 1, i - j);
        }
        return PACF;
    }
    Matrix<double> calculatePACF(int n)const
        {return calculatePACFHelper(n, gammaAll(n + 1));}
    Matrix<double> calculatePACFEmpirical(int n, Vector<double> const& x)const
        {return calculatePACFHelper(n, gammaAllEmpirical(n + 1, x));}
public:
    ARMA(): std(1){}//dummy constructor for fitting in future
    ARMA(Vector<double> const& thePhi, Vector<double> theTheta, double theStd):
        phi(thePhi), theta(theTheta), std(theStd){}//direct copy
    ARMA(int p, int q, Vector<double> const& x): phi(p, 0), theta(q, 0), std(1)
    {//fit the model; want more data than parameters
        assert(p >= 0 && q >= 0 && x.getSize() > p + q + 1);
        if(p + q > 0)
        {
            //uninformed initial solution - all 0
            Vector<double> x0(p + q, 0);
            LLFunctor ll = {x, p};
            //try Yule-Walker preliminary estimator for theta
            if(p > 0)
            {
                Vector<double> x1(p + q, 0);
                Matrix<double> pacfe = calculatePACFEmpirical(p, x);
                for(int i = 0; i < p; ++i) x1[i] = pacfe(p, i + 1);
                if(ll(x1) < ll(x0)) x0 = x1;
            }
            //maximum likelihood
            pair<Vector<double>, double> result =
                hybridLocalMinimize(x0, ll, 200);//don't need more
            pair<Vector<double>, Vector<double> > result2 =
                unpackParams(result.first, p);
            phi = result2.first;
            theta = result2.second;
            std = sqrt(LLHelper(x).second/x.getSize());
        }
        else
        {//noise model
            IncrementalStatistics s;
            for(int i = 0; i < x.getSize(); ++i) s.addValue(x[i]);
            std = s.stdev();
        }
    }
    ARMA(Vector<double> const& x, int maxP = 5, int maxQ = 5): std(1)
    {//select using smallest BIC based from 0 <= p, q <= maxOrder
        double bestBIC = numeric_limits<double>::infinity();
        for(int p = 0; p <= maxP; ++p)
            for(int q = 0; q <= maxQ; ++q)
            {//TODO - use prev model coeffs and 0 to init next model coeffs
                ARMA b(p, q, x);
                double bBIC = b.BIC(x);
                if(bBIC < bestBIC)
                {//executes at least once unless NaN data
                    *this = b;
                    bestBIC = bBIC;
                }
            }
    }
    Vector<pair<double, double> > operator()(//predictions and their variances
        Vector<double> const& xPast, int m)const
    {//based on truncated evaluation
        assert(m > 0);
        int p = phi.getSize(), q = theta.getSize();
        //the most recent value at index size - 1 in both
        //initialize with 0 unavailable past values
        Queue<double> xPastRolling(p), wPastRolling(q);
        for(int i = 0; i < p; ++i) xPastRolling.push(0);
        for(int i = 0; i < q; ++i) wPastRolling.push(0);
        //for prediction with last values
        double psiSum = 0;
        Vector<double> psi = calculatePsi(m);
        Vector<pair<double, double> > predictionsAndIntervals(m);
        //process relevant past data, then all future data
        for(int i = q > 0 ? 0 : max(0, xPast.getSize() - 1 - p);
            i < xPast.getSize() + m; ++i)
        {//both buildup and prediction
            double xPredicted = 0;//first evaluate prediction
            for(int j = 0; j < p; ++j)
                xPredicted += phi[p - 1 - j] * xPastRolling[j];
            for(int j = 0; j < q; ++j)
                xPredicted += theta[q - 1 - j] * wPastRolling[j];
            double xi = xPredicted;//prediction case
            if(i < xPast.getSize()) xi = xPast[i];//past data case
            else
            {//prediction case
                int predI = i - xPast.getSize();
                psiSum += psi[predI] * psi[predI];
                predictionsAndIntervals[predI] =
                    make_pair(xPredicted, std * std * psiSum);
            }
            double wPredicted = xi - xPredicted;//0 for prediction case
            //update queues
            if(p > 0)
            {
                xPastRolling.pop();
                xPastRolling.push(xi);
            }
            if(q > 0)
            {
                wPastRolling.pop();
                wPastRolling.push(wPredicted);
            }
        }
        return predictionsAndIntervals;
    }
    //xPast is in the usual order, least recent to most recent
    double operator()(Vector<double> const& xPast)const
        {return operator()(xPast, 1)[0].first;}//single truncated prediction
    Vector<double> rValues(int n)const
    {
        double r = calculateR0();
        Matrix<double> PACF = calculatePACF(min<int>(n, HISTORY_LIMIT));
        Vector<double> result;
        for(int i = 0; i < n; ++i)
        {
            if(i < HISTORY_LIMIT)
            {//check for numerical error - likely process not causal
                if(r < 0) return Vector<double>(n,
                    numeric_limits<double>::infinity());
                if(i > 0) r *= 1 - PACF(i, i) * PACF(i, i);
            }
            else r = 1;//limiting r value for a causal process
            result.append(r);
        }
        return result;
    }
    pair<double, double> LLHelper(Vector<double> const& x)const
    {
        int n = x.getSize();
        double rSum = 0, SSum = 0;
        Vector<double> innovations = calculateTruncatedInnovations(x),
            r = rValues(n);
        if(!isfinite(r[0]) || log10(r[0]) > 3)//heuristic for inf
        {//not a causal process if r0 diverged
            double inf = numeric_limits<double>::infinity();
            return make_pair(inf, inf);
        }
        for(int i = 0; i < n; ++i)
        {
            //r[i] = 1;//uncomment to get unconditional least squares
            rSum += log(r[i]);
            SSum += innovations[i] * innovations[i]/r[i];
        }
        return make_pair(rSum, SSum);
    }
    double concentratedLL(Vector<double> const& x)const
    {//for minimization
        int n = x.getSize();
        pair<double, double> sums = LLHelper(x);
        if(!isfinite(sums.first) || !isfinite(sums.second))
            return numeric_limits<double>::infinity();
        return log(sums.second/n) + sums.first/n;
    }
    double LL(Vector<double> const& x)const
    {//for BIC
        pair<double, double> sums = LLHelper(x);
        //not a causal process if r0 diverged
        if(!isfinite(sums.first) || !isfinite(sums.second))
            return -numeric_limits<double>::infinity();
        int n = x.getSize();
        double temp = 2 * std * std;
        return -(log(temp * PI()) * n/2 + sums.first/2 + sums.second/temp);
    }
    double BIC(Vector<double> const& x, int extraParams = 0)const
    {
        int n = x.getSize(),
            k = phi.getSize() + theta.getSize() + 1 + extraParams;//1 for std
        assert(n >= 0 && extraParams >= 0);
        return k * log(n) - 2 * LL(x);
    }
};

struct ARIMA
{
    ARMA a;
    mutable double mean;
    int d;
    Vector<double> transformData(Vector<double>& x, bool computeMean)const
    {//difference
        Vector<double> xRemoved;
        for(int j = 0; j < d; ++j)
        {
            xRemoved.append(x.lastItem());
            for(int i = 0; i + 1 < x.getSize(); ++i) x[i] = x[i + 1] - x[i];
            x.removeLast();
        }
        if(computeMean)
        {//estimate mean
            IncrementalStatistics s;
            for(int i = 0; i < x.getSize(); ++i) s.addValue(x[i]);
            mean = s.getMean();
        }
        //remove mean from data
        for(int i = 0; i < x.getSize(); ++i) x[i] -= mean;
        return xRemoved;
    }
public:
    ARIMA(int p, int theD, int q, Vector<double> x): d(theD)
    {//want more data than params
        assert(d >= 0 && p >= 0 && q >= 0 && x.getSize() > p + q + d + 3);
        transformData(x, true);
        //fit predictor
        a = ARMA(p, q, x);
    }
    ARIMA(Vector<double> const& x, int maxD = 1, int maxP = 5, int maxQ = 5)
    {//want more data than params
        assert(maxD >= 0 && maxP >= 0 && maxQ >= 0 &&
            x.getSize() > maxP + maxQ + maxD + 2);
        double bestBIC = numeric_limits<double>::infinity(), bestMean = 0;
        int bestD = 0;
        for(d = 0; d <= maxD; ++d)
        {//BIC across d comparable because MDL compression still holds
            Vector<double> xCopy = x;
            transformData(xCopy, true);
            //fit predictor
            ARMA aNew = ARMA(xCopy, maxP, maxQ);
            double bic = aNew.BIC(xCopy, d + 1);
            if(bic < bestBIC)
            {
                a = aNew;
                bestBIC = bic;
                bestMean = mean;
                bestD = d;
            }
        }
        mean = bestMean;
        d = bestD;
    }
    Vector<pair<double, double> > operator()(Vector<double> x, int m)const
    {//predictions and their variances
        assert(x.getSize() >= d + 1 && m > 0);//need at least 1 data point left
        Vector<double> xRemoved = transformData(x, false);
        Vector<pair<double, double> > predictions = a(x, m);
        //add back the mean
        for(int i = 0; i < m; ++i) predictions[i].first += mean;
        //integrate
        for(int i = 0; i < d; ++i)
            for(int j = 0; j < m; ++j)
            {
                predictions[j].first +=
                    j == 0 ? xRemoved[d - 1 - i] : predictions[j - 1].first;
                predictions[j].second += //variances add up
                    j == 0 ? 0 : predictions[j - 1].second;
            }
        return predictions;
    }
    double operator()(Vector<double> x)const{return operator()(x, 1)[0].first;}
    double BIC(Vector<double> x)const
    {
        transformData(x, false);
        return a.BIC(x, d + 1);
    }
};

}//end namespace
#endif
