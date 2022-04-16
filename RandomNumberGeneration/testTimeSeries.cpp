#include "TimeSeries.h"
#include "Statistics.h"
#include <cmath>
using namespace igmdk;

double rho(int h, double f, double t)
{
    return (1 + t * f)* (t + f)/(1 + 2 * t * f + t * t) * pow(f, h - 1);
}

int main(int argc, char *argv[])
{
    double std = 1;
    double f = 0.9, t = 0.5;
    Vector<double> phi(1, f), theta(1, t);
    ARMA a(phi, theta, std);
    DEBUG("psi");
    a.calculatePsi(100).debug();
    DEBUG("expect first 3 psi to be");
    DEBUG(1);
    DEBUG(1.4);
    DEBUG(1.4 * 0.9);
    DEBUG("gamma all");
    Vector<double> gammas = a.gammaAll(10);
    gammas.debug();
    double y0 = std * std * (1 + 2 * t * f + t * t)/(1 - f * f);
    double y1 = std * std * (1 + t * f)* (t + f)/(1 - f * f);
    double y2 = y1 * f;
    DEBUG(y0);
    DEBUG(y1);
    DEBUG(y2);
    DEBUG(sqrt(gammas[0]/a.calculateR0()));
    DEBUG("rho");
    for(int i = 0; i < 10; ++i) DEBUG(a.rho(i, gammas));
    DEBUG(1);
    for(int i = 1; i < 10; ++i) DEBUG(rho(i, f, t));
    DEBUG("pacf");
    a.calculatePACF(10).debug();

    Vector<double> x;
    double prevX = 0, prevNormal = 0;

    IncrementalStatistics s;
    int n = 1000, burnIn = 100;
    for(int i = 0; i < burnIn + n; ++i)
    {
        double variate = GlobalRNG().normal(0, std);
        //DEBUG(variate);
        prevX = prevX * f + prevNormal * t + variate;
        prevNormal = variate;
        if(i >= burnIn)
        {
            x.append(prevX);
            s.addValue(prevX);
        }
    }
    DEBUG("pacf empirical");
    a.calculatePACFEmpirical(10, x).debug();

    /*DEBUG("innovations");
    a.calculateInnovations(x).debug();
    DEBUG("innovationsTruncated");
    a.calculateTruncatedInnovations(x).debug();*/
    DEBUG("R for AR1");
    DEBUG(1/(1 - f * f));
    DEBUG(1);
    DEBUG("rest is 1");
    DEBUG("LL");
    DEBUG(a.LL(x));
    DEBUG(sqrt(a.LLHelper(x).second/x.getSize()));

    DEBUG("stats");
    DEBUG(s.getMean());
    DEBUG(s.stdev());

    DEBUG("rValues - must converge to std fast");
    a.rValues(10).debug();

    DEBUG("fitting ARMA from simulation");
    ARMA b(1, 1, x);
    //ARMA b(x);
    b.phi.debug();
    b.theta.debug();
    DEBUG(b.std);

    //test prediction on all data 10 steps
    int m = 1000;
    double l2 = 0;
    int confMissCount = 0;
    for(int i = 0; i < m; ++i)
    {
        double variate = GlobalRNG().normal(0, std);
        //DEBUG(variate);
        double xNew = prevX * f + prevNormal * t + variate;
        pair<double, double> result = b(x, 1)[0];
        double prediction = result.first;
        l2 += (xNew - prediction) * (xNew - prediction);
        double plusMinus = 2 * sqrt(result.second);
        if(i == 0) DEBUG(sqrt(result.second));
        if(xNew < prediction - plusMinus || xNew > prediction + plusMinus)
            ++confMissCount;
    }
    DEBUG("1 step error");
    DEBUG(sqrt(l2/m));
    DEBUG(1.0 * confMissCount/m);

    //test prediction on all data 1 step
    l2 = 0;
    confMissCount = 0;

    for(int i = 0; i < m; ++i)
    {
        double prevNormal2 = prevNormal;
        int k = 10;
        for(int j = 0; j < k; ++j)
        {
            double variate = GlobalRNG().normal(0, std);
            //DEBUG(variate);
            prevX = prevX * f + prevNormal2 * t + variate;
            prevNormal2 = variate;
        }


        double variate = GlobalRNG().normal(0, std);
        //DEBUG(variate);
        double xNew = prevX * f + prevNormal2 * t + variate;
        pair<double, double> result = b(x, k + 1)[k];
        double prediction = result.first;
        l2 += (xNew - prediction) * (xNew - prediction);
        double plusMinus = 2 * sqrt(result.second);
        if(i == 0) DEBUG(sqrt(result.second));
        if(xNew < prediction - plusMinus || xNew > prediction + plusMinus)
            ++confMissCount;
    }
    DEBUG("k + 1 step error");
    DEBUG(sqrt(l2/m));
    DEBUG(1.0 * confMissCount/m);



    DEBUG("fitting ARIMA from simulation - takes about a minute");
    ARIMA c(x);
    DEBUG(c.d);
    DEBUG(c.mean);
    c.a.phi.debug();
    c.a.theta.debug();
    DEBUG(c.a.std);
    DEBUG(c.BIC(x));

    return 0;
}
