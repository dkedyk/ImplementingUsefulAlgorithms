#ifndef OPTION_PRICING_H
#define OPTION_PRICING_H
#include "../RandomNumberGeneration/Statistics.h"
namespace igmdk{

double priceEuropeanCall(double price, double r, double q,
    double strikePrice, double strikeTime, double t0 = 0)
{//Black-Scholes formula
//r, q, and time have same period unit, usually annual
    assert(price > 0 && isfinite(r) && r > 0 && q > 0 && strikePrice > 0 &&
        strikeTime > 0 && t0 <= strikeTime);
    double t = strikeTime - t0, temp = q * sqrt(t),
        d1 = (log(price/strikePrice) + (r + q * q/2) * t)/temp,
        d2 = d1 - temp, discountedStrike = strikePrice * exp(-r * t);
    return price * approxNormalCDF(d1) - discountedStrike * approxNormalCDF(d2);
}

double priceEuropeanPut(double price, double riskFreeRate, double r, double q,
    double strikePrice, double strikeTime, double t0 = 0)
{//price call, then use put-call parity with continuous compounding discounting
    assert(price > 0 && isfinite(r) && r > 0 && q > 0 && strikePrice > 0 &&
        strikeTime > 0 && t0 <= strikeTime);
    double callPrice = priceEuropeanCall(price, r, q, strikePrice, strikeTime,
        t0), discountedStrike = strikePrice * exp(-riskFreeRate * strikeTime);
    return callPrice + discountedStrike - price;
}

}//end namespace
#endif
