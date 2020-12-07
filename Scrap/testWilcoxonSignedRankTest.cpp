#include "WilcoxonSignedRankTest.h"
#include "../Utils/Debug.h"
#include "../RandomNumberGeneration/Statistics.h"
using namespace igmdk;

void testWilcoxon()
{
    Vector<pair<double, double> > a;
    a.append(pair<double, double>(51.2, 45.8));
    a.append(pair<double, double>(46.5, 41.3));
    a.append(pair<double, double>(24.1, 15.8));
    a.append(pair<double, double>(10.2, 11.1));
    a.append(pair<double, double>(65.3, 58.5));
    a.append(pair<double, double>(92.1, 70.3));
    a.append(pair<double, double>(30.3, 31.6));
    a.append(pair<double, double>(49.2, 35.4));
    double z = signedRankZ(a);
    DEBUG(z);
    DEBUG(1 - approxNormal2SidedConf(z));
}

void testWilcoxon2()
{
    Vector<pair<double, double> > a;
    a.append(pair<double, double>(125, 110));
    a.append(pair<double, double>(115, 122));
    a.append(pair<double, double>(130, 125));
    a.append(pair<double, double>(140, 120));
    a.append(pair<double, double>(140, 140));
    a.append(pair<double, double>(115, 124));
    a.append(pair<double, double>(140, 123));
    a.append(pair<double, double>(125, 137));
    a.append(pair<double, double>(140, 135));
    a.append(pair<double, double>(135, 145));
    DEBUG(signedRankZ(a));
}

int main(int argc, char *argv[])
{
    testWilcoxon();
    return 0;
    testWilcoxon2();
    return 0;
}
