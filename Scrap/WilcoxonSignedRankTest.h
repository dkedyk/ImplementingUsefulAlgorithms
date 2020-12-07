#ifndef IGMDK_WILCOXON_SIGNED_RANK_TEST_H
#define IGMDK_WILCOXON_SIGNED_RANK_TEST_H
#include "../Sorting/Sort.h"
namespace igmdk{

struct SignedRankComparator
{
    typedef pair<double, double> P;
    int sign(P const& p)const{return p.first - p.second > 0 ? 1 : -1;}
    double diff(P const& p)const{return abs(p.first - p.second);}
    bool operator()(P const& lhs, P const& rhs)const
        {return diff(lhs) < diff(rhs);}
    bool isEqual(P const& lhs, P const& rhs)const
        {return diff(lhs) == diff(rhs);}
};
double signedRankZ(Vector<pair<double, double> > a)
{//same test, use Conover version!
    SignedRankComparator c;
    quickSort(a.getArray(), 0, a.getSize() - 1, c);
    int nP = a.getSize(), i = 0;
    //if odd number of 0's, drop first, distribute rest evenly
    while(i < a.getSize() && c.diff(a[i]) == 0) ++i;
    if(i % 2) --nP;
    double signedRankSum = 0, rank2Sum = 0;
    for(i = i % 2; i < a.getSize(); ++i)
    {//rank lookahead to scan for ties, then sum computation
        int j = i;
        while(i + 1 < a.getSize() && c.isEqual(a[i], a[i + 1])) ++i;
        double rank = (i + j)/2.0 + 1 + nP - a.getSize();
        while(j <= i)
        {
            signedRankSum += c.sign(a[j++]) * rank;
            rank2Sum += rank * rank;
        }
    }
    return rank2Sum == 0 ? 0 : abs(signedRankSum)/sqrt(rank2Sum);
}

}
#endif
