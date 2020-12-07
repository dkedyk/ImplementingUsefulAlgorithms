#ifndef IGMDK_SATISFIABILITY3_H
#define IGMDK_SATISFIABILITY3_H
#include "../RandomNumberGeneration/Random.h"
#include "SearchAlgorithms.h"
#include "MetaHeuristics.h"
namespace igmdk{

class Satisfiability3RandomInstance
{
    Vector<Vector<int> > clauses;
    int n;
public:
    Satisfiability3RandomInstance(int theN, int nc): n(theN),
        clauses(nc, Vector<int>(3))//3-sat
    {//use n clauses for n variables
        for(int i = 0; i < nc; ++i)
        {
            Vector<int>& clause = clauses[i];
            for(int j = 0; j < clause.getSize(); ++j)
            {//random variables and signs
                int variable = GlobalRNG().mod(n);
                clause[j] = variable * GlobalRNG().sign();
            }
        }
    }
    int getSize()const{return n;}
    int getClauseSize()const{return clauses.getSize();}
    Vector<int> const& getClause(int i)const
    {
        assert(i >= 0 && i < getClauseSize());
        return clauses[i];
    }
};

template<typename PROBLEM> class SatisfiabilityProxy
{
    PROBLEM const& p;
    Vector<Vector<int> > varToClauseMap;
    bool evaluateClause(Vector<bool> const& subset, Vector<int> const& clause,
        int overrideI = -1)const
    {
        bool value = false;
        for(int j = 0; j < clause.getSize(); ++j)
        {
            int i = abs(clause[j]);
            bool valVar = subset[i];
            if(i == overrideI) valVar = !valVar;//override means flip
            if(clause[j] < 0) valVar = !valVar;//negated variable
            value |= valVar;//sat is "and" of "or" clauses
        }
        return value;
    }
public:
    typedef int INCREMENTAL_STATE;//nSatisfied
    SatisfiabilityProxy(PROBLEM const& theP): p(theP),
        varToClauseMap(p.getSize())
    {//use n clauses for n variables
        for(int i = 0; i < p.getClauseSize(); ++i)
        {
            Vector<int> const& clause = p.getClause(i);
            for(int j = 0; j < clause.getSize(); ++j)
            {//variables have signs
                int variable = abs(clause[j]);
                varToClauseMap[variable].append(i);
            }
        }
    }
    INCREMENTAL_STATE updateIncrementalState(int iToFlip,
        Vector<bool> const& subset, INCREMENTAL_STATE const& is)const
    {
        Vector<int> const& affectedClauses = varToClauseMap[iToFlip];
        int nSatisfiedOld = 0, nSatisfied = 0;
        for(int i = 0; i < affectedClauses.getSize(); ++i)
        {
            Vector<int> const& clause = p.getClause(affectedClauses[i]);
            nSatisfiedOld += evaluateClause(subset, clause);
            nSatisfied += evaluateClause(subset, clause, iToFlip);
        }
        return is + (nSatisfied - nSatisfiedOld);
    }
    double evalStep(int iToFlip,//evaluate single step difference
        Vector<bool> const& subset, INCREMENTAL_STATE const& is)const
        {return -(updateIncrementalState(iToFlip, subset, is) - is);}
    INCREMENTAL_STATE getIncrementalState(Vector<bool> const& subset)const
    {
        assert(subset.getSize() == varToClauseMap.getSize());
        int nSatistied = 0;
        for(int i = 0; i < p.getClauseSize(); ++i)
            nSatistied += evaluateClause(subset, p.getClause(i));
        return nSatistied;
    }
    double operator()(Vector<bool> const& subset)const
        {return -getIncrementalState(subset);}
};

}//end namespace
#endif
