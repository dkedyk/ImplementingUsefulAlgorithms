#ifndef IGMDK_CONSTRAINT_PROCESSING_H
#define IGMDK_CONSTRAINT_PROCESSING_H
#include "../Utils/Queue.h"
#include "../Graphs/Graph.h"
#include "../Utils/Bitset.h"
#include "../HashTable/LinearProbingHashTable.h"
namespace igmdk{

template<typename CONSTRAINT> struct ConstraintGraph
{
    typedef GraphAA<CONSTRAINT> GRAPH;
    GRAPH g;
    Vector<Bitset<> > variables;
    void addVariable(int domain)
    {
        g.addVertex();
        variables.append(Bitset<>(domain));
    }
    void addConstraint(int v1, int v2, CONSTRAINT const& constraint)
    {
        assert(v1 != v2);
        g.addUndirectedEdge(v1, v2, constraint);
    }
    void disallow(int variable, int value)
        {variables[variable].set(value, false);}
    bool hasSolution(int variable){return !variables[variable].isZero();}

    bool isAllowed(int variable, int value, int otherVariable,
        CONSTRAINT const& constraint)
    {
        for(int i = 0; i < variables[otherVariable].getSize(); ++i)
            if(variables[otherVariable][i] && constraint.isAllowed(variable,
                value, otherVariable, i)) return true;
        return false;
    }
    bool revise(int variable, int otherVariable, CONSTRAINT const& constraint)
    {
        bool changed = false;
        for(int i = 0; i < variables[variable].getSize(); ++i)
            if(variables[variable][i] &&
                !isAllowed(variable, i, otherVariable, constraint))
            {
                disallow(variable, i);
                changed = true;
            }
        return changed;
    }
    bool SAC3Helper(int v, Queue<int>& q, Vector<bool>& onQ, bool isFirstPass)
    {
        onQ[v] = false;
        for(typename GRAPH::AdjacencyIterator i = g.begin(v);
            i != g.end(v); ++i)
        {//in the first pass revise the variables and in subsequent passes
            //the neighbors
            int revisee = i.to(), against = v;
            if(isFirstPass) swap(revisee, against);
            if(revise(revisee, against, i.data()))
            {
                if(!hasSolution(revisee)) return false;//problem unsatisfiable
                if(!onQ[revisee]) q.push(revisee);
                onQ[revisee] = true;
            }
        }
        return true;
    }
    bool SAC3()
    {
        Queue<int> q;
        Vector<bool> onQ(g.nVertices(), true);
        for(int j = 0; j < g.nVertices(); ++j)
            if(!SAC3Helper(j, q, onQ, true)) return false;
        while(!q.isEmpty())
            if(!SAC3Helper(q.pop(), q, onQ, false)) return false;
        return true;
    }
};

struct AllDifferent
{
    LinearProbingHashTable<int, bool> variables;
    void addVariable(int variable){variables.insert(variable, true);}
    struct Handle
    {
        LinearProbingHashTable<int, bool>& variables;
        bool isAllowed(int variable, int value, int variable2, int value2)
            const
        {
            if(variables.find(variable) && variables.find(variable2))
                return value != value2;
            return true;
        }
        Handle(LinearProbingHashTable<int, bool>& theVariables):
            variables(theVariables) {}
    } handle;
    AllDifferent(): handle(variables) {}
};

struct Sudoku
{
    AllDifferent ad[3][9];
    ConstraintGraph<AllDifferent::Handle> cg;
    Sudoku(int values[81])
    {
        for(int i = 0; i < 81; ++i)
        {
            cg.addVariable(9);
            if(values[i])
            {
                cg.variables[i].setAll(false);
                cg.variables[i].set(values[i] - 1, true);
            }
            else cg.variables[i].setAll(true);
        }
        for(int i = 0; i < 9; ++i)
        {
            int rowStart = i * 9, columnStart = i,
                boxStart = i/3 * 27 + (i % 3) * 3;
            for(int j = 0; j < 9; ++j)
            {
                int rowMember = rowStart + j;
                int columnMember = columnStart + j*9;
                int boxMember = boxStart + j/3 * 9 + j % 3;
                ad[0][i].addVariable(rowMember);
                ad[1][i].addVariable(columnMember);
                ad[2][i].addVariable(boxMember);
                if(j == 8) continue;
                for(int k = j+1; k < 9; ++k)
                {
                    int boxMember2 = boxStart + k/3 * 9 + k % 3;
                    cg.addConstraint(rowMember, rowStart + k, ad[0][i].handle);
                    cg.addConstraint(columnMember, columnStart + k * 9,
                        ad[1][i].handle);
                    cg.addConstraint(boxMember, boxMember2, ad[2][i].handle);
                }
            }
        }
        cg.SAC3();
    }
    void printSolution()const
    {
        for(int i = 0; i < 81; ++i)
        {
            if(i % 9 == 0) cout << '\n';
            int count = 0, value = -1;
            for(int j = 0; j < 9; ++j)
            {
                if(cg.variables[i][j])
                {
                    ++count;
                    value = j + 1;
                }
            }
            if(count > 1) cout << "x";
            else cout << value;

        }
        cout << endl;
    }
};

}//end namespace
#endif
