#ifndef IGMDK_REG_EX_H
#define IGMDK_REG_EX_H
#include "../Utils/Utils.h"
#include "../Utils/Vector.h"
#include "../Utils/Stack.h"
#include "../Utils/Bits.h"
#include "../Graphs/Graph.h"
namespace igmdk{

class RegularExpressionMatcher
{
    string re;
    int m;
    GraphAA<bool> g;//dummy edge data
    Vector<int> findActiveStates(Vector<int> const& sources)
    {
        Vector<bool> visited(g.nVertices(), false);
        DefaultDFSAction a;
        for(int i = 0; i < sources.getSize(); ++i) if(!visited[sources[i]])
            {
                visited[sources[i]] = true;
                DFSComponent(g, sources[i], visited, a);
            }
        Vector<int> activeStates;
        for(int i = 0; i < visited.getSize(); ++i)
            if(visited[i]) activeStates.append(i);
        return activeStates;
    }
public:
    RegularExpressionMatcher(string const& theRe): re(theRe), m(re.length()),
        g(m + 1)
    {
        Stack<int> clauses;
        for(int i = 0; i < m; ++i)
        {
            int clauseStart = i;
            if(re[i] == '(' || re[i] == '|') clauses.push(i);
            else if(re[i] == ')')
            {
                int clauseOp = clauses.pop();
                if(re[clauseOp] == '|')
                {
                    clauseStart = clauses.pop();
                    g.addEdge(clauseStart, clauseOp + 1);
                    g.addEdge(clauseOp, i);
                }
                else clauseStart = clauseOp;
            }
            if(i < m - 1 && re[i + 1]=='*')//to next start from clause start
                g.addUndirectedEdge(clauseStart, i + 1);
            if(re[i] == '(' || re[i] == '*' || re[i] == ')')
                g.addEdge(i, i + 1);//to next state from current
        }
    }

    bool matches(string const& text)
    {
        Vector<int> activeStates = findActiveStates(Vector<int>(1, 0));
        for(int i = 0; i < text.length() && activeStates.getSize() > 0; ++i)
        {//must be in >= 1 active state to keep going
            Vector<int> stillActive;
            for(int j = 0; j < activeStates.getSize(); ++j)
                if(activeStates[j] < m && re[activeStates[j]] == text[i])
                    stillActive.append(activeStates[j] + 1);
            activeStates = findActiveStates(stillActive);
        }
        for(int j = 0; j < activeStates.getSize(); ++j)
            if(activeStates[j] == m) return true;
        return false;
    }
};

class ShiftAndExtended
{//Joker handling omitted for simplicity
    enum{ALPHABET_SIZE = 1 << numeric_limits<unsigned char>::digits};
    unsigned char *pattern;
    int patternSize, position;//patternSize must be before masks
    unsigned long long charPos[ALPHABET_SIZE], O, P, L, R, state;
    unsigned long long makeMask(Vector<int> const& positions)const
    {
        unsigned long long mask = 0;
        for(int i = 0; i < positions.getSize(); ++i)
        {
            assert(positions[i] >= 0 && positions[i] < patternSize);
            Bits::set(mask, positions[i], true);
        }
        return mask;
    }
public:
    ShiftAndExtended(unsigned char* thePattern, int thePatternSize,
        Vector<int> const& repeatedPositions = Vector<int>(),
        Vector<int> const& optionalPositions = Vector<int>()): position(0),
        state(0), patternSize(thePatternSize), pattern(thePattern),
        R(makeMask(repeatedPositions)), O(makeMask(optionalPositions))
    {//first precompute character bit strings
        assert(patternSize <= numeric_limits<unsigned long long>::digits &&
            !Bits::get(O, 0));//position 0 can't be optional
        for(int i = 0; i < ALPHABET_SIZE; ++i) charPos[i] = 0;
        for(int i = 0; i < patternSize; ++i)
            Bits::set(charPos[pattern[i]], i, true);
        //then masks for optional characters
        unsigned long long sides = O ^ (O >> 1);
        P = (O >> 1) & sides;
        L = O & sides;
    }
    int findNext(unsigned char* text, int textSize)
    {
        while(position < textSize)
        {//first regular and repeatable update
            state = (((state << 1) | 1) | (state & R)) &
                charPos[text[position++]];
            //then optional character update
            unsigned long long sL = state | L;
            state |= O & (sL ^ ~(sL - P));
            if(Bits::get(state, patternSize - 1))return position - patternSize;
        }
        return -1;
    }
};

}//end namespace
#endif

