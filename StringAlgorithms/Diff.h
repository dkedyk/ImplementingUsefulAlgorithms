#ifndef IGMDK_DIFF_H
#define IGMDK_DIFF_H
#include "../Utils/Utils.h"
#include "../Utils/Vector.h"
#include "../Utils/GCFreeList.h"
namespace igmdk{

template<typename CHAR> class Diff
{
public:
    struct EditResult
    {
        CHAR c;
        int position;
        bool isInsert;
    };
private:
    struct Edit
    {
        Edit* prev;//used only for intermediate work and not the final result
        int position;
        bool isInsert;
    };
    static void extendDiagonal(int d, Vector<int>& frontierX, Vector<Edit*>&
        edits, Vector<CHAR> const& a, Vector<CHAR> const& b, Freelist<Edit>& f)
    {//pick next best edit
        int x = max(frontierX[d - 1] + 1, frontierX[d + 1]),
            y = x - (d - 1 - a.getSize());
        if(x != -1 || y != -1)
        {//apply it if not base case
            bool isInsert = x != frontierX[d + 1];
            edits[d] = new(f.allocate())Edit();
            edits[d]->isInsert = isInsert;
            edits[d]->prev = edits[d + (isInsert ? -1 : 1)];
            edits[d]->position = isInsert ? x : y;
        }//move diagonally as much as possible
        while(y + 1 < a.getSize() && x + 1 < b.getSize() &&
            a[y + 1] == b[x + 1])
        {
            ++y;
            ++x;
        }
        frontierX[d] = x;
    }
    static Vector<EditResult> DiffInternal(Vector<CHAR> const& a,
        Vector<CHAR> const& b, CHAR const& nullC)
    {
        int M = a.getSize(), N = b.getSize(), size = M + N + 3,
            mainDiagonal = N + 1;
        assert(M <= N);//a must be shorter then b for this helper
        Vector<int> frontierX(size, -2);
        Vector<Edit*> edits(size, 0);
        Freelist<Edit> f;
        for(int p = 0; frontierX[mainDiagonal] < N - 1; ++p)
        {//from lower left to main
            for(int d = M + 1 - p; d < mainDiagonal; ++d)
                extendDiagonal(d, frontierX, edits, a, b, f);
            //from upper right to main
            for(int d = mainDiagonal + p; d >= mainDiagonal; --d)
                extendDiagonal(d, frontierX, edits, a, b, f);
        }//retrieve the computed path in reverse order
        Vector<EditResult> result;
        for(Edit* link = edits[mainDiagonal]; link; link = link->prev)
        {
            EditResult er = {nullC, link->position, link->isInsert};
            result.append(er);
        }//fix the order
        result.reverse();
        return result;
    }
public:
    static Vector<EditResult> diff(Vector<CHAR> const& a, Vector<CHAR> const&b,
        CHAR const& nullC = CHAR())//null char used for delete action as dummy
    {//edits needed to get a into b - positions are with respect to b
        bool exchange = a.getSize() > b.getSize();
        Vector<EditResult> result = exchange ? DiffInternal(b, a, nullC) :
            DiffInternal(a, b, nullC);
        for(int i = 0, netInserted = 0; i < result.getSize(); ++i)
        {//exchange if needed, set characters, and adjust deletion positions
            if(exchange) result[i].isInsert = !result[i].isInsert;
            if(result[i].isInsert)
            {
                ++netInserted;
                result[i].c = b[result[i].position];
            }
            else result[i].position += netInserted--;
        }
        return result;
    }
    static Vector<CHAR> applyDiff(Vector<CHAR> const& a,
        Vector<EditResult> const& script)
    {
        Vector<CHAR> b;
        int nextA = 0;
        for(int i = 0; i < script.getSize(); ++i)
        {//take chars from a until next edit position
            while(b.getSize() < script[i].position)
            {//basic input check - must not run out of a before next position
                assert(nextA < a.getSize());
                b.append(a[nextA++]);
            }
            if(script[i].isInsert) b.append(script[i].c);
            else ++nextA;//skip one a char on delete
        }//done with script, append the rest from a
        while(nextA < a.getSize()) b.append(a[nextA++]);
        return b;
    }
};

}//end namespace
#endif

