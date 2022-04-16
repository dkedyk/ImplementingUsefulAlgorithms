#ifndef IGMDK_SEARCH_ALGORITHMS_H
#define IGMDK_SEARCH_ALGORITHMS_H
#include "../Utils/Vector.h"
#include "../Utils/Stack.h"
#include "../HashTable/ChainingHashTable.h"
#include "../HashTable/LinearProbingHashTable.h"
#include "../Heaps/Heap.h"
#include "../Heaps/IndexedHeap.h"
#include "../Sorting/Sort.h"
namespace igmdk{

template<typename PROBLEM> pair<typename PROBLEM::X, bool> branchAndBound(
    PROBLEM const& p, int maxLowerBounds = 1000000)
{//require one complete solution before stopping
    typename PROBLEM::X best;//assumed to be default-constructable
    double bestScore = numeric_limits<double>::infinity();
    typename PROBLEM::INCREMENTAL_STATE is = p.getInitialState();
    bool foundCompleteSolution = false;//will guarantee complete solution
    struct BBState
    {
        Vector<pair<double, typename PROBLEM::Move> > moves;
        int prevMove;
    };
    Stack<BBState> states;
    bool start = true;
    while(!states.isEmpty() || start)
    {
        bool goNextLevel = start, isSolutionComplete = false;
        if(start) start = false;
        else//act on current level
        {
            BBState& level = states.getTop();
            int& i = level.prevMove;
            //undo prev move if any
            if(i != -1) p.undoMove(is, level.moves[i].second);
            ++i;//check next move if any
            if(i < level.moves.getSize() && level.moves[i].first < bestScore)
            {//if not out of moves and not pruned do next move
                p.move(is, level.moves[i].second);
                if(p.isSolutionComplete(is))
                {//update best
                    double score = p.evaluate(is);
                    if(score < bestScore)
                    {
                        best = p.extractSolution(is);
                        bestScore = score;
                    }//update flags
                    isSolutionComplete = true;
                    foundCompleteSolution = true;
                }
                else goNextLevel = true;
            }
        }
        if(goNextLevel && (!foundCompleteSolution || maxLowerBounds > 0))
        {//setup next level
            BBState levelNext = {p.generateMoves(is), -1};
            int m = levelNext.moves.getSize();
            assert(m > 0);//else bad problem
            maxLowerBounds -= m;
            quickSort(levelNext.moves.getArray(), 0, m - 1,
                PairFirstComparator<double, typename PROBLEM::Move>());
            states.push(levelNext);
        }//for a complete solution first undo current move as generate no more
        else if(!isSolutionComplete) states.pop();//come back to prev level
    }
    return make_pair(best, maxLowerBounds > 0);
}

template<typename PROBLEM> typename PROBLEM::X realtimeAStar(PROBLEM const& p)
{//require one complete solution before stopping
    typename PROBLEM::INCREMENTAL_STATE is = p.getInitialState();
    do
    {
        Vector<pair<double, typename PROBLEM::Move> > moves =
            p.generateMoves(is);
        assert(moves.getSize() > 0);//else bad problem
        int best = argMin(moves.getArray(), moves.getSize(),
            PairFirstComparator<double, typename PROBLEM::Move>());
        p.move(is, moves[best].second);
    }while(!p.isSolutionComplete(is));
    return p.extractSolution(is);
}

template<typename PROBLEM> struct AStar
{//closed set paths stored as parent pointer tree
    typedef typename PROBLEM::STATE_ID STATE_ID;
    typedef typename PROBLEM::HASHER HASHER;
    typedef ChainingHashTable<STATE_ID, STATE_ID, HASHER> P;
    class PredVisitor
    {
        P& pred;
    public:
        PredVisitor(P& thePred): pred(thePred){}
        STATE_ID const* getPred(STATE_ID x)const{return pred.find(x);}
        Vector<STATE_ID> getPath(STATE_ID x)const
        {
            Vector<STATE_ID> path;
            for(;;)
            {
                path.append(x);
                STATE_ID* px = pred.find(x);
                if(px) x = *px;//form path
                else break;//no pred
            }
            path.reverse();//need path in travel not parent pointer order
            return path;
        }
    };
    static pair<Vector<STATE_ID>, bool> solve(PROBLEM const& p,
        int maxSetSize = 1000000)
    {//parent of current state as data
        typedef pair<double, STATE_ID> QNode;
        IndexedHeap<QNode,
            PairFirstComparator<double, STATE_ID>, STATE_ID, HASHER> pQ;
        P pred;
        PredVisitor v(pred);
        bool foundGoal = false;
        STATE_ID j = p.start();//start has no predecessor
        pQ.insert(QNode(p.remainderLowerBound(p.nullState(), j, v),
            p.nullState()), j);
        while(!pQ.isEmpty() && pred.getSize() + pQ.getSize() < maxSetSize)
        {
            pair<QNode, STATE_ID> step = pQ.deleteMin();
            j = step.second;
            if(j != p.start())
                pred.insert(j, step.first.second);//now know best predecessor
            if(p.isGoal(j, v))
            {
                foundGoal = true;
                break;
            }//subtract the last move's lower bound to get the exact distance
            double dj = step.first.first -
                p.remainderLowerBound(step.first.second, j, v);
            Vector<STATE_ID> next = p.nextStates(j, v);
            for(int i = 0; i < next.getSize(); ++i)
            {
                STATE_ID to = next[i];
                double newChildLowerBound = dj + p.distance(j, to, v) +
                    p.remainderLowerBound(j, to, v);
                QNode const* current = pQ.find(to);
                if((current && newChildLowerBound < current->first) ||
                   (!current && !pred.find(to)))//update if better or new
                   pQ.changeKey(QNode(newChildLowerBound, j), to);
            }
        }//form path to goal or best current partial solution
        return make_pair(v.getPath(j), foundGoal);
    }
};

template<typename PROBLEM> struct RecursiveBestFirstSearch
{
    typedef typename PROBLEM::STATE_ID STATE_ID;
    typedef Stack<STATE_ID> P;
    P pred;//path to the goal, which is top
    PROBLEM const& p;
    enum{SUCCESS = -1, FAILURE = -2};
    typedef pair<double, STATE_ID> INFO;//lower bound and state
    bool foundGoal;
    Vector<STATE_ID> bestPath;
    class PredVisitor
    {//assume top of stack always current node, so pred is next
        Stack<STATE_ID>& pred;
    public:
        PredVisitor(P& thePred): pred(thePred){}
        STATE_ID const* getPred(STATE_ID dummy)const
        {
            Vector<STATE_ID>& storage = pred.storage;
            return storage.getSize() > 1 ? &storage[storage.getSize() - 2] : 0;
        }
        Vector<STATE_ID> getPath(STATE_ID dummy)const
        {
            Vector<STATE_ID> path;
            Stack<STATE_ID> s = pred;
            while(!s.isEmpty()) path.append(s.pop());
            path.reverse();//need path in travel not parent pointer order
            return path;
        }
    };
    double work(INFO state, double alternative, double pathCost,
        int& maxLowerBounds)
    {//stop if found goal, of out of moves, or exceed computation budget
        PredVisitor v(pred);
        if(p.isGoal(state.second, v)) return SUCCESS;
        Vector<STATE_ID> next = p.nextStates(state.second, v);
        if(next.getSize() == 0) return numeric_limits<double>::infinity();
        if(maxLowerBounds < next.getSize()) return FAILURE;
        maxLowerBounds -= next.getSize();
        //sort children by lower bound
        Heap<INFO, PairFirstComparator<double, STATE_ID> > children;
        for(int i = 0; i < next.getSize(); ++i)
            children.insert(INFO(max(state.first, pathCost +
                p.distance(state.second, next[i], v) +
                p.remainderLowerBound(state.second, next[i], v)), next[i]));
        for(;;)
        {
            INFO best = children.deleteMin();
            //don't process remaining children if alternative better and
            //return the current best child value
            if(best.first > alternative) return best.first;
            //compute d before push best, else break visitor invariant
            double d = p.distance(state.second, best.second, v);
            pred.push(best.second);
            //as alternative use better of alternative and next best child
            best.first = work(best, children.isEmpty() ?
                alternative : min(children.getMin().first, alternative),
                pathCost + d, maxLowerBounds);
            if(best.first == SUCCESS) return SUCCESS;
            else if(best.first == FAILURE) return FAILURE;
            children.insert(best);//enqueue child with revised estimate
            pred.pop();//undo move
        }
    }
    RecursiveBestFirstSearch(PROBLEM const& theProblem,
        int maxLowerBounds = 10000000): p(theProblem), foundGoal(false)
    {
        pred.push(p.start());
        foundGoal = (work(INFO(0.0, pred.getTop()),
            numeric_limits<double>::infinity(), 0, maxLowerBounds) == SUCCESS);
        PredVisitor v(pred);
        bestPath = v.getPath(pred.getTop());
    }
};

}//end namespace
#endif
