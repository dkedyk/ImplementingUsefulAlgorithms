#include "ReinforcementLearning.h"
#include "../RandomNumberGeneration/Random.h"
using namespace igmdk;

struct GridWorld
{
    DiscreteValueFunction u;
    int state, nEpisodes, nextState;
    double reward()
    {
        if(state == 3) return 1;
        if(state == 7) return -1;
        return -0.04;
    }
    double discountRate(){return 1;}
    double goToNextState(){state = nextState;}
    double pickNextState()
    {
        int row = state % 4, column = state / 4;
        int rows[4] = {row+1,row,row,row-1};
        int columns[4] = {column,column+1,column-1,column};
        bool set = false;
        for(int i = 0; i < 4; ++i)
        {
            if(rows[i] >= 0 && rows[i] <= 3 && columns[i] >= 0 && columns[i] <= 2 && !(rows[i] == 1 && columns[i] == 1))
            {
                int newState = rows[i] + columns[i] * 4;
                assert(state !=newState);
                if(!set || u.values[newState].first > u.values[nextState].first) {nextState = newState; set = true;}
            }
        }
        assert(set);
        return u.values[nextState].first;
    }
    bool isInFinalState(){return state == 3 || state == 7;}
    double learningRate(){return u.learningRate(state);}
    bool hasMoreEpisodes(){return nEpisodes;}
    double startEpisode()
    {
        do{state = GlobalRNG().mod(12);} while(state == 5);
        --nEpisodes;
        return u.values[state].first;}
    void updateCurrentStateValue(double delta){u.updateValue(state, delta);}
    GridWorld():nEpisodes(100), u(12){}
    void debug()
    {
        for(int i = 0; i < 3; ++i)
        {
            for(int j = 0; j < 4; ++j)
            {
                cout << " " << u.values[j + i * 4].first;
            }
            cout << endl;
        }
    }
};

void testReinforcement()
{
    GridWorld g;
    TDLearning(g);
    g.debug();
}

int main(int argc, char *argv[])
{
    testReinforcement();
}


