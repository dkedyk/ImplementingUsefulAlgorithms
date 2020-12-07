#include "EMLinearHashing.h"
#include "../Utils/Debug.h"
#include <cmath>
using namespace std;

using namespace igmdk;

template<typename VERTEX_DATA, typename EDGE_DATA>
struct EMGraph
{
    struct Vertex
    {
        long long upto;
        VERTEX_DATA data;
    };
    EMVector<Vertex> vertices;
    struct Edge
    {
        long long to;
        EDGE_DATA data;
    };
    EMVector<Edge> edges;
    void addVertex(VERTEX_DATA const& data)
    {
        Vertex entry = vertices.lastItem();
        entry.data = data;
        vertices.append(entry);
    }
    void addEdgeToLastVertex(long long to, EDGE_DATA data)
    {
        Edge entry = {to, data};
        edges.append(entry);
        Vertex lastVertex = vertices.lastItem();
        ++lastVertex.upto;
        vertices[vertices.getLast()] = lastVertex;
    }
    long long edgeStart(long long v)
    {
        assert(v >= 0 && v < vertices.getSize());
        return v ? vertices[v-1].upto : 0;
    }
};

int main()
{
    EMLinearHashTable<int, int> trie;
    //EMLinearProbingHashTable<int, int> trie;
    int N = 150000;
    for(int i = 0; i < N; ++i)
    {
        trie.insert(i, i);
    }
    DEBUG("Done inserting");
    for(int j = 0; j < 1; ++j)
    {
        for(int i = 0; i < N; ++i)
        {
            bool status;
            int item = trie.find(i, status);
            assert(status);
            //trie.remove(i);
            //if(status) DEBUG(item);
        }
    }
    /*DEBUG((EMLinearHashTable<int, int>::N));
    DEBUG(trie.ioFindCount);
    DEBUG(trie.listCount);
    DEBUG(trie.maxLength);
    DEBUG(trie.bitSize);
    DEBUG(trie.table.getSize());*/
    /*for(int i = 0; i < 30; ++i)
    {
        DEBUG(i);
        DEBUG(trie.counts[i]);
    }*/



	return 0;
}
