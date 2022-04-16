#ifndef IGMDK_HUFFMAN_TREE_H
#define IGMDK_HUFFMAN_TREE_H
#include "../Heaps/Heap.h"
#include "Stream.h"
#include "StaticCodes.h"
#include "../Utils/GCFreeList.h"
#include <cstdlib>
namespace igmdk{

struct HuffmanTree
{
    enum{W = numeric_limits<unsigned char>::digits, N = 1 << W};
    struct Node
    {
        unsigned char letter;
        int count;
        Node *left, *right;
        Node(int theCount, Node* theLeft, Node* theRight,
             unsigned char theLetter): left(theLeft), right(theRight),
             count(theCount), letter(theLetter) {}
        bool operator<(Node const& rhs)const{return count < rhs.count;}
        void traverse(Bitset<unsigned char>* codebook,
            Bitset<unsigned char>& currentCode)
        {
            if(left)//internal node
            {
                currentCode.append(false);//went left
                left->traverse(codebook, currentCode);
                currentCode.removeLast();
                currentCode.append(true);//went right
                right->traverse(codebook, currentCode);
                currentCode.removeLast();
            }
            else codebook[letter] = currentCode;//leaf
        }
        void append(Bitset<unsigned char>& result)
        {
            result.append(!left);//0 for nonleaf, 1 for leaf
            if(left)
            {
                left->append(result);
                right->append(result);
            }
            else result.appendValue(letter, W);
        }
    }* root;
    Freelist<Node> f;

    HuffmanTree(Vector<unsigned char> const& byteArray)
    {//calculate frequencies
        int counts[N];
        for(int i = 0; i < N; ++i) counts[i] = 0;
        for(int i = 0; i < byteArray.getSize(); ++i) ++counts[byteArray[i]];
        //create leaf nodes
        Heap<Node*, PointerComparator<Node> > queue;
        for(int i = 0; i < N; ++i) if(counts[i] > 0) queue.insert(
            new(f.allocate())Node(counts[i], 0, 0, i));
        //merge leaf nodes to create the tree
        while(queue.getSize() > 1)//until forest merged
        {
            Node *first = queue.deleteMin(), *second = queue.deleteMin();
            queue.insert(new(f.allocate())
                Node(first->count + second->count, first, second, 0));
        }
        root = queue.getMin();
    }

    Node* readHuffmanTree(BitStream& text)
    {
        Node *left = 0, *right = 0;
        unsigned char letter;
        if(text.readBit()) letter = text.readValue(W);//got to a leaf
        else
        {//process internal nodes recursively
            left = readHuffmanTree(text);
            right = readHuffmanTree(text);
        }
        return new(f.allocate())Node(0, left, right, letter);
    }
    HuffmanTree(BitStream& text){root = readHuffmanTree(text);}

    void writeTree(Bitset<unsigned char>& result){root->append(result);}
    void populateCodebook(Bitset<unsigned char>* codebook)
    {
        Bitset<unsigned char> temp;
        root->traverse(codebook, temp);
    }

    Vector<unsigned char> decode(BitStream& text)
    {//wrong bits will give wrong result, but not a crash
        Vector<unsigned char> result;
        for(Node* current = root;;
            current = text.readBit() ? current->right : current->left)
        {
            if(!current->left)
            {
                result.append(current->letter);
                current = root;
            }
            if(!text.bitsLeft()) break;
        }
        return result;
    }
};

Vector<unsigned char> HuffmanCompress(Vector<unsigned char> const& byteArray)
{
    HuffmanTree tree(byteArray);
    Bitset<unsigned char> codebook[HuffmanTree::N], result;
    tree.populateCodebook(codebook);
    tree.writeTree(result);
    for(int i = 0; i < byteArray.getSize(); ++i)
        result.appendBitset(codebook[byteArray[i]]);
    return ExtraBitsCompress(result);
}

Vector<unsigned char> HuffmanUncompress(Vector<unsigned char> const& byteArray)
{
    BitStream text(ExtraBitsUncompress(byteArray));
    HuffmanTree tree(text);
    return tree.decode(text);
}

}//end namespace
#endif
