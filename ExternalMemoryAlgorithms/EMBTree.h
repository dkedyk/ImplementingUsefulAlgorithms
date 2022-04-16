#ifndef IGMDK_EMBT_H
#define IGMDK_EMBT_H
#include "EMVector.h"
#include "EMFreelist.h"
namespace igmdk{

template<typename KEY, typename VALUE, typename KEY_SERIALIZER =
    CastSerializer<KEY>, typename VALUE_SERIALIZER = CastSerializer<VALUE> >
    class EMBPlusTree
{
    typedef KVPair<KEY, long long> Key;
    typedef KVPair<KEY, VALUE> Record;
    //satisfy the constraints on M and L, but first from internal node make
    //space for size and from leaf also the next pointer
    enum{NULL_IO_POINTER = -1, NODE_SIZE_BYTES = 4, POINTER_SIZE = 8,
        KEY_SIZE = KEY_SERIALIZER::byteSize() + POINTER_SIZE, RECORD_SIZE =
        KEY_SERIALIZER::byteSize() + VALUE_SERIALIZER::byteSize(), B =
        BlockFile::targetBlockSize() - NODE_SIZE_BYTES, M = 2 * min<int>(2,
        B/2/KEY_SIZE), L = 2 * min<int>(1, (B - POINTER_SIZE)/2/RECORD_SIZE)
    };
    struct Node
    {
        int size;
        Key next[M];
        Node(): size(1) {next[0].value = NULL_IO_POINTER;}
        int findChild(KEY const& key)
        {//last child not larger than the key
            int i = 0;
            while(i < size - 1 && key >= next[i].key) ++i;
            return i;
        }
        struct Serializer
        {
            KEY_SERIALIZER ks;
            constexpr static int byteSize()
                {return NODE_SIZE_BYTES + int(M) * int(KEY_SIZE);}
            Node operator()(Vector<unsigned char> const& bytes)
            {
                assert(bytes.getSize() == byteSize());//basic file check
                Node node;
                BitStream bs(bytes);//first decode size
                node.size = ReinterpretDecode(bs.readBytes(NODE_SIZE_BYTES));
                for(int i = 0; i < M; ++i)//then the next pointers
                {
                    node.next[i].key =
                        ks(bs.readBytes(KEY_SERIALIZER::byteSize()));
                    node.next[i].value =
                        ReinterpretDecode(bs.readBytes(POINTER_SIZE));
                }
                return node;
            }
            Vector<unsigned char> operator()(Node const& node)
            {
                BitStream bs;//first encode size
                bs.writeBytes(ReinterpretEncode(node.size, NODE_SIZE_BYTES));
                for(int i = 0; i < M; ++i)//then the next pointers
                {
                    bs.writeBytes(ks(node.next[i].key));
                    bs.writeBytes(ReinterpretEncode(node.next[i].value,
                        POINTER_SIZE));
                }
                return bs.bitset.getStorage();
            }
        };
    };
    struct Leaf
    {
        int size;
        long long next;
        Record records[L];
        Leaf(): size(0), next(NULL_IO_POINTER) {}
        int inclusiveSuccessorRecord(KEY const& key)
        {//last record smaller than the key
            int i = 0;
            while(i < size && key > records[i].key) ++i;
            return i;
        }
        struct Serializer
        {
            KEY_SERIALIZER ks;
            VALUE_SERIALIZER vs;
            constexpr static int byteSize()
            {
                return int(NODE_SIZE_BYTES) + int(POINTER_SIZE) +
                    int(L) * int(RECORD_SIZE);
            }
            Leaf operator()(Vector<unsigned char> const& bytes)
            {
                assert(bytes.getSize() == byteSize());//basic file check
                Leaf leaf;
                BitStream bs(bytes);//first decode size
                leaf.size = ReinterpretDecode(bs.readBytes(NODE_SIZE_BYTES));
                //then next pointer
                leaf.next = ReinterpretDecode(bs.readBytes(POINTER_SIZE));
                for(int i = 0; i < L; ++i)//then the records
                {
                    leaf.records[i].key =
                        ks(bs.readBytes(KEY_SERIALIZER::byteSize()));
                    leaf.records[i].value =
                        vs(bs.readBytes(VALUE_SERIALIZER::byteSize()));
                }
                return leaf;
            }
            Vector<unsigned char> operator()(Leaf const& leaf)
            {
                BitStream bs;//first encode size
                bs.writeBytes(ReinterpretEncode(leaf.size, NODE_SIZE_BYTES));
                //the next pointer
                bs.writeBytes(ReinterpretEncode(leaf.next, POINTER_SIZE));
                for(int i = 0; i < L; ++i)//then the records
                {
                    bs.writeBytes(ks(leaf.records[i].key));
                    bs.writeBytes(vs(leaf.records[i].value));
                }
                return bs.bitset.getStorage();
            }
        };
    };//leaf indices start at -2
    long long leafIndex(long long index){return -(index + 2);}
    long long inverseLeafIndex(long long lIndex){return -lIndex - 2;}
    long long root;
    File header;//for the root pointer, uncached access
    EMVector<Node, typename Node::Serializer> nodes;
    EMFreelist<Leaf, typename Leaf::Serializer> leaves;

    pair<long long, long long> findLeaf(KEY const& key)
    {
        long long current = root, parent = NULL_IO_POINTER;
        while(current >= 0)//stop at leaf or null
        {
            parent = current;
            Node node = nodes[current];
            current = node.next[node.findChild(key)].value;
        }
        return make_pair(current, parent);
    }

    void splitInternal(long long index, int child)
    {
        Node parent = nodes[index];
        long long childIndex = parent.next[child].value;
        Node left = nodes[childIndex], right;
        //copy middle item key into the parent, shifting the latter's other
        //keys
        for(int i = parent.size++; i > child; --i)
            parent.next[i] = parent.next[i - 1];
        parent.next[child].key = left.next[M/2 - 1].key;
        parent.next[child + 1].value = nodes.getSize();
        //move items starting from middle into right
        right.size = M/2 + 1;
        for(int i = 0; i < right.size; ++i)
            right.next[i] = left.next[i + M/2 - 1];
        left.size = M/2;
        nodes.append(right);//write the nodes
        nodes.set(left, childIndex);
        nodes.set(parent, index);
    }
    void splitLeaf(long long index, int child)
    {
        Node parent = nodes[index];
        long long childIndex = parent.next[child].value,
            newChildIndex = inverseLeafIndex(leaves.allocate());
        Leaf left = leaves[leafIndex(childIndex)], right;
        //copy middle item key into the parent internal node, shifting the
        //latter's other keys
        for(int i = parent.size++; i > child; --i)
            parent.next[i] = parent.next[i - 1];
        parent.next[child].key = left.records[L/2].key;
        parent.next[child + 1].value = newChildIndex;
        //move items starting from middle into right
        left.size = right.size = L/2;
        for(int i = 0; i < right.size; ++i)
            right.records[i] = left.records[i + L/2];
        right.next = left.next;
        left.next = newChildIndex;
        leaves.set(right, leafIndex(newChildIndex));//write the nodes
        leaves.set(left, leafIndex(childIndex));
        nodes.set(parent, index);
    }
    EMBPlusTree(EMBPlusTree const&);//no copying allowed
    EMBPlusTree& operator=(EMBPlusTree const&);
public:
    EMBPlusTree(string const& filenameSuffix): root(NULL_IO_POINTER),
        header(("Header" + filenameSuffix).c_str(), false),
        nodes("Keys" + filenameSuffix, 8), leaves("Records" + filenameSuffix)
    {
        Vector<unsigned char> temp(POINTER_SIZE);
        if(header.getSize() > 0)
        {
            header.read(temp.getArray(), POINTER_SIZE);
            root = ReinterpretDecode(temp);
        }
        else header.append(temp.getArray(), POINTER_SIZE);//make space for root
    }
    ~EMBPlusTree()
    {//write root to header
        header.setPosition(0);
        header.write(ReinterpretEncode(root, POINTER_SIZE).getArray(),
            POINTER_SIZE);
    }
    VALUE find(KEY const& key, bool& status)
    {
        status = true;
        long long current = findLeaf(key).first;
        if(current != NULL_IO_POINTER)
        {
            Leaf leaf = leaves[leafIndex(current)];
            int i = leaf.inclusiveSuccessorRecord(key);
            if(i < leaf.size && key == leaf.records[i].key)
                return leaf.records[i].value;
        }
        status = false;
        return VALUE();
    }
    bool shouldSplit(long long node)
    {
        return node < NULL_IO_POINTER ?
            leaves[leafIndex(node)].size == L : nodes[node].size == M;
    }
    void insert(KEY const& key, VALUE const& value)
    {//the first node is root as leaf
        if(root == NULL_IO_POINTER) root = inverseLeafIndex(leaves.allocate());
        else if(shouldSplit(root))
        {//split the root if needed
            Node newRoot;
            newRoot.next[0].value = root;
            bool wasLeaf = root < NULL_IO_POINTER;
            root = nodes.getSize();
            nodes.append(newRoot);
            wasLeaf ? splitLeaf(root, 0) : splitInternal(root, 0);
        }
        long long index = root;
        while(index > NULL_IO_POINTER)//internal node work
        {//go down, inserting and splitting
            Node node = nodes[index];
            int childI = node.findChild(key),child = node.next[childI].value;
            if(shouldSplit(child))
            {//split child if needed
                child < NULL_IO_POINTER ? splitLeaf(index, childI) :
                    splitInternal(index, childI);
                if(key > nodes[index].next[childI].key)//go to next child
                    child = nodes[index].next[childI + 1].value;
            }
            index = child;
        }
        //insert the item into the leaf
        Leaf leaf = leaves[leafIndex(index)];
        int i = leaf.inclusiveSuccessorRecord(key);
        if(i < leaf.size && key == leaf.records[i].key)
            leaf.records[i].value = value;
        else
        {
            for(int j = leaf.size++; j > i; --j)
                leaf.records[j] = leaf.records[j - 1];
            leaf.records[i] = Record(key, value);
        }
        leaves.set(leaf, leafIndex(index));
    }
    void remove(KEY const& key)
    {
        pair<long long, long long> pointerAndParent = findLeaf(key);
        long long pointer = pointerAndParent.first;
        if(pointer != NULL_IO_POINTER)
        {
            Leaf leaf = leaves[leafIndex(pointer)];
            int i = leaf.inclusiveSuccessorRecord(key);
            if(i < leaf.size && key == leaf.records[i].key)
            {
                --leaf.size;
                for(int j = i; j < leaf.size; ++j)
                    leaf.records[j] = leaf.records[j + 1];
            }
            if(leaf.size > 0) leaves.set(leaf, leafIndex(pointer));
            else
            {//remove leaf
                leaves.deallocate(leafIndex(pointer));
                Node parent = nodes[pointerAndParent.second];
                parent.next[parent.findChild(key)].value = NULL_IO_POINTER;
                nodes.set(parent, pointerAndParent.second);
            }
        }
    }
};

}//end namespace
#endif
