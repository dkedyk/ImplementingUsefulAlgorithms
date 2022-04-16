#ifndef IGMDK_VECTOR_H
#define IGMDK_VECTOR_H
#include <cmath>
#include "Utils.h"
namespace igmdk{

template<typename ITEM> class Vector: public ArithmeticType<Vector<ITEM> >
{
    enum{MIN_CAPACITY = 8};
    int capacity, size;
    ITEM* items;
    void resize()
    {//allocate
        ITEM* oldItems = items;
        capacity = max(2 * size, int(MIN_CAPACITY));
        items = rawMemory<ITEM>(capacity);
        for(int i = 0; i < size; ++i) new(&items[i])ITEM(oldItems[i]);//copy
        rawDestruct(oldItems, size);//deallocate
    }
public:
    ITEM* getArray(){return items;}
    ITEM* const getArray()const{return items;}
    int getSize()const{return size;}
    ITEM& operator[](int i)
    {
        assert(i >= 0 && i < size);
        return items[i];
    }
    ITEM const& operator[](int i)const
    {
        assert(i >= 0 && i < size);
        return items[i];
    }

    explicit Vector(): capacity(MIN_CAPACITY), size(0),
        items(rawMemory<ITEM>(capacity)) {}//no default ITEM constructor needed
    explicit Vector(int initialSize, ITEM const& value = ITEM()): size(0),
        capacity(max(initialSize, int(MIN_CAPACITY))),
        items(rawMemory<ITEM>(capacity))
        {for(int i = 0; i < initialSize; ++i) append(value);}
    Vector(Vector const& rhs): capacity(max(rhs.size, int(MIN_CAPACITY))),
        size(rhs.size), items(rawMemory<ITEM>(capacity))
        {for(int i = 0; i < size; ++i) new(&items[i])ITEM(rhs.items[i]);}
    Vector& operator=(Vector const& rhs){return genericAssign(*this, rhs);}
    ~Vector(){rawDestruct(items, size);}

    void append(ITEM const& item)
    {
        if(size >= capacity) resize();
        new(&items[size++])ITEM(item);
    }
    void removeLast()
    {
        assert(size > 0);
        items[--size].~ITEM();
        if(capacity > MIN_CAPACITY && size * 4 < capacity) resize();
    }
    void swapWith(Vector& other)
    {
        swap(items, other.items);
        swap(size, other.size);
        swap(capacity, other.capacity);
    }

    ITEM const& lastItem()const{return items[size - 1];}
    ITEM& lastItem(){return items[size - 1];}
    void reverse(int left, int right)
        {while(left < right) swap(items[left++], items[right--]);}
    void reverse(){reverse(0, size - 1);}
    void appendVector(Vector const& rhs)
        {for(int i = 0; i < rhs.getSize(); ++i) append(rhs[i]);}
    Vector& operator+=(Vector const& rhs)
    {
        assert(size == rhs.size);
        for(int i = 0; i < size; ++i) items[i] += rhs.items[i];
        return *this;
    }
    Vector& operator-=(Vector const& rhs)
    {
        assert(size == rhs.size);
        for(int i = 0; i < size; ++i) items[i] -= rhs.items[i];
        return *this;
    }
    template<typename SCALAR> Vector& operator*=(SCALAR const& scalar)
    {
        for(int i = 0; i < size; ++i) items[i] *= scalar;
        return *this;
    }
    friend Vector operator*(Vector const& a, ITEM const& scalar)
    {
        Vector result(a);
        return result *= scalar;
    }
    friend ITEM dotProduct(Vector const& a, Vector const& b)
    {
        assert(a.size == b.size);
        ITEM result(0);
        for(int i = 0; i < a.size; ++i) result += a[i] * b[i];
        return result;
    }
    Vector operator-()const{return *this * -1;}
    bool operator==(Vector const& rhs)const
    {
        if(size == rhs.size)
        {
            for(int i = 0; i < size; ++i)
                if(items[i] != rhs[i]) return false;
            return true;
        }
        return false;
    }

    void debug()const
    {
        cout << "[";
        for(int i = 0; i < getSize(); ++i)
        {
            cout << items[i];
            if(i + 1 < getSize()) cout  << ", ";
        }
        cout << "]" << endl;
    }
	void insert(int i, ITEM const& item)
	{
		assert(i >= 0 && i <= size);
		if(size >= capacity) resize();
		for(int j = size++; j > i; --j) items[j] = items[j-1];
		items[i] = item;
	}
};

double norm(Vector<double> const& x){return sqrt(dotProduct(x, x));}

}//end namespace
#endif
