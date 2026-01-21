#ifndef _tvector_h
#define _tvector_h

#if defined(DEBUG) || defined(ARRAY_BOUNDS_CHECK)
#include <iostream>
#include <cstdlib> /* abort() */
#ifndef _COREDEBUG_STR
#define _COREDEBUG_STR "To debug the core, use this command:\n"\
"    gdb -c core /path/to/executable             (gdb) bt"
#endif
#endif

/**
 * Array class with debugging support and std::vector-compatible methods.
 * Extra method: erase(T);
 *
 * @version 0.3, 02/05/2003
 * @author Peter Csizmadia (1993-95,99-2001)
 * @author Sen Cheng (2001)
 */
template <class T>
class tvector
{
public:
    /** Iterator. */
    typedef T* iterator;

    /** Constant iterator. */
    typedef const T* const_iterator;

private:
    T* x;
    unsigned max;
    unsigned n;

public:
    /** Create an empty array. */
    tvector() {
	x = 0;
	max = 0;
       	n = 0;
    }

    /** Create an array with n elements. */
    tvector(unsigned n) {
//	if(n < 0)
//	    n = 0;
	this->n = max = n;
	if(n > 0) {
	    x = new T[n];
	    fill(begin(), begin() + n, T());
	} else {
	    x = 0;
	}
    }

    /** Create an array with n initialized elements. */
    tvector(unsigned n, const T& value);

    /** Copy constructor. */
    tvector(const tvector& a) {
	n = max = a.n;
	if(n > 0) {
	    x = new T[n];
	    a.copyInto(x);
	} else {
	    x = 0;
	}
    }

    ~tvector() {
       if(x != 0) {
	    delete [] x;
       }
    }

    /** Gets the number of elements. */
    unsigned size() const {
	return n;
    }

    /** Gets the capacity. */
    unsigned capacity() const {
	return max;
    }

    /** Gets the underlying array. */
    operator T*() {
	return x;
    }

    /** Gets the underlying array. */
    operator const T*() const {
	return x;
    }

    /** Copy. */
    tvector& operator =(const tvector& a) {
	if(n != a.n) {
	    resize_nofill(a.n);
	}
	a.copyInto(x);
       	return *this;
    }

    /** Remove all elements and delete the underlying array. */
    void clear() {
	delete[] x;
	x = 0;
	n = max = 0;
    }

    /** Erase an element. */
    void erase(iterator first) {
	erase(first, first + 1);
    }

    /** Erase elements. */
    void erase(iterator first, iterator last);

    /** Erase elements with the specified value. */
    void erase(const T& x);

    /** Resize array and initialize the new elements. */
    void resize(unsigned n) {
	resize(n, T());
    }

    /** Resize array and initialize the new elements. */
    void resize(unsigned n, const T& value) {
	unsigned prev = this->n;
	resize_nofill(n);
	if(this->n > prev) {
	    fill(begin() + prev, begin() + this->n, value);
	}
    }

    /** Insert a new element. */
    void insert(iterator p, const T& value);

    /** Insert new elements. */
    void insert(iterator p, const T* first, const T* last);

    /** Add a new element. */
    void push_back(const T& value) {
	if(n < max) {
	    x[n++] = value;
	} else {
	    resize_nofill(n + 1);
	    x[n - 1] = value;
	}
    }

    /** Remove the last element. */
    void pop_back() {
#if defined(DEBUG) || defined(ARRAY_BOUNDS_CHECK)
	if(n == 0) {
	    std::cerr << "Fatal error: pop_back() called for empty tvector\n"
		 _COREDEBUG_STR << std::endl;
	    ::abort();
	}
#endif
	--n;
    }

    /** Gets the kth element. */
    T& operator [](int k) {
#if defined(DEBUG) || defined(ARRAY_BOUNDS_CHECK)
	if(k < 0 || k >= (int)n) {
	    std::cerr << "Fatal error: tvector index " << k
		 << " out of range [0, "
		 << n << "[\n" _COREDEBUG_STR << std::endl;
	    ::abort();
	}
#endif
	return x[k];
    }

// method commented out because of gcc 3.3.3 compilation error:
// In member function `virtual gromit::TwoBodyInter*
// gromit_interactions_impl::CollisionTableImp::getInteraction(int, int) const':
// error: ISO C++ says that `const T& tvector<T>::operator[](int) const [with
// T = gromit::TwoBodyInter*]' and `operator[]' are ambiguous even though the
// worst conversion for the former is better than the worst conversion for the
// latter (06/22/2004 cspeter)
//    /** Gets the kth element. */
//    const T& operator [](int k) const {
//#if defined(DEBUG) || defined(ARRAY_BOUNDS_CHECK)
//	if(k < 0 || k >= (int)n) {
//	    std::cerr << "Fatal error: tvector index " << k
//		 << " out of range [0, "
//		 << n << "[\n" << _COREDEBUG_STR << std::endl;
//	    ::abort();
//	}
//#endif
//	return x[k];
//    }

    /** Gets the first element. */
    T& front() {
	return *begin();
    }

    /** Gets the first element as a constant. */
    const T& front() const {
	return *begin();
    }

    /** Gets the last element. */
    T& back() {
	return *(end() - 1);
    }

    /** Gets the last element as a constant. */
    const T& back() const {
	return *(end() - 1);
    }

    /** Returns an iterator to the beginning. */
    iterator begin() {
	return x;
    }

    /** Returns an iterator to the beginning. */
    const_iterator begin() const {
	return x;
    }

    /** Returns an iterator to end. */
    iterator end() {
	return x+n;
    }

    /** Returns an iterator to end. */
    const_iterator end() const {
	return x+n;
    }

    /** quicksort the contend of the vector, T must define operator < */
    void sort();

private:
    /** Resize array without initializing the new elements. */
    void resize_nofill(unsigned _n);

    /** Copy tvector elements into a normal array. */
    void copyInto(T* y) const {
	for(unsigned i = 0; i < n; ++i) {
	    y[i] = x[i];
	}
    }

    /** Initialize elements in the specified range. */
    void fill(iterator begin, iterator end, const T& value);

    /** helper function for quicksort */
    void sort(int lo, int hi);
};

template <class T> tvector<T>::tvector(unsigned _n, const T& value)
{
//    if(_n < 0)
//	_n = 0;
    n = max = _n;
    if(n > 0) {
	x = new T[n];
	fill(begin(), begin() + n, value);
    } else {
	x = 0;
    }
}

template <class T> void tvector<T>::erase(iterator first, iterator last)
{
    int d = last - first;
    if(d <= 0) {
	return;
    }
    for(unsigned i = last - x; i < n; ++i) {
	x[i - d] = x[i];
    }
    n -= d;
}

template <class T> void tvector<T>::erase(const T& v)
{
    for(int i = n - 1; i >= 0; --i) {
	if(x[i] == v) {
	    erase(begin() + i);
	}
    }
}
 
template <class T> void tvector<T>::insert(iterator p, const T& value)
{
    unsigned k = p - x;
    if(n < max) {
	for(unsigned i = n; i > k; --i) {
	    x[i] = x[i - 1];
	}
	x[k] = value;
    } else {
	max = 2*(n + 1);
	T* y = new T[max];
	for(unsigned i = 0; i < k; ++i) {
	    y[i] = x[i];
	}
	for(unsigned i = k; i < n; ++i) {
	    y[i + 1] = x[i];
	}
	y[k] = value;
	delete [] x;
	x = y;
    }
    ++n;
}
 
template <class T> void tvector<T>::insert(iterator p,
					   const T* first, const T* last)
{
    int delta = last - first + 1;
#if defined(DEBUG) || defined(ARRAY_BOUNDS_CHECK)
	if(delta < 1) {
	    std::cerr << "Fatal error: cannot insert zero or less elements "
		    "into tvector\n"
		 _COREDEBUG_STR << std::endl;
	    ::abort();
	}
#endif
    int k = p - x;
    int k_plus_delta = k + delta;
    if(n + delta <= max) {
	for(int i = n - 1; i >= k; --i) {
	    x[i + delta] = x[i];
	}
	for(int i = k; i < k_plus_delta; ++i) {
	    x[i] = *(first++);
	}
    } else {
	max = 2*(n + delta);
	T* y = new T[max];
	for(int i = 0; i < k; ++i) {
	    y[i] = x[i];
	}
	for(int i = n - 1; i >= k; --i) {
	    y[i + delta] = x[i];
	}
	for(int i = k; i < k_plus_delta; ++i) {
	    y[i] = *(first++);
	}
	delete [] x;
	x = y;
    }
    n += delta;
}

template <class T> void tvector<T>::resize_nofill(unsigned _n)
{
//    if(_n < 0) {
//	n = 0;
//    } else
    if(_n <= max) {
	n = _n;
    } else {
	for(max = 1; max < _n; max <<= 1);
	T* y = new T[max];
	copyInto(y);
	delete [] x;
	x = y;
	n = _n;
    }
}

template <class T> void tvector<T>::fill(iterator p, iterator end, const T& val)
{
    while(p != end) {
	*(p++) = val;
    }
}

// quicksort the contend of the vector, T must define operator <
template <class T> void tvector<T>::sort()
{
    sort(0, n-1);
}

// quicksort the contend of the vector, T must define operator <
template <class T> void tvector<T>::sort(int lo, int hi)
{
    if(lo>=hi) return;
    T ref= x[(lo+hi)/2];
    int left= lo, right= hi;
    while(left <= right) {
	while(left < hi && x[left] < ref) left++;
	while(right > lo && ref < x[right]) right--;
	if(left < right) {
	    T tmp= x[left];
	    x[left]= x[right];
	    x[right]= tmp;
	    left++; right--;
	} else if(left== right) {
	    left++; right--;
	}
    }
    sort(lo, right);
    sort(left, hi);
}

#ifdef DEBUG
template <class T>
inline std::ostream& operator <<(std::ostream& out, const tvector<T>& x)
{
    for(unsigned i = 0; i < x.size(); ++i) {
	out << x[i];
	if(i+1 < x.size() ) {
	    out << ", ";
	}
    }
    return out;
}
#endif /* DEBUG */

#endif /* _tvector_h */
