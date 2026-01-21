#ifndef FFT_h
#define FFT_h

#include <cmath>

enum fft_type {FFT_FORWARD, FFT_BACKWARD, FFT_INVERSE};

/**
 * Fast Fourier Transform.
 * Forward Fourier transform is
 *         A_k + i B_k := sum((A_n + i B_n)*exp(-i 2 pi k n/N), n = 0..N-1)
 * For the type parameter values FFT_FORWARD, FFT_BACKWARD and FFT_INVERSE,
 * the fft function is equivalent to the GSL functions
 * gsl_fft_complex_radix2_forward, gsl_fft_complex_radix2_backward and
 * gsl_fft_complex_radix2_inverse, respectively.
 * @param A     data array, real part
 * @param B     data array, imaginary part
 * @param N     size of data array
 * @param type  type of transform, FFT_FORWARD, FFT_BACKWARD or FFT_INVERSE
 * @since 06/26/1992
 * @author Peter Csizmadia
 */
template <class T>
int fft(T A[], T B[], int N, enum fft_type type)
{
    int nstages; // = log2(N)
    for(nstages = 8*sizeof(int)-1; nstages >= 0 && !(N & (1 << nstages)); --nstages);

    if(N != (1<<nstages)) {
	return -1;
    }
    for(int i = 1; i < N; ++i) {
	int n = 0;
	for(int j = 1, k = N>>1; k != 0; j <<= 1, k >>=1 ) {
	    if(i & j) {
		n |= k;
	    }
	}
	if(n > i) {
	    T w;
	    w = A[n]; A[n] = A[i]; A[i] = w;
	    w = B[n]; B[n] = B[i]; B[i] = w;
	}
    }
    if(type == FFT_BACKWARD || type == FFT_INVERSE) {
	if(type == FFT_INVERSE) {
	    for(int i = 0; i < N; ++i) {
		A[i] /= N;
	    }
	    for(int i = 0; i < N; ++i) {
		B[i] /= -N;
	    }
	} else {
	    for(int i = 0; i < N; ++i) {
		B[i] /= -B[i];
	    }
	}
    }
    int n = 1;
    for(int i = 0; i < nstages; ++i) {
	int prev_n = n;
	T phi = (T)M_PI/n;
	T vr = std::cos(phi);
	T vi = -std::sin(phi);
	T wr = 1;
	T wi = 0;
	n <<= 1;
	for(int j = 0; j < prev_n; ++j) {
	    for(int g = j; g < N; g += n) {
		int h = g + prev_n;
		T tmpr = A[h]*wr - B[h]*wi;
		T tmpi = A[h]*wi + B[h]*wr;
		A[h] = A[g] - tmpr;
		A[g] += tmpr;
		B[h] = B[g] - tmpi;
		B[g] += tmpi;
	    }
	    T tmp = wr*vr - wi*vi;
	    wi = wr*vi + wi*vr;
	    wr = tmp;
	}
    }
    return 0;
}

#endif /* FFT_h */
