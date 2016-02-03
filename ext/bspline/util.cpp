#include <string.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
#include "basis/util.h"

class pascal
{
	int *pas;
public:
	pascal(int K) : pas(new int[K]) {}
	~pascal() { delete pas; }
	int* operator[](int j)
	{
		int a, b;
		pas[0] = a = b = 1;
		for (int i = 1; i <= j; ++i) {
			a *= j - i + 1;
			b *= i;
			pas[i] = a / b;
		}
		return pas;
	}
};
/*******************************************************************************
	微分係数(differential coefficient)
	Ｎ個の関数の積のＪ階微分係数を求める。
	data : 微分係数の配列 [[f(t),D1f(t),...,DJf(t)],...]
	   s : 配列のストライド [st1,...,stN,K] (K = J+1)
	   N : 関数の積の項数
*******************************************************************************/
template<typename T> T *diff(T **data, const int *s, int N)
{
	int K = s[N];
	T *W = T_ALLOC(T,K);
	if (N == 1)
		for (int i = 0; i < K; ++i) W[i] = (*data)[(*s)*i];
	else {
		const T *b = *data;
			  T *c = diff(&data[1], &s[1], N-1);
		pascal pas(K);
		for (int j = 0; j < K; ++j) {
			int *a = pas[j];
			W[j] = 0;
			for (int k = 0; k <= j; ++k)
				W[j] += a[k] * b[(*s)*k] * c[j-k];
		}
		FREE(c);
	}
	return W;
}
template<typename T> T coeff(T **data, const int *s, int N, int jbn)
{
	T *W = diff(data, s, N);
	T result = W[jbn];
	FREE(W);
	return result;
}
/*******************************************************************************
	全微分(Total derivative)
	Ｎ個の関数の積のＪ階全微分を求める。
	data : 微分係数の配列 [[f(t),D1f(t),...,DJf(t)],...]
	   s : 配列のストライド [st1,...,stN,K] (K = J+1)
	   N : 関数の積の項数
	  ds : 方向ベクトルの成分の配列 [d1,...dN]
*******************************************************************************/
template<typename T>
static T *sub_derive(T **data, const int *s, int N, const T *ds)
{
	int k, K = s[N];
	T dt, *W = T_ALLOC(T,K);
	if (N == 1)
		for (dt = 1.0, k = 0; k < K; ++k, dt *= (*ds)) W[k] = (*data)[(*s)*k] * dt;
	else {
		const T *b = *data;
			  T *c = sub_derive(&data[1], &s[1], N-1, &ds[1]);
		pascal pas(K);
		for (int j = 0; j < K; ++j) {
			int *a = pas[j];
			W[j] = 0;
			for (dt = 1.0, k = 0; k <= j; ++k, dt *= (*ds))
				W[j] += a[k] * b[(*s)*k] * c[j-k] * dt;
		}
		FREE(c);
	}
	return W;
}
template<typename T> T total_derive(T **data, const int *s, int N, int jbn, const T *ds)
{
	T *W = sub_derive(data, s, N, ds);
	T result = W[jbn];
	FREE(W);
	return result;
}

/*******************************************************************************
	Matrix operation
*******************************************************************************/
/*
  Doolittle : LU分解 ドゥーリトル法

	Li0 = ai0;						i =  0,...,n-1
	Lij = aij - ΣLik*Ukj;			i >= j, k = 0,...,j-1
	Uij =(aij - ΣLik*Ukj)/Lii; 	i <  j, k = 0,...,i-1
*/
template <class T> void lud_decomp(T * a, size_t n, size_t * p, int & s)
{
	size_t i, j, k, l, L;

	s = 1;
	for (i = 0; i < n; i++) {
		T *ai = &a[i*n];
		T aii = ai[i]; L = i;
		for (j = 1; j < n; j++) {
			l = i < j ? i : j;
			T aij = ai[j];
			for (k = 0; k < l; k++) aij -= ai[k] * a[k*n+j];
			ai[j] = aij;
			// ピボット選択
			if ((j == i) || ((j > i) && (gabs(aii) < gabs(aij)))) { L = j; aii = aij; }
		}
		// ピボット列交換
		if (L != i) {
			for (k = 0; k < n; ++k) {
				T *au = &a[k*n];
				T  av = au[i]; au[i] = au[L]; au[L] = av;	// a[*,i] <=> a[*,L]
			}	s *= -1;
		}	p[i] = L;
		if (i < n-1) for (j = i+1; j < n; j++) ai[j] /= aii;
	}
}

template <class T, class S> void lud_subst(T * a, size_t n, size_t * p, S * b)
{
	size_t i, j, k, js = n; S sum;
	// 前進代入
	for (i = 0; i < n; ++i) {
		T *lu = &a[i*n];
		sum = b[i];
		if (js < n)
			for (j = js; j < i; ++j) sum -= lu[j] * b[j];
		else if (sum != 0.0) js = i;
		b[i] = sum / lu[i];
	}
	// 後退代入
	for (k = n-1; k > 0; --k) {
		i = k - 1;
		T *lu = &a[i*n];
		sum = b[i];
		for (j = n-1; j > i; --j) sum -= lu[j] * b[j];
		b[i] = sum;
	}
	// 解の保存
	for (k = n-1; k > 0; --k) {
		i = k - 1; j = p[i];
		if (i != j) { sum = b[j]; b[j] = b[i]; b[i] = sum; }
	}
}

template <class T, class S> void lud_subst(T * a, size_t n, size_t * p, S * b, int K)
{
	size_t i, j, k, l, js = n; S sum;
	// 前進代入
	for (i = 0; i < n; ++i) {
		T *lu = &a[i*n];
		for (l = 0; l < size_t(K); ++l) {
			sum = b[i*K+l];
			if (js < n)
				for (j = js; j < i; ++j) sum -= lu[j] * b[j*K+l];
			else if (sum != 0.0) js = i;
			b[i*K+l] = sum / lu[i];
		}
	}
	// 後退代入
	for (k = n-1; k > 0; --k) {
		i = k - 1;
		T *lu = &a[i*n];
		for (l = 0; l < size_t(K); ++l) {
			sum = b[i*K+l];
			for (j = n-1; j > i; --j) sum -= lu[j] * b[j*K+l];
			b[i*K+l] = sum;
		}
	}
	// 解の保存
	for (k = n-1; k > 0; --k) {
		i = k - 1; j = p[i];
		if (i != j) for (l = 0; l < size_t(K); ++l) {
			sum = b[j*K+l]; b[j*K+l] = b[i*K+l]; b[i*K+l] = sum;
		}
	}
}
	
template<class T> void lud_decomp(T ** a, size_t n, size_t * p, int & s)
{
	size_t i, j, k, L;
	double big, tmp, *v = new double[n];

	s = 1;
	for (i = 0; i < n; ++i) {
		big = 0.0;
		for (j = 0; j < n; ++j)
			if ((tmp = gabs(a[i][j])) > big) big = tmp;
		if (big == 0.0)
			throw "Singular matrix in routine lud_decomp";
		v[i] = 1.0 / big;
	}
	for (k = 0; k < n-1; ++k) {
		// 陰的ピボット選択
		L = k; big = 0.0;
		for (i = k+1; i < n; ++i)
			if ((tmp = gabs(a[i][k]) * v[i]) > big ) {
				big = tmp; L = i;
			}
		// ピボット行交換
		if (L != k) {
			T *w = a[k]; a[k] = a[L]; a[L] = w;	// a[k,*] <=> a[L,*]
			v[L] = v[k];
			  s *= -1;
		}	p[k] = L;
		// 前進消去
		T akk = a[k][k];
		for (j = k+1; j < n; ++j) {
			T akj = a[k][j] / akk;
			for (i = k+1; i < n; ++i) a[i][j] -= a[i][k] * akj;
			a[k][j] = akj;
		}
	}
	delete[] v;
}

template<class T, class S> void lud_subst(T ** a, size_t n, size_t * p, S * b)
{
	size_t i, j, k, js = n; S sum;
	// 前進代入
	for (i=0;i<n;++i) {
		T *lu = a[i];
		k = (i < n-1) ? p[i] : i;
		sum = b[k]; if (i != k) b[k] = b[i];
		if (js < n)
			for (j=js;j<i;++j) sum -= lu[j] * b[j];
		else if (sum != 0.0) js = i;
		b[i] = sum / lu[i];
	}
	// 後退代入
	for (k=n-1;k>0;--k) {
		i = k - 1;
		T *lu = a[i];
		sum = b[i];
		for (j=n-1;j>i;--j) sum -= lu[j] * b[j];
		b[i] = sum;
	}
}

template <class T, class S> void lud_subst(T ** a, size_t n, size_t * p, S * x, int K)
{
	S sum, *su, *sv, **b = create_marray_view(n, K, x);
	size_t i, j, k, l, js = n;
	for (i=0;i<n;++i)
	{
		T *lu = a[i];
		k = (i<n-1) ? p[i] : i;
		su = b[k]; sv = b[i];
		for (l=0;l<size_t(K);++l) {
			sum = su[l]; if (i != k) su[l] = sv[l];
			if (js < n)
				for (j=js;j<i;++j) sum -= lu[j] * b[j][l];
			else if (sum != 0.0) js = i;
			sv[l] = sum / lu[i];
		}
	}
	for (k=n;k>0;--k)
	{
		i = k - 1;
		T *lu = a[i];
		su =  b[i];
		for (l=0;l<size_t(K);++l) {
			sum = su[l];
			for (j=n-1;j>i;--j) sum -= lu[j] * b[j][l];
			su[l] = sum;
		}
	}
	FREE(b);
}
/*
  Crout : LU分解 クラウト法

	U0j = a0j;						j =  0,...,n-1
	Uij = aij - ΣLik*Ukj;			i <= j, k = 0,...,i-1
	Lij =(aij - ΣLik*Ukj)/Ujj; 	i >  j, k = 0,...,j-1
*/
template<class T> void luc_decomp(T * a, size_t n, size_t * p, int & s)
{
	size_t i, j, k, L;

	s = 1;
	for (k = 0; k < n-1; ++k) {
		T *ak = &a[k*n];
		T akk = ak[k]; L = k;
		// ピボット選択
		for (j = k+1; j < n; ++j)
			if (gabs(akk) < gabs(ak[j])) { L = j; akk = ak[j]; }
		// ピボット列交換
		if (L != k) {
			for (i = 0; i < n; ++i) {
			T *au = &a[i*n];
			T  av = au[k]; au[k] = au[L]; au[L] = av;	// a[*,k] <=> a[*,L]
			}  s *= -1;
		}	 p[k] = L;
		// 前進消去
		for (i = k+1; i < n; ++i) {
			T *ai = &a[i*n];
			T aik = ai[k] / akk;
			for (j = k+1; j < n; ++j) ai[j] -= aik * a[k*n+j];
			ai[k] = aik;
		}
	}
}

template<class T, class S> void luc_subst(T * a, size_t n, size_t * p, S * b)
{
	size_t i, j, js = n, k; S sum;
	// 前進代入
	for (i=0;i<n;++i) {
		T* lu = &a[i*n];
		sum = b[i];
		if (js < n)
			for (j=js;j<i;++j) sum -= lu[j] * b[j];
		else if (sum != 0.0) js = i;
		b[i] = sum;
	}
	// 後退代入
	for (k=n;k>0;--k) {
		i = k - 1;
		T* lu = &a[i*n];
		sum = b[i];
		for (j=n-1;j>i;--j) sum -= lu[j] * b[j];
		b[i] = sum / lu[i];
	}
	// 解の保存
	for (k=n-1;k>0;--k) {
		i = k - 1; j = p[i];
		if (i != j) { sum = b[j]; b[j] = b[i]; b[i] = sum; }
	}
}

template<class T, class S> void luc_subst(T * a, size_t n, size_t * p, S * b, int K)
{
	size_t i, j, k, l, js = n; S sum;
	// 前進代入
	for (i=0;i<n;++i) {
		T* lu = &a[i*n];
		for (l=0;l<size_t(K);++l) {
			sum = b[i*K+l];
			if (js < n)
				for (j=js;j<i;++j) sum -= lu[j] * b[j*K+l];
			else if (sum != 0.0) js = i;
			b[i*K+l] = sum;
		}
	}
	// 後退代入
	for (k=n;k>0;--k) {
		i = k - 1;
		T* lu = &a[i*n];
		for (l=0;l<size_t(K);++l) {
			sum = b[i*K+l];
			for (j=n-1;j>i;--j) sum -= lu[j] * b[j*K+l];
			b[i*K+l] = sum / lu[i];
		}
	}
	// 解の保存
	for (k = n-1; k > 0; --k) {
		i = k - 1; j = p[i];
		if (i != j) for (l=0;l<size_t(K);++l) {
			sum = b[j*K+l]; b[j*K+l] = b[i*K+l]; b[i*K+l] = sum;
		}
	}
}

template <class T> void luc_decomp(T ** a, size_t n, size_t * p, int & s)
{
	size_t i, j, k, l, L;
	double big, tmp, *v = new double[n];

	s = 1;
	for (i = 0; i < n; ++i) {
		big = 0.0;
		for (j = 0; j < n; ++j)
			if ((tmp = gabs(a[i][j])) > big) big = tmp;
		if (big == 0.0)
			throw "Singular matrix in routine luc_decomp";
		v[i] = 1.0 / big;
	}
	for (j = 0; j < n; j++) {
		L = j; big = 0.0;
		for (i = 1; i < n; i++) {
			l = j < i ? j : i;
			T aij = a[i][j];
			for (k = 0; k < l; k++) aij -= a[i][k] * a[k][j];
			a[i][j] = aij;
			// ピボット選択
			if (i >= j) 
			if ((tmp = gabs(aij) * v[i]) > big) { big = tmp; L = i; }
		}	p[j] = L;
		// ピボット行交換
		if (L != j) {
			T *w = a[j]; a[j] = a[L]; a[L] = w;	// a[j,*] <=> a[L,*]
			v[L] = v[j];
			  s *= -1;
		}  T ajj = a[j][j];
		if (ajj == 0.0) ajj = DBL_MIN;
		if (j < n-1) for (i = j+1; i < n; i++) a[i][j] /= ajj;
	}
	delete[] v;
}

template <class T, class S> void luc_subst(T ** a, size_t n, size_t * p, S * b)
{
	size_t i, j, k, js = n; S sum;
	// 前進代入
	for (i=0;i<n;++i) {
		T *lu = a[i]; k = p[i];
		sum = b[k]; if (i != k) b[k] = b[i];
		if (js < n) 
			for (j=js;j<i;++j) sum -= lu[j] * b[j];
		else if (sum != 0.0) js = i;
		b[i] = sum;
	}
	// 後退代入
	for (k=n;k>0;--k) {
		i = k - 1;
		T *lu = a[i];
		sum = b[i];
		for (j=n-1;j>i;--j) sum -= lu[j] * b[j];
		b[i] = sum / lu[i];
	}
}

template <class T, class S> void luc_subst(T ** a, size_t n, size_t * p, S * x, int K)
{
	S sum, *su, *sv, **b = create_marray_view(n, K, x);
	size_t i, j, k, l, js = n;
	// 前進代入
	for (i=0;i<n;++i) {
		T *lu = a[i]; k = p[i];
		su = b[k]; sv = b[i];
		for (l=0;l<size_t(K);++l) {
			sum = su[l]; if (i != k) su[l] = sv[l];
			if (js < n) 
				for (j=js;j<i;++j) sum -= lu[j] * b[j][l];
			else if (sum != 0.0) js = i;
			sv[l] = sum;
		}
	}
	// 後退代入
	for (k=n;k>0;--k) {
		i = k - 1;
		T *lu = a[i];
		su = b[i];
		for (l=0;l<size_t(K);++l) {
			sum = su[l];
			for (j=n-1;j>i;--j) sum -= lu[j] * b[j][l];
			su[l] = sum / lu[i];
		}
	}
	free(b);
}

template double coeff(double**, const int*, int, int);
template double total_derive(double**, const int*, int, int, const double*);
template void luc_decomp(double *, size_t, size_t*, int&);
template void luc_decomp(double**, size_t, size_t*, int&);
template void lud_decomp(double *, size_t, size_t*, int&);
template void lud_decomp(double**, size_t, size_t*, int&);
template void luc_subst(double *, size_t, size_t*, double*);
template void luc_subst(double**, size_t, size_t*, double*);
template void lud_subst(double *, size_t, size_t*, double*);
template void lud_subst(double**, size_t, size_t*, double*);
template void luc_subst(double *, size_t, size_t*, double*, int);
template void luc_subst(double**, size_t, size_t*, double*, int);
template void lud_subst(double *, size_t, size_t*, double*, int);
template void lud_subst(double**, size_t, size_t*, double*, int);
