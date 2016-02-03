#ifndef _UTIL_H_INCLUDED_
#define _UTIL_H_INCLUDED_

#include <complex.h>

//   Ｔ配列	[T0,...,Ti-1] : T*
#define T_ALLOC(T,i) (T*)malloc((i)*sizeof(T))
//   Ｔ行列	[[T00,...],...] : T**
#define T_MALLOC(T,i,j) create_marray<T>((i),(j))
// destructor
#define FREE(m) free((char*)(m))

template <typename T> inline
T pow(T a, int n)
{
	return n == 0 ? 1 : a * pow(a, n-1);
}

template <typename T> T coeff(T**, const int*, int, int = 0);
template <typename T> T total_derive(T**, const int*, int, int, const T*);

/*******************************************************************************
	define print utility.
*******************************************************************************/
template <class T> inline int sprintg(char* &c, T o);

template <> inline int sprintg(char* &c, int a)    {return sprintf(c, "%d", a);}
template <> inline int sprintg(char* &c, size_t a) {return sprintf(c, "%ld",a);}
template <> inline int sprintg(char* &c, double a) {return sprintf(c, "%g", a);}

template <class T> inline double gabs(T o);

template <> inline double gabs(double a) {return fabs(a);}
template <> inline double gabs(std::complex<double> c) {return abs(c);}

template <typename T> void poly_to_s(char* &c, int N, const int *size, const T *data)
{
	int n = *size;
	if (N == 0) {
		if (n > 1) sprintf(c, "("); while(*c)++c;
		for (int i = 0; i < n; ++i) {
			sprintg<T>(c, data[i]); while(*c)++c;
			if (i < n-1) sprintf(c, ","); while(*c)++c;
		}
		if (n > 1) sprintf(c, ")"); while(*c)++c;
	}
	else {
		sprintf(c, "["); while(*c)++c;
		int m = 1; for (int i = 1; i <= N; ++i) m *= size[i];
		for (int i = 0; i < n; ++i) {
			poly_to_s(c, N-1, size+1, &data[i*m]);
			if (i < n-1) sprintf(c, ","); while(*c)++c;
		}
		sprintf(c, "]"); while(*c)++c;
	}
}

/*******************************************************************************
	Matrix operation
*******************************************************************************/
template <class T>
void luc_decomp(T *, size_t, size_t*, int&);

template <class T>
void lud_decomp(T *, size_t, size_t*, int&);

template <class T>
void luc_decomp(T**, size_t, size_t*, int&);

template <class T>
void lud_decomp(T**, size_t, size_t*, int&);

template <class T, class S = T>
void luc_subst(T *, size_t, size_t*, S*);

template <class T, class S = T>
void lud_subst(T *, size_t, size_t*, S*);

template <class T, class S = T>
void luc_subst(T**, size_t, size_t*, S*);

template <class T, class S = T>
void lud_subst(T**, size_t, size_t*, S*);

template <class T, class S = T>
void luc_subst(T *, size_t, size_t*, S*, int);

template <class T, class S = T>
void lud_subst(T *, size_t, size_t*, S*, int);

template <class T, class S = T>
void luc_subst(T**, size_t, size_t*, S*, int);

template <class T, class S = T>
void lud_subst(T**, size_t, size_t*, S*, int);

template <class T> T luc_solve(T * a, size_t n, T * b, int K = 1)
{
	size_t *p = new size_t[n];
	int s;	luc_decomp(a, n, p, s);
	if (K>1) luc_subst(a, n, p, b, K); else luc_subst(a, n, p, b);
	T det = s;
	for (size_t k = 0; k < n; ++k) det *= a[k*n+k];
	delete[] p;
	return det;
}

template <class T> T luc_solve(T ** a, size_t n, T * b, int K = 1)
{
	size_t *p = new size_t[n];
	int s;	luc_decomp(a, n, p, s);
	if (K>1) luc_subst(a, n, p, b, K); else luc_subst(a, n, p, b);
	T det = s;
	for (size_t k=0;k<n;++k) det *= a[k][k];
	delete[] p;
	return det;
}

template <class T> T lud_solve(T * a, size_t n, T * b, int K = 1)
{
	size_t *p = new size_t[n];
	int s;	lud_decomp(a, n, p, s);
	if (K>1) lud_subst(a, n, p, b, K); else lud_subst(a, n, p, b);
	T det = s;
	for (size_t k = 0; k < n; ++k) det *= a[k*n+k];
	delete[] p;
	return det;
}

template <class T> T lud_solve(T ** a, size_t n, T * b, int K = 1)
{
	size_t *p = new size_t[n];
	int s;	lud_decomp(a, n, p, s);
	if (K>1) lud_subst(a, n, p, b, K); else lud_subst(a, n, p, b);
	T det = s;
	for (size_t k = 0; k < n; ++k) det *= a[k][k];
	delete[] p;
	return det;
}

#define lu_decomp lud_decomp
#define lu_subst  lud_subst
#define lu_solve  lud_solve

template <class T>
T ** marray_alloc(char *mm, int nr, size_t sr, size_t sc)
{
	T **m = (T **)mm;
	 m[0] = (T *)(mm = mm + sr);
	for (int i = 1; i < nr; ++i) m[i] = (T *)(mm = mm + sc);
	return m;
}

template <typename T>
T ** create_marray(int nr, int nc)
{
	size_t sr = nr * sizeof(T *);
	size_t sc = nc * sizeof(T);
	char * mm = (char *)malloc(sr + nr * sc);
	if (mm == NULL) throw "allocate error, create_marray";
	return marray_alloc<T>(mm, nr, sr, sc);
}

template <typename T>
T ** marray_view_alloc(T **m, int nr, int nc, T *v)
{
	m[0] = v;
	for (int i = 1; i < nr; i++) m[i] = m[i-1] + nc;
	m[nr] = v;
	return m;
}

template <typename T>
T ** create_marray_view(int nr, int nc, T *v)
{
	T **m = (T**)malloc((nr+1) * sizeof(T*));
	return marray_view_alloc(m, nr, nc, v);
}

template <typename T>
T ** carray_alloc(char *mm, size_t sr, int n, int *s)
{
	T **m = (T**)mm;
	 m[0] = (T*)(mm = mm + sr);
	for (int i = 1; i < n; ++i) m[i] = (T*)(mm = mm + s[i-1] * sizeof(T));
	return m;
}

template <typename T>
T ** create_carray(int n, int *s)
{
	int c = 0;
	for (int i = 0; i < n; ++i) c += s[i];
	size_t sr = n * sizeof(T *);
	size_t sc = c * sizeof(T);
	char * mm = (char *)malloc(sr + sc);
	if (mm == NULL) throw "allocate error, create_carray";
	return carray_alloc<T>(mm, sr, n, s);
}

template <typename T>
T ** carray_view_alloc(T **m, int n, int *s, T *a)
{
	m[0] = a;
	for (int i = 1; i < n; ++i) m[i] = m[i-1] + s[i-1];
	m[n] = a;
	return m;
}

template <typename T>
T ** create_carray_view(int n, int *s, T *a)
{
	T **m = (T**)malloc((n+1) * sizeof(T*));
	return carray_view_alloc<T>(m, n, s, a);
}

template <typename T> class varray
{
	T *data; int s;
  public:
	varray(int n) : s(n)
	{
		data = (T*)malloc(s * sizeof(T));
	}
	varray(const varray& v) : s(v.s)
	{
		data = (T*)malloc(s * sizeof(T));
		memcpy(data, v.data, s * sizeof(T));
	}
	~varray() { free((void*)data); }
	varray& operator = (const varray& v)
	{
		if (this != &v) {
			free((void*)data);
			s = v.s;
			data = (T*)malloc(s * sizeof(T));
			memcpy(data, v.data, s * sizeof(T));
		}
		return *this;
	}
	operator T* () { return data; }
T	operator * () { return *data; }
T & operator [] (int i) { return data[i]; }
	int size() const { return s; }
};

template <typename T> class varray_view
{
	T *data; int s, t;
  public:
	varray_view(T *v, int n, int st = 1) : data(v), s(n), t(st) {}
	~varray_view() {}
	varray_view& operator = (const varray_view& v)
	{
		if (this != &v) {
			assert(s < v.s);
			if (s > v.s)
			{
				int i, j, k(s - v.s), m((k + 1) / 2), o(k - m);
				for (i=m-1,  k=s-m;i>=0;--i) { j=i+k; data[j*t] = v.data[i*v.t]; }
				for (i=v.s-1,k=o;  i>=0;--i) { j=i+k; data[j*t] = v.data[i*v.t]; }
				for (i=o-1,  k=v.s;i>=0;--i) { j=i+k; data[i*t] =   data[j*t];   }
			}
			else if (data != v.data)
				for (int i = 0; i < s; ++i) data[i*t] = v.data[i*v.t];
		}
		return *this;
	}
	operator T* () { return data; }
T	operator * () { return *data; }
T & operator [] (int i) const { return data[i*t]; }
	int size() const { return s; }
};

template <typename T> class marray
{
	T **index; int r, c;
  public:
	marray(int nr, int nc) : r(nr), c(nc)
	{
		index = create_marray<T>(nr, nc);
	}
	marray(const marray& m) : r(m.r), c(m.c)
	{
		index = create_marray<T>(r, c);
		memcpy(*index, *m, r * c * sizeof(T));
	}
	~marray() { free((void*)index); }
	marray& operator = (const marray& m)
	{
		if (this != &m) {
			free((void*)index);
			r = m.r; c = m.c;
			index = create_marray<T>(r, c);
			memcpy(*index, *m, r * c * sizeof(T));
		}
		return *this;
	}
    operator T** () const { return index; }
T * operator * () const { return (T*)(index + r); }
T * operator [] (int i) { return index[i]; }
	void swap(int i, int j) { T *w = index[i]; index[i] = index[j]; index[j] = w; }
	int rows() const { return r; }
	int cols() const { return c; }
	varray_view<T> row(int i) const { return varray_view<T>(index[i], c); }
	varray_view<T> col(int i) const { return varray_view<T>((T*)(index + r) + i, r, c); }
};

template <typename T> class marray_view
{
	T **index; int r, c;
  public:
	marray_view(T *m, int nr, int nc) : r(nr), c(nc)
	{
		index = create_marray_view(r, c, m);
	}
	marray_view(const marray<T>& m) : r(m.rows()), c(m.cols())
	{
		index = create_marray_view(r, c, *m);
	}
	marray_view(const marray_view& m) : r(m.r), c(m.c)
	{
		index = create_marray_view<T>(r, c, *m);
	}
	~marray_view() { free((void*)index); }
	marray_view& operator = (const marray_view& m)
	{
		if (this != &m) {
			free((void*)index);
			r = m.r; c = m.c;
			index = create_marray_view<T>(r, c, *m);
		}
		return *this;
	}
	operator T** () { return index; }
T * operator * () { return index[r]; }
T * operator [] (int i) const { return index[i]; }
	void swap(int i, int j) { T *w = index[i]; index[i] = index[j]; index[j] = w; }
	int rows() const { return r; }
	int cols() const { return c; }
	varray_view<T> row(int i) const { return varray_view<T>(index[i], c); }
	varray_view<T> col(int i) const { return varray_view<T>(index[r] + i, r, c); }
};

template <typename T> class marray_lu
{
	int  n;
	T ** m;
	size_t *p;
	int  s;
  public:
	marray_lu(int N, const T *dat) : n(N) {
		m = create_marray<T>(N, N);
		p = (size_t*)malloc(N * sizeof(size_t));
		memcpy((T*)(m+n), dat, n * n * sizeof(T));
		lu_decomp<T>(m, size_t(N), p, s);
	}
	marray_lu(const marray_lu& lu) : n(lu.n) {
		m = create_marray<T>(n, n);
		p = (size_t*)malloc(n * sizeof(size_t));
		memcpy((T*)(m+n), *lu, n * n * sizeof(T));
		memcpy(p, lu.p, n * sizeof(size_t));
		s = lu.s;
	}
	~marray_lu() { FREE(p); FREE(m); }
	marray_lu& operator= (const marray_lu& lu) {
		if (this != &lu) {
			FREE(p); FREE(m);
			n = lu.n;
			m = create_marray<T>(n, n);
			p = (size_t*)malloc(n * sizeof(size_t));
			memcpy((T*)(m+n), *lu, n * n * sizeof(T));
			memcpy(p, lu.p, n * sizeof(size_t));
			s = lu.s;
		}
		return *this;
	}
	T * operator[](int i) { return m[i]; }
	T * operator * () { return (T*)(m + n); }
	template <class S>
	void solve(S * b, S * x, int K = 1) {
		for (int i = 0; i < n*K; ++i) x[i] = b[i];
		if (K > 1)
			lu_subst<T,S>(m, size_t(n), p, x, K);
		else
			lu_subst<T,S>(m, size_t(n), p, x);
	}
	template <class S>
	S *solve(S *b, int K = 1) {
		S *result = new S[n*K];
		for (int i = 0; i < n*K; ++i) result[i] = b[i];
		if (K > 1)
			lu_subst<T,S>(m, size_t(n), p, result, K);
		else
			lu_subst<T,S>(m, size_t(n), p, result);
		return result;
	}
// end of marray_lu;
};

template <typename T>
class marray_lu_view
{
	int  n;
	T ** m;
	int  s;
	size_t *p;
  public:
	marray_lu_view(int N, T *dat) : n(N) {
		m = create_marray_view<T>(n, n, dat);
		p = (size_t*)malloc(n * sizeof(size_t));
		lu_decomp<T>(m, size_t(n), p, s);
	}
	marray_lu_view(const marray<T>& mm) : n(mm.rows()) {
		m = create_marray_view(n, n, *mm);
		p = (size_t*)malloc(n * sizeof(size_t));
		lu_decomp<T>(m, size_t(n), p, s);
	}
	marray_lu_view(const marray_lu<T>& lu) : n(lu.n) {
		m = create_marray_view(n, n, *lu);
		p = (size_t*)malloc(n * sizeof(size_t));
		memcpy(p, lu.p, n * sizeof(size_t));
		s = lu.s;
	}
	marray_lu_view(const marray_lu_view& lu) : n(lu.n) {
		m = create_marray_view(n, n, *lu);
		p = (size_t*)malloc(n * sizeof(size_t));
		memcpy(p, lu.p, n * sizeof(size_t));
		s = lu.s;
	}
	~marray_lu_view() { FREE(p); FREE(m); }
	marray_lu_view& operator= (const marray_lu_view& lu) {
		if (this != &lu) {
			FREE(p); FREE(m);
			n = lu.n;
			m = create_marray_view(n, n, *lu);
			p = (size_t*)malloc(n * sizeof(size_t));
			memcpy(p, lu.p, n * sizeof(size_t));
			s = lu.s;
		}
		return *this;
	}
	T * operator[](int i) { return m[i]; }
	T * operator * () { return m[n]; }
	template <class S>
	void solve(S * b, S * x, int K = 1) {
		for (int i = 0; i < n*K; ++i)
			x[i] = b[i];
		if (K > 1)
			lu_subst<T,S>(m, n, p, x, K);
		else
			lu_subst<T,S>(m, n, p, x);
	}
	template <class S>
	S *solve(S *b, int K = 1) {
		S *result = new S[n*K];
		for (int i = 0; i < n*K; ++i)
			result[i] = b[i];
		if (K > 1)
			lu_subst<T,S>(m, n, p, result, K);
		else
			lu_subst<T,S>(m, n, p, result);
		return result;
	}
// end of marray_lu_view;
};

template <typename T>
void kronecker(const varray_view<T>& p, const varray_view<T>& q, const varray_view<T>& r)
{
	int s = p.size(), t = q.size();
	for (int i = 0; i < s; i++)
		for (int j = 0; j < t; j++) 
			r[i*t+j] = p[i] * q[j];
}

template <typename T>
void kronecker(const marray<T>& p, const marray<T>& q, const marray<T>& r)
{
	int s = p.rows(), t = q.rows();
	for (int i = 0; i < s; i++)
		for (int j = 0; j < t; j++)
			kronecker(p.row(i), q.row(j), r.row(i*t+j));
}

template <typename T>
marray<T> operator^(const marray<T>& p, const marray<T>& q)
{
	int r = p.rows() * q.rows(), c = p.cols() * q.cols();
	marray<T> result(r, c);
	kronecker(p, q, result);
	return result;
}

template <typename T> class parray
{
	T **data; int o, s, k;
  public:
	parray(int no, int ns, int nk) : o(no), s(ns), k(nk)
	{
		data = create_marray<T>(k, s);
	}
	~parray() { free((void*)data); }
    operator T** () const { return data; }
T * operator * () const { return (T*)(data + k); }
T * operator [] (int i) { return data[i]; }
int operator () (int i) { return i == 0 ? o : i == 1 ? s : k; }
	int offset() const { return o; }
	int size() const { return s; }
	int rows() const { return k; }
	varray_view<T> row(int i) const { return varray_view<T>(data[i], s); }
	varray_view<T> col(int i) const { return varray_view<T>((T*)(data + k) + i, k, s); }
	T saylor(int i, T h)
	{
		T q = 0.0;
		for (int j=k;j>0;--j) q = (q + data[j-1][i]) * h / (T)j;
		return q;
	}
	T taylor(int i, T h)
	{
		T q = 0.0;
		for (int j=k;j>0;--j) q = q * h / (T)j + data[j-1][i];
		return q;
	}
};

typedef int * narray;
typedef varray<double> vector;
typedef marray<double> matrix;

// 整数配列	[I0,...,Ii-1] : narray
#define NALLOC(i) (int*)malloc((i)*sizeof(int))
// 実数配列	[X0,...,Xi-1] : varray
#define VALLOC(i) (double*)malloc((i)*sizeof(double))
//  dis配列	[D0,...,Di-1] : di_array
#define DI_ALLOC(i) (dis*)malloc((i)*sizeof(dis))
//  bis配列	[B0,...,Bi-1] : bi_array
#define BI_ALLOC(i) (bis*)malloc((i)*sizeof(bis))
// tris配列	[R0,...,Ri-1] : tri_array
#define TRI_ALLOC(i) (tris*)malloc((i)*sizeof(tris))
//tetra配列	[E0,...,Ei-1] : tetra_array
#define TETRA_ALLOC(i) (tetris*)malloc((i)*sizeof(tetris))
// 実数行列	[[X00,...,X0(c-1)],...,[X(r-1)0,...,X(r-1)(c-1)]] : marray
#define MALLOC(r,c) create_marray<double>(r,c)
//  dis行列 [[D00,...,D0(c-1)],...,[D(r-1)0,...,D(r-1)(c-1)]] : di_marray
#define DI_MALLOC(r,c) create_marray<dis>(r,c)
//  bis行列 [[B00,...,B0(c-1)],...,[B(r-1)0,...,B(r-1)(c-1)]] : bi_marray
#define BI_MALLOC(r,c) create_marray<bis>(r,c)
// tris行列 [[R00,...,R0(c-1)],...,[R(r-1)0,...,R(r-1)(c-1)]] : tri_marray
#define TRI_MALLOC(r,c) create_marray<tris>(r,c)
//tetra行列 [[E00,...,E0(c-1)],...,[E(r-1)0,...,E(r-1)(c-1)]] : tetra_marray
#define TETRA_MALLOC(r,c) create_marray<tetris>(r,c)
// 実数行列	[[X00,...],...] : double**
#define S_MALLOC(d) create_marray<double>(d[0],d[1])
//  dis行列 [[D00,...],...] : dis**
#define DIS_MALLOC(d) create_marray<dis>(d[0],d[1])
//  bis行列 [[B00,...],...] : bis**
#define BIS_MALLOC(d) create_marray<bis>(d[0],d[1])
// tris行列 [[R00,...],...] : tris**
#define TRIS_MALLOC(d) create_marray<tris>(d[0],d[1])
//tetra行列 [[E00,...],...] : tetris**
#define TETRIS_MALLOC(d) create_marray<tetris>(d[0],d[1])
// ポインター演算
#define CDR(T,N,p) (poly<T,N>*)(&(p)+1)
#define CAR(T,N,p) *(poly<T,N>*)&(p)

template<class T, class Object, T (Object::*getter) () const> class Rproperty
{
Object *orner;
public:
void operator()(Object* obj) { orner = obj; }
operator T () { return (orner->*getter)(); }
};

#endif
