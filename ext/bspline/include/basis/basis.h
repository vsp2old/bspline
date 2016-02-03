#ifndef _BASIS_H_INCLUDED_
#define _BASIS_H_INCLUDED_
#include <stdlib.h>
#include <string.h>

template <typename T>
int deboor(int, T*, int, T, T*, int&, int = 0);
template<typename T>
void deboor_cox(int, T*, int, T, T*, int&, int, int = 0);

/*******************************************************************************
	Ｂスプライン基底（節点列によるＢスプライン関数計算）
*******************************************************************************/
template <typename S = double> class base_spline
{
  protected:
	int imax, icox, rank, jisu, shuki, maxq;
	S *knots;
  public:
	// constructor
	base_spline() : knots(NULL) {}
	base_spline(int, const S*, int, int = 0, int* = NULL);
	base_spline(const base_spline&);
	// destructor
	~base_spline() { delete[] knots; }
	// operator
	base_spline& operator=(const base_spline&);
	int ranks(int b=0) const {return rank-b;}
	int Imax() const { return imax; }
	int Icox() const { return icox; }
	int Rank() const { return rank; }
	int Jisu() const { return jisu; }
	int Shuki()const { return shuki;}
	int Maxq() const { return maxq; }
	int Kset(S x) const;
	S x_max() const {return knots[icox];}
	S x_min() const {return knots[jisu];}
	S operator[] (int i) const {return knots[i];}
	int basic(S, int, S*, int = 0) const;
	void   bibun_keisu(int, int, S**) const;
	void sekibun_keisu(int, int, S**) const;
	void gyoretu(S**, int, const S* = NULL, const int* = NULL) const;
	parray<S> basis(S, int = 0) const;
	parray<S> bases(S, int = 0) const;
/*******************************************************************************
	線形結合係数の適用 Σa[i]*b[i]; i = kset-Kai+1,...,kset;
	  d: d > 0 ? 微分階数 : d == 0 ? 補間値
	  t: データ点の座標値
	alp: 線形結合係数の配列
*******************************************************************************/
template <typename T> T apply(S t, T *alp, int d = 0) const
{
	int k = rank - d;
	S *b = T_ALLOC(S, icox);
	int kset = basic(t, icox, b, k);
	k = kset - k + 1;
	if (k < 0) k = 0;
	T sum = 0.0;
	for (int i = k; i <= kset; i++) sum += alp[i] * b[i];
	FREE(b);
	return sum;
}
/*******************************************************************************
	定積分の計算
*******************************************************************************/
template <typename T> T sekibun(S x, T *alp) const
{
	S  *q = knots;
	int K = rank, J = K - 1;
	int Imax = Kset(x);
	T sum = 0.0;
	for (int ql = J; ql <= Imax; ++ql) {
		S t = q[ql];
		parray<S> fl = bases(t, J);
		int ks = fl(0), kt = fl(1), kset = ks + J;
		S *ft = *fl, Qh = (ql < Imax ? q[ql+1] : x) - t;
		for (int i = ks; i <= kset; ++i) {
			S Qs = 0.0; int j = i - ks + K * kt;
			for (int r = K; r > 0; --r) {
				S Bik = ft[j -= kt];
				Qs = (Qs + Bik) * Qh / (S)r;
			}
			sum += alp[i] * Qs;
		}
	}
	return sum;
}
/*******************************************************************************
	高階積分の計算
*******************************************************************************/
template <typename T> T sekibun(S t, T *alp, int Jsk) const
{
	int ks, kset, k = rank + Jsk;
	S b[k];
	ks = deboor<S>(rank, knots, icox, t, b, kset, k);
	T **a = MALLOC(2, kset+1);
	for (int i = 0; i <= kset; ++i) a[Jsk%2][i] = alp[i];
	while (Jsk) {
		T *a0 = a[Jsk%2], *a1 = a[(Jsk+1)%2];
		for (int i = 0; i <= kset; ++i)
			a1[i] = (i > 0 ? a1[i-1] : 0.0) + a0[i] * (knots[i+(k-Jsk)] - knots[i]) / (k-Jsk);
		Jsk--;
	}
	T sek = 0.0;
	for (int i = 0; i <= kset; ++i)
		if (ks <= i && i <= kset) sek += a[0][i] * b[i-ks];
	FREE(a);
	return sek;
}
/*******************************************************************************
	多次元積分の計算
*******************************************************************************/
template <typename T> void sekibun(S x, T *alp, int n, T *p) const
{
	 S *q = knots;
	int K = rank, J = K - 1;
	int Imax = Kset(x);
	for (int ql = J; ql <= Imax; ++ql) {
		S t = q[ql];
		parray<S> fl = bases(t, J);
		int ks = fl(0), kt = fl(1), kset = ks + J;
		S *ft = *fl, Qh = (ql < Imax ? q[ql+1] : x) - t;
		for (int i = ks; i <= kset; ++i) {
			S Qs = 0.0; int j = i - ks + K * kt;
			for (int r = K; r > 0; --r) {
				S Bik = ft[j -= kt];
				Qs = (Qs + Bik) * Qh / (S)r;
			}
			for (int k = 0; k < n; ++k) p[k] += alp[n*i+k] * Qs;
		}
	}
}

};	// end of base_spline;

#endif
