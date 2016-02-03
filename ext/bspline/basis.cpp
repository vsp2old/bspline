#include <cstring>
#include <cstdlib>
#include <cstddef>
#include <cstdio>
#include <cassert>
#include "bspline_Config.h"
#include "basis/util.h"
#include "basis/basis.h"

static const double eps = 1.0e-10;

/*******************************************************************************
	de Boor-Cox の漸化式によるＢ−スプライン関数値の計算
	  k: Ｂスプラインの階数（j+1）
	  q: 節点の配列
	cox: Ｂスプラインの基底数
	  t: データ点の座標（q[j] <= t <= q[cox]）
	  d: d > 0 ? スプライン階数 : 微分階数
*******************************************************************************/
template <typename T>
int deboor(int k, T *q, int c, T t, T *b, int& kset, int d)
{
	int h = c-1, i = k-1, j, K = d > 0 ? d : k, L = k + d;
	T a, q1, q2, dt;
	
	while (i < h) { if(q[i] <= t && t < q[i+1]) break; i++; }
	kset = i; if (q[h] <= t && t < q[c] + eps) kset = h;
	b[0] = 1;
	for (j = 1; j < K; ++j) {
		if (L) --L;
		b[j] = 0;
		for (i = 0; i < j; ++i) {
			h = j - i;
			if (L) {
				q1 = q[kset + h] - t;
				q2 = t - q[kset - i];
				dt = q1 + q2;
				a = 0; if (dt > 0.0) a = b[h-1] / dt;
				b[ h ] += a * q2;	// 決定項
				b[h-1]  = a * q1;	// 準備項
			} else {
				dt = q[kset + h] - q[kset - i];
				a = 0; if (dt > 0.0) a = b[h-1] * (T)j / dt;
				b[ h ] += a;	// 決定項
				b[h-1] = -a;	// 準備項
			}
		}
	}
	return kset - j + 1;
}
/*******************************************************************************
	de Boor-Cox の漸化式によるＢ−スプライン関数値の計算
	  K: Ｂスプラインの階数（j+1）
	  q: 節点の配列
	cox: Ｂスプラインの基底数
	  t: データ点の座標（q[j] <= t <= q[cox]）
	  J: 漸化式の適用回数
	  L: (0 >= L >= -j) ? 微分階数
*******************************************************************************/
template<typename T>
void deboor_cox(int K, T *q, int cox, T t, T *b, int &kset, int J, int L)
{
	int M = J + L;
	int H = M + K + L;
	if (J == 0 && M >= 0) {
		int j = K - 1, c = cox - 1;
		for(int i=j;i<c;++i)if(q[i]<=t && t<q[i+1]){kset=i;break;}
		if(q[c]<=t && t<=q[cox]+eps){kset=c;}
		b[0] = 1;
	} else if (M >= 0 || H < 0) {
		int j = J - (M < 0 ? M + 1 : 0);
		T *b0 = &b[j+1];
		deboor_cox(K, q, cox, t, b0, kset, J-1, L+1);
		b[j] = 0.0;
		for(int i = 0; i < j; ++i) {
			int h = j - i;
			T a, q1, q2, dt;
			if (L >= 0) {
				q1 = q[kset + h] - t;
				q2 = t - q[kset - i];
				dt = q1 + q2;
				a = 0; if (dt > 0.0) a = b0[h-1] / dt;
				b[ h ] += a * q2;
				b[h-1]  = a * q1;
			} else {
				dt = q[kset + h] - q[kset - i];
				a = 0; if (dt > 0.0) a = b0[h-1] * (T)j / dt;
				b[ h ] += a;
				b[h-1] = -a;
			}
		}
	}
}
/*******************************************************************************
	シェーンバーグ・ホイットニー条件を満たす節点の設定
    q0, ..., qj, qk, ..., qi, ..., qc-1, qc, ..., qc+j
    x0, ..., x0, | , ..., | , ..., | , xn-1, ..., xn-1
            (x0+xk)/2 (xi-k+xi)/2 (xn-k-1+xn-1)/2
	 q : 節点の配列
	 m : 節点配列の長さ (c + k)
	 c : スプライン基底の個数
	 n : データ点の個数
	 k : スプライン階数 (k <= n)
	 d : 境界条件数 (d < k)
	 s : 周期
	 x : データ点の配列（x == NULL: パラメトリックスプラインを示す）
	ms : 
*******************************************************************************/
template<typename T>
void setten(T *q, int m, int c, int n, int k, int d = 0, int s = 0, const T *x = NULL, const int *ms = NULL)
{
	if (s == 0)
		for (int i = 0; i < m; i++)
			if (x)		/* シェーンバーグ・ホイットニー条件 */
				q[i] = i < k ? x[0] : i < c ? (x[i-k] + x[i-d]) / 2.0 : x[n-1];
			else		/* パラメトリック */
				q[i] = i < k ? 0 : i < c ? i - T(k+d) / 2.0 : n-1;
	else
	if (s > 0)			/* 周期境界条件による節点の設定 */
		for (int i = 0; i < m; i++) {
			 int j = i - k + 1;
			if (x)
				q[i] = j < 0 ? x[0] + x[j+s] - x[s] : j <= s ? x[j] : x[s] + x[j-s] - x[0];
			else
				q[i] = j;
		}
	else
	if (ms) {
		int k0 = k, n1, k1;
		while (s++) {
			n1 = *ms++; k1 = *ms++;
			setten(q, k0+n1-1, k0+n1-1, n1, k0, 0, 0, x);
			q += k0+n1-1; m -= k0+n1-1; c -= n1-1; n -= n1-1; k0 = k1; x += n1-1;
		}	setten(q, m, c, n, k0, 0, 0, x);
	}
}
/*******************************************************************************
	初期化（d == 0 ? 境界条件のないスプライン補間）
	n: データ点の個数
	x: データ点の座標	x[0,...,n-1]
	j: スプラインの次数（スプライン階数 K = j + 1）
	d: d >= 0 ? データ点に追加される境界条件数 : 周期境界条件（d = T - n - 1)
	T: 周期（境界条件数 = d + K）
	d > j : 多重節点の付加、m = {n1, k1, ... }
*******************************************************************************/
template <typename S>
base_spline<S>::base_spline(int n, const S *x, int j, int d, int *m)
{
	int s = 0;
	if (d > j) { s = j - d; d = 0; }
//	周期境界条件フラグ
	shuki = d < 0 ? 1 : 0;
//	スプライン階数
	rank = j + 1;
//	スプライン次数
	jisu = j;
//	周期（データ点の個数）
	imax = n + d + shuki;
//	Ｂスプライン関数基底の個数
	icox = imax + jisu * shuki;
//	節点の個数
	maxq = icox + rank;
//	節点の配列
	knots = new S[maxq];
//	節点の配列を初期化
	setten(knots, maxq, icox, n, rank, d, s + imax * shuki, x, m);
}

//	コピーコンストラクター
template <typename S>
base_spline<S>::base_spline(const base_spline& b)
 : imax(b.imax), icox(b.icox), rank(b.rank), jisu(b.jisu), shuki(b.shuki), maxq(b.maxq)
{
	knots = new S[maxq];
	memcpy(knots, b.knots, maxq * sizeof(S));
}
//	代入オペレーター
template <typename S>
base_spline<S>& base_spline<S>::operator=(const base_spline& b)
{
	if (this != &b) {
		delete[] knots;
		imax = b.imax; icox = b.icox; rank = b.rank;
		jisu = b.jisu; shuki = b.shuki; maxq = b.maxq;
		knots = new S[maxq];
		memcpy(knots, b.knots, maxq * sizeof(S));
	}
	return *this;
}
/*******************************************************************************
    Ｂスプライン関数値の計算 Ｂi,K(t): i = kset-jisu,...,kset
	  t: データ点の座標値
	c,b: 結果を入れる配列（c:サイズ，b:配列へのポインタ）
	  k: k > 0 ? Ｂi,k(x) : k == 0 ? Ｂi,rank(x) : k < 0 ? Ｄ-kＢi,rank(x)
*******************************************************************************/
template <typename S>
int base_spline<S>::basic(S t, int c, S* b, int k) const
{
	int kset = 0, j = k > 0 ? k : rank;
	S bc[j];
	int ks = deboor<S>(rank, knots, icox, t, bc, kset, k);
	if (c <= j)
		for (int i=0;i<c;++i) b[i] = bc[ks < 0 ? i - ks : i];
	else
		for (int i=0;i<c;++i) b[i] = (ks <= i) && (i <= kset) ? bc[i - ks] : 0;
	return kset;
}
/*******************************************************************************
    係数行列の構成
*******************************************************************************/
static void shuki_spline(int k, double *b, int n, double *bc)
{
	int j = k - 1;
	int m = k / 2;	/* m = Mjis = (Jisu+1)/2 = kai/2 */
	int m_1 = j - m;
	int i, o = n - m + 1;
	for (i = 0; i < m; i++) bc[i] = b[i+m_1] + b[i+m_1+n];
	for (i = m; i < o; i++) bc[i] = b[i+m_1];
	for (i = o; i < n; i++) bc[i] = b[i+m_1] + b[i+m_1-n];
}

template <typename S>
void base_spline<S>::gyoretu(S **bd, int n, const S *x, const int *d) const
{
	if (shuki == 0)
		for (int j = 0; j < n; j++)
			basic((x ? x[j] : j), imax, bd[j], (d ? -d[j] : 0));
	else {
		S *b = new S[icox];	/* 周期スプライン係数行列の計算 */
		for (int j = 0; j < n; j++) {
			basic((x ? x[j] : j), icox, b, (d ? -d[j] : 0));
			shuki_spline(rank, b, imax, bd[j]);
		}
		delete[] b;
	}
}

template <typename S>
parray<S> base_spline<S>::basis(S t, int jbn) const
{
	S data[rank];
	int kset, offset = deboor<S>(rank, knots, icox, t, data, kset, -jbn);
	parray<S> result(offset, rank, 1);
	S *p = *result; for (int i = 0; i < rank; ++i) p[i] = data[i];
	return result;
}
	
template <typename S>
parray<S> base_spline<S>::bases(S t, int jbn) const
{
	int N = 1; for (int i = 2; i <= rank; i++) N += i;
	S *bv = new S[N];
	int kset; deboor_cox<S>(rank, knots, icox, t, bv, kset, jisu);
	int offset = kset - jisu, k = jbn + 1;
	parray<S> data(offset, rank, k);
	S *mp = *data;
	for (int j = 0; j <= jbn; ++j) {
		if (j > 0) deboor_cox<S>(rank, knots, icox, t, bv, kset, rank-j, -rank);
		for(int i = 0; i < rank; ++i) *mp++ = bv[i];
	}
	delete[] bv;
	return data;
}	
/*******************************************************************************
	微分係数の計算
*******************************************************************************/
template <typename S> void base_spline<S>::bibun_keisu(int Jbn, int n, S** ad) const
{
	S zk;
	S *a0, *a1, *q = knots;

	for (int j = 0; j < Jbn; j++) {
		a0 = ad[j];
		a1 = ad[j+1];
		a1[0] = 0.0;
		for (int i = 1; i < n; i++) {
			if ((zk = q[i+Jbn-j] - q[i]) != 0.0)
				a1[i] = (Jbn-j)*(a0[i] - a0[i-1])/zk;
			else
				a1[i] = 0.0;
		}
	}
}
/*******************************************************************************
	積分係数の計算
*******************************************************************************/
static double qvalue(int m, double *q, int i, int j, int c)
{
	int k = i - c + j;
	return ((k < m) ? q[k] : qvalue(m, q, k, j, c)) + q[c] - q[j];
}

#define QVALUE(i) ((i)<m ? q[i] : shuki ? qvalue(m,q,i,jisu,icox) : q[m-1])

template <typename S> void base_spline<S>::sekibun_keisu(int Jsk, int n, S** ai) const
{
	int ki = rank, m = Maxq();
	S *a0, *a1, *q = knots;

	for (int js = 1; js <= Jsk; js++) {
		a0 = ai[(js-1)%2];
		a1 = ai[js%2];
		a1[0] = a0[0] * (QVALUE(ki)-q[0])/ki;
		for (int j = 1; j < n; j++)
			a1[j] = a1[j-1] + a0[j]*(QVALUE(j+ki)-q[j])/ki;
		ki++;
	}
}

template <typename S> int base_spline<S>::Kset(S x) const
{
	int L = 0; S *q = knots;
	for(int i = jisu; i < icox; ++i)
		if (q[i] <= x && x <  q[i+1]){ L = i; break;}
	if(q[icox-1] <= x && x <= q[icox]+eps){ L = icox-1;}
	return L;
}




template class base_spline<double>;

