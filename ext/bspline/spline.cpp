#include <math.h>
#include <stdarg.h>
#include <string.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <assert.h>
#include "bspline_Config.h"
#include "basis/util.h"
#include "basis/basis.h"
#include "basis/spline.h"

/*******************************************************************************
	Ｂスプラインの初期化
	 n: データ点の数
	Xp: データ点（＋条件点）のＸ座標の配列 : Xp = X[n+d];
	Yp: データ点（＋条件点）のＹ座標の配列 : Yp = Y[n+d];
	 j: 次数
	 d: 条件数(j >= d > 0)
	Dp: 条件点の微分階数
*******************************************************************************/
spline::spline(int n, double *Xp, double *Yp, int j, int d, int *Dp)
: base_spline(n, Xp, j, d, Dp)
{
	alp = VALLOC(icox);
	int nc = imax < n ? imax : n;
	marray<double> bd(imax, imax);
	gyoretu(bd, nc, Xp);	/* 係数行列を作成 */
	if (j >= d && d > 0) for (int i = 0; i < d; i++) {
		int l = i + nc;
		if (Dp == NULL || Dp[i] == 0)
			basic(Xp[l], imax, bd[l]);
		else
			basic(Xp[l], imax, bd[l], -Dp[i]);
	} /* d == 0 ? Spline without the boundary condition : d < 0 ? Period Spline */
	marray_lu_view<double> Lu(imax, bd[0]);
	Lu.solve(Yp, alp);	/* 線形結合係数を作成 */
}
//	コピーコンストラクター
spline::spline(const spline& s) : base_spline(s)
{
	alp = VALLOC(icox);
	for (int i = 0; i < icox; i++) alp[i] = s.alp[i];
}
//	代入オペレーター
spline& spline::operator=(const spline& s)
{
	if (this != &s) {
		base_spline::operator=(s);
		FREE(alp);
		alp = VALLOC(icox);
		for (int i = 0; i < icox; i++) alp[i] = s.alp[i];
	}
	return *this;
}
//	Ｂスプライン関数の微分
double spline::bibun(double t, int Jbn) const
{
	if (t < x_min() || t > x_max() || Jbn > Jisu())
		throw "out of range, spline.bibun";
	varray<double> b(rank);
	int kset = basic(t, rank, b, -Jbn);
	int ks = kset - jisu;
	double sum = 0.0;
	for ( int i = 0; i < rank; ++i) sum += alp[ks + i] * b[i];
	return sum;
}
//	境界条件のないＢスプライン補間
Bspline::Bspline(int n, double(*XY)[2], int j, int d, int *m)
{
	double *Wx = new double[n];
	double *Wy = new double[n];
	for (int i = 0; i < n; i++) { Wx[i] = XY[i][0]; Wy[i] = XY[i][1]; }
	sp = new spline(n, Wx, Wy, j, d, m);
	alp = MALLOC(Rank(), Icox());
	alp[0] = *sp;					/* 線形結合係数 */
	bibun_keisu(j, Icox(), alp);	/* 微分線形結合係数 */
	delete[] Wx;
	delete[] Wy;
}
// 	関数値境界条件によるＢスプライン補間
Bspline::Bspline(int n, double(*XY)[2], int j, int d, double(*Cp)[2])
{
	double *Wx = new double[n+d];
	double *Wy = new double[n+d];
	for (int i = 0; i < n; i++) { Wx[i]  =  XY[i][0]; Wy[i]  =  XY[i][1]; }
	for (int i = 0; i < d; i++) { Wx[n+i] = Cp[i][0]; Wy[n+i] = Cp[i][1]; }
	sp = new spline(n, Wx, Wy, j, d);
	alp = MALLOC(Rank(), Icox());
	alp[0] = *sp;					/* 線形結合係数 */
	bibun_keisu(j, Icox(), alp);	/* 微分線形結合係数 */
	delete[] Wx;
	delete[] Wy;
}
// 	微分境界条件によるＢスプライン補間
Bspline::Bspline(int n, double(*XY)[2], int j, int d, int *Dp, double(*Cp)[2])
{
	double *Wx = new double[n+d];
	double *Wy = new double[n+d];
	for (int i = 0; i < n; i++) { Wx[i]  =  XY[i][0]; Wy[i]  =  XY[i][1]; }
	for (int i = 0; i < d; i++) { Wx[n+i] = Cp[i][0]; Wy[n+i] = Cp[i][1]; }
	sp = new spline(n, Wx, Wy, j, d, Dp);
	alp = MALLOC(Rank(), Icox());
	alp[0] = *sp;					/* 線形結合係数 */
	bibun_keisu(j, Icox(), alp);	/* 微分線形結合係数 */
	delete[] Wx;
	delete[] Wy;
}
//	周期境界条件によるスプライン補間
Bspline::Bspline(int n, double(*XY)[2], int j, int t)
{
	double *Wx = new double[n];
	double *Wy = new double[n];
	for (int i = 0; i < n; i++) { Wx[i] = XY[i][0]; Wy[i] = XY[i][1]; }
	sp = new spline(n, Wx, Wy, j, t-n-1);
	alp = MALLOC(Rank(), Icox());	/* Icox（基底数） > Imax（周期） */
	varray_view<double> alv(alp[0], Icox()), spv((double*)(*sp), Imax());
	alv = spv;
	bibun_keisu(j, Icox(), alp);	/* 微分線形結合係数 */
	delete[] Wx;
	delete[] Wy;
}
//	コピーコンストラクター
Bspline::Bspline(const Bspline& b)
{
	sp = new spline(*b.sp);
	int r = Rank(), c = Icox();
	alp = MALLOC(r, c);
	for (int i = 0; i < r; i++)
		for (int j = 0; j < c; j++) alp[i][j] = b.alp[i][j];
}
//	代入オペレーター
Bspline& Bspline::operator=(const Bspline& b)
{
	if (this != &b) {
		delete sp; FREE(alp);
		sp = new spline(*b.sp);
		int r = Rank(), c = Icox();
		alp = MALLOC(r, c);
		for (int i = 0; i < r; i++)
			for (int j = 0; j < c; j++) alp[i][j] = b.alp[i][j];
	}
	return *this;
}
/*******************************************************************************
//	関数値の補間計算
//		B : spline Object
//		x : 区間の座標
//		Xp: プロットの位置
//		Yp: プロットの値
//		Dp: 補間の区分数
*******************************************************************************/
int plots(const Bspline& B, int n, double *x, double **Xp, double **Yp, int Dp, int b)
{

	int Npx, L, N0, Nd, Np;
	double Tp, DT;

	Npx = (n - 1) * Dp + 1;
	*Xp = VALLOC(Npx);
	*Yp = VALLOC(Npx);
	for (L = 1; L < n; L++) {
		DT = (x[L] - x[L-1]) / Dp;
		N0 = (L == 1 ? 0 : 1);
		for (Nd = N0; Nd <= Dp; Nd++) {
			Np = Nd + (L-1) * Dp;
			Tp = x[L-1] + Nd * DT;
			(*Xp)[Np] = Tp;
			(*Yp)[Np] = B(Tp, b);
		}
	}
	return Npx;
}

