#include <string.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include "basis/util.h"
#include "basis/gspline.h"

/*******************************************************************************
    係数 Bi の計算
*******************************************************************************/
void Gspline::heikatu(GARRAY cb) const
{
	int i, j, k, jst = 0, jsp = 0;

	for (i = 0; i < Kosu; i++)
	{
		for (j = 0; j < Kosu; j++) cb[i][j] = 0.0;
		if (0 <= i && i <= Mjis - 1)
		{
			jst = 0;
			jsp = i;
		}
		if (Mjis <= i && i <= Kosu - Mjis - 1)
		{
			jst = i - Mjis;
			jsp = i;
		}
		if (Kosu - Mjis <= i && i <= Kosu - 1)
		{
			jst = i - Mjis;
			jsp = Kosu - Mjis - 1;
		}
		for (j = jst; j <= jsp; j++)
		{
			double psk = 1.0;
			for (k = j; k <= j + Mjis; k++)
				if (i != k) psk *= (dpx[i] - dpx[k]);
			cb[i][j] = 1.0 / psk;
		}
	}
}
/*******************************************************************************
    xおよびtについての差分商による係数行列の構成
*******************************************************************************/
void Gspline::sabunsho(GARRAY Bd, GARRAY pnw, double *sp) const
{
	int i, j, k;

	double gw = 1.0; for (j = 1; j <= Jisu; j++) gw *= j;
	double fu = 1.0; for (j = 1; j <= Mjis; j++) fu = -fu;
	gw *= fu * G;
	for (i = 0; i < Kosu; i++)
		for (j = 0; j < Kosu; j++)
			Bd[i][j] = (i > j) ? pow((dpx[i] - dpx[j]), Jisu) : 0.0;
	for (k = 1; k <= Mjis; k++)
		for (i = 0; i < Kosu; i++)
			for (j = 0; j+k < Kosu; j++)
				Bd[i][j] = (Bd[i][j+1] - Bd[i][j]) / (dpx[j+k] - dpx[j]);
	matrix cb(Kosu, Kosu);
	heikatu(cb);
	for (i = 0; i < Kosu; i++) {
		double gwi = (dpw[i] != 0) ? gw / dpw[i] : 1.0E+20;
		for (j = 0; j < Imax; j++)
			Bd[i][j]   += gwi * cb[i][j];
			Bd[i][Imax] = dpy[i];
	}	for (j = 0; j < Imax; j++) pnw[0][j] = Bd[0][j];
		sp[0] = dpy[0];
	for (k = 1; k <= Mjis; k++) {
		for (i = 0; i+k < Kosu; i++)
			for (j = 0; j <= Imax; j++)
				Bd[i][j] = (Bd[i+1][j] - Bd[i][j]) / (dpx[i+k] - dpx[i]);
		if (k == Mjis) break;
		for (j = 0; j < Imax; j++) pnw[k][j] = Bd[0][j];
		sp[k] = Bd[0][Imax];
	}
}
/*******************************************************************************
    K-1次の多項式の微分
*******************************************************************************/
void Gspline::takoushiki(double t, int Jbn, double* sbn) const
{
	int j, jj, i, m, ist, mmx;
	narray lm = NALLOC(Jisu);
//	set r to Jbn! = Jbn*(Jbn-1)*...*1
	double r = 1.0; for (j = 1; j <= Jbn; j++) r *= j;
	for (m = 0; m <= Mjis - 1; m++) {
		if      (m <  Jbn) sbn[m] = 0.0;
		else if (m == Jbn) sbn[m] = 1.0;
		else {
			mmx = m - Jbn;
			ist = 0;
			sbn[m] = 0.0;
			for (i = 0; i < mmx; i++) lm[i] = -1;
			for (jj = 1; jj <= 1000; jj++) {
				if (-1 < lm[ist] && lm[ist] <  (m - mmx + ist)) ist = mmx - 1;
				if (                lm[ist] == (m - mmx + ist)) ist = ist - 1;
				if (ist <= -1) break;
				lm[ist] = lm[ist] + 1;
				if (mmx - 1 > ist)
					for (i = ist; i < mmx - 1; i++) lm[i+1] = lm[i] + 1;
				double h = 1.0; for (i = 0; i < mmx; i++) h *= (t - dpx[lm[i]]);
				sbn[m] += h;
			}
		}
		sbn[m] *= r;
	}
	FREE(lm);
}
/*******************************************************************************
    B-スプラインの線形結合係数の計算
*******************************************************************************/
void Gspline::keisu_smo(GARRAY bd, double *eo) const
{
	matrix ld(Imax, Imax);
	for (int i = 0; i < Imax; i++)
	{
		memcpy(ld[i], bd[i], Imax * sizeof(double));
		eo[i] = bd[i][Imax];
	}
	int s; size_t p[Imax];
	lu_decomp<double>(ld, Imax, p, s);
	lu_subst<double>(ld, Imax, p, eo);
}
/*******************************************************************************
    剰余項の関数値計算
*******************************************************************************/
void Gspline::zyoyo(double* fp, double *eo, GARRAY pnw, double* sp) const
{
	for (int k = 0; k < Mjis; k++)
	{
		double app = 0.0;
		for (int j = 0; j < Imax; j++) app += eo[j] * pnw[k][j];
		fp[k] = sp[k] - app;
	}
}
/*******************************************************************************
    Ｇスプラインコンストラクタ
*******************************************************************************/
Gspline::Gspline(int n, double *x, double *y, double *w, int j, double g)
 : Kosu(n), Jisu(j), Kai(j+1)
{
	Mjis = Kai / 2;
	Imax = Kosu - Mjis;

	dpx = VALLOC(4 * n);
	dpy = dpx + n;
	dpw = dpy + n;
	alp = dpw + n;		//線形結合係数
	gfp = alp + Imax;	//剰余項
	memcpy(dpx, x, n * sizeof(double));
	memcpy(dpy, y, n * sizeof(double));
	memcpy(dpw, w, n * sizeof(double));
	G = g;
	matrix Bd (Kosu,   Kosu);
	matrix pnw(Mjis+1, Imax);
	vector sp (Mjis+1      );
	sabunsho (Bd,       pnw, sp);
	keisu_smo(Bd,  alp         );
	zyoyo    (gfp, alp, pnw, sp);
}
Gspline::Gspline(const Gspline& gs) : Kosu(gs.Kosu), Jisu(gs.Jisu), Kai(gs.Kai)
{
	Mjis = gs.Mjis; Imax = gs.Imax; G = gs.G;
	dpx = VALLOC(4 * Kosu);
	dpy = dpx + Kosu;
	dpw = dpy + Kosu;
	alp = dpw + Kosu;	//線形結合係数
	gfp = alp + Imax;	//剰余項
	memcpy(dpx, gs.dpx, 4 * Kosu * sizeof(double));
}
Gspline& Gspline::operator=(const Gspline& gs)
{
	if (this != &gs) {
		FREE(dpx);
		Kosu = gs.Kosu; Jisu = gs.Jisu; Kai = gs.Kai;
		Mjis = gs.Mjis; Imax = gs.Imax;   G = gs.G;
		dpx = VALLOC(4 * Kosu);
		dpy = dpx + Kosu;
		dpw = dpy + Kosu;
		alp = dpw + Kosu;	//線形結合係数
		gfp = alp + Imax;	//剰余項
		memcpy(dpx, gs.dpx, 4 * Kosu * sizeof(double));
	}
	return *this;
}
/*******************************************************************************
    平滑化スプライン関数値計算
*******************************************************************************/
double Gspline::operator()(double t, int Jbn) const
{
	vector ty(Kosu);
	int i, j, k, imx, jmx;
//	set yks to J*(J-1)*...*(J-Jbn+1)
	double yks = 1.0;
	if (Jbn > 0) for (j = 0; j < Jbn; j++) yks *= (Jisu - j);
//	initialize imx : x[imx] <= t < x[imx+1]
	imx = Kosu - 1; while (t < dpx[imx] && imx > 0) imx--;
//	initialize ty[]
	for (i = 0; i < Kosu; i++)
		if (i <= imx) ty[i] = yks * pow(t - dpx[i], Jisu - Jbn);
		else          ty[i] = 0.0;
//	calculate  ty[]
	for (k = 1; k <= Mjis; k++) {
		jmx = Kosu - k;
		for (j = 0; j < jmx; j++) ty[j] = (ty[j+1] - ty[j]) / (dpx[j+k] - dpx[j]);
	}
	double y1 = 0.0;
		for (j = 0; j < Imax; j++)  y1 += alp[j] * ty[j];
	double y2 = 0.0;
	if (Jbn > 0)
	{
		vector sbn(Jisu);	//差分商
		takoushiki(t, Jbn, sbn);
		for (k = 0; k <= Mjis; k++) y2 += gfp[k] * sbn[k];
	}
	else
	{
		double sek = 1.0;
		for (k = 0; k <= Mjis; k++) {
			if (k > 0) sek *= (t - dpx[k-1]);
			y2 += gfp[k] * sek;
		}
	}
	return y1 + y2;
}

static void intde(Gspline *f, double a, double b, double eps, double *i, double *err)
{
    /* ---- adjustable parameter ---- */
    int mmax = 256;
    double efs = 0.1, hoff = 8.5;
    /* ------------------------------ */
    int m;
    double pi2, epsln, epsh, h0, ehp, ehm, epst, ba, ir, h, iback, 
        irback, t, ep, em, xw, xa, wg, fa, fb, errt, errh, errd;
    
    pi2 = 2 * atan(1.0);
    epsln = 1 - log(efs * eps);
    epsh = sqrt(efs * eps);
    h0 = hoff / epsln;
    ehp = exp(h0);
    ehm = 1 / ehp;
    epst = exp(-ehm * epsln);
    ba = b - a;
    ir = (*f)((a + b) * 0.5) * (ba * 0.25);
    *i = ir * (2 * pi2);
    *err = fabs(*i) * epst;
    h = 2 * h0;
    m = 1;
    do {
        iback = *i;
        irback = ir;
        t = h * 0.5;
        do {
            em = exp(t);
            ep = pi2 * em;
            em = pi2 / em;
            do {
                xw = 1 / (1 + exp(ep - em));
                xa = ba * xw;
                wg = xa * (1 - xw);
                fa = (*f)(a + xa) * wg;
                fb = (*f)(b - xa) * wg;
                ir += fa + fb;
                *i += (fa + fb) * (ep + em);
                errt = (fabs(fa) + fabs(fb)) * (ep + em);
                if (m == 1) *err += errt * epst;
                ep *= ehp;
                em *= ehm;
            } while (errt > *err || xw > epsh);
            t += h;
        } while (t < h0);
        if (m == 1) {
            errh = (*err / epst) * epsh * h0;
            errd = 1 + 2 * errh;
        } else {
            errd = h * (fabs(*i - 2 * iback) + 4 * fabs(ir - 2 * irback));
        }
        h *= 0.5;
        m *= 2;
    } while (errd > errh && m < mmax);
    *i *= h;
    if (errd > errh) {
        *err = -errd * m;
    } else {
        *err = errh * epsh * m / (2 * efs);
    }
}

double Gspline::sekibun(double t)
{
	double t0 = dpx[0];
	double err, eps = 1e-10, i = 0.0;
	if (t > t0) intde(this, t0, t, eps, &i, &err);
	return i;
}

double Gspline::sekibun(double ta, double tb)
{
	double t0 = dpx[0];
	if (ta <= t0) ta = t0;
	double err, eps = 1e-10, i = 0.0;
	if (tb > ta) intde(this, ta, tb, eps, &i, &err);
	return i;
}
