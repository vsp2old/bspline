#ifndef _SPLINE_H_INCLUDED_
#define _SPLINE_H_INCLUDED_

//
//	Ｂスプライン関数
//
class spline : public base_spline<>
{
	double *alp;
  public:
	// constructor
	spline() : alp(NULL) {}
	spline(int n, double *x, double *y, int j, int d = 0, int *dp = NULL);
	spline(const spline& s);
	// destructor
	~spline() { FREE(alp); }
	// operator
	spline& operator=(const spline& s);
	double operator[](double t) const {
		if (t < x_min() || t > x_max())
			throw "out of range, spline[]";
		return apply(t, alp);
	}
	double& operator[](int i) const {
		if (i < 0 || i >= Icox())
			throw "out of range, spline.alp[]";
		return alp[i];
	}
	double operator()(double t, double *a, int b=0) const {
		if (t < x_min() || t > x_max() || b > Jisu())
			throw "out of range, spline()";
		return b < 0 ? sekibun(t, -b) : apply(t, a, b);
	}
	operator double* () const {return alp;}
	double bibun(double t, int d = 1) const;
	double sekibun(double t, int s) const { return base_spline::sekibun(t, alp, s); }
	double sekibun(double x) const { return base_spline::sekibun(x, alp); }
};
//
//	Ｂスプライン補間
//
class Bspline
{
	spline *sp;
	double **alp;
  public:
	// constructor
	Bspline() : sp(NULL), alp(NULL) {}
	// bspline : 境界条件のないＢスプライン補間
	Bspline(int n, double(*XY)[2], int j, int d, int *m);
	// cspline : 関数値境界条件によるＢスプライン補間
	Bspline(int n, double(*XY)[2], int j, int d, double(*Cp)[2]);
	// dspline : 微分境界条件によるＢスプライン補間
	Bspline(int n, double(*XY)[2], int j, int d, int *Dp, double(*Cp)[2]);
	// tspline : 周期境界条件によるスプライン補間
	Bspline(int n, double(*XY)[2], int j, int t);
	// copy constructor
	Bspline(const Bspline& b);
	// destructor
	~Bspline() { delete sp; FREE(alp); }
	// operator
	Bspline& operator=(const Bspline& b);
	int Imax() const { return sp->Imax(); }
	int Icox() const { return sp->Icox(); }
	int Maxq() const { return sp->Maxq(); }
	int Rank() const { return sp->Rank(); }
	int Jisu() const { return sp->Jisu(); }
	int Shuki()const { return sp->Shuki();}
	double operator[](double t) const { return (*sp)[t]; }
	double operator()(double t, int b=0) const
	{
		if (b > Jisu())
			throw "out of range, Bspline()";
		return (*sp)(t, (b >= 0) ? alp[b] : NULL, b);
	}
	double x_max() const { return sp->x_max(); }
	double x_min() const { return sp->x_min(); }
	void bibun_keisu(int j, int i, double **a) const { sp->bibun_keisu(j, i, a); }
	double bibun(double t, int jbn = 1)  const { return sp->bibun(t, jbn); }
	double sekibun(double t, int jsk)    const { return sp->sekibun(t, jsk); }
	double sekibun(double ta, double tb) const { return sp->sekibun(tb, 1) - sp->sekibun(ta, 1); }
	double sekibun(double x) const { return sp->sekibun(x); }
};

//
//	misc
//
int plots(const Bspline&, int, double *, double **, double **, int, int = 0);

#endif
