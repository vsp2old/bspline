#ifndef __GSPLINE_H__
#define __GSPLINE_H__

#ifdef ARRAY_OF_GSL

template <typename T> class block
{
	 int * cp;
  protected:
	 T * data;
  public:
	block() : cp(NULL),data(NULL) {}
	virtual ~block() {if(cp && --*cp == 0){delete cp;delete[] data;}}
	block(const block& ba) : cp(ba.cp),data(ba.data) {++*cp;}
	block& operator=(const block& ba)
	{
		if(this != &ba){
			if(cp && --*cp == 0){delete cp;delete[] data;}
			cp = ba.cp; data = ba.data; ++*cp;
		}
		return *this;
	}
	block(int l) : cp(new int(1)),data(new T[l]) {}
	operator T * () { return data; }
};

class narray : public block<int>
{
  public:
	int size;

	narray(){}
	narray(int i) : block<int>(i),size(i) {}
	narray(int i, int *a) : block<int>(i),size(i)
	{
		memcpy(data,a,i*sizeof(int));
	}
	int & operator[](int i) { return data[i]; }
	operator int * () { return data; }
};
#endif

typedef double** GARRAY;
typedef double * VARRAY;

class Gspline
{
	int Kosu, Jisu, Kai, Mjis, Imax; double G;
	double *dpx, *dpy, *dpw, *alp, *gfp;
  public:
	// constructor
	Gspline(): dpx(NULL) {}
	Gspline(int, double *, double *, double *, int, double = 0.0);
	Gspline(const Gspline& gs);
	// destructor
	~Gspline() { FREE(dpx); }
	// operator
	Gspline& operator=(const Gspline& gs);
	double operator[](double t) const { return (*this)(t); }
	double operator()(double t, int b = 0) const;
	double   getx  (int i)     { return dpx[i]; }
	double   gety  (int i)     { return dpy[i]; }
	double   getw  (int i)     { return dpw[i]; }
	   int   size  (int i = 0) { return Kosu-i; }
	double & domain(int i = 0) { return dpx[i]; }
	double & values(int i = 0) { return dpy[i]; }
	double & weight(int i = 0) { return dpw[i]; }
	double  sekibun(double);
	double  sekibun(double, double);
  private:
	void sabunsho  (GARRAY,         GARRAY, VARRAY) const;
	void heikatu   (GARRAY                        ) const;
	void keisu_smo (GARRAY, VARRAY                ) const;
	void zyoyo     (VARRAY, VARRAY, GARRAY, VARRAY) const;
	void takoushiki(double, int,    VARRAY        ) const;
};

#endif
