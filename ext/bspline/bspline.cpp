/* Extention for ruby */

#include "ruby.h"

#ifdef __cplusplus
# define VALUEFUNC(f) ((VALUE (*)(...))f)
# define VOIDFUNC(f) ((void (*)(...))f)
#else
# define VALUEFUNC(f) (f)
# define VOIDFUNC(f) (f)
#endif

static VALUE mBspline;

#include "bspline.h"

static VALUE cBspline;

static void free_Bspline(Bspline *arg0) { delete arg0; }

#define Wrap_Bspline(klass, ptr) (\
	(ptr) ? Data_Wrap_Struct(klass, 0, VOIDFUNC(free_Bspline), ptr) : Qnil )

#define Get_Bspline(val, ptr) {\
	if (NIL_P(val)) ptr = NULL;\
	else {\
	    if (!rb_obj_is_kind_of(val, cBspline))\
	        rb_raise(rb_eTypeError, "wrong argument type (expected Bspline)");\
	    Data_Get_Struct(val, Bspline, ptr);\
	    if (!ptr) rb_raise(rb_eRuntimeError, "This Bspline already released");\
	}\
}

//	Bspline#new()
static VALUE
_wrap_new_Bspline(int argc, VALUE *argv, VALUE klass)
{
	VALUE self = Data_Wrap_Struct(klass, 0, VOIDFUNC(free_Bspline), 0);
	rb_obj_call_init(self, argc, argv);
	return self;
}
/*******************************************************************************
	境界条件のないＢスプライン補間
	Bspline.new([[X0,Y0], ..., [Xn,Yn]], J)
	  X0, ..., Xn:点のX座標
	  Y0, ..., Yn:点のY座標
	  J:スプライン次数（階数K=J+1）
	Result:Ｂスプラインオブジェクト
*******************************************************************************/
//	Bspline#initialize
static VALUE
_wrap_init_Bspline(VALUE self, VALUE varg, VALUE vargj)
{
	int argn, argj;
	double *xy;
	VALUE vx, vy, va;
	
	Check_Type(varg,T_ARRAY);
	argn = RARRAY_LEN(varg);
	argj = NUM2INT(vargj);
	if (argn < argj)
		rb_raise(rb_eArgError, "This Array is too small");
	xy = ALLOC_N(double, argn * 2);
	for (int i = 0; i < argn; i++) {
		va = RARRAY_PTR(varg)[i];
		Check_Type(va,T_ARRAY);
		if (RARRAY_LEN(va) < 2)
			rb_raise(rb_eArgError, "Need both x,y value");
		vx = RARRAY_PTR(va)[0];
		vy = RARRAY_PTR(va)[1];
		xy[i*2]   = NUM2DBL(vx);
		xy[i*2+1] = NUM2DBL(vy);
	}
	try {
		int M;
		DATA_PTR(self) = new Bspline(argn, (double(*)[2])xy, argj, 0, &M);
	} catch (const char *c) {
		rb_raise(rb_eRuntimeError, "%s", c);
	}
	free(xy);
	return self;
}
/*
// 関数値境界条件によるＢスプライン補間
*/
static VALUE
_wrap_init_Cspline(VALUE self, VALUE varg, VALUE vargj, VALUE vargc) {
	int argn, argj;
	double *xy, *c;
	VALUE va, vc;

	Check_Type(varg,T_ARRAY);
	argn = RARRAY_LEN(varg);
	argj = NUM2INT(vargj);
	xy = ALLOC_N(double, argn * 2);
	c = ALLOC_N(double, (argj-1)*2);
	for (int i = 0; i < argn; i++) {
		va = RARRAY_PTR(varg)[i];
		Check_Type(va,T_ARRAY);
		xy[i*2]   = NUM2DBL(RARRAY_PTR(va)[0]);
		xy[i*2+1] = NUM2DBL(RARRAY_PTR(va)[1]);
	}
	Check_Type(vargc,T_ARRAY);
	for (int i = 0; i < argj-1; i++) {
		vc = RARRAY_PTR(vargc)[i];
		Check_Type(vc,T_ARRAY);
		c[i*2]   = NUM2DBL(RARRAY_PTR(vc)[0]);
		c[i*2+1] = NUM2DBL(RARRAY_PTR(vc)[1]);
	}
	try {
	DATA_PTR(self) = new Bspline(argn, (double(*)[2])xy, argj, argj-1, (double(*)[2])c);
	} catch (const char *c) {
		rb_raise(rb_eRuntimeError, "%s", c);
	}
	free(c);
	free(xy);
	return self;
}
/*
// 微分境界条件によるＢスプライン補間
*/
static VALUE
_wrap_init_Dspline(VALUE self, VALUE varg, VALUE vargj, VALUE vargc, VALUE vargd)
{
	int argn, argj;
	double *xy, *c; int *d;
	VALUE va, vc, vd;

	Check_Type(varg,T_ARRAY);
	argn = RARRAY_LEN(varg);
	argj = NUM2INT(vargj);
	xy = ALLOC_N(double, argn * 2);
	c = ALLOC_N(double, (argj-1)*2);
	d = ALLOC_N(int, argj-1);
	for (int i = 0; i < argn; i++) {
		va = RARRAY_PTR(varg)[i];
		Check_Type(va,T_ARRAY);
		xy[i*2]   = NUM2DBL(RARRAY_PTR(va)[0]);
		xy[i*2+1] = NUM2DBL(RARRAY_PTR(va)[1]);
	}
	Check_Type(vargc,T_ARRAY);
	Check_Type(vargd,T_ARRAY);
	for (int i = 0; i < argj-1; i++) {
		vc = RARRAY_PTR(vargc)[i];
		Check_Type(vc,T_ARRAY);
		vd = RARRAY_PTR(vargd)[i];
		c[i*2]   = NUM2DBL(RARRAY_PTR(vc)[0]);
		c[i*2+1] = NUM2DBL(RARRAY_PTR(vc)[1]);
		d[i] = NUM2INT(vd);
	}
	try {
	DATA_PTR(self) = new Bspline(argn, (double(*)[2])xy, argj, argj-1, d, (double(*)[2])c);
	} catch (const char *c) {
		rb_raise(rb_eRuntimeError, "%s", c);
	}
	free(d);
	free(c);
	free(xy);
	return self;
}
static VALUE
_wrap_init_Mspline(int argc, VALUE *argv, VALUE self)
{
	VALUE varg0, varg1, varg2, varg3;
	if (argc >= 1) { varg0 = argv[0];
	if (argc >= 2) { varg1 = argv[1];
	if (argc >= 3) { varg2 = argv[2];
	if (argc >= 4) { varg3 = argv[3];
	if (argc >= 5) { rb_raise(rb_eArgError, "This is too many args");
	} else return _wrap_init_Dspline(self, varg0, varg1, varg2, varg3);
	} else return _wrap_init_Cspline(self, varg0, varg1, varg2);
	} else return _wrap_init_Bspline(self, varg0, varg1);
	}
	}	rb_raise(rb_eArgError, "This is few args");
	return Qnil;
}
static VALUE
_wrap_ary_flatten(VALUE arg)
{
	VALUE arg1;
	VALUE result = rb_ary_new();
	if (TYPE(arg) == T_ARRAY) {
		if(RARRAY_LEN(arg) == 0)
			return result;
		else {
			arg1 = rb_ary_pop(arg);
			return rb_ary_concat(_wrap_ary_flatten(arg), _wrap_ary_flatten(arg1));
		}
	} else 
		return rb_ary_push(result, arg);
}
/*
	関数値

	Bspline#[x[, ...]]
	  x:点のX座標（数値，または数値の配列）

	Result:補間値
*/
static VALUE
_wrap_Bspline_bracket(int argc, VALUE *argv, VALUE self)
{
	Bspline *bsp;
	double arg;
	VALUE val, vargs = rb_ary_new();
	VALUE vresult = rb_ary_new();
	Get_Bspline(self, bsp);
	for (int i = 0; i < argc; i++) {
		if (TYPE(argv[i]) == T_ARRAY) {
			val = rb_obj_dup(argv[i]);
			vargs = rb_ary_concat(vargs, _wrap_ary_flatten(val));
		}
		else
			rb_ary_push(vargs, argv[i]);
	}
	argc = RARRAY_LEN(vargs);
	try {
		if (argc == 1) {
			arg = NUM2DBL(RARRAY_PTR(vargs)[0]);
			vresult = rb_float_new((*bsp)[arg]);
		}
		else
			for (int i = 0; i < argc; i++) {
				arg = NUM2DBL(RARRAY_PTR(vargs)[i]);
				val = rb_float_new((*bsp)[arg]);
				rb_ary_push(vresult, val);
			}
	} catch (const char *c) {
		rb_raise(rb_eRuntimeError, "%s", c);
	}
	return vresult;
}
/*
	関数値
	Bspline#value(x, d = 0)
	  x:点のX座標
	  d:微分階数（省略可）
	Result:補間値，または微分値
*/
static VALUE
_wrap_Bspline_value(int argc, VALUE *argv, VALUE self)
{
	VALUE varg1, varg2 ;
	Bspline *arg0 ;
	double arg1 ;
	int arg2 = 0 ;
	double result ;
	VALUE vresult = Qnil;
	
	rb_scan_args(argc, argv, "11", &varg1, &varg2);
	Get_Bspline(self, arg0);
	arg1 = NUM2DBL(varg1);
	if (argc > 1) arg2 = NUM2INT(varg2);
	try {
		result = (*arg0)(arg1, arg2);
	} catch (const char *c) {
		rb_raise(rb_eRuntimeError, "%s", c);
	}
	vresult = rb_float_new(result);
	return vresult;
}
/*
	定積分
	Bspline#sekibun(ta, tb)
	  ta:区間の左端座標
	  tb:区間の右端座標
	Result:積分値
*/
static VALUE
_wrap_Bspline_sekibun(int argc, VALUE *argv, VALUE self)
{
	VALUE varg1, varg2;
	Bspline *arg0;
	double arg1, arg2, result;
	VALUE vresult = Qnil;
	
	rb_scan_args(argc, argv, "11", &varg1, &varg2);
	Get_Bspline(self, arg0);
	arg1 = NUM2DBL(varg1);
	try {
		result = (*arg0)(arg1, -1);
		if (argc > 1) {
			arg2 = NUM2DBL(varg2);
			result = arg0->sekibun(arg2) - result;
		}
	} catch (const char *c) {
		rb_raise(rb_eRuntimeError, "%s", c);
	}
	vresult = rb_float_new(result);
	return vresult;
}
/*
	Bspline#plot([], D, d = 0)
	 []:空配列
	  D:分割数
	  d:微分階数（省略可）
	Result:配列の要素数
*/
static VALUE
_wrap_Bspline_plot(int argc, VALUE *argv, VALUE self)
{
	Bspline *arg0;
	VALUE varg1; double *arg1, *arg2;
	VALUE varg3; int arg3;
	VALUE varg4; int arg4 = 0;
	rb_scan_args(argc, argv, "21", &varg1, &varg3, &varg4);
	Check_Type(varg1, T_ARRAY);
	arg3 = NUM2INT(varg3);
	if (argc > 2) arg4 = NUM2INT(varg4);
	Get_Bspline(self, arg0);
	int argn = RARRAY_LEN(varg1);
	double *argx = ALLOC_N(double, argn);
	for (int i = 0; i < argn; ++i) {
		VALUE x = rb_ary_shift(varg1);
		argx[i] = NUM2DBL(x);
	}
	int result;
	try {
		result = plots(*arg0, argn, argx, &arg1, &arg2, arg3, arg4);
	} catch (const char *c) {
		rb_raise(rb_eRuntimeError, "%s", c);
	}
	for (int i = 0; i < result; i++) {
		VALUE v1 = rb_float_new(arg1[i]);
		VALUE v2 = rb_float_new(arg2[i]);
		VALUE vv = rb_ary_new();
		rb_ary_push(vv, v1);
		rb_ary_push(vv, v2);
		rb_ary_push(varg1, vv);
		if (rb_block_given_p()) rb_yield(vv);
	}
	free(arg1);
	free(arg2);
	free(argx);
	return INT2NUM(result);
}
/*
// 周期境界条件によるＢスプライン補間
*/

static VALUE cTspline;

static void free_Tspline(Bspline *arg0) { delete arg0; }

#define Wrap_Tspline(klass, ptr) (\
	(ptr) ? Data_Wrap_Struct(klass, 0, VOIDFUNC(free_Tspline), ptr) : Qnil )

#define Get_Tspline(val, ptr) {\
	if (NIL_P(val)) ptr = NULL;\
	else {\
	    if (!rb_obj_is_kind_of(val, cTspline))\
	        rb_raise(rb_eTypeError, "wrong argument type (expected Tspline)");\
	    Data_Get_Struct(val, Bspline, ptr);\
	    if (!ptr) rb_raise(rb_eRuntimeError, "This Tspline already released");\
	}\
}

//  Tspline#new
static VALUE
_wrap_new_Tspline(int argc, VALUE *argv, VALUE klass) {
	VALUE self = Data_Wrap_Struct(klass, 0, VOIDFUNC(free_Tspline), 0);
	rb_obj_call_init(self, argc, argv);
	return self;
}
//	Tspline#initialize
static VALUE
_wrap_init_Tspline(VALUE self, VALUE varg, VALUE vargj)
{
	int argn, argj;
	double *xy;
	VALUE vx, vy, va;
	
	Check_Type(varg,T_ARRAY);
	argn = RARRAY_LEN(varg);
	argj = NUM2INT(vargj);
	if (argn < argj)
		rb_raise(rb_eArgError, "This Array is too small");
	xy = ALLOC_N(double, argn * 2);
	for (int i = 0; i < argn; i++) {
		va = RARRAY_PTR(varg)[i];
		Check_Type(va,T_ARRAY);
		if (RARRAY_LEN(va) != 2)
			rb_raise(rb_eArgError, "Need both x,y value");
		vx = RARRAY_PTR(va)[0];
		vy = RARRAY_PTR(va)[1];
		xy[i*2]   = NUM2DBL(vx);
		xy[i*2+1] = NUM2DBL(vy);
	}
	try {
		DATA_PTR(self) = new Bspline(argn, (double(*)[2])xy, argj, argn-1);
	} catch (const char *c) {
		rb_raise(rb_eRuntimeError, "%s", c);
	}
	free(xy);
	return self;
}
static VALUE
_wrap_Tspline_bracket(int argc, VALUE *argv, VALUE self)
{
	Bspline *bsp;
	double arg;
	VALUE val, vargs = rb_ary_new();
	VALUE vresult = rb_ary_new();
	Get_Tspline(self, bsp);
	for (int i = 0; i < argc; i++) {
		if (TYPE(argv[i]) == T_ARRAY) {
			val = rb_obj_dup(argv[i]);
			vargs = rb_ary_concat(vargs, _wrap_ary_flatten(val));
		}
		else
			rb_ary_push(vargs, argv[i]);
	}
	argc = RARRAY_LEN(vargs);
	try {
		if (argc == 1) {
			arg = NUM2DBL(RARRAY_PTR(vargs)[0]);
			vresult = rb_float_new((*bsp)[arg]);
		}
		else
			for (int i = 0; i < argc; i++) {
				arg = NUM2DBL(RARRAY_PTR(vargs)[i]);
				val = rb_float_new((*bsp)[arg]);
				rb_ary_push(vresult, val);
			}
	} catch (const char *c) {
		rb_raise(rb_eRuntimeError, "%s", c);
	}
	return vresult;
}
static VALUE
_wrap_Tspline_value(int argc, VALUE *argv, VALUE self)
{
	VALUE varg1 ;
	VALUE varg2 ;
	Bspline *arg0 ;
	double arg1 ;
	int arg2 = 0 ;
	double result ;
	VALUE vresult = Qnil;
	
	rb_scan_args(argc, argv, "11", &varg1, &varg2);
	Get_Tspline(self, arg0);
	arg1 = NUM2DBL(varg1);
	if (argc > 1) {
	    arg2 = NUM2INT(varg2);
	}
	try {
		result = (*arg0)(arg1, arg2);
	} catch (const char *c) {
		rb_raise(rb_eRuntimeError, "%s", c);
	}
	vresult = rb_float_new(result);
	return vresult;
}
static VALUE
_wrap_Tspline_sekibun(int argc, VALUE *argv, VALUE self)
{
	VALUE varg1, varg2;
	Bspline *arg0;
	double arg1, arg2, result;
	VALUE vresult = Qnil;
	
	rb_scan_args(argc, argv, "11", &varg1, &varg2);
	Get_Tspline(self, arg0);
	arg1 = NUM2DBL(varg1);
	try {
		result = (*arg0)(arg1, -1);
		if (argc > 1) {
			arg2 = NUM2DBL(varg2);
			result = arg0->sekibun(arg2) - result;
		}
	} catch (const char *c) {
		rb_raise(rb_eRuntimeError, "%s", c);
	}
	vresult = rb_float_new(result);
	return vresult;
}
static VALUE
_wrap_Tspline_plot(int argc, VALUE *argv, VALUE self)
{
	Bspline *arg0;
	VALUE varg1; double *arg1, *arg2;
	VALUE varg3; int arg3;
	VALUE varg4; int arg4 = 0;
	rb_scan_args(argc, argv, "21", &varg1, &varg3, &varg4);
	Check_Type(varg1, T_ARRAY);
	arg3 = NUM2INT(varg3);
	if (argc > 2) arg4 = NUM2INT(varg4);
	Get_Tspline(self, arg0);
	int argn = RARRAY_LEN(varg1);
	double *argx = ALLOC_N(double, argn);
	for (int i = 0; i < argn; ++i) {
		VALUE x = rb_ary_shift(varg1);
		argx[i] = NUM2DBL(x);
	}
	int result;
	try {
		result = plots(*arg0, argn, argx, &arg1, &arg2, arg3, arg4);
	} catch (const char *c) {
		rb_raise(rb_eRuntimeError, "%s", c);
	}
	for (int i = 0; i < result; i++) {
		VALUE v1 = rb_float_new(arg1[i]);
		VALUE v2 = rb_float_new(arg2[i]);
		VALUE vv = rb_ary_new();
		rb_ary_push(vv, v1);
		rb_ary_push(vv, v2);
		rb_ary_push(varg1, vv);
		if (rb_block_given_p()) rb_yield(vv);
	}
	free(arg1);
	free(arg2);
	free(argx);
	return INT2NUM(result);
}

static void
_wrap_Init_Bspline(void)
{
	cBspline = rb_define_class_under(mBspline, "Bspline", rb_cObject);
	rb_include_module(cBspline, rb_mEnumerable);
	rb_define_singleton_method(cBspline, "new", VALUEFUNC(_wrap_new_Bspline), -1);
	rb_define_method(cBspline, "initialize", VALUEFUNC(_wrap_init_Mspline), -1);
	rb_define_method(cBspline, "[]", VALUEFUNC(_wrap_Bspline_bracket), -1);
	rb_define_method(cBspline, "value", VALUEFUNC(_wrap_Bspline_value), -1);
	rb_define_method(cBspline, "sekibun", VALUEFUNC(_wrap_Bspline_sekibun), -1);
	rb_define_method(cBspline, "plot", VALUEFUNC(_wrap_Bspline_plot), -1);

	cTspline = rb_define_class_under(mBspline, "Tspline", rb_cObject);
    rb_include_module(cTspline, rb_mEnumerable);
	rb_define_singleton_method(cTspline, "new", VALUEFUNC(_wrap_new_Tspline), -1);
	rb_define_method(cTspline, "initialize", VALUEFUNC(_wrap_init_Tspline), 2);
	rb_define_method(cTspline, "[]", VALUEFUNC(_wrap_Tspline_bracket), -1);
	rb_define_method(cTspline, "value", VALUEFUNC(_wrap_Tspline_value), -1);
	rb_define_method(cTspline, "sekibun", VALUEFUNC(_wrap_Tspline_sekibun), -1);
	rb_define_method(cTspline, "plot", VALUEFUNC(_wrap_Tspline_plot), -1);
}

static VALUE cGspline;

#include "basis/gspline.h"

static void Wrap_free_Gspline(Gspline *arg0) { delete arg0; }

#define Wrap_Gspline(klass, ptr) (\
    (ptr) ? Data_Wrap_Struct(klass, 0, VOIDFUNC(Wrap_free_Gspline), ptr) : Qnil )

#define Get_Gspline(val, ptr) {\
    if (NIL_P(val)) ptr = NULL;\
    else {\
        if (!rb_obj_is_kind_of(val, cGspline))\
            rb_raise(rb_eTypeError, "wrong argument type (expected Gspline)");\
        Data_Get_Struct(val, Gspline, ptr);\
        if (!ptr) rb_raise(rb_eRuntimeError, "This Gspline already released");\
    }\
}

//	Gspline#new
static VALUE
_wrap_new_Gspline(int argc, VALUE *argv, VALUE klass)
{
	VALUE self = Data_Wrap_Struct(klass, 0, VOIDFUNC(Wrap_free_Gspline), 0);
	rb_obj_call_init(self, argc, argv);
	return self;
}
//	Gspline#initialize
static VALUE
_wrap_init_Gspline(VALUE self, VALUE varg, VALUE vargj, VALUE vargg)
{
    int argn, argj; double argg;
    VALUE va, vx, vy, vw; int vl;
    double *x, *y, *w;

    Check_Type(varg,T_ARRAY);
	argn = RARRAY_LEN(varg);
	argj = NUM2INT(vargj);
    if (argn < argj)
    	rb_raise(rb_eArgError, "This Array is too small");
	argg = NUM2DBL(vargg);
	x = ALLOC_N(double, argn);
	y = ALLOC_N(double, argn);
	w = ALLOC_N(double, argn);
	for (int i = 0; i < argn; i++) {
		va = RARRAY_PTR(varg)[i];
		Check_Type(va,T_ARRAY);
		vl = RARRAY_LEN(va);
		if ( vl < 2 || vl > 3 )
			rb_raise(rb_eArgError, "Need both x,y value");
		vx = RARRAY_PTR(va)[0];
		vy = RARRAY_PTR(va)[1];
		x[i] = NUM2DBL(vx);
		y[i] = NUM2DBL(vy);
		if (vl < 3)
			w[i] = 1.0;
		else {
			vw = RARRAY_PTR(va)[2];
			w[i] = NUM2DBL(vw);
		}
	}
	try {
		DATA_PTR(self) = new Gspline(argn, x, y, w, argj, argg);
	} catch (const char *c) {
		rb_raise(rb_eRuntimeError, "%s", c);
	}
	free(w);
	free(y);
	free(x);
    return self;
}
//	Gspline#[]
static VALUE
_wrap_Gspline_bracket(int argc, VALUE *argv, VALUE self)
{
	Gspline *gsp;
	double arg;
	VALUE val, vargs = rb_ary_new();
	VALUE vresult = rb_ary_new();
	Get_Gspline(self, gsp);
	for (int i = 0; i < argc; i++) {
		if (TYPE(argv[i]) == T_ARRAY) {
			val = rb_obj_dup(argv[i]);
			vargs = rb_ary_concat(vargs, _wrap_ary_flatten(val));
		}
		else
			rb_ary_push(vargs, argv[i]);
	}
	argc = RARRAY_LEN(vargs);
	try {
		if (argc == 1) {
			arg = NUM2DBL(RARRAY_PTR(vargs)[0]);
			vresult = rb_float_new((*gsp)(arg));
		}
		else
			for (int i = 0; i < argc; i++) {
				arg = NUM2DBL(RARRAY_PTR(vargs)[i]);
				val = rb_float_new((*gsp)(arg));
				rb_ary_push(vresult, val);
			}
	} catch (const char *c) {
		rb_raise(rb_eRuntimeError, "%s", c);
	}
	return vresult;
}
//	Gspline#()
static VALUE
_wrap_Gspline_value(int argc, VALUE *argv, VALUE self)
{
    VALUE varg1 ;
    VALUE varg2 ;
    Gspline *arg0 ;
    double arg1 ;
    int arg2 = 0 ;
    double result ;
    VALUE vresult = Qnil;
    
    rb_scan_args(argc, argv, "11", &varg1, &varg2);
    Get_Gspline(self, arg0);
    arg1 = NUM2DBL(varg1);
    if (argc > 1) {
        arg2 = NUM2INT(varg2);
    }
    try {
 	   result = (*arg0)(arg1, arg2);
	} catch (const char *c) {
		rb_raise(rb_eRuntimeError, "%s", c);
	}
    vresult = rb_float_new(result);
    return vresult;
}
//	Gspline#sekibun
static VALUE
_wrap_Gspline_sekibun(VALUE self, VALUE varg1)
{
	Gspline *arg0;
	double arg1, result;
	VALUE vresult = Qnil;
	
	Get_Gspline(self, arg0);
	arg1 = NUM2DBL(varg1);
	try {
		result = arg0->sekibun(arg1);
	} catch (const char *c) {
		rb_raise(rb_eRuntimeError, "%s", c);
	}
	vresult = rb_float_new(result);
	return vresult;
}
static VALUE
_wrap_Gspline_size(VALUE self)
{
    Gspline *arg0 ;
    int result ;
    VALUE vresult = Qnil;
    
    Get_Gspline(self, arg0);
    try {
    	result = arg0->size();
	} catch (const char *c) {
		rb_raise(rb_eRuntimeError, "%s", c);
	}
    vresult = INT2NUM(result);
    return vresult;
}
static VALUE
_wrap_Gspline_getx(VALUE self, VALUE varg1)
{
    Gspline *arg0 ;
    int arg1 ;
    double result ;
    VALUE vresult = Qnil;
    
    Get_Gspline(self, arg0);
    arg1 = NUM2INT(varg1);
    try {
    	result = arg0->getx(arg1);
	} catch (const char *c) {
		rb_raise(rb_eRuntimeError, "%s", c);
	}
    vresult = rb_float_new(result);
    return vresult;
}
static VALUE
_wrap_Gspline_gety(VALUE self, VALUE varg1)
{
    Gspline *arg0 ;
    int arg1 ;
    double result ;
    VALUE vresult = Qnil;
    
    Get_Gspline(self, arg0);
    arg1 = NUM2INT(varg1);
    try {
		result = arg0->gety(arg1);
	} catch (const char *c) {
		rb_raise(rb_eRuntimeError, "%s", c);
	}
    vresult = rb_float_new(result);
    return vresult;
}
static VALUE
_wrap_Gspline_getw(VALUE self, VALUE varg1)
{
    Gspline *arg0 ;
    int arg1 ;
    double result ;
    VALUE vresult = Qnil;
    
    Get_Gspline(self, arg0);
    arg1 = NUM2INT(varg1);
    try {
	    result = arg0->getw(arg1);
	} catch (const char *c) {
		rb_raise(rb_eRuntimeError, "%s", c);
	}
    vresult = rb_float_new(result);
    return vresult;
}
/* iterator: each() */
static VALUE
_wrap_Gspline_each(VALUE self)
{
	Gspline *gsp;

	Get_Gspline(self, gsp);
	try {
		int size = gsp->size();
		for (int i = 0; i < size; i++) {
			VALUE vv = rb_ary_new();
			rb_ary_push(vv, rb_float_new(gsp->getx(i)));
			rb_ary_push(vv, rb_float_new(gsp->gety(i)));
			rb_yield(vv);
		}
	} catch (const char *c) {
		rb_raise(rb_eRuntimeError, "%s", c);
	}
	return Qnil;
}

static void
_wrap_Init_Gspline(void)
{
    cGspline = rb_define_class_under(mBspline, "Gspline", rb_cObject);
    rb_include_module(cGspline, rb_mEnumerable);
	rb_define_singleton_method(cGspline, "new", VALUEFUNC(_wrap_new_Gspline), -1);
	rb_define_method(cGspline, "initialize", VALUEFUNC(_wrap_init_Gspline), 3);
    rb_define_method(cGspline, "[]", VALUEFUNC(_wrap_Gspline_bracket), -1);
    rb_define_method(cGspline, "value", VALUEFUNC(_wrap_Gspline_value), -1);
	rb_define_method(cGspline, "sekibun", VALUEFUNC(_wrap_Gspline_sekibun), 1);
    rb_define_method(cGspline, "size", VALUEFUNC(_wrap_Gspline_size), 0);
    rb_define_method(cGspline, "getx", VALUEFUNC(_wrap_Gspline_getx), 1);
    rb_define_method(cGspline, "gety", VALUEFUNC(_wrap_Gspline_gety), 1);
    rb_define_method(cGspline, "getw", VALUEFUNC(_wrap_Gspline_getw), 1);
    rb_define_method(cGspline, "each", VALUEFUNC(_wrap_Gspline_each), 0);
}

#ifdef __cplusplus
extern "C" {
#endif

void Init_bspline(void)
{
	mBspline = rb_define_module("BSPLINE");

	_wrap_Init_Bspline();
	_wrap_Init_Gspline();
}

#ifdef __cplusplus
}
#endif

