#! /usr/local/bin/ruby
require 'bspline'
include BSPLINE

=begin

Interpolation with boundary condition by additional data points

【Module】		BSPLINE
【Class】 		Bspline
【Method】
（1）new			Initialize
	obj = Bspline.new([[x1,y1],...,[xn,yn]], j, [[xn+1,yn+1],...,[xn+d,yn+d]])
	:1 list of data points.
	:2 dimension
	:3 list of additional data points
（2）[]			Calculate interpolation
	obj[x]
	obj[x0,...,xi]
（3）value		Calculate interpolation with differential value
	obj.value(x, b = 0)
	b:order of differential value （optional）
（4）sekibun		Calculate interpolation with integrated value
	obj.sekibun(a)		:result is a indefinite integral of point a
	obj.sekibun(a,b)	:result is a definite integral of section [a,b]
（5）plot
	obj.plot([x0,...,xn], d, b = 0) { |x,y| ... }
	:1 list of data points
	:2 number of the division
	:3 order of differential value （optional）
=end

def plotsub(bp, vp, dp, b = 0)
#
#	bp: Bspline Object
#	vp: list of data points
#	dp: number of the division
#   b : order of differential value （optional）
#
	n = vp.size
	xp = []
	for i in 0...n
		x = vp.shift
		xp.push x
	end
	s = (n - 1) * dp + 1
	for l in 1...n
		dt = (xp[l] - xp[l-1]) / dp
		n0 = (l == 1 ? 0 : 1)
		for nd in n0..dp
			np = nd + (l-1) * dp;
			tp = xp[l-1] + nd * dt;
			yp = (b == 0) ? bp[tp] : bp.value(tp, b);
			vp.push [tp, yp]
			yield tp, yp
		end
	end
	return s;
end
#
puts "# Interpolation of the Bessel function"
#
ad = [[0.0, 1.0    ],[0.8, 0.84629],[1.6, 0.45540],[2.0, 0.22389],
      [2.4, 0.00251],[3.2,-0.32019],[4.0,-0.39715]]
ac = [[0.4, 0.96040],[1.2, 0.67113],[2.8,-0.18504],[3.6,-0.39177]]
Jbn = ARGV[0].to_i
Dp = 10
bp = ad.bspline order:5, condition:ac
vv = []
ad.each do |p, q|
	printf "% .2f, % f\n", p, q
	vv.push p
end
printf "# value of interpolation points, Dp = %d", Dp
if Jbn == 0
	print "\n"
else
	printf ", Jbn = %d\n", Jbn
end
s = plotsub(bp, vv, Dp, Jbn) do |a, b|
	printf "% .2f, % f\n", a, b
end
# STDERR.puts s
require "gnuplot"
 
Gnuplot.open do |gp|
	Gnuplot::Plot.new( gp ) do |plot|
		plot.title  'Bessel'
		plot.ylabel 'Y'
		plot.xlabel 'X'
 
		x = vv.map {|v| v[0]}
		y = vv.map {|v| v[1]}
 
		plot.data << Gnuplot::DataSet.new( [x, y] ) do |ds|
			ds.with = "lines"
			ds.linewidth = 2
			ds.notitle
		end

		y = x.map {|v| bp.sekibun(v)}
 
		plot.data << Gnuplot::DataSet.new( [x, y] ) do |ds|
			ds.with = "lines"
			ds.title = "Integral"
		end
	end
end
