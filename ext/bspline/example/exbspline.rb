#! /usr/local/bin/ruby
require 'bspline'
include BSPLINE

=begin

Basic interpolation without boundary condition

【Module】		BSPLINE
【Class】 		Bspline
【Method】
（1）new			Initialize
	obj = Bspline.new([[x1,y1],...,[xn,yn]], j)
	:1 list of data points.
	:2 dimension
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
#
puts "# Interpolation of the sec2(x) function"
#
ad = [
	[-4.0, 1.34095e-3],[-3.6, 2.98189e-3],[-3.2, 6.62420e-3],[-2.8, 1.46827e-2],
	[-2.4, 3.23838e-2],[-2.0, 7.06508e-2],[-1.6, 1.50527e-1],[-1.2, 3.05020e-1],
	[-0.8, 5.59055e-1],[-0.4, 8.55639e-1],[ 0.4, 8.55639e-1],[ 0.8, 5.59055e-1],
	[ 1.2, 3.05020e-1],[ 1.6, 1.50527e-1],[ 2.0, 7.06508e-2],[ 2.4, 3.23838e-2],
	[ 2.8, 1.46827e-2],[ 3.2, 6.62420e-3],[ 3.6, 2.98189e-3],[ 4.0, 1.34095e-3]]
Jbn = ARGV[0].to_i
Dp = 10
bp = Bspline.new(ad, 5)
vv = []
ad.each do |p|
	printf "% .2f, %f\r\n", p[0], p[1]
	vv.push p[0];
end
printf "# value of interpolation points, Dp = %d", Dp
if Jbn == 0
	print "\n"
else
	printf ", Jbn = %d\n", Jbn
end
s = bp.plot(vv, Dp, Jbn) do |a,b|
	printf "% .2f, % f", a, b
	printf ", % f\n", bp.sekibun(a)
end
# STDERR.puts s
require "gnuplot"
 
Gnuplot.open do |gp|
	Gnuplot::Plot.new( gp ) do |plot|
		plot.title  'sec2(x)'
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
