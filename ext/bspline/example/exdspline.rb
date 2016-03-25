#! /usr/local/bin/ruby
require 'bspline'
include BSPLINE

=begin

# Interpolation with boundary condition by differential value

【Module】		BSPLINE
【Class】 		Bspline
【Method】
（1）new		Initialize
	obj = Bspline.new([[x1,y1],...,[xn,yn]],j,[[xn+1,yn+1],...,[xn+d,yn+d],[b1,...,bd])
	:1 list of data points.
	:2 dimension
	:3 list of additional data points
	:4 list of order of differential value
（2）[]			Calculate interpolation
	obj[x]
	obj[x0,...,xi]
（3）value		Calculate interpolation with differential value
	obj.value(x, b = 0)
	b: order of differential value （optional）
（4）sekibun		Calculate interpolation with integrated value
	obj.sekibun(a)		: result is a indefinite integral of point a
	obj.sekibun(a,b)	: result is a definite integral of section [a,b]
（5）plot
	obj.plot([x0,...,xn], d, b = 0) { |x,y| ... }
	:1 list of data points
	:2 number of the division
	:3 order of differential value （optional）
=end

puts "# Interpolation of the sin(x) function"
xy = [[0.0],[1.0],[2.0],[3.0],[4.0],[5.0],[6.283185]]
xy.each do |a|
	a.push Math.sin(a[0])
end
C = [[0.0, 1.0],[0.0, 0.0],[6.283185, 1.0],[6.283185, 0.0]]
D = [1, 2, 1, 2]
Jbn = ARGV[0].to_i
Dp = 10
dp = xy.bspline order:5, condition:C, differential_order:D
vv = []
xy.each do |p|
	printf "%f, %f\n", p[0], p[1]
	vv.push p[0]
end
printf "# value of interpolation points, Dp = %d", Dp
if Jbn == 0
	print "\n"
else
	printf ", Jbn = %d\n", Jbn
end
s = dp.plot(vv, Dp, Jbn) do |a,b|
	printf "%f, %f\n", a, b
end
# Draw Graph
require "gnuplot"

Gnuplot.open do |gp|
	Gnuplot::Plot.new( gp ) do |plot|
		plot.title  'sin(x)'
		plot.ylabel 'Y'
		plot.xlabel 'X'
		x = vv.map {|v| v[0]}
		y = vv.map {|v| v[1]}
		plot.data << Gnuplot::DataSet.new( [x, y] ) do |ds|
			ds.with = "lines"
			ds.linewidth = 2
			ds.notitle
		end
		y = x.map {|v| dp.sekibun(v) - 1}
		plot.data << Gnuplot::DataSet.new( [x, y] ) do |ds|
			ds.with = "lines"
			ds.title = "Integral"
		end
	end
end
