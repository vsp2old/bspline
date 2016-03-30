#! /usr/local/bin/ruby
require 'bspline'
include BSPLINE

=begin

# Interpolation with periodic boundary condition

【Module】		BSPLINE
【Class】 		Tspline
【Method】
（1）new			Initialize
	obj = Tspline.new([[x0,y0],...,[xn,yn]], j)
	:1 list of data points.
	:2 dimension
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

puts "# Interpolation of the Jacobi function"
xy = [
	[ 0.0,    0.0    ],[ 0.5474, 0.5    ],[ 1.2837, 0.86603],[ 2.7681, 1.0    ],
	[ 4.2525, 0.86603],[ 4.9888, 0.5    ],[ 5.5362, 0.0    ],[ 6.0836,-0.5    ],
	[ 6.8199,-0.86603],[ 8.3043,-1.0    ],[ 9.7887,-0.86603],[10.5250,-0.5    ],
	[11.0724, 0.0    ]]
Jbn = ARGV[0].to_i
Dp = 10
tp = xy.tspline order:5
vv = []
xy.each do |p|
	printf "%8.5f, % f\n", p[0], p[1]
	vv.push p[0]
end
printf "# value of interpolation points, Dp = %d", Dp
if Jbn == 0
	print "\n"
else
	printf ", Jbn = %d\n", Jbn
end
s = tp.plot(vv, 10, Jbn) do |u, v|
	printf "%8.5f, % f\n", u, v
end
# Draw Graph
require "gnuplot"
 
File.open( "extspline.plt", "w") do |gp|
	Gnuplot::Plot.new( gp ) do |plot|
		plot.title  'Jacobi'
		plot.ylabel 'Y'
		plot.xlabel 'X'
		x = vv.map {|v| v[0]}
		y = vv.map {|v| v[1]}
		plot.data << Gnuplot::DataSet.new( [x, y] ) do |ds|
			ds.with = "lines"
			ds.linewidth = 2
			ds.notitle
		end
		y = x.map {|v| tp.value(v, 1)}
		plot.data << Gnuplot::DataSet.new( [x, y] ) do |ds|
			ds.with = "lines"
			ds.title = "Differential"
		end
	end
end
