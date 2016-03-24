#! /usr/local/bin/ruby
require 'bspline'
include BSPLINE

fdata = [[0.0],[0.2],[0.4],[0.6],[0.8],[1.0]]
fdata.each do |a|
	x = a[0]
	a[1] = (1-x)*(1-x)*(1-x)*x
end
Dp = 10

bp = Bspline.new(fdata, 4)

vv = []
fdata.each do |p|
#	printf "%f, %f\n", p[0], p[1]
	vv.push p[0]
end
#printf "# value of interpolation points, Dp = %d\n", Dp
s = bp.plot(vv, Dp, -1) do |a,b|
	printf "% .10f % .10f\n", a, b
end

#STDERR.puts bp.sekibun(1)
require "gnuplot"
 
Gnuplot.open do |gp|
	Gnuplot::Plot.new( gp ) do |plot|
		plot.title  'Sample'
		plot.ylabel 'Y'
		plot.xlabel 'X'
 
		x = vv.map {|v| v[0]}
		y = vv.map {|v| v[1]}
 
		plot.data << Gnuplot::DataSet.new( [x, y] ) do |ds|
			ds.with = "lines"
			ds.linewidth = 2
			ds.notitle
		end
	end
end
