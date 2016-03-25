require "bspline/version"
require "bspline.so"

class Array
	def bspline order:, condition:nil, differential_order:nil
		if differential_order then
			BSPLINE::Bspline.new self, order, condition, differential_order
		elsif condition then
			BSPLINE::Bspline.new self, order, condition
		else
			BSPLINE::Bspline.new self, order
		end
	end
	def tspline order:
		BSPLINE::Tspline.new self, order
	end
	def gspline order:, smoothing_factor:0
		BSPLINE::Gspline.new self, order, smoothing_factor
	end
end
