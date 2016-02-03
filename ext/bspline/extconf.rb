require 'mkmf'

dir_config('BSPLINE','.')
$CPPFLAGS += " -DNDEBUG -std=c++11"
if have_header('bspline.h')
	create_makefile('bspline')
end
