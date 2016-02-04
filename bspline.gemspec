# coding: utf-8
lib = File.expand_path('../lib', __FILE__)
$LOAD_PATH.unshift(lib) unless $LOAD_PATH.include?(lib)
require 'bspline/version'

Gem::Specification.new do |spec|
  spec.name          = "bspline"
  spec.version       = BSPLINE::VERSION
  spec.authors       = ["vsp2old"]
  spec.email         = ["yatcho@hotmail.co.jp"]

  spec.description   = %q{Bspline interpolation libraly implemented in C-Extensions}
  spec.summary       = %q{Bspline interpolation libraly implemented in C-Extensions}
  spec.homepage      = ""

  spec.extensions    = %w[ext/bspline/extconf.rb]

  spec.files         = `git ls-files -z`.split("\x0").reject { |f| f.match(%r{^(test|spec|features)/}) }
  spec.bindir        = "exe"
  spec.executables   = spec.files.grep(%r{^exe/}) { |f| File.basename(f) }
  spec.require_paths = ["ext","lib"]

  spec.add_development_dependency "bundler", "~> 1.11"
  spec.add_development_dependency "rake", "~> 10.0"
end
