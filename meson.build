project(
	'digest',
	'cpp',
	version: '0.1.0',
	default_options : [
		'cpp_std=c++17',
		'optimization=3',
	]
)

nthash = subproject('ntHash')
nthash_dep = nthash.get_variable('lib_dep')
include_dirs = [include_directories('include'), nthash.get_variable('include_dirs')]
sources = ['src/digester.cpp', 'src/mod_minimizer.cpp', 'src/syncmer.cpp', 'src/window_minimizer.cpp']

digest_lib = static_library(
	'digest',
	sources,
	include_directories: include_dirs,
	dependencies: nthash_dep,
	install: true,
    install_dir: 'lib',
)
install_headers(
	'include/digest/digester.hpp', 'include/digest/mod_minimizer.hpp',
	'include/digest/syncmer.hpp', 'include/digest/window_minimizer.hpp',
	install_dir: 'include/digest'
)

digest_dep = declare_dependency(
	link_with: digest_lib,
	include_directories: include_dirs,
)

### test ###
catch2 = dependency('catch2-with-main')
tests = executable(
	'tests',
	'test/tests.cpp',
	dependencies : [catch2, digest_dep],
)
test('tests', tests)
