project(
	'digest',
	'cpp',
	default_options : [
		'cpp_std=c++17',
		'optimization=3',
	]
)

nthash_dep = subproject('ntHash').get_variable('lib_dep')

include_dirs = [include_directories('include/Digester'), subproject('ntHash').get_variable('include_dirs')]
sources = ['src/digester.cpp', 'src/mod_minimizer.cpp',
			'src/syncmer.cpp', 'src/window_minimizer.cpp']

digest_lib = static_library(
	'digest',
	sources,
	include_directories: include_dirs,
	dependencies: nthash_dep,
)

digest_dep = declare_dependency(
	link_with: digest_lib,
	include_directories: include_dirs,
)

catch2 = dependency('catch2-with-main')
executable(
	'digest',
	'test.cpp',
	dependencies : [catch2, digest_dep],
	install : true
)
