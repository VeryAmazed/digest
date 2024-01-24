project(
	'digest',
	'cpp',
	version: '0.2.0',
	default_options : [
		'cpp_std=c++17',
		'optimization=3',
        'werror=true',
        'warning_level=3',
        # 'b_sanitize=thread',
	],
)

nthash = subproject('ntHash')
nthash_dep = nthash.get_variable('lib_dep')
include_dirs = [include_directories('include'), nthash.get_variable('include_dirs')]
sources = ['src/digester.cpp', 'src/mod_minimizer.cpp']

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
	'include/digest/data_structure.hpp',
  'include/digest/syncmer.tpp', 'include/digest/window_minimizer.tpp',
	'include/digest/thread_out.hpp', 'include/digest/thread_out.tpp',
	install_dir: 'include/digest'
)
install_headers(
	'include/nthash/kmer.hpp', 'include/nthash/internal.hpp',
	install_dir: 'include/nthash'
)

digest_dep = declare_dependency(
	link_with: digest_lib,
	include_directories: include_dirs,
)

if get_option('buildtype') != 'release'
  ### test ###
  catch2 = dependency('catch2-with-main')
  executable(
  	'tests',
  	'tests/test/test.cpp',
  	dependencies : [catch2, digest_dep],
  )

  ### benchmark ###
  bench = dependency('benchmark')
  executable(
    'bench',
    'tests/bench/benchmark.cpp',
    dependencies : [bench, digest_dep],
  )

  ### benchmark data structures ###
  # executable(
  #   'bench_ds',
  #   'tests/data_structure/bench_ds.cpp',
  #   dependencies : [bench, digest_dep],
  # )

  ### Expected Density ACTG ###
  executable(
    'expected_ACTG',
  	'tests/density/ACTG.cpp',
  	dependencies : [digest_dep]
  )
  
  ### Expected Density non-ACTG ###
  executable(
  	'expected_non-ACTG',
  	'tests/density/non-ACTG.cpp',
  	dependencies : [digest_dep]
  )
  
  ### test thread functions ###
  executable(
  	'test_thread',
  	'tests/test/test_thread.cpp',
  	dependencies : [catch2, digest_dep],
  )
endif
