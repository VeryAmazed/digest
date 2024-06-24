project(
	'digest',
	'cpp',
	version: '0.2.0',
	default_options : [
		'cpp_std=c++17',
		'optimization=3',
        'werror=true',
        'warning_level=3',
        'buildtype=release'
        # 'b_sanitize=thread',
	],
)

nthash = subproject('ntHash')
nthash_dep = nthash.get_variable('lib_dep')
include_dirs = [include_directories('include'), nthash.get_variable('include_dirs')]

install_headers(
	'include/digest/digester.hpp', 'include/digest/mod_minimizer.hpp',
	'include/digest/syncmer.hpp', 'include/digest/window_minimizer.hpp',
    'include/digest/thread_out.hpp',
	'include/digest/data_structure.hpp',
	install_dir: 'include/digest'
)

digest_dep = declare_dependency(
	include_directories: include_dirs,
	dependencies: nthash_dep,
)

if get_option('buildtype') != 'release'	
  ### test ###
  catch2 = dependency('catch2-with-main')
  executable(
   'tests',
   'tests/test/test.cpp',
    dependencies : [catch2, digest_dep],
  )
  
  thread_dep = dependency('threads')

  ### benchmark ###
  bench = dependency('benchmark')
  executable(
   'bench',
   'tests/bench/benchmark.cpp',
   dependencies : [bench, digest_dep, thread_dep],
  )

  ### benchmark data structures ###
  executable(
    'bench_ds',
    'tests/data_structure/bench_ds.cpp',
    dependencies : [bench, digest_dep],
  )

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
   dependencies : [catch2, digest_dep, thread_dep],
  )

  doxygen = find_program('doxygen', required: false)
  if doxygen.found()
      run_target('docs', command: [meson.project_source_root() + '/docs/build.sh'])
  else
      warning('Documentation disabled without doxygen')
  endif
endif