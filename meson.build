project('digest', 'cpp')

add_project_arguments('-std=c++17', language : 'cpp')
add_project_arguments('-O3', language : 'cpp')


nthash_dep = subproject('ntHash').get_variable('lib_dep')

# catch2 = cmake.subproject('catch2')
# cm_lib = catch2.dependency('Catch2WithMain')
cm_lib = subproject('catch2').get_variable('catch2_with_main_dep')


include_dirs = [include_directories('digester/include/Digester'), 'subprojects/ntHash/include/nthash']
sources = ['digester/src/Digester/digester.cpp', 'digester/src/Digester/mod_minimizer.cpp',
            'digester/src/Digester/syncmer.cpp', 'digester/src/Digester/window_minimizer.cpp']

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

executable('digest',
           'test.cpp',
           dependencies : [cm_lib, digest_dep],
           install : true)
