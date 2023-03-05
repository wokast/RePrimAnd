
project('RePrimAnd', ['cpp','c'], default_options : ['cpp_std=c++11'],
        version : '1.6', license : 'CC BY-NC-SA 4.0')

project_headers_dest = 'reprimand'

dep_boost = dependency('boost')
dep_gsl   = dependency('gsl', version : '>=2.0')
dep_h5    = dependency('hdf5')

subdir('EOS')

subdir('library')

subdir('tests')

if get_option('build_documentation')
  subdir('docsrc')
endif


clangtidy = find_program('clang-tidy', required: false)
if clangtidy.found()
    run_target(
        'tidy',
        command: [
            clangtidy,
            '-checks=cppcoreguidelines-*,bugprone-*,hicpp-*,misc-*,modernize-*,readability-*,-cppcoreguidelines-pro-bounds-array-to-pointer-decay,-hicpp-no-array-decay,-hicpp-special-member-functions,-cppcoreguidelines-pro-bounds-constant-array-index,-cppcoreguidelines-special-member-functions,-misc-non-private-member-variables-in-classes,-cppcoreguidelines-non-private-member-variables-in-classes,-hicpp-braces-around-statements,-readability-braces-around-statements,-misc-unused-parameters,-misc-non-private-member-variables-in-classes,-modernize-use-trailing-return-type',
            '-header-filter=.*', 
            '-p', meson.build_root()
        ] + sources_lib)
endif
