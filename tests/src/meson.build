

configure_file(input : 'test_config.h.in',
               output : 'test_config.h',
               configuration : conf_tests)

dep_utf = dependency('boost', static : false, modules : ['unit_test_framework'])

src_tst_c2p = ['test_con2prim_mhd.cc', 'test_utils.cc']

exe_tst_c2p = executable('test_c2p', sources : src_tst_c2p, 
                      dependencies : [dep_reprim, dep_utf], 
                      cpp_args : '-DBOOST_TEST_DYN_LINK')

test('Con2Prim', exe_tst_c2p, timeout : 90)


src_tst_eos = ['test_eos.cc', 'test_utils.cc']

exe_tst_eos = executable('test_eos', sources : src_tst_eos, 
                      dependencies : [dep_reprim, dep_utf], 
                      cpp_args : '-DBOOST_TEST_DYN_LINK')

test('EOS', exe_tst_eos)
