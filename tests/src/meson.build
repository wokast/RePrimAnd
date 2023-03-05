

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

src_tst_interp = ['test_interpol.cc', 'test_utils.cc']

exe_tst_interp = executable('test_interp', sources : src_tst_interp, 
                      dependencies : [dep_reprim, dep_utf], 
                      cpp_args : '-DBOOST_TEST_DYN_LINK')

test('INTERP', exe_tst_interp)

src_tst_tovsol = ['test_tovsol.cc', 'test_utils.cc']

exe_tst_tovsol = executable('test_tovsol', sources : src_tst_tovsol, 
                      dependencies : [dep_reprim, dep_utf], 
                      cpp_args : '-DBOOST_TEST_DYN_LINK')

test('TOVSOL', exe_tst_tovsol)

src_tst_nsseqs = ['test_nsseqs.cc', 'test_utils.cc']

exe_tst_nsseqs = executable('test_nsseqs', sources : src_tst_nsseqs, 
                      dependencies : [dep_reprim, dep_utf], 
                      cpp_args : '-DBOOST_TEST_DYN_LINK')

test('NSSEQS', exe_tst_nsseqs)


src_gen_tov = ['gen_tov_ref_data.cc']

exe_gen_tov = executable('gen_tov_ref_data', sources : src_gen_tov, 
                      dependencies : [dep_reprim])

