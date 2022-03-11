

configure_file(input : 'bench_config.h.in',
               output : 'bench_config.h',
               configuration : conf_tests)


sources_bench = ['benchmark_con2prim_mhd.cc']

exe_bench = executable('benchmark_c2p', sources : sources_bench, 
                       dependencies : [dep_reprim])

sources_acc = ['accuracy_con2prim_mhd.cc']

exe_acc = executable('accuracy_c2p', sources : sources_acc, 
                      dependencies : [dep_reprim])


sources_bench_tov = ['benchmark_tovsolver.cc']

exe_bench_tov = executable('bench_tov', sources : sources_bench_tov, 
                           dependencies : [dep_reprim])

