

configure_file(input : 'sroot_config.h.in',
               output : 'sroot_config.h',
               configuration : conf_tests)

sources_sroot = ['froot_samp.cc']

exe_sroot = executable('froot_samp', sources : sources_sroot, 
                       dependencies : [dep_reprim])
