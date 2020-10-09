
data_eos = eos_dir / 'MS1_PP.eos.h5'
data_eos_hyb = eos_dir / 'HYB1.80_MS1_PP.eos.h5'

conf_tests = configuration_data()
conf_tests.set('PATH_EOS', data_eos)
conf_tests.set('PATH_EOS_HYB', data_eos_hyb)

if get_option('build_tests')
  subdir('src')
endif

if get_option('build_benchmarks')
  subdir('benchmarks')
  subdir('sample_root')
endif 
