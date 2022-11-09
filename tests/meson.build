
subdir('data')

data_eos = eos_dir / 'MS1_PP.eos.h5'
data_eos_pp = eos_dir / 'MS1_Read_PP.eos.h5'
data_eos_hyb = eos_dir / 'HYB1.80_MS1_PP.eos.h5'

conf_tests = configuration_data()
conf_tests.set('PATH_EOS', data_eos)
conf_tests.set('PATH_EOS_PP', data_eos_pp)
conf_tests.set('PATH_EOS_HYB', data_eos_hyb)
conf_tests.set('PATH_TOV_REF', test_data_dir)
conf_tests.set('PATH_TOV_EOS', eos_dir)

if get_option('build_tests')
  subdir('src')
endif

if get_option('build_benchmarks')
  subdir('benchmarks')
  subdir('sample_root')
endif 
