

subdir('src')


data_sroot = ['cold_vzero_bzero.dat',   \
              'cold_vlarge_bzero.dat',  \
              'cold_vzero_blarge.dat',  \
              'cold_vlarge_blarge.dat', \
              'hot_vzero_bzero.dat',    \
              'hot_vzero_blarge.dat',   \
              'hot_vlarge_bzero.dat',   \
              'hot_vlarge_blarge.dat']

gen_sroot = custom_target('gen_sroot_', 
                          command : [exe_sroot, '@OUTDIR@'], \
                          output : data_sroot)

alias_target('srootdata', gen_sroot)

