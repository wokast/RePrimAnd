

subdir('src')
subdir('scripts')

data_bench = ['perf_eosig_z_eps_Bzero.dat', \
              'perf_eosig_z_eps_Blarge.dat', \
              'perf_eosig_z_b_cold.dat', 'perf_eosig_z_b_hot.dat', \
              'perf_eoshyb_z_eps_Bzero.dat', \
              'perf_eoshyb_z_eps_Blarge.dat', \
              'perf_eoshyb_z_b_cold.dat', 'perf_eoshyb_z_b_hot.dat']
              
gen_bench = custom_target('gen_bench_', 
                          command : [exe_bench, '@OUTDIR@'], \
                          output : data_bench)

alias_target('benchdata', gen_bench)



data_acc = ['acc_eosig_z_eps_Bzero.dat', \
            'acc_eosig_z_eps_Blarge.dat', \
            'acc_eosig_z_b_cold.dat', 'acc_eosig_z_b_hot.dat', \
            'acc_eoshyb_z_eps_Bzero.dat', \
            'acc_eoshyb_z_eps_Blarge.dat', \
            'acc_eoshyb_z_b_cold.dat', 'acc_eoshyb_z_b_hot.dat']
              
gen_acc = custom_target('gen_acc_', 
                          command : [exe_acc, '@OUTDIR@'], \
                          output : data_acc)

alias_target('accuracydata', gen_acc)


figs_bench = ['perf_eosig_z_eps.pdf', \
              'perf_eosig_z_b.pdf', \
              'perf_eoshyb_z_eps.pdf', \
              'perf_eoshyb_z_b.pdf']

plt_bench  = custom_target('plt_bench_', 
                          command : [plot_bench, '@OUTDIR@'], \
                          depend_files : dep_plot_bench, \
                          input : gen_bench, \
                          output : figs_bench)

alias_target('benchplots', plt_bench)


figs_acc = ['acc_eosig_z_eps.pdf', \
              'acc_eosig_z_b.pdf', \
              'acc_eoshyb_z_eps.pdf', \
              'acc_eoshyb_z_b.pdf']

plt_acc  = custom_target('plt_acc_', 
                          command : [plot_acc, '@OUTDIR@'], \
                          depend_files : dep_plot_acc, \
                          input : gen_acc, \
                          output : figs_acc)

alias_target('accuracyplots', plt_acc)



