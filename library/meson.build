

subdir('BasicStuff')
subdir('EOS_Thermal')
subdir('EOS_Thermal_Idealgas')
subdir('EOS_Thermal_Hybrid')
subdir('EOS_Barotropic')
subdir('Con2Prim_IMHD')
subdir('NeutronStar')


src_dir_lib = meson.current_source_dir()

sources_lib = [sources_basic_stuff, sources_eos_thermal, \
               sources_eos_idealgas, sources_eos_hybrid, \
               sources_eos_barotr, sources_c2p_imhd, \
               sources_tovsolver]
include_lib = [include_basic_stuff, include_eos_thermal, \
               include_eos_idealgas, include_eos_hybrid, \
               include_eos_barotr, include_c2p_imhd, \
               include_tovsolver]

headers_lib = [headers_basic_stuff, headers_eos_thermal, \
               headers_eos_idealgas, headers_eos_hybrid, \
               headers_eos_barotr, headers_c2p_imhd, \
               headers_tovsolver]

dep_extern  = [dep_boost, dep_gsl, dep_h5]


lib_reprim  = library('RePrimAnd', sources_lib, \
                include_directories : include_lib, \
                implicit_include_directories : true, \
                dependencies : dep_extern, \
                install : true)

dep_reprim = declare_dependency( \
    dependencies: dep_extern, \
    include_directories : include_lib, \
    link_with : lib_reprim)
