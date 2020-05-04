

doxygen  = find_program('doxygen')

bdir_doxy = meson.current_build_dir()

conf_doxy = configuration_data()
conf_doxy.set('DOXYGEN_OUT_DIR', bdir_doxy)
conf_doxy.set('DOXYGEN_SOURCE_DIR', src_dir_lib)

doxyfile = configure_file(input : 'Doxyfile.in',
                          output : 'Doxyfile',
                          configuration : conf_doxy)

doxy_index   = 'index.xml'
doxy_out_dir = meson.current_build_dir()

gen_doxy = custom_target('doxygendoc', 
                command : [doxygen, doxyfile], \
                input : doxyfile, \
                output : doxy_index,
                depend_files : [sources_lib, headers_lib])

alias_target('doxygendoc', gen_doxy)
