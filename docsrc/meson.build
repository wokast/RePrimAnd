

bdir_doc = meson.current_build_dir() 

subdir('doxygen')
subdir('sphinx')

install_subdir(sphinx_out_dir, 
               install_dir : 'share/doc/libreprimand/html', \
               strip_directory : true)
