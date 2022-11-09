
sphinx = find_program('sphinx-build')

sphinx_src_dir = meson.current_source_dir() 
sphinx_out_dir = meson.current_build_dir() 

conf_sphinx = configuration_data()
conf_sphinx.set('SPHINX_DOXYGEN', doxy_out_dir)
conf_sphinx.set('SPHINX_AUTODOC', eoscreate_dir)

sphinxfile = configure_file(input : 'conf.py.in',
                            output : 'conf.py',
                            configuration : conf_sphinx)

sphinx_sources = files('index.rst', 'notation.rst', \
  'eos_barotr.rst', 'eos_barotr_feat.rst', 'eos_barotr_ref.rst', \
  'eos_barotr_interf.rst', 'eos_barotr_available.rst', \
  'eos_barotr_custom.rst', 'eos_barotr_files.rst', \
  'eos_thermal.rst', 'eos_thermal_feat.rst', \
  'eos_thermal_ref.rst', 'eos_thermal_interf.rst', \
  'eos_thermal_available.rst', 'eos_thermal_custom.rst', \
  'eos_thermal_files.rst', \
  'c2p_imhd.rst', 'c2p_imhd_feat.rst', 'c2p_imhd_interf.rst', \
  'c2p_imhd_ref.rst', 'little_helpers.rst', 'installing.rst', \
  'tov_solver_feat.rst', 'tov_solver_interf.rst', \
  'tov_solver_ref.rst', 'tov_solver.rst', \
  'ns_seqs_feat.rst', 'ns_seqs_interf.rst', \
  'ns_seqs_files.rst', 'ns_seqs_ref.rst', 'ns_seqs.rst', \
  'python_interf.rst', 'refs.bib', \
  'minimal_tov.cc', 'minimal_seq.cc')

gen_sphinx = custom_target('sphinxdoc', \
  command : [sphinx, '-b', 'html', \
             '-c', sphinx_out_dir, \
             sphinx_src_dir, sphinx_out_dir], \
  input : [sphinx_sources, sphinxfile, eoscreate_src], \
  output : 'index.html', depends : gen_doxy)

alias_target('documentation', gen_sphinx)



