py3_mod = import('python')
py3_inst = py3_mod.find_installation('python3', required:false,\
                                     modules:['pybind11'])

pyreprimand_path = meson.current_build_dir()

if py3_inst.found()
  dep_py3 = py3_inst.dependency()
  
  fnd_pyb11 = run_command(py3_inst.path(), '-m', 'pybind11', 
                           '--cmakedir')
  dep_pyb11 = dependency('pybind11', 
                         cmake_module_path : fnd_pyb11.stdout().strip(),
                         fallback : ['pybind11', 'pybind11_dep'])
  
  #pybind11_proj = subproject('pybind11')
  #dep_pyb11 = pybind11_proj.get_variable('pybind11_dep')
  
  #fnd_pyb11 = run_command(py3_inst.path(), '-m', 'pybind11', 
  #                          '--includes')  
  #dep_pyb11 = declare_dependency(
  #  include_directories: [fnd_pyb11.stdout().split('-I')[-1].strip()],
  #)
  
  
  py3_binding = py3_inst.extension_module('pyreprimand',
      sources:['pyreprimand.cpp', sources_lib], 
      include_directories : include_lib,
      implicit_include_directories : true, 
      dependencies: [dep_extern, dep_py3, dep_pyb11],
      install: true,
      subdir: py3_inst.get_path('platlib')
  )

endif
