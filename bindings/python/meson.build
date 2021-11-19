py3_mod = import('python')
py3_inst = py3_mod.find_installation('python3', required:false,\
                                     modules:['pybind11'])

if py3_inst.found()
  dep_py3 = py3_inst.dependency()
  
  fnd_pyb11 = run_command(py3_inst.path(), '-m', 'pybind11', 
                            '--cmakedir')
  dep_pyb11 = dependency('pybind11', 
                         cmake_module_path : fnd_pyb11.stdout().strip())
  
  py3_binding = py3_inst.extension_module('pyreprimand',
      sources:['pyreprimand.cpp', sources_lib], 
      include_directories : include_lib,
      implicit_include_directories : true, 
      dependencies: [dep_extern, dep_py3, dep_pyb11],
      install: true,
      subdir: py3_inst.get_path('purelib')
  )

endif
