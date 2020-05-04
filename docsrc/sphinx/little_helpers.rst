

Little Helpers
--------------

Tensors
^^^^^^^

The library contains classes representing fixed-size tensors (not
tensor fields). The functionality is limited mostly to the needs
of the primitive recovery algorithm, but might also be useful 
elsewhere. The objects used by the library itself are:

* :cpp:type:`~EOS_Toolkit::sm_vec3u`: a 3-dimensional 
  vector (upper index)
* :cpp:type:`~EOS_Toolkit::sm_vec3l`: a 3-dimensional 
  vector (lower index)
* :cpp:type:`~EOS_Toolkit::sm_symt3u`: a 3-dimensional rank-2 
  symmetric tensor (upper index)
* :cpp:type:`~EOS_Toolkit::sm_symt3l`: a 3-dimensional rank-2 
  symmetric tensor (lower index)
* :cpp:type:`~EOS_Toolkit::sm_metric3`: a 3-metric containing 
  upper+lower components and determinant

Standard arithmetic operations are supported via operator overloading. 
Contraction operations between tensors are supported via the 
multiplication operator. Indices can be raised and lowered via the
metric. Consistent use of co- and contra-variant 
indices is enforced at compile time, e.g. one cannot add an upper index 
and a lower index vector.




.. doxygenclass:: EOS_Toolkit::sm_tensor1
   :project: RePrimAnd
   :members:

.. doxygenclass:: EOS_Toolkit::sm_tensor2_sym
   :project: RePrimAnd
   :members:

.. doxygenclass:: EOS_Toolkit::sm_metric
   :project: RePrimAnd
   :members:

.. doxygentypedef:: sm_vec3u

.. doxygentypedef:: sm_vec3l

.. doxygentypedef:: sm_symt3u

.. doxygentypedef:: sm_symt3l

.. doxygentypedef:: sm_metric3

Units
^^^^^

.. doxygenclass:: EOS_Toolkit::units
   :project: RePrimAnd
   :members:


Other
^^^^^

.. doxygenclass:: EOS_Toolkit::interval
   :project: RePrimAnd
   :members:
