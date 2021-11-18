.. _c2p_imhd_interface:

Interface
---------


The primitive recovery is provided in form of a function object
:cpp:class:`~EOS_Toolkit::con2prim_mhd`
which stores all fixed parameters such as EOS, required accuracy, and 
error policies. Once created, it can be used repeatetly as a function 
to recover primitives. The evolved variables are passed as a 
structure :cpp:class:`~EOS_Toolkit::cons_vars_mhd` and the 
recovered primitives are stored in a structure 
:cpp:class:`~EOS_Toolkit::prim_vars_mhd`. 


The recovery also 
requires the metric, which is passed in another object. The framework
provides objects representing tensors (*not* tensor fields),
used here to represent co- and contravariant 3-vectors and the 
3-metric. We use the convention that the 3-metric has positive
signature.


If the evolved variables are 
invalid, the code checks whether the given error policy allows to apply 
corrections to the conserved variables to make them valid. If yes, the 
structure with the evolved varibles is modified to be consistent with 
the corrected primitives. 
The conserved variables and metric are checked to ensure all values
are finite (not NAN or INF) and that the metric determinent is 
positive. 

Whether or not the evolved variables are valid after allowed corrections 
is reported via another object 
:cpp:class:`~EOS_Toolkit::c2p_mhd_report`. 
If the evolved variables violate the error policy, this also contains
the type of violation. In addition, it reports if atmosphere
was enforced, and if the evolved variables were corrected.
The object also provides a method to assemble a string with a 
diagnostic message.


If the recovery fails for any reason the conserved and primitive 
variables are all set to NAN. Otherwise, they are consistent
up to the root solving accuracy. The quantities :math:`\rho,W, v^i, D`
are always consistent to machiene precision.

A minimal example can be found in `example/minimal.cc`.


.. note::

   Since the function operates pointwise and does not need to know 
   anything about numerical grids, it is usable in any evolution code 
   based on (or able to link to) C++.

.. note::

   It is save to use the same function object in multiple threads. 
   Copying the function object inside an innermost multithreaded loop 
   should be avoided for performance reasons, since copying the EOS is 
   an atomic operation.
