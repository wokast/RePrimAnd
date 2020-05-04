
RePrimAnd Documentation
=======================

Recovery of Primitives And EOS framework


Overview
--------

`RePrimAnd` is a support library for numerical simulations of general 
relativistic magnetohydrodynamics. 
If provides methods for recovering primitive variables like pressure
and velocity from the variables evolved in quasi-conservative 
formulations. Further, it provides a general framework for handling
matter equations of state. 

This version of the library was made public together with the article 
describing the recovery scheme:

Wolfgang Kastaun, Jay Vijay Kalinani, Riccardo Ciolfi, 
"Robust Recovery of Primitive Variables in Relativistic Ideal 
Magnetohydrodynamics" (2020).

Prospective users should be sure to check for updated versions
that will likely be made available in a public repository.


Developers
^^^^^^^^^^
`RePrimAnd` was written by Dr. Wolfgang Kastaun 
<physik@fangwolg.de>.
Thanks for additional testing go to Jay Vijay Kalinani.

License
^^^^^^^

.. raw:: html

   <a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License</a>.


Contents
--------

.. toctree::
   :maxdepth: 2

   installing
   notation
   c2p_imhd
   eos_thermal
   eos_barotr
   little_helpers


Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
