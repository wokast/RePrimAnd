.. _thermal_interface:

Interface
---------


The EOS with thermal+composition degrees of freedom are represented
by objects of type :cpp:class:`~EOS_Toolkit::eos_thermal`. 
These objects 
represent the generic EOS interface and do not implement any EOS. 
Instead, they forward any request to an actual implementation of a 
particular EOS. To obtain an 
:cpp:class:`~EOS_Toolkit::eos_thermal`
object representing a particular EOS, one uses functions such as 
:cpp:func:`~EOS_Toolkit::make_eos_hybrid` 
provided for each type of EOS, or reads from an universal EOS file 
using :cpp:func:`~EOS_Toolkit::load_eos_thermal`.

Calling the EOS is a two stage process: first, one needs to specify
the state of matter. This can be done either in terms of mass density,
specific internal energy, and electron fraction, or in terms of 
mass density, temperature, and electron fraction. For this, use the
methods 
:cpp:func:`~eos_thermal::at_rho_eps_ye` and
:cpp:func:`~eos_thermal::at_rho_temp_ye`. 
Both return an object of type
:cpp:class:`~EOS_Toolkit::eos_thermal::state`. 
All EOS provide the first form, which is required for numerical
simulations. The second form based on temperature is optional and 
might not be supported by all EOS implementations.

Second, one can use methods of this object to obtain pressure, 
soundspeed, internal specific energy, temperature, specific entropy, 
derivatives of pressure with respect to mass density and internal 
specific energy. Further, one can use the state object like a boolean 
variable to check if the state is in the valid regime of the EOS. 
Trying to compute any of the above from an invalid state will throw 
an exception. The following snipped demonstrates the EOS use:

.. code-block:: cpp

   auto s = eos.at_rho_eps_ye(rho, eps, ye);
   if (s) {
     auto pressure   = s.press()
     auto soundspeed = s.csnd()
   }
   
The functionality above is also provided by another set of EOS methods 
with names such as :cpp:func:`~eos_thermal::press_at_rho_eps_ye`, 
:cpp:func:`~eos_thermal::press_at_rho_temp_ye`, etc. Those do not throw 
exceptions when called outside the valid ranges, but instead return
NANs. When computing several quantities for the same state, the 
alternative syntax might be less efficient because the validity has to 
be checked and the thermal degree of freedom mapped onto the 
EOS-internal representation each time instead of only once as for the
primary syntax.


The :cpp:class:`~EOS_Toolkit::eos_thermal` class
also provides methods to query the valid ranges of the independent 
variables. For mass density and electron fraction, the validity region
is a simple range. For specific internal energy or temperature,
the valid range depends on mass density and electron fraction.
For convenience, there are methods to check if given parameter
combinations are in the valid region.

.. tip::

   The interface objects are designed to be used as ordinary variables,
   which can be copied around cheaply and without worrying about ownership 
   of any resource allocations. Internally, they just contain a shared 
   pointer to an actual (immutable) implementation.

.. note::

   The EOS interface is thread-safe in the sense that the same EOS 
   oject can be used in multiple threads. Since the EOS implementations
   are immutable, data-races cannot occur. Copying an EOS object inside
   a thread is save but might affect performance since the copying is 
   an atomic operation. Passing the EOS inside an innermost loop should
   better be done by reference instead of value. 

