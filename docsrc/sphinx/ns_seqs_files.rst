.. _ns_seqs_files:

Sequence Files
--------------

The library defines a hdf5-based file format for storing star sequences and branches.
There is no need for users to deal with those files directly since there is are C++ 
and Python interfaces for that.


To avoid unit confusion, the files use SI units internally. When saving a sequence,
conversion happens automatically using the unit system stored within the sequence.
When loading, one has to specify the unit system, which is assumed to be geometric.


For saving a sequence to a hdf5 file, use
:cpp:func:`~EOS_Toolkit::save_star_seq`, specifying the path 
and the sequence. For saving a branch, use
:cpp:func:`~EOS_Toolkit::save_star_branch` instead.

To load sequences or branches, use
:cpp:func:`~EOS_Toolkit::load_star_seq` or
:cpp:func:`~EOS_Toolkit::load_star_branch`, respectively.

The following example demonstrates basic use of stable TOV branches.


.. literalinclude:: minimal_seq_file.cc
   :language: cpp


