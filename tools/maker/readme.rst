Running with MPI
================

By default, Maker will only run on 1 cpu. If you want to parallelize on any number
of cpus and nodes, you can enable MPI by setting the ``MAKER_MPI`` variable to 1
in the job destination, and setting ``GALAXY_SLOTS`` to the desired number of cores
(see https://galaxyproject.org/admin/config/galaxy_slots/#galaxy_slots-for-server-admins).
