Subroutine CB_system_energy_GPU(state)

   Use set_precision
   Use Cell
   Use configuration
   Use potential
   Use rundata
   Use utils
   Use properties

   Integer      ::  i, j, k, l, indx, ngpus
   Integer      ::  nmol,state
   Real(wp)     ::  energy(Natoms), Ener_tot


   call computeenergygpu(bool_tor, state, Ndim, Nat_cell_max, Ntot_Cell, Npart_types, &
                     &   Nrb_sites_max, Ntor_max, Nrb_sites, patch, ref_tor_vec, &
                     &   sigma_jon, sigma_tor_jon, sigma_LJ, rangeP, VL0, xop, &
                     &   Ntor_angle, tor_angle, Vpot_Matrix, &
                     &   h, Sphere_Cell, List_Cell, Map_Cell, Nl_Cell, Nat_Cell, DR_Cell, Ener_tot)
!!
   En_tot = Ener_tot

End Subroutine CB_system_energy_GPU

