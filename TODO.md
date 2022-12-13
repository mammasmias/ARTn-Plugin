
## TO-DO List

- nsteppos in ARTn doesn't have the same meaning for QE and LAMMPS in FIRE algo

- Add the output filename custom

- Create option to read a configuration as a reference (energy and position) for the rest of the computation. Kind of `ref_config = file`. Usefule for the refine saddle

- `nperp`: Follows antoine method: progressive increase of nperp after the inflection line. Or maybe to be proportional to the fperp magnitude because happen when the magnitude is too high the perp-relax lead the lost of saddle point. :ok:

  - Done by the routine `nperp_limitation_*()` routines:

  ```fortran
    module artn_param_mod
        logical :: lnperp_limitation
        integer, allocatable :: nperp_limitation(:)
        integer :: def_nperp_limitation(5) = [4,8,12,16,0]
        integer :: noperp, nperp_step

    subroutine initialize_artn()
        [...]
        !! After read artn_params namelist
        call nperp_limitation_init( lnperp_limitation )

    subroutine check_force()
        [...]
        call nperp_limitation_step( 0 ) !! stay in same nperp_step
        [...]
        call nperp_limitation_step( 1 ) !! go to the next nperp_step
        [...]
        call nperp_limitation_step( -1 ) !! return to the first nperp to the list

    ```

- **RESTART**: Fast Restart procedure for lammps and binary - Write the restart file with lammps take too mush time

- Do **time profiler** for the ARTn library.
- **NEXTMIN**: Verify the Threshold to load the new minimum. HARD-CODED!!
- **displacement_validation**: Nico proposed a new way but has to be tested
- Nico: I Have a doubt on the routine fpbc() in pbc.f90
