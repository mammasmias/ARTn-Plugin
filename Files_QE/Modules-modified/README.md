# A summary of changes to Modules files of QE  
  - `input_parameters.f90`: [**fire**] 

    **line 1035**: add `fire` keywords to `ion_dynamics_allowed`

    ```fortran
    CHARACTER(len=80) :: ion_dynamics_allowed(11)
     DATA ion_dynamics_allowed / 'none', 'sd', 'cg', 'langevin', &
                                'damp', 'verlet', 'bfgs', 'beeman',&
                                'langevin-smc', 'ipi', 'fire' 
    ```
     
    **lines 1185-1203**: added the default values of FIRE global variables and to the ions namelist

    ```fortran
         INTEGER  :: fire_nmin = 5 ! minimum number of steps for time step increase 
         REAL(DP) :: fire_f_inc = 1.1_DP ! factor for time step increase  
         REAL(DP) :: fire_f_dec = 0.5_DP ! factor for time step decrease
         REAL(DP) :: fire_alpha_init = 0.1_DP ! initial value of mixing factor
         REAL(DP) :: fire_falpha = 0.99_DP ! modify the mixing factor
         REAL(DP) :: fire_dtmax = 10.0_DP ! maximum time step; calculated as dtmax = fire_dtmax*dt 
         !
         NAMELIST / ions / ion_dynamics, iesr, ion_radius, ion_damping,         &
                           ion_positions, ion_velocities, ion_temperature,      &
                           tempw, fnosep, nhgrp, fnhscl, nhpcl, nhptyp, ndega, tranp,   &
                           amprp, greasp, tolp, ion_nstepe, ion_maxstep,        &
                           refold_pos, upscale, delta_t, pot_extrapolation,     &
                           wfc_extrapolation, nraise, remove_rigid_rot,         &
                           trust_radius_max, trust_radius_min,                  &
                           trust_radius_ini, w_1, w_2, bfgs_ndim,l_mplathe,     &
                           n_muller,np_muller,l_exit_muller,                    &
                           fire_nmin, fire_f_inc, fire_f_dec, fire_alpha_init,  &
                           fire_falpha, fire_dtmax 

    ```
    

  - `io_files.f90`: [**fire**] 

     **line 219**: add the test on keyword `fire` (maybe a file extension):

     ```fortran
     CALL delete_if_present( trim( file_path ) // '.fire' )
     ```
  - `read_namelists.f90`: [**fire**]
     **lines 547-552**: added the defaults for fire global variables
     ```fortran
        fire_nmin = 5 ! minimum number of steps P > 0 before dt incread
        fire_f_inc = 1.1_DP ! factor for time step increase
        fire_f_dec = 0.5_DP ! factor for time step decrease 
        fire_alpha_init = 0.1_DP ! initial value of mixing factor
        fire_falpha = 0.99_DP ! modification of the mixing factor
        fire_dtmax = 10.0_DP ! factor for calculating dtmax 

     ```
