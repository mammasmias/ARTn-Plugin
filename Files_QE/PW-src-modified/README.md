# A summary of changes to the PW/src/  files of QE: 

  - `dynamics_modules.f90`: [**fire**]

    **line 47, lines 76-85**:  global variables of the fire minimization algorithm

    ```fortran
    PUBLIC :: fire_nmin, fire_f_inc, fire_f_dec, fire_alpha_init, fire_falpha, fire_dtmax
    ```

    and
    ```fortran
    !! 
    INTEGER ::  fire_nmin
    !! FIRE: minimum number of steps for time step increase 
    REAL(DP) :: fire_f_inc
    !! FIRE: factor for time step increase  
    REAL(DP) :: fire_f_dec
    !! FIRE: factor for time step decrease
    REAL(DP) :: fire_alpha_init
    !! FIRE: initial value of mixing factor
    REAL(DP) :: fire_falpha
    !! FIRE: modify the mixing factor
    REAL(DP) :: fire_dtmax
    ```
    **lines 1009-1277**: add the subroutine fire

    ```fortran
    SUBROUTINE fire( conv_ion )
      LOGICAL, INTENT( INOUT ) :: conv_ion
    END SUBROUTINE fire
    ```

  - `move_ions.f90`: [**fire**]

    **line 48**: add access to module variable `fire` in `dynamics_module`

    ```fortran
    USE dynamics_module,        ONLY : verlet, terminate_verlet, proj_verlet, fire
    ```
    **lines 50-51**: add access to global variables of `fire` in `dynamics_module`

    ```fortran
    USE dynamics_module,        ONLY : fire_nmin, fire_f_inc, fire_f_dec, &
                                    fire_alpha_init, fire_falpha, fire_dtmax
    ```

     **lines 277:286**: add the test on `calc` keyword for fire algorithm 

    ```fortran
    ELSEIF ( calc == 'fi' ) THEN
       !
       CALL fire( conv_ions)
       ! 
       IF ( .NOT. conv_ions .AND. idone >= nstep ) THEN
          WRITE( UNIT = stdout, FMT =  &
            '(/,5X,"The maximum number of steps has been reached.")' )
          WRITE( UNIT = stdout, &
             FMT = '(/,5X,"End of FIRE minimization")' )
       ENDIF
    ```

    

  - `input.f90`: [**fire**]
  
    **lines 53-63**: add fire global variables to the USE dynamics_module statement;
    ```fortran
      USE dynamics_module, ONLY : control_temp, temperature, thermostat, &
                              dt_         => dt, &
                              delta_t_    => delta_t, &
                              nraise_     => nraise, &
                              refold_pos_ => refold_pos, &
                              fire_nmin_ => fire_nmin, &
                              fire_f_inc_ => fire_f_inc, &
                              fire_f_dec_ => fire_f_dec,  &
                              fire_alpha_init_ => fire_alpha_init, &  
                              fire_falpha_ => fire_falpha, &
                              fire_dtmax_ => fire_dtmax 
    ```
    **lines 294-301**:  added fire global variables to USE
      input_parameters, so that they are read correctly from the input
      namelist

    ```fortran
     USE input_parameters, ONLY : ion_dynamics, ion_positions, tolp, &
                               tempw, delta_t, nraise, ion_temperature,        &
                               refold_pos, remove_rigid_rot, upscale,          &
                               pot_extrapolation,  wfc_extrapolation,          &
                               w_1, w_2, trust_radius_max, trust_radius_min,   &
                               trust_radius_ini, bfgs_ndim, &
                               fire_nmin, fire_f_inc, fire_f_dec, &
                               fire_alpha_init, fire_falpha, fire_dtmax
    ```
     **lines 407-417** : overwrite the fire variables of the dynamics_module with values from input_parameters

     ```fortran
        CASE ( 'fire' )
        !
        lmd     = .true.
        calc    = 'fi'
        ! set fire variables
        fire_nmin_  = fire_nmin
        fire_f_inc_  = fire_f_inc
        fire_f_dec_  = fire_f_dec
        fire_alpha_init_  = fire_alpha_init
        fire_falpha_  = fire_falpha
        fire_dtmax_  = fire_dtmax
        !
        ntcheck = nstep + 1
     ```
   
  - `plugin_ext_forces.f90`: [**ARTn**] added a call to the artn subroutine   

    ```fortran
     USE ions_base,     ONLY : nat, tau, if_pos, ityp, atm, amass
     USE cell_base,     ONLY : alat, at
     USE force_mod,     ONLY : force
     USE ener,          ONLY : etot 
     USE relax,         ONLY : epsf, epse
     USE control_flags, ONLY : istep
     USE dynamics_module, ONLY : vel, dt, fire_alpha_init
     USE io_files,      ONLY : prefix,tmp_dir
     !
     IMPLICIT NONE
     ! 
     LOGICAL :: lconv
     !
     ! ARTn convergence flag 
     ! 
     lconv = .false. 
     !
     IF ( ionode ) THEN
        CALL artn(force,etot,epsf,nat,ityp,atm,tau,at,alat,istep,if_pos,vel,dt,fire_alpha_init,lconv,prefix,tmp_dir) 
     ENDIF
     IF ( lconv ) THEN
        WRITE (*,*) "ARTn calculation converged, stopping" 
        STOP 1
     END IF
     
    ```
