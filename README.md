# ARTn-plugin-QE

This is a working repo for the current version of the ARTn plugin  

## Contain:

- `Modules-modified/`: Modification on QE's Modules

  - `input_parameters.f90`: [**fire**] 

    **line 1035**: add `fire` keywords to `ion_dynamics_allowed`

    ```
    CHARACTER(len=80) :: ion_dynamics_allowed(9)
            DATA ion_dynamics_allowed / 'none', 'sd', 'cg', 'langevin', &
                                        'damp', 'verlet', 'bfgs', 'beeman',&
                                        'langevin-smc' /
    ```

    to 

    ```
    CHARACTER(len=80) :: ion_dynamics_allowed(11)
    DATA ion_dynamics_allowed / 'none', 'sd', 'cg', 'langevin', &
                                'damp', 'verlet', 'bfgs', 'beeman',&
                                'langevin-smc', 'ipi', 'fire' /
    ```

  - `io_files.f90`: [**fire**] 

    **line 219**: add the test on keyword `fire` (maybe a file extension):

    ```
    CALL delete_if_present( trim( file_path ) // '.fire' )
    ```

     

- `PW-src-modified/`: Routine relative to plugin_ext_forces routine of QE

  - `dynamics_modules.f90`: [**fire**] 

    **lines 1009:1277**: add the subroutine fire

    ```
    SUBROUTINE fire( conv_ion )
      LOGICAL, INTENT( INOUT ) :: conv_ion
    END SUBROUTINE fire
    ```

  - `move_ion.f90`: [**fire**]

    **line 48**: add access to module variable `fire` in `dynamics_module`

    ```
    USE dynamics_module,        ONLY : verlet, terminate_verlet, proj_verlet, fire
    ```

     **lines 277:286**: add the on `calc` keyword

    ```
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

    

  - `input.f90`: **???**

  - `forces.f90`: [**ARTn**] change the call interface of plugin_ext_forces()

    ```
    call plugin_ext_forces()
    ```

    to 

    ```
    call plugin_ext_forces( forces(:,:) ) !! forces in input/output
    ```

    **Has to be deleted**

  - `plugin_ext_forces.f90`: [**ARTn**] Change the interface of plugin_ext_forces subroutine

    ```
    SUBROUTINE plugin_ext_forces()
    END SUBROUTINE plugin_ext_forces
    ```

    to 

    ```
    SUBROUTINE plugin_ext_forces( forces )
      REAL(DP), INTENT( INOUT ) :: forces(3,nat)
    END SUBROUTINE plugin_ext_forces
    ```

    **Has to be deleted**: forces is global variable in `MODULE force_mod` in  `QE/PW/src/pwcom.f90` Don't need to give it in argument

- `examples/`: Contains `Al-vacancy.d` and `H2+H.d` 

- `qe-6.6/`: Qantum ESPRESSO software

- `patch-plugin-ext-forces.sh`: Update script

- `patch.sh`: Install script

- `README.md`: What you read!!

## How to use:

**first installation:**  Execute the installation script `patch.sh`

**Udpate ARTn-plugin**: Execute the update script `patch-plugin-ext-forces.sh`