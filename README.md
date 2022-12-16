# plugin-ARTn

This is a working repository of the current version of the plugin-ARTn; currently it can be used with Quantum ESPRESSO and LAMMPS.
This code has been developped in collaboration by Matic Poberznic, Miha Gunde, Nicolas Salles and Antoine Jay.

The repository is developped on [GiLab](https://gitlab.com/mammasmias/artn-plugin) and a copy of the `master` branch is on [GitHub](https://github.com/mammasmias/ARTn-Plugin). Please put your issue on [GiLab](https://gitlab.com/mammasmias/artn-plugin).

<img src="/home/mammasmias-1/Nico/Project_pARTn/artn-plugin/.extra/ARTn_workflow-1.png" alt="ARTn-Plugin Work Flow" width="400" size="auto" />


## Contains:


- `examples/`: Contains many example including `Al-vacancy.d` and `H2+H.d` 
- `Files_LAMMPS/`: Contains the fix of lammps to interface LAMMPS/ARTn
- `Files_QE/`: Contains the file plugin_ext_forces.f90 which call the ARTn library
- `README.md`: The file you are reading
- `src/`: ARTn plugin subroutines 
- `Makefile`: Command to patch the engine and compile the library. Use the variables defined in file `environment_variables`
- `environment_variables`: User costum file in which it should be define the fortran compiler to compile the library in the variable `F90` and the path where the engine can be found in variables `LAMMPS+PATH` or `QE_PATH`

## Interface with engine

[QE Interface](./Files_QE/README.md)

[LAMMPS-old Interface](./Files_LAMMPS/README-old.md)

[LAMMPS-plugin interface](./Files_LAMMPS/README.md)

## Example

The list of [example](./examples/README.md)


## Issues, bugs, requests

Use the [issue](https://gitlab.com/mammasmias/artn-plugin/-/issues) traker to report the bugs.

## License

[Terms of use](./TERMS_OF_USE)
