# plugin-ARTn

This is a working repository of the current version of the plugin-ARTn; currently it can be used with Quantum ESPRESSO and LAMMPS.
This code has been developed in collaboration by Matic Poberznik, Miha Gunde, Nicolas Salles and Antoine Jay.

The repository is developed on [GitLab](https://gitlab.com/mammasmias/artn-plugin) and a copy of the `master` branch is on [GitHub](https://github.com/mammasmias/ARTn-Plugin).
Please post your issue(s) on [GitLab](https://gitlab.com/mammasmias/artn-plugin).

<img src="./.extra/ARTn_workflow-1.png" alt="ARTn-Plugin Work Flow" width="400" size="auto" />


## Contains:


- `examples/`: Contains many examples, from molecules to surfaces;
- `Files_LAMMPS/`: Contains the lammps fix, for the LAMMPS/ARTn interface;
- `Files_QE/`: Contains the file `plugin_ext_forces.f90`, for the QE/ARTn interface;
- `README.md`: This file;
- `src/`: ARTn plugin subroutines;
- `Makefile`: Compilation commands, uses the environment variables defined in `environment_variables`;
- `environment_variables`: custom file defining:
    - compilers `F90`, `CXX`/`CC`;
    - the paths to ARTn (current directory), BLAS and FORTRAN libraries, and chosen engine(s) QE/LAMMPS;
    - `PARA_PREFIX` prefix for launching examples via provided run scripts.

## Interface with engine

Two interfaces has been developed for the moment:

- One for **Quantum ESPRESSO**. To use it read the [manual](./Files_QE/README.md).
- One for **LAMMPS**. Two version exist, one using the class [Plugin](https://docs.lammps.org/plugin.html) of LAMMPS, for this version please read the [manual](./Files_LAMMPS/README.md); The second one does not use the class plugin of LAMMPS because this class exist only since 2022. If you use a version older than 2022 please read the [manual](./Files_LAMMPS/README-old.md)

## Examples

The list of [examples](./examples/README.md) using both interfaces.


## Using ARTn

To customise the input of ARTn please read the [manual](./MANUAL.md)

## Issues, bugs, requests

Use the [issue](https://gitlab.com/mammasmias/artn-plugin/-/issues) tracker to report the bugs.

## License

[Terms of use](./TERMS_OF_USE)
