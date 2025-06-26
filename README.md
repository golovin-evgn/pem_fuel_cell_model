# pem_fuel_cell_model
Along-the-channel PEM fuel cell model (Uni Magdeburg)

This repository contains the reference implementation of the fuel cell model.
The model was created by Uni Magdeburg, Ievgen Golovin and Christian Kunde, as
part of the publicly funded project "KI Embedded".

## Getting started

Open MATLAB in this root directory. The `startup.m` script should have 
automatically added the necessary paths. Try starting a basic sample run via
```
> script_PEMFC
```
On modern hardware this shouldn't take longer than a couple of minutes.
After that a couple of plots will open to display some characteristic curves
of the PEM fuel cell model.

## Directory structure

`src/` contains all MATLAB source code.
`doc/` contains documentation.
`res/` contains MAT-files, e.g., containing exogeneous inputs.
`bin/` contains binaries, e.g. MEX-files from code generation.

## Authors

- [Ievgen Golovin | ievgen.golovin@ovgu.de](mailto:ievgen.golovin@ovgu.de)
- [Christian Kunde | christian.kunde@ovgu.de](mailto:christian.kunde@ovgu.de)

## Contributors

- [Michael Tiemann | michael.tiemann@de.bosch.com](mailto:michael.tiemann@de.bosch.com)
- [Kilian Knoppe | kilian.knoppe@ovgu.de](mailto:kilian.knoppe@ovgu.de)
