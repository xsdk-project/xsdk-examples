# deal.II-preCICE example

A simple example on how to use deal.II with preCICE: a Laplace equation implemented in deal.II is coupled to a stand-alone C++-only boundary condition by calling preCICE in both codes. This example is taken from the [deal.II code gallery](https://github.com/dealii/code-gallery/tree/master/coupled_laplace_problem) and was contributed there by @davidscn.

- `boundary_condition`: Boundary condition calling preCICE for coupling
- `laplace_problem.cc`: Laplace equation implemented in deal.II and calling preCICE for coupling
- `precice-config.xml`: preCICE configuration file, read by both codes
- `clean`: Script to clean temporary, result, and log files

To run the example, run both codes from the level where `precice-config.xml` is.

```txt
./laplace_problem & ./boundary_condition
```
