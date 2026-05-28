# Robust Control of a Hydropower Turbine with a Long Penstock

This repository contains MATLAB files for modelling, simulation, and robust
control of a hydropower plant with a long penstock, developed as part of a
Master’s thesis.

## Main Simulation Scripts

The scripts with the prefix `HP_` are the main simulation scripts.

* `HP_hpp_cycle_nonlin.m`
  Simulates the full operating cycle of the hydropower plant using the nonlinear model with a constant cross-section penstock.

* `HP_hpp_cycle_geom.m`
  Compares the full operating cycle for the constant cross-section penstock model and the segmented penstock model.

* `HP_lin_nonlin_compare.m`
  Compares the linearized model with the nonlinear model for the constant cross-section penstock.

* `HP_weight_functions_moc.mlx`
  Computes model uncertainties and the weighting functions (W_1) and (W_3).

## Controller Synthesis Scripts

The scripts with the prefix `tds_` are used for controller synthesis. These scripts require the TDS-CONTROL Toolbox.

* `tds_design_island_moc.mlx`
  Performs robust PID controller synthesis for island operation.

* `tds_design_sync_moc.mlx`
  Performs robust PID controller synthesis for synchronous operation.

## Supporting Functions

The remaining files are supporting functions and parameter files used by the main simulation and controller synthesis scripts.

* `MOC_step.m`
  Method of Characteristics step for the constant cross-section penstock model.

* `MOC_step_geom.m`
  Method of Characteristics step for the segmented penstock model.

* `mechanics_step_island.m`
  Mechanical dynamics of the turbine-generator unit in island operation.

* `mechanics_step_sync.m`
  Mechanics of the turbine-generator in synchronous operation.

* `PID_parallel.m`
  Discrete-time parallel PID controller implementation.

* `params.m`
  Main model and simulation parameters.

* `params_pid.m`
  PID controller parameters used in the simulations.

* `steady_equations.m`
  Steady-state equations used for model initialization.

## Model Variants

The constant cross-section penstock model is used as the main model for linearization, controller design, and numerical simulations.

The segmented penstock model is used to assess the influence of penstock geometry on the closed-loop response.

