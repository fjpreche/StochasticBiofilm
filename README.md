# StochasticBiofilm
This is a code written in Fortran to simulate the stochastic growth of a biofilm. The code was used to run the simulations presented in 

*Paulina A. Dzianach, Gary A. Dykes, Norval J. C. Strachan, Ken J. Forbes and Francisco J. Pérez-Reche "Unveiling the Mechanisms for Campylobacter jejuni Biofilm Formation Using a Stochastic Mathematical Model" Hygiene (2024).*

The parameters of the model can be set in the code file.

On a Linux machine, the code can be compiled as follows:

$ gfortran -o StochasticBiofilm StochasticBiofilm.f95

and run as

$ ./StochasticBiofilm

The code produces two output files (appropriate names for the files can be specified in the code before compiling):

* **File with information on events** (e.g. data_v3_1_s1_c_6_5_Aero_tmax_12_rls_100_diff_E11_rE_0_5_betty_1000_Kg_6_5.csv)

  Contains the following columns:
    - rls_bcellcount(rls,tstep): Number of cells as a funtion of time
    - rls: Realization number
    - rls_eps(rls,tstep): Number of patches in EPS state
    - rls_events(rls,tstep): Event type (Gro: Growth, Dbs: Death by starvation, - Dbl: Death by lysis, Eps: 
    - rls_ev_loc(rls,tstep): Index i of the patch where the event occurs. A patch -with coordinates (x,y) has an index i = x+(y-1)L.
    - rls_t(rls,tstep): Time of the event. 
    - mean_uptk_c(rls,tstep): Mean uptake of carbon from live cells before the event.
    - mean_uptk_o(rls,tstep): Mean uptake of oxygen from live cells before the event.
    - i_uptk_c(rls,tstep): Carbon uptake rate at the event patch.
    - i_uptk_o(rls,tstep): Oxygen uptake rate at the event patch.
    - i_growth(rls,tstep): Growth rate at the event patch.
    - i_lysis(rls,tstep): Lysis rate at the event patch.
    - i_starvation(rls,tstep): Starvation rate at the event patch.
    - i_conc(rls,tstep): Carbon concentration at the event patch.
    - i_cono(rls,tstep): Oxygen concentration at the event patch.
    - mean_c(rls,tstep): Mean carbon concentration per live cell.
    - mean_o(rls,tstep): Mean oxygen concentration per live cell.
    - Diff(tstep): Effective diffusion coefficient at the time of the event.
    - mean_growth(rls,tstep): Mean growth rate per live cell before an event occurs.
    - mean_lysis(rls,tstep): Mean lysis rate per live cell before an event occurs.
    - mean_starvation(rls,tstep): Mean starvation rate per live cell before an event occurs.
  - waiting_t(rls,tstep): Waiting time for the event to occur.

* **File with final configuration for several stochastic realizations** (e.g. fbiopic_v3_1_s1_c_6_5_Aero_tmax_12_rls_100_diff_E11_rE_0_5_betty_1000_Kg_6_5.csv).

  This may be useful for plotting the state of the biofilm at the final time. It contains two columns:

  - final patch state of each patch (Cell: 1, EPS: 0.5, Deactivated: 0.25).
  - stochastic realization number.
 
The *Examples_Results* folder contains the files obtained for 100 stochastic realizations for a system of size 20x50, $S_c^H = 6.5$, $S_o^H = 0.26$ mmol L$-1$ (aerobic conditions), effective diffusion coefficient in the biofilm $D_B = 9 \times 10^{-11}$ m$^2$ h$^-1$.
