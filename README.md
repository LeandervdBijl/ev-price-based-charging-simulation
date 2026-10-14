# Energy Price Coordination

This repository contains the simulation code used for the research presented in [our paper](https://ieeexplore.ieee.org/document/11305459).

The code simulates the charging behaviour of individual electric vehicle (EV) owners, assuming that they make rational charging decisions based on electricity prices, subject to vehicle availability and charging requirements.

## Usage

The main simulation script is `run_sim.py`. It can be used to simulate EV charging behaviour based on vehicle availability and charging demand.

The dataset used for the simulations presented in the paper is publicly available through [4TU.ResearchData](https://data.4tu.nl/datasets/80ef3824-3f5d-4e45-8794-3b8791efbd13/1).

The FOCS and Price files that are necessary are obtained from: https://github.com/lwinschermann/FlowbasedOfflineChargingScheduler

## Notes

This repository contains the research code as developed for the paper. It is provided primarily for transparency and reproducibility, rather than as a general-purpose software package. The code has not been extensively refactored or documented for external use.