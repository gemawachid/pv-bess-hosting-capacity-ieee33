# PV & BESS Hosting Capacity Sizing — IEEE 33-Bus Distribution System

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue?logo=python)](https://www.python.org/)
[![Pyomo](https://img.shields.io/badge/Pyomo-6.x-orange)](http://www.pyomo.org/)
[![Solver](https://img.shields.io/badge/Solver-HiGHS%20%7C%20CBC%20%7C%20Gurobi-green)](https://highs.dev/)
[![License](https://img.shields.io/badge/License-MIT-lightgrey)](LICENSE)
[![Status](https://img.shields.io/badge/Status-Active%20Research-brightgreen)]()

> **Thesis project** — Optimal siting and sizing of rooftop PV and battery energy storage (BESS) on the IEEE 33-bus radial distribution feeder using Mixed-Integer Linear Programming (MILP) with full LinDistFlow (P+Q) and 7-day stochastic load & PV profiles.

---

## Table of Contents

- [Overview](#overview)
- [System Description](#system-description)
- [7-Day Stochastic Profile Model](#7-day-stochastic-profile-model)
- [Mathematical Formulation](#mathematical-formulation)
- [Repository Structure](#repository-structure)
- [Requirements](#requirements)
- [How to Run](#how-to-run)
- [Key Results](#key-results)
- [Debugging Journey](#debugging-journey)
- [Parameters Reference](#parameters-reference)
- [Citation](#citation)

---

## Overview

This repository solves the **PV and BESS hosting capacity problem** on the IEEE 33-bus benchmark distribution network. The optimizer simultaneously decides:

- **Where** to place PV and BESS (bus selection via binary variables)
- **How large** each unit should be (continuous sizing variables)
- **How to dispatch** BESS charge/discharge across 7 days x 24 hours

The objective minimizes **Energy Not Served (ENS)** while keeping all bus voltages and line flows within limits across a full representative week.

---

## System Description

| Parameter | Value |
|-----------|-------|
| Network | IEEE 33-bus radial feeder (Baran & Wu, 1989) |
| Buses | 33 |
| Branches | 32 |
| Nominal voltage | 12.66 kV |
| System base | 10 MVA |
| Impedance base | 16.034 Ohm |
| Peak active load | 3,715 kW |
| Peak reactive load | 2,300 kVAr |
| Time horizon | 7 days x 24 hours (hourly dispatch) |
| Solver | HiGHS 1.x (auto-detects CBC, Gurobi, GLPK as fallback) |

The feeder has three laterals branching from the main feeder:

```
Bus 1 -- 2 -- 3 -- 4 -- 5 -- 6 -- 7 -- ... -- 18   (main feeder)
          |    |         |
         19   23        26
         20   24        27
         21   25        28
         22             29 -- 30 -- 31 -- 32 -- 33
```

> **Weakest bus:** Bus 18 — V = 0.9115 pu at peak load (no DER)
> **Worst power factor:** Bus 30 — Q = 600 kVAr, P = 200 kW, PF = 0.316

---

## 7-Day Stochastic Profile Model

Instead of a single 24-hour profile, this study uses a **7-day stochastic model** with three layers of randomness — all seeded for full reproducibility (SEED_BASE = 42).

### Load profile — 3-layer model

**Layer 1 — Day-type scaling:**
```
lambda_scaled[d,t] = gamma[d] * lambda_base[t - delta[d]]
```

| Day | gamma | Type |
|-----|-------|------|
| Monday | 1.00 | Weekday |
| Tuesday | 1.03 | Weekday |
| Wednesday | 1.05 | Weekday (peak) |
| Thursday | 1.02 | Weekday |
| Friday | 0.98 | Weekday |
| Saturday | 0.86 | Weekend |
| Sunday | 0.80 | Weekend |

**Layer 2 — AR(1) slow drift** (temperature-correlated, tau ~= 6 h):
```
eps_slow[d,t] = 0.75 * eps_slow[d,t-1] + 0.04 * N(0,1)
```

**Layer 3 — Fast appliance noise:**
```
eps_fast[d,t] ~ N(0, 0.013^2)
```

Final profile smoothed with Gaussian kernel (sigma = 1.2 h), clipped to [0.40, 1.15].

### PV profile — 3-state weather model

| State | Label | Probability | CF scale | Noise sigma |
|-------|-------|:-----------:|:--------:|:-----------:|
| s=0 | Clear sky | 0.50 | 1.00 | 0.02 |
| s=1 | Partly cloudy | 0.35 | 0.72 | 0.09 |
| s=2 | Overcast | 0.15 | 0.38 | 0.15 |

---

## Mathematical Formulation

### LinDistFlow voltage drop (Baran & Wu 1989)

```
V[j,d,t] = V[i,d,t] - R_l * P[l,d,t] - X_l * Q[l,d,t]
```

### Two-stage MILP structure

**Stage 1 — planning (day-independent):**
```
xPV[n], xB[n], PPVmax[n], Pchmax[n], Pdismax[n], Emax[n]
```

**Stage 2 — dispatch (per day d, hour t):**
```
Ppv[n,d,t], Pcurt[n,d,t], Pch[n,d,t], Pdis[n,d,t], E[n,d,t]
Pflow[l,d,t], Qflow[l,d,t], V[n,d,t]
```

### BESS state of charge

```
E[n,d,t] = E[n,d,t-1] + (eta_ch * Pch[n,d,t] - Pdis[n,d,t]/eta_dis) * dt
0 <= E[n,d,t] <= Emax[n]       (usable energy window = 20%-90% SOC)
E[n,d,23] = E_init[n,d]        (daily cycle constraint)
```

> Physical rated capacity = Emax[n] / 0.70
> Physical SOC for plotting = 0.20 + (E[n,d,t] / Emax[n]) * 0.70

### Objective
```
min  ENS  =  sum over n,d,t  of  LS[n,d,t] * dt
```

---

## Repository Structure

```
pv-bess-hosting-capacity-ieee33/
├── Distribution-PVBESS-7day.ipynb       <- Main MILP notebook (7-day model)
├── loadprofile.ipynb                    <- Load profile exploration notebook
├── week_profiles.py                     <- Standalone 7-day profile generator
├── pv_bess_hosting_capacity_ieee33.csv  <- Output results after solving
├── figures/
│   ├── fig_basecase_voltage.png         <- Base-case voltage (no DER, peak t=19h)
│   ├── fig_voltage_7day.png             <- Post-solve voltage (all 7 days)
│   ├── fig_soc_7day.png                 <- BESS SOC profiles (all 7 days)
│   └── fig1_weekly_profiles.png         <- 7-day load & PV profiles
└── README.md                            <- This file
```

---

## Requirements

```bash
pip install pyomo numpy pandas matplotlib scipy highspy
```

**Solver — choose one:**

```bash
# Option 1: HiGHS (recommended, free)
pip install highspy

# Option 2: CBC (free)
conda install -c conda-forge coincbc

# Option 3: Gurobi (commercial, free academic licence)
pip install gurobipy
```

**Tested on:**
- Python 3.9+, Pyomo 6.7, HiGHS 1.11.0, Anaconda on Windows 11

---

## How to Run

**1. Clone the repository**

```bash
git clone https://github.com/<your-username>/pv-bess-hosting-capacity-ieee33.git
cd pv-bess-hosting-capacity-ieee33
```

**2. Install dependencies**

```bash
pip install pyomo numpy pandas matplotlib scipy highspy
```

**3. Open and run the main notebook**

```bash
jupyter notebook Distribution-PVBESS-7day.ipynb
```

> Always use **Kernel -> Restart & Run All**.
> Never run only the solver cell — Jupyter caches old model objects in memory.

**4. Standalone profile generator (no MILP needed)**

```bash
python week_profiles.py
```

Generates `fig1_weekly_profiles.png` and `fig2_statistics.png` in the current directory.

---

## Key Results

After solving, the notebook reports:

- **HC_PV** — total PV hosting capacity (kW) across selected buses
- **HC_BESS_E** — usable BESS energy (kWh); multiply by 1/0.70 for physical rated
- **HC_BESS_P** — BESS power capacity (kW)
- **Voltage profiles** — worst-case across all 7 days x 24 hours
- **SOC daily cycles** — for each BESS site across all 7 days

Results exported to `pv_bess_hosting_capacity_ieee33.csv`.

---

## Debugging Journey

Eight critical bugs were identified and fixed. Full report: `project_summary.docx`.

| # | Bug | Symptom | Fix |
|---|-----|---------|-----|
| 1 | Wrong Sbase (1 MVA) | Min voltage ~0.20 pu (should be 0.91) | Sbase_kW = 10000, Zbase = 16.034 Ohm |
| 2 | Bilinear SOC constraints | Presolve: Infeasible at 0 nodes | Redefine Emax as usable energy; E in [0, Emax] |
| 3 | Slack bus in Q-balance | Qflow[branch 1] forced to 0 | Constraint.Skip for bus 1 in P/Q balance |
| 4 | Vmin = 0.95 pu (too tight) | Infeasible without any DER | Vmin = 0.90 pu |
| 5 | No reactive flow modelled | Voltage drops underestimated | Added Qflow vars + full LinDistFlow |
| 6 | No PV curtailment variable | Excess PV forced into network | Added Pcurt[n,d,t] |
| 7 | Wrong Vbase (20 kV) | Zbase = 40 Ohm instead of 16.034 | Vbase_kV = 12.66 |
| 8 | SOC plot formula (E/Emax) | SOC appears to violate bounds | SOC = 0.20 + (E/Emax) x 0.70 |

---

## Parameters Reference

| Parameter | Value | Description |
|-----------|-------|-------------|
| `Sbase_kW` | 10,000 | System MVA base (kW) |
| `Vbase_kV` | 12.66 | Nominal line voltage (kV) |
| `Zbase_ohm` | 16.034 | Impedance base (Ohm) |
| `Vmin` | 0.90 pu | Minimum allowable bus voltage |
| `Vmax` | 1.05 pu | Maximum allowable bus voltage |
| `Npv_max` | 5 | Maximum number of PV sites |
| `Nbess_max` | 3 | Maximum number of BESS sites |
| `PV_cap_max_kW` | 5,000 | Per-site PV capacity ceiling (kW) |
| `BESS_p_max_kW` | 2,000 | Per-site BESS power ceiling (kW) |
| `BESS_e_max_kWh` | 5,000 | Per-site BESS usable energy ceiling (kWh) |
| `SOC_min_frac` | 0.20 | Minimum state of charge |
| `SOC_max_frac` | 0.90 | Maximum state of charge |
| `eta_ch` | 0.95 | Charging efficiency |
| `eta_dis` | 0.95 | Discharging efficiency |
| `N_DAYS` | 7 | Number of representative days |
| `SEED_BASE` | 42 | Random seed (full reproducibility) |
| `dt` | 1 h | Time step |

---

## Citation

```bibtex
@misc{pvbess2026,
  author       = {Aryasatya, Gema Wachid},
  title        = {PV and BESS Optimal Hosting Capacity Sizing on IEEE 33-Bus Distribution System},
  year         = {2026},
  howpublished = {\url{https://github.com/gemawachid/pv-bess-hosting-capacity-ieee33}},
  note         = {project -- LinDistFlow MILP, Pyomo/HiGHS, 7-day stochastic profiles}
}
```

**Reference network:**
> Baran, M. E., & Wu, F. F. (1989). *Network reconfiguration in distribution systems for loss reduction and load balancing*. IEEE Transactions on Power Delivery, 4(2), 1401-1407.

---

## License

MIT License — see [LICENSE](LICENSE) for details.

---

<div align="center">
  <sub>Built with Pyomo · HiGHS · IEEE 33-bus · LinDistFlow · 7-day stochastic profiles · March 2026</sub>
</div>
