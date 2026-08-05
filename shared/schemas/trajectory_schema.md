# Trajectory Schema (Module 1 output) — DRAFT for Week 1 discussion

Format: Parquet file

Columns:
- time (float, simulation time units)
- species_1_count, species_2_count, ... (int, molecule counts per species per timestep)

Open questions for Student B:
- What timestep resolution does the population simulator need?
- Do we need per-cell trajectories or per-species-averaged?
- Should species names be fixed columns or a long-format table (time, species, count)?
