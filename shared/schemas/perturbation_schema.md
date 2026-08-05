# Perturbation Schema (Module 2 interface) — DRAFT for Week 1 discussion

Perturbation.apply(circuit_model) -> modified_circuit_model

Subclasses:
- GeneKnockout(gene_name)
- InducerLevelChange(inducer_name, new_level)
- CarbonSourceShift(from_source, to_source)

Open questions for Student B:
- Does Module 4 expect Perturbation objects directly, or a serialized config (e.g. YAML)?
- How does Module 4 pass a perturbation to hundreds/thousands of simulated cells at once — one Perturbation per population run, or per-cell variation?
