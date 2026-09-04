# QUARANTINE — NOT_FOR_CLAIM (adversarial-review fix 3, 2026-09-04)

The following M3-01 outputs are QUARANTINED and must not be cited in any claim:

- `M3-01_tRNA_supply.tsv`
- `M3-01_tRNA_charging.tsv` (empty — periodate filename-regex gap; superseded by L1c)
- `M3-01_tRNA_qc.tsv` (where it relates to supply/charging)

## Reasons (adversarial review B3)

1. **Reference namespace mismatch**: GSE226632 tRNA count features use an `HS_`-prefixed
   namespace (human-style tRNA gene models), not a P. falciparum-native tRNA reference.
   Family-level (AA-anticodon) aggregation is plausible but not validated against a
   Pf-native tRNA annotation.
2. **Charging table is empty**: the periodate-derived charging measurement never produced
   data in M3-01 (regex gap), so "charging" claims from this directory are vacuous.

## Status of the science without these tables

- L1c (`data/derived/WP4/L1c_trna_charging/`) re-derived family-level abundance and
  stage-resolved charging from the same raw tar with correct parsing; it carries the same
  namespace caveat and is labeled QUARANTINED-REF / descriptive-only in its claim_impact.
- No project claim (CLM registry) currently depends on the quarantined tables.
