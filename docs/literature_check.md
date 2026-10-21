# Literature check: is a split-half or noise-ceiling benchmark already used for virtual control groups?

Searched 2026-10-21, before the Contributions section of the README was written.
Scope: ALTEX, Frontiers in Pharmacology, Frontiers in Toxicology, Toxicologic
Pathology, Regulatory Toxicology and Pharmacology, Toxicological Sciences, the
eTRANSAFE and VICT3R outputs, and the EMA draft qualification opinion. Sources:
PubMed (`"virtual control group" toxicology`, 8 records), Crossref, and the full
text of the three documents that define the field's evaluation practice (the EMA
opinion, Andaya et al. 2024, Golden et al. 2024), searched for *split-half*,
*two halves*, *noise ceiling*, *upper bound*, *attainable*, *test-retest*,
*disjoint* and *reliability*.

## What the field uses to evaluate a virtual control group

| Study | Evaluation against the concurrent control | Ceiling or attainability benchmark |
|---|---|---|
| Wright et al. 2023, *Regul Toxicol Pharmacol* 138:105309 | balanced accuracy of treatment-relatedness calls against the original toxicologist's designation; covariate matching varied | none |
| Andaya et al. 2024, *ALTEX* 41(4):647-659 | percentage agreement of statistical calls, 100 VCGs per study | none |
| Gurjanov et al. 2024, *ALTEX* 41(4):660-665 | reproduction of body-weight outcomes under three sampling schemes | none |
| Golden et al. 2024, *ALTEX* 41(2):282-301 (workshop roadmap) | discusses what evidence regulators would need | none proposed |
| Steger-Hartmann et al. 2025, *ALTEX*, doi 10.14573/altex.2503061 | replicability framing across data sources | none |
| Duchateau-Nguyen et al. 2026, *Front Pharmacol* 17:1704002 | share of VCGs inside the concurrent effect-size CI; share with matching significance status; 1,000 draws | none; "high" agreement is set at 90% by convention |
| Kim et al. 2026, *Front Toxicol* 8:1888158 | five-category reproducibility by direction, magnitude and significance; 100 iterations, majority vote | none; states that "completely eliminating inter-study variance is practically unattainable" without quantifying it |
| EMA, *Draft qualification opinion for VCG to replace CCG in rat non-GLP DRF studies*, 31 March 2026 | context of use, matching criteria, expert review per study | none |

Also checked and not using one: Adedeji & Naor 2025 (*Toxicol Pathol*), Lotfi et
al. 2026 (*Regul Toxicol Pharmacol*, re-analysis of 20 studies), Chandra et al.
2026 (*Toxicol Sci*, generative synthetic controls).

**Finding: no study in this literature splits the concurrent control group to
measure how much agreement is attainable, and none reports any other ceiling.**
Agreement is reported against an implicit target of 100%, or against a
conventional threshold such as 90%.

## The closest prior work

* **Kim et al. 2026** come nearest in words: inter-study variance cannot be
  eliminated, so perfect reproduction is not expected. They do not quantify it.
* **Dertinger et al. 2023**, the IWGT report on historical negative control data
  (*Environ Mol Mutagen* 66:144-165, doi 10.1002/em.22541), decomposes historical
  control variability into animal and study components and recommends using
  historical controls only where animal rather than study factors dominate. This
  is the same decomposition this repository applies, in a different setting
  (genetic toxicology control limits, not virtual control groups) and without a
  link to an error rate.
* Outside toxicology the two ideas are standard. Split-half reliability goes back
  to Spearman and Brown (*Br J Psychol* 3, 1910); the noise ceiling as an upper
  bound on model-data agreement is routine in computational neuroscience, for
  example Nili et al. 2014 (*PLoS Comput Biol* 10:e1003553) and the comparison of
  estimators in Lage-Castellanos et al. 2019 (*PLoS Comput Biol* 15:e1006397).

## Conclusion carried into the README

The benchmark is not new as a method; it is, as far as this search goes, new to
the evaluation of virtual control groups. The README states it that way and does
not claim the method itself. Anyone who finds prior use of it in this field is
asked to open an issue so the claim can be corrected.
