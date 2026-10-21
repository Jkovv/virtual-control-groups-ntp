# Analysis plan: three improvements to the virtual control group

Written and committed before the comparison was run. The numeric settings are
mirrored in the "method comparison" block of `src/vcg/config.py`. Any change made
after results were seen is listed under *Deviations* at the end, with the reason.

## Question

Do three literature-based changes to how the virtual control group (VCG) is built
or analysed bring its conclusions closer to those of the concurrent control
group? "Closer" is always measured against the concurrent-control conclusion,
never against a ground truth.

## Data split

* Unit of analysis: a target study x sex; its concurrent control group is held out.
* **Development set**: target studies that started before **1990** (1,174 of the
  1,558 studies that contribute any endpoint).
* **Held-out set**: target studies that started in **1990 or later** (384 studies,
  24.6%). This is the final evaluation. Its results are reported as they come out.
* The split partitions *targets* only. Historical pools are built exactly as in
  the current method: any other study, excluding the same test article, matched
  at level L3 (species, strain, sex, design, route, vehicle, laboratory). A
  held-out target can therefore borrow controls from studies that started later.
  This is unchanged from the current method, to keep the comparison paired.
* No setting is tuned on either set. The development set is reported to show
  whether conclusions are stable across periods.

## Endpoint change that applies to every method, including the current one

The pre-dosing body weight (in-life "Body weight, initial") becomes a matching
covariate. It is removed from the evaluated endpoints for all methods. Body
weight *gain* from pre-dosing to week 13 / week 104 stays an endpoint.

## Methods

All methods use matching level L3, 20 draws where sampling is involved, alpha
0.05, and the seed in `config.py`.

| Id | What changes | Applies to |
|---|---|---|
| `baseline` | Current method: individual historical animals sampled uniformly; Dunnett, Fisher, poly-3 | all domains |
| `bw_matched` | Pool restricted to animals whose pre-dosing weight is within +/- 2 SD of the target study/sex mean; sampling weights proportional to the normal density of the target's pre-dosing weights (Gurjanov et al. 2024, weighted sampling). Same tests as baseline. | all domains with pre-dosing weights |
| `bw_ancova` | Same VCG as `bw_matched`; analysis is ANCOVA with pre-dosing weight as covariate, Dunnett-adjusted. The reference conclusion is the same ANCOVA with the concurrent control. | terminal body weight, in-life body weight, organ weights |
| `whole_group` | Whole historical control groups (one study x sex control group at a time, in random order) are added until the VCG has at least as many animals as it replaces; all animals of the chosen groups are used. Same tests as baseline. | all domains |
| `hierarchical` | Continuous endpoints. Full VCG: random-effects model of the matched historical control group means (DerSimonian-Laird); each treated group is compared with the meta-analytic-predictive (MAP) distribution of a new study's control mean, N(mu, tau^2 + se(mu)^2); Dunnett-type adjustment using the implied product correlation. Partial VCG: robust MAP prior, 0.8 x MAP + 0.2 x unit-information normal (Schmidli et al. 2014), updated with the kept half of the concurrent controls; treated groups are compared with the posterior of the control mean. No sampling for the full VCG (one deterministic result per endpoint). | continuous domains |
| `baseline_trend` | Neoplasms. Cochran-Armitage trend test (one-sided, increase) with the baseline VCG as control. | neoplasms |
| `tarone` | Neoplasms. Tarone (1982) trend test with the historical control rates as a beta-binomial prior, fitted by maximum likelihood to the matched historical control groups. Full VCG: prior only. Partial: prior plus the kept concurrent half. | neoplasms |

Reference conclusion for each method:

* `baseline`, `bw_matched`, `whole_group`, `hierarchical`: Dunnett / Fisher / poly-3
  with the concurrent control, as now.
* `bw_ancova`: ANCOVA with the concurrent control.
* `baseline_trend`, `tarone`: Cochran-Armitage with the concurrent control
  (Tarone's test with no historical information).

Requirements for a method to be evaluable on a unit x endpoint (otherwise it
is counted as not covered, never silently replaced):

* sampling methods: pool of at least twice the control size (full VCG) or at
  least the virtual part (partial), as now;
* `bw_matched`, `bw_ancova`: at least 5 pre-dosing weights in the target study/sex,
  pre-dosing weights for the pool animals used, and the size rule above after
  the +/- 2 SD window;
* `whole_group`, `hierarchical`, `tarone`: at least 3 historical control groups.

## Metrics

Exactly as in the current analysis, at endpoint level (an endpoint is significant
if any dose group is; trend tests are endpoint-level by construction): agreement,
Cohen's kappa, the "never significant" baseline, sensitivity, specificity, false
positive and false negative shares, and the concurrent-versus-VCG check:

* sampling methods: Welch test of the concurrent control against each drawn VCG;
* `bw_ancova`: the same comparison as a two-group ANCOVA;
* `hierarchical`: z test of the concurrent control mean against the MAP
  predictive distribution (is the real control inside the predicted range?).

## Comparison rule ("helped", "hurt", "no clear difference")

For each method, set and endpoint domain, the method is compared with `baseline`
(full with full, partial with partial; `tarone` with `baseline_trend`) **on the
same unit x endpoint combinations**. The kappa difference gets a 95% interval
from a study-cluster bootstrap with 1,000 resamples.

* **helped**: kappa gain of at least 0.05 and the interval excludes 0;
* **hurt**: kappa loss and the interval excludes 0;
* **no clear difference**: otherwise.

Changes in sensitivity, specificity and the false positive share are reported
next to each verdict, because a kappa gain can come from fewer false positives at
the cost of missed findings.

## Deviations

None so far.
