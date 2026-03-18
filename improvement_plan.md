# Improvement Plan: "High Secret Key Rate QKD Protocols"

## Current state of the manuscript

- Well-defined problem: optimize worst-case SKR $R_\epsilon$ over prepare-and-measure
  QKD protocols under a Choi-fidelity noise model.
- Clean SDP formulation for $\mathcal{P}_\epsilon(S)$, $\mathcal{P}_\epsilon(B=A|S)$,
  $\mathcal{P}_\epsilon(E=P|S)$ (Appendix B).
- Numerical results for $N \in \{1,2,3,4\}$ via Artificial Bee Colony heuristic.
- Key finding: rate-vs-robustness tradeoff (N=1 peaks high, N=4 is noise-robust).

### Gaps common to both target journals

| Gap | Severity |
|-----|----------|
| Empty References section | Critical |
| Section B.4 ($R^*_\epsilon$) is literally "?" | Critical |
| AI-generated disclaimer in Sec. 6 | Must be removed; rewrite those paragraphs |
| No comparison with known analytical SKR bounds (Devetak-Winter, Shor-Preskill, Scarani-Renner) | Critical |
| No discussion relating Choi-fidelity model to standard QBER or device-independent settings | Major |
| No convergence/validation analysis of the ABC heuristic | Major |
| Upper bound curve (brown line in Fig. 1 of `preliminary_results.pdf`) is plotted but never defined or proved | Major |
| Writing quality needs significant tightening throughout | Moderate |

---

## Route 1: Quantum (quantum-journal.org)

**Bar:** Clear, significant advance; rigorous (not necessarily purely analytical);
strong motivation for the broad quantum information community; polished writing.

### Minimum viable paper

1. **Analytical result for N=1**
   - The code (`generate_three_case_curves.jl`, `build_case1_protocol`) already
     parametrizes N=1 as a 1D family: $|\psi_{1,0}\rangle = |0\rangle$,
     $|\psi_{1,1}\rangle \propto (1, x)^T$, with $q_1$ free.
   - Derive a closed-form expression for $R^*_\epsilon$ in the N=1 case
     (or at least tight analytical bounds).
   - This would be the central theorem of the paper.

2. **Prove or tightly characterize the upper bound**
   - The brown "upper" curve in `preliminary_results.pdf` must be defined,
     proved, and cited. If it comes from a known result (e.g., Takeoka-Guha-Wilde
     or Pirandola-Laurenza-Ottaviani-Banchi), state so explicitly.
   - If it is a new bound, prove it rigorously.
   - Show where your numerical optima sit relative to this bound; any tightness
     result (even for a single $\epsilon$ value) is very strong.

3. **Relate the noise model to standard settings**
   - Show that $\Phi \sim \epsilon$ (Choi fidelity $\geq 1-\epsilon$) implies a QBER
     bound and vice versa for specific protocols.
   - Discuss how this compares to device-independent / semi-device-independent
     assumptions. The model constrains the *channel* rather than the *state*,
     which is a meaningful distinction worth highlighting.

4. **Structural insight into optimal protocols**
   - Explain *why* optimal N=1 protocols use nearly orthogonal states at low
     noise (connection to accessible information / Holevo bound).
   - Explain the structure of optimal N=4 protocols: do the pairs approximate
     a SIC-POVM or MUB structure? Are the $q_i$ values interpretable?

5. **Thorough literature review**
   - BB84, B92, six-state: Shor-Preskill, Lo-Chau, Renner's thesis.
   - Key rate optimization: Coles-Metodiev-Lutkenhaus framework,
     Winick-Lutkenhaus-Hoi-Kwong numerical methods.
   - Choi-fidelity / channel fidelity in QKD: any prior use of this model.
   - Metaheuristic optimization in QKD protocol design (if any prior work).

6. **Heuristic validation**
   - Multiple independent restarts with different seeds; report variance.
   - For N=1, compare ABC output against the analytical optimum (from item 1).
   - For small $\epsilon$, verify against known BB84/B92/six-state analytical rates.

7. **Writing**
   - Remove AI disclaimer; rewrite Sec. 6 from scratch.
   - Consolidate Secs. 6.5-6.7 (parameter definitions) into Sec. 1 where the
     protocol is defined.
   - Add an Introduction with motivation (why this noise model? what gap in
     the literature?) and a Conclusion.

### Stretch goals (would strengthen the paper significantly)

- Finite-key analysis or at least a discussion of finite-key corrections.
- Extend to $d > 2$ (qudits) to show the framework generalizes.
- SDP relaxation hierarchy for the outer optimization over protocols
  (replacing the heuristic with something provably convergent).

### Estimated effort: High (2-4 months of focused work, mainly on analytical N=1 result and upper bound proof)

---

## Route 2: IEEE Transactions on Information Theory

**Bar:** Provably optimal or near-optimal constructions; information-theoretic
proofs (converse bounds, capacity-style arguments); rigorous security proofs.

Everything from Route 1 is a prerequisite, plus the following additional items:

### Additional requirements beyond Route 1

1. **Converse bound (critical)**
   - Prove an upper bound on $R^*_\epsilon = \sup_{\text{QKD}} R_\epsilon$ that matches
     (or nearly matches) the numerical achievability results.
   - This is the single hardest and most important item. Without it, the paper
     is fundamentally not an IT paper.
   - Possible approaches:
     - Relate $R^*_\epsilon$ to a quantum channel capacity
       (private capacity of the worst-case channel in the $\Phi \sim \epsilon$ class).
     - Use entropic uncertainty relations to derive an upper bound on
       $H(A|E) - H(A|B)$ as a function of $\epsilon$.
     - Meta-converse techniques (Tomamichel-Hayashi style).

2. **Composable security proof**
   - IT referees will expect security in the composable framework
     (Renner, Portmann-Renner).
   - At minimum, prove that the SKR formula (Eq. 3) is achievable in the
     composable sense under the stated noise model.

3. **Characterize optimal protocol structure analytically**
   - For IT, numerical examples are supporting evidence, not results.
   - Ideally, prove that the optimal N=1 protocol has a specific form
     (e.g., parametrized by a single angle).
   - Prove or disprove: does $R^*_\epsilon$ increase monotonically with $N$?
     Is there a finite $N^*(\epsilon)$ that is optimal?

4. **Capacity interpretation**
   - Frame $R^*_\epsilon$ as the secret-key capacity of the class of channels
     with Choi fidelity $\geq 1-\epsilon$.
   - Compare with known private capacity bounds (Pirandola-Laurenza-Ottaviani-Banchi
     PLOB bound, etc.).

5. **Coding-theoretic framing**
   - IT expects information-theoretic language: achievability + converse.
   - Reframe the ABC heuristic as evidence for achievability; the converse
     bound closes the gap.
   - Discuss error correction explicitly: LDPC codes, polar codes for the
     reconciliation step; how close to Shannon limit?

6. **Writing style**
   - IT papers follow a specific structure: Introduction, System Model,
     Main Results (stated as Theorems), Proofs, Numerical Illustrations,
     Conclusion.
   - All results must be stated as formal Theorems/Propositions with proofs.

### Estimated effort: Very high (6-12 months; the converse bound alone is a research problem)

---

## Comparison summary

| Aspect | Quantum | IEEE IT |
|--------|---------|---------|
| Analytical N=1 result | Required | Required |
| Upper/converse bound | Strongly recommended | Required (tight) |
| Composable security | Not required | Required |
| Numerical results | Central contribution | Supporting evidence only |
| Noise model justification | Needed | Needed + capacity interpretation |
| Heuristic validation | Needed | Needed + coding-theoretic context |
| Writing effort | High | Very high (formal theorem/proof style) |
| Timeline estimate | 2-4 months | 6-12 months |
| Probability of acceptance | Moderate (if N=1 analytical result is strong) | Low-moderate (converse bound is hard) |

---

## Recommended strategy

1. **Start with the N=1 analytical result** -- this is the lynchpin for both routes.
   The 1D parametrization is already in `generate_three_case_curves.jl`; the SDP
   for N=1 may admit a closed-form dual solution.
2. **Characterize the upper bound** -- define and prove whatever the brown curve is.
3. **Write a proper introduction and literature review.**
4. **Submit to Quantum** once items 1-3 are done.
5. **If the converse bound falls out naturally**, pivot to IT. Otherwise, Quantum
   is the better target for the current scope of results.
