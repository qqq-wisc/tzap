# tzap-lean

`tzap-lean` is a formally verified Lean 4 port of [tzap](https://github.com/qqq-wisc/tzap), the Rust optimizer for Clifford+T/Rz circuits. It formalizes the circuit representation and channel semantics, then implements deterministic passes with unconditional correctness proofs and randomized phase folding with a proved failure-probability bound.

## Installation

Install [Lean through `elan`](https://lean-lang.org/install/), then clone and build the project:

```sh
git clone https://github.com/qqq-wisc/tzap.git
cd tzap/lean
lake exe cache get
lake build
```

## Running

Optimize an OpenQASM 2.0 circuit and write the result to a file:

```sh
lake exe tzap-lean input.qasm optimized.qasm
```

With `--verbose`, the input gate set and the exact MURM synthesis basis are shown in canonical
gate order:

```console
$ lake exe tzap-lean ../benchmarks/feynman/tof_5.qasm --verbose -O2
⚡️ tzap-lean
  Parsed ../benchmarks/feynman/tof_5.qasm in 0.000s
	└─ 9 qubits · 133 gates
	   └─ Circuit gates: {h, t, tdg, cx}

  Loaded MURM (549,456 unitaries) in 0.070s
    └─ Synthesis basis: {h, x, z, s, sdg, t, tdg, cx}

  Result  (0.012s)
    Gates     133 → 90  (32.3% reduction)
    2q gates  42 → 36  (14.2% reduction)
    T/Tdg     49 → 31  (36.7% reduction)
    Depth     86 → 69  (19.7% reduction)
```

CCX, CCZ, CZ, and Rz remain native by default. Use `--decompose-ccx` to decompose CCX and
CCZ, or `--decompose-cz` to decompose CZ. The Lean implementation does not decompose Rz.
Requested exact decompositions run after input optimization, followed by another optimization
stage.

SuperOpt uses `--superopt-gates auto` by default. `base` selects the fixed Clifford+T basis;
an exact comma-separated list such as `h,x,z,s,sdg,t,tdg,cx,cz` restricts what its MURM may
emit.

The supported OpenQASM gates are `h`, `x`, `z`, `s`, `sdg`, `t`, `tdg`, `rz`, `cx`, `cz`,
`ccx`, `ccz`, `measure`, and `reset`.

## The obligation

Deterministic passes return a checked circuit that is provably equivalent to their input:

```lean
structure Pass where
  name : String
  run : ∀ {n m}, Circuit n m → Circuit n m
  correct : ∀ {n m} (c : Circuit n m),
    (run c).Equivalent c
```

Randomized phase folding uses fresh 128-bit tags from the operating system. Its executable CCX rule matches Rust: evaluate `target + control₁·control₂` in packed `GF(2^128)`, propagate a formal-degree upper bound, and replace the target with a fresh opaque value when the next product would exceed `2^32`. Reset writes the field zero, and rotations on known zero/one values are removed. The formal polynomial analysis proves that this nonlinear CCX transfer simulates Toffoli on Boolean basis states, including the conservative cutoff fallback. `GF128Bridge.lean` proves that packed arithmetic implements the abstract field and that uniform 128-bit samples induce uniform field samples. The verified pipeline now uses the same nonlinear sampled transformation as the executable, with a whole-circuit failure bound of `comparison_pairs · 2^-96`. The bound assumes the operating system supplies independent uniform samples.

## What is trusted

| Component | Role in the trust boundary | Consequence if wrong |
|---|---|---|
| Formal circuit semantics and specification | Define what “correct” means: equality of the modeled quantum channels. | The development could prove preservation of an unintended model. |
| OpenQASM output boundary | `serializeChecked` reparses emitted text and refuses to write it unless it reconstructs the optimized circuit. | A failed check stops output rather than emitting a changed circuit. |
| Operating-system randomness | Supplies the independent uniform 128-bit tags assumed by the `PhaseFoldRand` probability theorem. | Biased, correlated, or adversarial entropy can invalidate the numerical failure bound. |
| Lean compiler and runtime | Execute the proved definitions and IO shell. | Outside the logic’s proof guarantee, as for any compiled verification artifact. |
