# tzap-lean

`tzap-lean` is a formally verified Lean 4 port of [tzap](https://github.com/qqq-wisc/tzap), the Rust quantum-circuit optimizer. It formalizes the circuit representation and channel semantics, then implements optimizer passes together with proofs that they preserve those semantics.

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

## The obligation

The optimizer is composed of a number of passes (`Pass`), each of which returns a new, optimized circuit that is provably equivalent (`correct`) to the input circuit. 

```lean
structure Pass where
  name : String
  run : ∀ {n m}, Circuit.Checked n m → Circuit.Checked n m
  correct : ∀ {n m} (c : Circuit.Checked n m),
    (run c).Equivalent c
```

## What is trusted

| Component | Role in the trust boundary | Consequence if wrong |
|---|---|---|
| Formal circuit semantics and specification | Define what “correct” means: equality of the modeled quantum channels. | The development could prove preservation of an unintended model. |
| OpenQASM output boundary | `serializeChecked` reparses emitted text and refuses to write it unless it reconstructs the optimized circuit. | A failed check stops output rather than emitting a changed circuit. |
| Lean compiler and runtime | Execute the proved definitions and IO shell. | Outside the logic’s proof guarantee, as for any compiled verification artifact. |
