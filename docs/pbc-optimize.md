# PBC rotation optimization

`PbcCircuit::optimize_rotations` lowers the T count of a converted PBC
circuit. It implements the design in `docs/pbc-phase-folding.tex`
(generalized MCR and a bounded left-to-right sweep), including its Clifford
normalization stage, with the changes listed below. On the command line:

```bash
tzap input.qasm --to-pbc --pbc-opt -o output.pbc
```

The CLI reports the change, e.g. `PBC rotation optimization: T 215 → 179
(18 merges, 0 MCR swaps, 13 Cliffords to frame, 0.000s)`. Without
`--pbc-opt`, it reports the PBC T count.

## What it does

The T count is the number of odd-angle rotations (`r` lines with odd `k`).
The pass preserves the exact quantum–classical channel. Output rotations
reuse existing Pauli handles where their axis is unchanged. New axes get
fresh handles, built from shared single-qubit leaves.

1. **Pack.** Every axis is evaluated once into canonical x/z bit planes by a
   dedicated evaluator over the shared Pauli DAG. That covers rotations,
   measurements, conditional rotations, and the output frame's images. The
   evaluator frees each node's value after its last use. Signs are absorbed
   into angles (`R_{-P}(k) = R_P(-k)`), and identical axes are interned.
2. **Stream: merge, and move Cliffords to the frame.** One pass walks all
   operations in order, carrying the Clifford F that is being moved to the
   end. At each operation it does the following:
   - **Conjugate:** replace the operation's axis P by `F† P F`. This applies
     to measurements and conditional rotations too, so outcomes and
     conditions are unchanged.
   - **Absorb:** absorb any unconditional Clifford rotation (`k = ±2, 4`)
     into F. A rotation merges into the latest earlier rotation with the same
     axis when every live rotation between them commutes with it. If the
     merged angle is even, that rotation also joins F. It commutes with
     everything after it (which is what allowed the merge), so moving it to
     the end is sound.

   Removing Cliffords this way unblocks merges later in the same pass.
   Measurements and conditional rotations end merging but not F. At the end,
   F is composed into the output frame (`f` records).
3. **MCR swaps.** Each segment between barriers is partitioned greedily into
   internally commuting runs of at most `window` rotations. For neighboring
   runs A|B|C, and up to `candidates` variants with A shortened to a suffix
   or C to a prefix, the pass tries the A,B swap, then the B,C swap. A swap
   is accepted only if two conditions hold:
   - **Certificate:** the generators commute, `[H_A, H_B] = 0`. This is
     checked with exact packed Pauli products and phases, as in the note.
   - **Improvement:** merging then lowers the T count. With `eager_swaps`
     (the default), it is enough that some rotations merge without raising
     T.

   Either way each accepted swap shrinks the circuit, so the search
   terminates.

Stream and swap rounds alternate until a round removes no T, up to `rounds`.
A final stream moves any Clifford rotations left by swaps into the frame.
So no unconditional Clifford rotation remains in the output.

## Differences from the design note

- **Global merge instead of a W-buffer.** The note merges only within a
  streaming buffer of W rotations. Merging across long distances matters in
  practice: on feynman after O3, a 4096-rotation lookback finds 0.08% T
  reduction and an unlimited lookback 0.66%. So merging is global, bounded
  by `lookback` (default 2^20). Two structures make it fast:
  - **Candidate lookup:** a per-axis table gives the latest rotation with the
    same axis in O(1), so nothing is scanned unless a merge is possible.
  - **Blocker check:** each rotation carries a 64-bit support signature (bit
    `q / L` for each qubit `q`, with `L = ceil(n/64)`). Per signature bit,
    the pass keeps the positions of the rotations that touch it. A blocker
    must share a qubit with the axis, so only those lists are checked.
- **Clifford normalization is interleaved with merging.** The note runs
  normalization as a separate stage after the sweep. Here each Clifford
  moves into F as soon as it appears. F is kept as packed images of every
  qubit's X, Y and Z; qubits it has not touched map to themselves, so
  conjugating an axis costs one product per touched qubit in its support.
  A first version alternated whole merge rounds with separate normalization
  passes. The streaming form reached a lower T count (cobble-t: 248,645 vs
  248,377) in about half the time.
- **Runs instead of re-optimizing the buffer per emitted rotation.** The
  note re-examines its buffer every time it emits one rotation. Here each run
  triple is examined once per round, and accepted rewrites are applied in one
  rebuild.
- **Cheap swap filter.** A and C each commute internally, so after either
  swap every axis they share merges. Triples without a shared axis (with
  strict swaps: a shared odd-angle axis) are skipped before any certificate
  is computed. Certificates are memoized across variants.
- **Packing is dense.** Each axis uses `2 * ceil(n/64)` words. There is no
  sparse-block representation yet. `max_packed_words` (default 64 Mi words)
  bounds the storage, and the pass fails, leaving the circuit unchanged,
  rather than exceeding it.

## Options

| field | default | meaning |
|---|---|---|
| `lookback` | 2^20 | how many rotations a merge may reach back |
| `window` | 8 | maximum MCR group length; 0 disables swaps |
| `candidates` | 16 | shorter A/C variants tried per run triple |
| `rounds` | 8 | maximum stream/swap rounds |
| `clifford_to_frame` | true | move Clifford rotations into the output frame |
| `eager_swaps` | true | accept swaps that merge without lowering T |
| `lazy_cliffords` | false | keep merged Cliffords as merge targets until an anticommuting rotation, barrier, or the end forces them into the frame |
| `max_packed_words` | 64 Mi | bound on packed axis storage |

## Results

Measured with `examples/pbc_opt_eval.rs`, which prints CSV and totals:

```bash
cargo run --release --example pbc_opt_eval -- [--level none|O1|O2|O3] [--decompose-ccx] \
    [--no-frame] [--no-eager] <dirs>
```

All figures are PBC T counts (odd-angle rotations after conversion), summed
over the circuits that parse and convert in every configuration. Circuits
with Rz or resets are skipped, as are two feynman circuits whose raw form
repeats a gate operand. qqq-wisc refers to
`quantum-compiler-benchmark-circuits`. Times are on an Apple-silicon
laptop.

### The pass compared with tzap O1, O2, and O3

"raw" is the converted input. "pass only" runs the pass on it directly.
"Ox" is tzap at that level (CCX kept native, the default), then conversion,
and "Ox+pass" runs the pass on top.

| corpus | n | raw | pass only | O1 | O1+pass | O2 | O2+pass | O3 | O3+pass |
|---|---|---|---|---|---|---|---|---|---|
| feynman | 41 | 913,262 | 504,676 | 508,172 | 504,676 | 505,762 | 504,676 | 505,054 | 504,676 |
| cobble-t | 6 | 587,731 | 248,377 | 336,151 | 248,377 | 284,789 | 248,377 | 248,403 | 248,377 |
| qft | 4 | 1,449,188 | 1,055,028 | 1,058,260 | 1,055,028 | 1,055,594 | 1,055,028 | 1,055,098 | 1,055,028 |
| arithmetic | 18 | 4,459 | 2,575 | 4,459 | 2,575 | 4,459 | 2,575 | 4,459 | 2,575 |
| mctoffoli | 8 | 616 | 324 | 616 | 324 | 616 | 324 | 616 | 324 |
| grover | 2 | 9,372 | 6,936 | 6,982 | 6,936 | 6,942 | 6,936 | 6,942 | 6,936 |
| shor | 1 | 16,670 | 12,214 | 13,216 | 12,214 | 13,212 | 12,214 | 13,212 | 12,214 |
| jku_suite | 152 | 1,265,915 | 590,041 | 590,917 | 590,041 | 590,583 | 590,039 | 590,403 | 590,041 |
| large_circuits | 58 | 2,761,421 | 1,685,399 | 1,690,659 | 1,685,399 | 1,687,513 | 1,685,397 | 1,686,817 | 1,685,399 |

The pass alone reaches a lower T count than O3 on every corpus, and the
same count as O3 followed by the pass. Running it after any level converges
to essentially the same T; the largest spread is 2 T (jku_suite,
large_circuits). With CCX native, tzap's phase folding cannot merge across
Toffolis, so the pass removes 42–47% more on arithmetic and mctoffoli.
Individual circuits improve by up to 50% even after O3 (adder_8, mod5_4).

These are T counts only. tzap's levels also reduce gate and two-qubit
counts and emit QASM; the pass works on the PBC form.

Time for tzap at each level, and for the pass alone or after each level:

| corpus | pass only | O1 | O2 | O3 | pass after O1 | pass after O2 | pass after O3 |
|---|---|---|---|---|---|---|---|
| feynman | 0.77 s | 0.15 s | 2.63 s | 3.35 s | 0.23 s | 0.23 s | 0.22 s |
| cobble-t | 0.17 s | 0.10 s | 0.99 s | 2.46 s | 0.13 s | 0.10 s | 0.06 s |
| qft | 0.58 s | 0.21 s | 2.20 s | 4.32 s | 0.36 s | 0.36 s | 0.27 s |
| arithmetic | 0.00 s | 0.00 s | 0.02 s | 0.02 s | 0.00 s | 0.00 s | 0.00 s |
| mctoffoli | 0.00 s | 0.00 s | 0.02 s | 0.02 s | 0.00 s | 0.00 s | 0.00 s |
| grover | 0.00 s | 0.00 s | 0.06 s | 0.07 s | 0.00 s | 0.00 s | 0.00 s |
| shor | 0.01 s | 0.00 s | 0.07 s | 0.07 s | 0.01 s | 0.01 s | 0.01 s |
| jku_suite | 0.49 s | 0.18 s | 2.27 s | 4.40 s | 0.28 s | 0.26 s | 0.25 s |
| large_circuits | 1.10 s | 0.39 s | 4.51 s | 8.86 s | 0.67 s | 0.64 s | 0.53 s |

The pass alone takes 4–14× less time than O3. On raw circuits it runs at
roughly 0.1–1 µs per rotation. Wide circuits cost the most: gf2^256_mult
(768 qubits, 459k rotations) takes 0.55 s, mostly moving 98k Cliffords
into a 768-qubit frame. Conversion itself takes milliseconds.

### Contribution of each stage

Raw circuits, pass only:

| corpus | merge + MCR (strict) | + Clifford to frame | + eager swaps |
|---|---|---|---|
| feynman | 507,794 | 504,676 | 504,676 |
| cobble-t | 335,285 | 248,377 | 248,377 |
| qft | 1,068,124 | 1,055,028 | 1,055,028 |

Moving Cliffords into the frame is the largest gain after merging (cobble-t:
−26%). Eager swaps changed nothing on these corpora. **MCR swaps are rare on
real circuits:** across every corpus, the certificate accepted one swap in
total. They do fire on synthetic circuits with planted patterns, including
the note's 6 → 2 and 12 → 4 examples.

### Pauli weight, and comparison with O1 and O3

Weight is the total number of non-identity single-qubit factors over all
rotation axes (`OptimizeStats::weight_before`/`weight_after`). Entries are T
count / weight, summed over circuits that convert in every configuration.
"O1"/"O3" are tzap at that level (CCX native), then conversion.

Per corpus: T count (equal to the number of rotations), total, average,
median and maximum weight per rotation, and time (tzap and/or the pass;
measured in one run, under some system load, so absolute times are about
1.5–2× those in the earlier tables). Weights come from
`PbcCircuit::rotation_weights`, merged over the corpus.

##### feynman (41 circuits)

| configuration | T (= rotations) | total weight | avg weight | median | max | time |
|---|---|---|---|---|---|---|
| raw | 913,262 | 9,674,956 | 10.59 | 2 | 135 | — |
| O1 | 508,172 | 7,426,956 | 14.62 | 6 | 135 | 0.28 s |
| O3 | 505,054 | 7,432,552 | 14.72 | 6 | 135 | 5.85 s |
| pass only | 504,676 | 7,667,212 | 15.19 | 6 | 135 | 1.36 s |
| O1 + pass | 504,676 | 7,666,314 | 15.19 | 6 | 135 | 0.68 s |
| O3 + pass | 504,676 | 7,669,174 | 15.20 | 6 | 135 | 6.24 s |

##### cobble-t (6 circuits)

| configuration | T (= rotations) | total weight | avg weight | median | max | time |
|---|---|---|---|---|---|---|
| raw | 587,731 | 1,003,092 | 1.71 | 1 | 7 | — |
| O1 | 336,151 | 1,162,424 | 3.46 | 2 | 16 | 0.18 s |
| O3 | 248,403 | 1,146,072 | 4.61 | 1 | 20 | 4.25 s |
| pass only | 248,377 | 785,677 | 3.16 | 1 | 18 | 0.29 s |
| O1 + pass | 248,377 | 795,843 | 3.20 | 1 | 20 | 0.40 s |
| O3 + pass | 248,377 | 1,165,024 | 4.69 | 1 | 19 | 4.35 s |

##### qft (4 circuits)

| configuration | T (= rotations) | total weight | avg weight | median | max | time |
|---|---|---|---|---|---|---|
| raw | 1,449,188 | 6,655,376 | 4.59 | 2 | 50 | — |
| O1 | 1,058,260 | 16,948,053 | 16.02 | 15 | 42 | 0.37 s |
| O3 | 1,055,098 | 18,042,219 | 17.10 | 16 | 45 | 7.65 s |
| pass only | 1,055,028 | 14,157,784 | 13.42 | 12 | 43 | 1.02 s |
| O1 + pass | 1,055,028 | 16,479,107 | 15.62 | 15 | 41 | 1.01 s |
| O3 + pass | 1,055,028 | 18,040,499 | 17.10 | 16 | 45 | 8.11 s |

##### arithmetic (18 circuits)

| configuration | T (= rotations) | total weight | avg weight | median | max | time |
|---|---|---|---|---|---|---|
| raw | 4,459 | 10,048 | 2.25 | 2 | 10 | — |
| O1 | 4,459 | 10,048 | 2.25 | 2 | 10 | 0.00 s |
| O3 | 4,459 | 10,200 | 2.29 | 2 | 10 | 0.04 s |
| pass only | 2,575 | 7,602 | 2.95 | 2 | 11 | 0.00 s |
| O1 + pass | 2,575 | 7,602 | 2.95 | 2 | 11 | 0.00 s |
| O3 + pass | 2,575 | 7,850 | 3.05 | 3 | 15 | 0.04 s |

##### shor (1 circuits)

| configuration | T (= rotations) | total weight | avg weight | median | max | time |
|---|---|---|---|---|---|---|
| raw | 16,670 | 37,914 | 2.27 | 2 | 8 | — |
| O1 | 13,216 | 33,584 | 2.54 | 3 | 12 | 0.00 s |
| O3 | 13,212 | 33,574 | 2.54 | 3 | 12 | 0.12 s |
| pass only | 12,214 | 37,052 | 3.03 | 3 | 10 | 0.02 s |
| O1 + pass | 12,214 | 37,438 | 3.07 | 3 | 12 | 0.02 s |
| O3 + pass | 12,214 | 35,446 | 2.90 | 3 | 12 | 0.14 s |

##### jku_suite (152 circuits)

| configuration | T (= rotations) | total weight | avg weight | median | max | time |
|---|---|---|---|---|---|---|
| raw | 1,265,915 | 2,690,062 | 2.12 | 2 | 11 | — |
| O1 | 590,917 | 3,687,784 | 6.24 | 7 | 15 | 0.33 s |
| O3 | 590,403 | 3,683,970 | 6.24 | 7 | 15 | 7.85 s |
| pass only | 590,041 | 3,640,718 | 6.17 | 6 | 15 | 0.83 s |
| O1 + pass | 590,041 | 3,683,454 | 6.24 | 7 | 15 | 0.81 s |
| O3 + pass | 590,041 | 3,682,026 | 6.24 | 7 | 15 | 8.29 s |

##### large_circuits (58 circuits)

| configuration | T (= rotations) | total weight | avg weight | median | max | time |
|---|---|---|---|---|---|---|
| raw | 2,761,421 | 9,439,587 | 3.42 | 2 | 50 | — |
| O1 | 1,690,659 | 20,808,216 | 12.31 | 9 | 42 | 0.71 s |
| O3 | 1,686,817 | 21,873,133 | 12.97 | 10 | 45 | 15.51 s |
| pass only | 1,685,399 | 17,944,030 | 10.65 | 8 | 43 | 1.89 s |
| O1 + pass | 1,685,399 | 20,322,396 | 12.06 | 9 | 41 | 1.86 s |
| O3 + pass | 1,685,399 | 21,871,281 | 12.98 | 10 | 45 | 16.43 s |

- **Median:** raw circuits have narrow axes (median 1–2). O1 and O3 raise
  the median sharply on qft (15–16), jku_suite (7) and large_circuits
  (9–10). The pass raises it less (qft: 12).
- **Max:** the widest axis is set by the circuit, not the optimizer, on
  feynman (135, from gf2^256_mult) and on qft/large_circuits (50 qubits,
  raw max 50). Frame moves narrow the widest qft axes slightly (43–45),
  while making typical ones much wider. On cobble-t they raise the max from
  7 to 16–20.

Time for tzap and for the pass:

| corpus | O1 | O3 | pass only | O3 + pass |
|---|---|---|---|---|
| feynman | 0.15 s | 3.36 s | 0.76 s | 3.59 s |
| cobble-t | 0.10 s | 2.47 s | 0.17 s | 2.53 s |
| qft | 0.21 s | 4.35 s | 0.59 s | 4.63 s |
| arithmetic | 0.00 s | 0.03 s | 0.00 s | 0.03 s |
| shor | 0.00 s | 0.07 s | 0.01 s | 0.08 s |
| jku_suite | 0.18 s | 4.42 s | 0.49 s | 4.67 s |
| large_circuits | 0.40 s | 8.98 s | 1.10 s | 9.51 s |

- **T count:** the pass alone beats both O1 and O3 everywhere; any level
  followed by the pass reaches the same T.
- **Weight:** the pass alone has lower weight than O1 and O3 except on
  feynman (+3%) and shor (+10%).
- **Time:** O1 is 1.7–5× faster than the pass; O3 is 4–15× slower.

#### A minimal case: O1 then conversion is wider than conversion then the pass

```
cx q[0],q[1]; t q[1]; t q[1]; tdg q[1]; cx q[0],q[1]; h q[0]; cx q[1],q[0]; t q[0];
```

Both pipelines end with T = 2. tzap -O1 combines T·T·T† on the parity
Z0Z1 into one T, so conversion gives `r 1 Z0 Z1` and `r 1 X0 Z1` (weight
4). The pass (defaults) merges the first two into an S, moves it into the
frame, and keeps the T† as `r -1 Z0 Z1`. The last T's axis X0Z1 is then
conjugated by that S, and the shared Z1 cancels: X0Z1 · Z0Z1 ~ Y0. The
result is `r -1 Z0 Z1` and `r -1 Y0` (weight 3). Conjugation by a frame
Clifford widens an axis when the supports are disjoint, and can narrow it
when they overlap; this is also why `lazy_cliffords` (which avoids that S
here) is not uniformly better for weight: it raises weight 11–14% on qft
and large_circuits at the same T.

#### Why tzap's optimizations raise weight

A rotation's axis is the pullback of its T gate's Z through every Clifford
gate before it, so weight depends on how the Clifford gates are arranged, not
only on the unitary. The raw qft circuits are ideal for narrow axes. Their
Clifford gates are only H and CX (no S, Sdg, or Z), 47% of axes have weight
1, and the average is 2.7. Measured on `qft_q020_d32421` (20 qubits):

| circuit | S / Sdg / Z gates | H gates | T | average T-axis weight | same, phase gates kept as local rotations |
|---|---|---|---|---|---|
| raw | 0 / 0 / 0 | 141,858 | 167,567 | 2.70 | 2.70 |
| PhaseFoldRand only | 10,599 / 8,822 / 2,057 | 141,858 | 123,853 | 8.82 | 2.66 |
| CancelGates only | 17,403 / 17,334 / 0 | 123,064 | 130,409 | 9.31 | 6.52 |
| O1 | 6,217 / 234 / 6,252 | 123,064 | 122,695 | 8.87 | 6.50 |

- **Phase folding** merges T gates on equal parities, producing new phase
  gates (T·T = S, and so on). Conversion folds each new Clifford into the
  frame, so every later rotation whose axis Q anticommutes with its axis P
  becomes (up to phase) PQ. Repeated thousands of times, this spreads the
  axes. Converting those phase gates as local `r 2`/`r 4` rotations instead
  brings the average back to 2.66, below the raw circuit: the merges alone
  do not widen anything.
- **Gate cancellation's Hadamard reduction** (e.g. H S H → S† H S†)
  removes and moves H gates. That changes which basis later T gates sit in
  relative to the CX gates, and it widens axes even when its phase gates
  stay local (6.52).
- **The pass's own frame moves** are the same mechanism as the first
  effect.

A converter option that keeps phase gates (and the pass's merged Cliffords)
as local rotations would recover narrow axes at the cost of extra Clifford
rotations in the output; a weight-aware policy could choose per gate.

### FTCircuitBench suite

`benchmarks/ftcircuitbench` holds 63 Clifford+T circuits synthesized with
Gridsynth (precision 5) from FTCircuitBench's inputs (see its README). T counts
per family, with the average weight per rotation:

| family | raw T | O1 T | O3 T | pass only T | O3 + pass T | avg weight: raw / O3 / pass only |
|---|---|---|---|---|---|---|
| adder | 624 | 360 | 360 | 360 | 360 | 4.6 / 6.5 / 6.5 |
| hamiltonians | 2,070,200 | 2,070,200 | 2,041,572 | 2,041,412 | 2,041,412 | 26.2 / 31.6 / 21.1 |
| hamiltonians_5trotter | 6,747,130 | 6,710,368 | 6,658,504 | 6,648,664 | 6,648,664 | 56.7 / 49.4 / 50.8 |
| hhl | 128,332 | 126,732 | 125,928 | 125,786 | 125,806 | 6.7 / 6.4 / 6.9 |
| qft | 156,614 | 156,614 | 156,134 | 156,134 | 156,134 | 20.6 / 19.0 / 20.5 |
| qpe | 466,229 | 460,617 | 456,523 | 456,163 | 456,163 | 5.4 / 5.4 / 5.5 |
| qsvt | 506,632 | 500,352 | 494,972 | 494,290 | 494,290 | 4.8 / 5.4 / 5.0 |

| whole suite | T | vs raw | time |
|---|---|---|---|
| raw | 10,075,761 | | |
| O1 | 10,025,243 | −0.5% | 3.6 s |
| O3 | 9,933,993 | −1.4% | 86.7 s |
| pass only | 9,922,809 | −1.5% | 11.5 s |
| O3 + pass | 9,922,829 | −1.5% | 93.1 s |

- **Little redundancy:** Gridsynth sequences are close to T-optimal per
  rotation, and consecutive T gates within a sequence are separated by
  Hadamards, so their axes anticommute and cannot merge. Every optimizer
  removes only 0.5–1.5% of T. The pass alone removes the most, 7.5×
  faster than O3.
- **Weight:** on the full-length Hamiltonians the pass lowers average weight
  (26.2 → 21.1) while O3 raises it (31.6); on the Trotter-5 set both lower it.
  Where frame moves cancel overlapping factors, they narrow axes (see
  `docs/pbc-weight-example.md`).
- **Adders:** FTCircuitBench's pipeline already decomposed their Toffolis
  into Clifford+T; every optimizer takes them from 624 to the same 360 T.

### FTCircuitBench large tier

`benchmarks/ftcircuitbench-large`: the 32 FTCircuitBench inputs with more than
60,000 rotations (0.5–24.6M gates each), same synthesis.

Gates:

| family | raw | O1 | O2 | O3 |
|---|---|---|---|---|
| hamiltonians | 51,068,156 | 42,613,252 (−16.6%) | 34,629,545 (−32.2%) | 34,207,580 (−33.0%) |
| hhl | 2,949,323 | 2,862,635 (−2.9%) | 2,323,835 (−21.2%) | 2,305,815 (−21.8%) |
| qpe | 101,414,705 | 98,454,136 (−2.9%) | 79,902,669 (−21.2%) | 79,294,783 (−21.8%) |
| qsvt | 54,541,232 | 53,424,217 (−2.0%) | 43,904,347 (−19.5%) | 43,651,390 (−20.0%) |
| **total** | 209,973,416 | 197,354,240 (−6.0%) | 160,760,396 (−23.4%) | 159,459,568 (−24.1%) |

T (PBC rotations):

| family | raw | O1 | O2 | O3 | pass only | O1 + pass | O2 + pass | O3 + pass |
|---|---|---|---|---|---|---|---|---|
| hamiltonians | 15,820,480 | 15,820,480 | 15,796,896 | 15,791,298 | 15,776,698 | 15,776,698 | 15,776,698 | 15,776,698 |
| hhl | 1,133,342 | 1,120,882 | 1,114,574 | 1,112,138 | 1,110,784 | 1,111,540 | 1,110,896 | 1,110,930 |
| qpe | 38,995,149 | 38,545,197 | 38,335,109 | 38,257,297 | 38,178,591 | 38,224,221 | 38,183,181 | 38,183,181 |
| qsvt | 21,283,980 | 21,003,022 | 20,870,738 | 20,852,690 | 20,838,122 | 20,839,362 | 20,838,866 | 20,840,602 |
| **total** | 77,232,951 | 76,489,581 | 76,117,317 | 76,013,423 | 75,904,195 | 75,951,821 | 75,909,641 | 75,911,411 |

Weight (whole tier): total / avg / median / max

| configuration | total | avg | median | max |
|---|---|---|---|---|
| raw | 1,164,986,406 | 15.08 | 8 | 116 |
| O1 | 1,388,217,064 | 18.15 | 9 | 114 |
| O2 | 1,389,234,848 | 18.25 | 9 | 116 |
| O3 | 1,439,148,971 | 18.93 | 9 | 118 |
| pass only | 1,152,452,505 | 15.18 | 9 | 115 |
| O1 + pass | 1,383,148,272 | 18.21 | 9 | 115 |
| O2 + pass | 1,386,436,165 | 18.26 | 9 | 115 |
| O3 + pass | 1,429,551,431 | 18.83 | 9 | 119 |

Time (whole tier):

| | tzap | pass | total |
|---|---|---|---|
| pass only | | 44.3 s | 44.3 s |
| O1 + pass | 19.5 s | 40.0 s | 59.5 s |
| O2 + pass | 145.6 s | 36.8 s | 182.4 s |
| O3 + pass | 598.8 s | 35.3 s | 634.0 s |

- **T:** O1 −0.96%, O2 −1.44%, O3 −1.58%, pass only −1.72%. As on the small
  tier, the pass alone removes the most T, and adding it after a tzap level
  converges to within 0.01% of that.
- **Gates:** O2 and O3 remove 23–24% of gates, nearly all of it Clifford
  overhead from Gridsynth; O3 adds only 0.7 points over O2 at 4× the time.
- **Weight:** the pass alone keeps average weight at the raw level (15.2 vs
  15.1) and lowers the total, because it removes rotations; every tzap level
  raises average weight by 20–25%.
- **Time:** the pass takes 44 s on 210M gates; O3 takes 10 minutes.

### MCR window sizes and Litinski's algorithm

Raw circuits, pass only, with Clifford moves on. "window 0" disables MCR:
rotations only move past rotations they commute with, one at a time.
"Litinski" is `Strategy::Litinski`, the greedy layering of "A Game of
Surface Codes" (Sec. 1). That algorithm partitions rotations into layers of
mutually commuting rotations. It then repeatedly moves each rotation of
layer i+1 into layer i if it commutes with all of layer i. Equal rotations
that meet in a layer combine, and a combined Clifford is commuted to the
end. T-depth is the number of such layers, filled as early as possible and
measured the same way for every configuration.

| corpus | raw T / T-depth | Litinski T / T-depth (time) | window 0 | window 8 (default) | window 32 |
|---|---|---|---|---|---|
| feynman | 913,262 / 29,884 | 504,678 / 29,205 (204 s) | 504,678 / 29,205 (0.68 s) | 504,676 / 29,206 (0.77 s) | 504,678 / 29,205 (0.95 s) |
| cobble-t | 587,731 / 157,730 | 248,377 / 148,996 (511 s) | 248,377 / 148,996 (0.13 s) | 248,377 / 148,996 (0.17 s) | 248,377 / 148,996 (0.18 s) |
| jku_suite | 1,265,915 / 137,070 | 590,049 / 136,970 (24 s) | 590,049 / 136,970 (0.33 s) | 590,041 / 136,967 (0.49 s) | 590,033 / 136,965 (0.84 s) |
| shor | 16,670 / 1,906 | 12,214 / 1,875 (0.29 s) | 12,214 / 1,875 | 12,214 / 1,875 | 12,214 / 1,875 |
| arithmetic | 4,459 / 87 | 2,575 / 87 | 2,575 / 87 | 2,575 / 87 | 2,575 / 87 |
| qft | 1,449,188 / 203,264 | did not finish in 50 min | 1,055,028 / 173,864 (0.40 s) | 1,055,028 / 173,864 (0.59 s) | 1,055,028 / 173,864 (0.59 s) |
| large_circuits | 2,761,421 / 356,599 | not run (contains the qft circuits) | 1,685,407 / 324,779 (0.74 s) | 1,685,399 / 324,776 (1.10 s) | 1,685,395 / 324,775 (1.45 s) |

- **Litinski's algorithm reaches exactly the result of window 0.** Both
  give the same T count and T-depth on every corpus; its combine step is
  the same commute-and-merge. As literally specified, it is superlinear: each
  repetition moves rotations only one layer, so it repeats about as many
  times as there are layers. It takes 70–4000× longer; gf2^256_mult alone
  takes 167 s.
- **Larger windows barely matter.** MCR adds 8–16 T of reduction on
  jku_suite with windows of 8 or more, 8–12 T on large_circuits, 2 T on
  feynman, and nothing elsewhere. Eager swaps on cobble-t (51) and qft (248) merge rotations
  without changing T.
- **T-depth.** The pass lowers T-depth only as a side effect of removing
  rotations (feynman: 29,884 → 29,205). Litinski's layering targets depth,
  yet ends at the same depth, because equal commutation freedom yields the
  same as-soon-as-possible layers.

#### Window sweep over every corpus

The pass alone (no tzap optimization first) on all 11 corpora, 385 circuits
and 94.3M T. "w=0" disables MCR.

T after the pass (MCR swaps in parentheses), by window size:

| corpus | raw T | w=0 | w=1 | w=2 | w=4 | w=8 | w=16 | w=32 | w=64 |
|---|---|---|---|---|---|---|---|---|---|
| feynman (41) | 913,262 | 504,678 (0) | 504,678 (0) | 504,678 (0) | 504,678 (0) | 504,676 (1) | 504,678 (0) | 504,678 (0) | 504,678 (0) |
| cobble-t (6) | 587,731 | 248,377 (0) | 248,377 (0) | 248,377 (52) | 248,377 (51) | 248,377 (51) | 248,377 (51) | 248,377 (51) | 248,377 (51) |
| qft (4) | 1,449,188 | 1,055,028 (0) | 1,055,028 (0) | 1,055,028 (248) | 1,055,028 (248) | 1,055,028 (248) | 1,055,028 (248) | 1,055,028 (248) | 1,055,028 (248) |
| arithmetic (18) | 4,459 | 2,575 (0) | 2,575 (0) | 2,575 (0) | 2,575 (0) | 2,575 (0) | 2,575 (0) | 2,575 (0) | 2,575 (0) |
| mctoffoli (8) | 616 | 324 (0) | 324 (0) | 324 (0) | 324 (0) | 324 (0) | 324 (0) | 324 (0) | 324 (0) |
| grover (2) | 9,372 | 6,936 (0) | 6,936 (0) | 6,936 (0) | 6,936 (0) | 6,936 (0) | 6,936 (0) | 6,936 (0) | 6,936 (0) |
| shor (1) | 16,670 | 12,214 (0) | 12,214 (0) | 12,214 (0) | 12,214 (0) | 12,214 (0) | 12,214 (0) | 12,214 (0) | 12,214 (0) |
| jku_suite (152) | 1,265,915 | 590,049 (0) | 590,049 (0) | 590,049 (0) | 590,049 (1) | 590,041 (3) | 590,033 (5) | 590,033 (5) | 590,033 (5) |
| large_circuits (58) | 2,761,421 | 1,685,407 (0) | 1,685,407 (0) | 1,685,407 (260) | 1,685,407 (261) | 1,685,399 (263) | 1,685,395 (264) | 1,685,395 (264) | 1,685,395 (264) |
| ftcircuitbench (63) | 10,075,761 | 9,923,101 (0) | 9,923,101 (0) | 9,923,085 (3376) | 9,923,049 (3386) | 9,922,809 (3446) | 9,922,809 (3446) | 9,922,809 (3446) | 9,922,809 (3446) |
| ftcircuitbench-large (32) | 77,232,951 | 75,956,755 (0) | 75,956,755 (0) | 75,950,077 (6866) | 75,946,439 (6878) | 75,904,195 (14644) | 75,904,195 (14644) | 75,904,195 (14644) | 75,904,195 (14644) |
| **all** | 94,317,346 | 89,985,444 (0) | 89,985,444 (0) | 89,978,750 (10802) | 89,975,076 (10825) | 89,932,574 (18656) | 89,932,564 (18658) | 89,932,564 (18658) | 89,932,564 (18658) |

Time of the pass (s):

| corpus | w=0 | w=1 | w=2 | w=4 | w=8 | w=16 | w=32 | w=64 |
|---|---|---|---|---|---|---|---|---|
| feynman | 0.70 | 0.73 | 0.72 | 0.75 | 0.80 | 0.83 | 0.97 | 1.21 |
| cobble-t | 0.13 | 0.17 | 0.17 | 0.17 | 0.17 | 0.17 | 0.17 | 0.17 |
| qft | 0.40 | 0.58 | 0.59 | 0.59 | 0.60 | 0.58 | 0.59 | 0.58 |
| arithmetic | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| mctoffoli | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| grover | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| shor | 0.00 | 0.00 | 0.00 | 0.00 | 0.01 | 0.01 | 0.01 | 0.01 |
| jku_suite | 0.32 | 0.37 | 0.36 | 0.43 | 0.48 | 0.67 | 0.83 | 0.83 |
| large_circuits | 0.71 | 0.94 | 0.95 | 1.02 | 1.07 | 1.26 | 1.43 | 1.44 |
| ftcircuitbench | 5.16 | 6.46 | 6.48 | 6.49 | 6.45 | 6.45 | 6.47 | 6.46 |
| ftcircuitbench-large | 29.46 | 43.49 | 42.87 | 43.04 | 44.13 | 44.51 | 44.85 | 44.44 |
| **all** | 36.9 | 52.7 | 52.1 | 52.5 | 53.7 | 54.5 | 55.3 | 55.1 |

Total weight after the pass:

| | w=0 | w=1 | w=2 | w=4 | w=8 | w=16 | w=32 | w=64 |
|---|---|---|---|---|---|---|---|---|
| **all** | 1,582,936,139 | 1,582,936,139 | 1,586,557,790 | 1,582,675,096 | 1,586,349,973 | 1,586,349,353 | 1,586,349,353 | 1,586,349,353 |

- **MCR's total contribution:** 52,880 T over window 0, 0.06% of the
  post-pass T count. Almost all of it (52,560 T) comes from the
  FTCircuitBench large tier, where swaps jump from 6,878 to 14,644 at
  window 8. It is the only corpus where MCR matters noticeably.
- **Saturation:** results stop changing at window 16. Windows 8 and 16
  differ by 10 T in total.
- **Window 1 is pure overhead:** single-rotation groups can only swap
  rotations that commute, which ordinary merging already handles. So it gives
  window 0's result and costs as much as window 64.
- **Cost:** enabling MCR at all costs about 45% more time (36.9 s → 52–55 s).
  Beyond that, window size barely matters: the swap scan is dominated by
  building runs, not by certificate size.
- **Weight:** unchanged to within 0.2%.

### SIMD

There is no hand-written SIMD. LLVM auto-vectorizes the multi-word kernels
(`anticommutes`, `product`) into unrolled NEON loops: four 128-bit
registers per iteration, AND/XOR, then a vector population count. That is
the design note's SIMD plan, generated by the compiler. Circuits with up to
about 500 qubits (fewer than 8 words per plane) use the scalar path, where
a check is a handful of instructions. On most corpora (64 qubits or fewer)
an axis is one word per plane, so SIMD cannot help. The frame's image scan
(`Clifford::absorb`) is not vectorized; it dominates only on very wide
circuits (gf2^256_mult).

## Validation

- **Kernels.** Commutation and ordered products with phases are checked
  against the single-qubit Pauli table: exhaustively for two qubits, and on
  random words of 1–130 qubits across word boundaries.
- **Packing.** Packed axes are checked against `expand()`, including signs
  inside the DAG and canonical handles.
- **Worked examples.** The note's positive example (6 → 2), its negative
  variant (rejected), and the four-versus-four example (12 → 4) are tests.
- **Exact fuzzing.** Exact-unitary fuzzing covers random rotation sequences
  under several option settings, converted gate circuits with CCX/CCZ, and
  planted MCR patterns (letter-permuted, embedded, with spectators and
  noise).
- **Channels.** The mid-circuit measurement fuzzer checks the pass against
  exact channels. Measurements are merge barriers, and their axes are
  conjugated as Cliffords move past them. This includes the extended mode,
  which optimizes at every tzap level first.
- **Clifford moves and eager swaps.** Dedicated tests cover a Clifford move
  that unblocks a merge (Z T; X S; Y T → no T), measurement conjugation, and
  eager versus strict swaps. The exact-unitary oracle applies the output
  frame, so every fuzz case checks the rewritten frame too.

## Future work

- **Two-way frame tableau.** Absorbing a Clifford visits every image of
  F (3n of them). Keeping F⁻¹ as well would give the exact set of
  anticommuting images, which matters for very wide circuits such as
  gf2^256_mult.
- **Commuting past measurements.** Let rotations that commute with a
  measurement move across it, instead of treating every measurement as a
  barrier.
- **Sparse packing.** Use sparse-block packing for very wide circuits.
- **Richer MCR candidates.** Search larger or non-greedy partitions, since
  the greedy runs rarely satisfy the certificate on real circuits.
