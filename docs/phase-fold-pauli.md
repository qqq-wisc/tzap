# Implementing `PhaseFoldPauli`

For a formal view of Pauli-axis propagation as an abstract interpretation,
including its abstract domain, gate transformers, and soundness theorem, see
[Pauli-axis propagation as abstract interpretation](pauli-axis-abstract-interpretation.md).

`PhaseFoldPauli` is a circuit-level Pauli-rotation folding pass. It merges T,
T†, and Rz rotations that have the same Pauli axis after the intervening
Clifford circuit is taken into account.

The implementation uses random 128-bit fingerprints to find possible equal
axes quickly. Fingerprints never decide whether a rewrite is valid. Before a
rewrite, the pass compares the complete Pauli axes and checks every live
intervening rotation exactly. Randomness can therefore affect running time and,
if a work budget is exhausted, how many optimizations are found. It cannot make
the pass accept an invalid rewrite.

Every optimization level (`-O1` to `-Osuper`) runs it right after
`PhaseFoldRand`. To run it alone:

```bash
tzap input.qasm --passes PhaseFoldPauli -o output.qasm
```

It folds T, T†, and finite-angle Rz rotations across H, X, Z, S, S†, CX, and
CZ, and handles every other gate:

- **CCX and CCZ** stay in the circuit. Like `to_pbc`, the pass models each as
  seven Pauli rotations about products of its operands' frame rows (the Z rows
  of the controls and the X row of a CCX's target, or its Z row for CCZ).
  Those rotations are never folded, but they block any fold they anticommute
  with. A T on a Toffoli's control therefore folds across it; one on a CCX's
  target does not.
- **Measurements and resets** are barriers: no rotation folds across one.
  Folds on either side still happen.

A malformed circuit (a qubit out of range, a repeated operand, or a non-finite
Rz angle) is returned unchanged.

## 1. Algebra behind the rewrite

Write a Pauli rotation as

$$
R_P(\theta)=\exp(-i\theta P/2).
$$

If two rotations have the same axis and that axis commutes with every rotation
between them, the earlier rotation can move to the later site:

$$
R_P(\alpha)\;V\;R_P(\beta)
=V\;R_P(\alpha+\beta),
$$

provided every Pauli rotation in $V$ commutes with $P$. Clifford gates do not
appear in this condition because the pass continuously changes the Pauli frame
as it walks through them.

Axes are unsigned for matching: $P$ and $-P$ have the same packed axis. Their
signs are stored separately because

$$
R_{-P}(\theta)=R_P(-\theta).
$$

Suppose the earlier and current frame rows are $s_iW$ and $s_jW$, where
$s_i,s_j\in\{-1,+1\}$. At the current physical Z-rotation site, the replacement
angle is

$$
\theta'=s_i s_j\theta_i+\theta_j.
$$

In the code, exact angles are counted in eighths of a turn (multiples of
$\pi/4$), the `eighths` field of `Angle`.

### Commuting example

The Rz on qubit 1 commutes with Z on qubit 0, so the two T gates merge into S:

```qasm
t q[0];
rz(0.3) q[1];
t q[0];
```

becomes

```qasm
rz(0.3) q[1];
s q[0];
```

### Blocking example

In the following circuit, the middle Rz has frame axis X on qubit 0 because it
is surrounded by H gates. X anticommutes with the Z axis of the outer T gates,
so the fold is rejected:

```qasm
t q[0];
h q[0];
rz(0.3) q[0];
h q[0];
t q[0];
```

### Opposite-sign example

X conjugates Z to -Z. The two T gates below therefore cancel:

```qasm
t q[0];
x q[0];
t q[0];
```

and the result is just `x q[0]`.

## 2. State maintained during the scan

The implementation makes one left-to-right pass over the gates and maintains
the following state.

### Exact Clifford frame

For each qubit, store the image of its X and Z generator under the Clifford
prefix seen so far. A row represents

$$
i^p X^x Z^z,
$$

where `x` and `z` are packed `u64` bit vectors and `p` is in
$\{0,1,2,3\}$. With

```text
L = max(1, ceil(number_of_qubits / 64))
```

one row occupies `2 * L` words: all X words followed by all Z words. Keeping
the phase is necessary to distinguish $P$ from $-P$.

Ordered row multiplication is:

```text
phase(a * b) = phase(a) + phase(b) + 2 parity(z_a & x_b)  (mod 4)
x(a * b)     = x_a xor x_b
z(a * b)     = z_a xor z_b
```

Update the generator rows for each Clifford gate:

```text
H(q):       swap X[q], Z[q]
X(q):       negate Z[q]
Z(q):       negate X[q]
S(q):       X[q] = X[q] * Z[q], with the convention's phase correction
Sdg(q):     X[q] = X[q] * Z[q], with the inverse phase correction
CX(c,t):    X[c] = X[c] * X[t]
             Z[t] = Z[t] * Z[c]
CZ(c,t):    X[c] = X[c] * Z[t]
             X[t] = X[t] * Z[c]
```

The concrete code uses phase corrections 3 for S and 1 for S† after row
multiplication. An independent implementation can instead use any standard
signed stabilizer-tableau convention, as long as multiplication, generator
updates, and sign extraction all use the same convention.

At a T, T†, or non-Clifford Rz on qubit `q`, copy `Z[q]`. Its packed X/Z bits
are the exact unsigned rotation axis; its tableau phase determines the sign.

### Candidate fingerprints

Assign independent random 128-bit labels to the initial X and Z generator of
each qubit. The fingerprint of an unsigned Pauli is the XOR of the labels for
all set X and Z bits. Maintain fingerprints for the current frame generators
with the unsigned versions of the same Clifford updates:

```text
H(q):       swap x_hash[q], z_hash[q]
S/Sdg(q):   x_hash[q] ^= z_hash[q]
CX(c,t):    x_hash[c] ^= x_hash[t]
             z_hash[t] ^= z_hash[c]
CZ(c,t):    x_hash[c] ^= z_hash[t]
             x_hash[t] ^= z_hash[c]
X/Z:        no unsigned-axis change
```

The candidate fingerprint at a rotation on `q` is `z_hash[q]`. Equal exact
axes always have equal fingerprints because the map is linear. Consequently,
fingerprinting cannot hide a legal fold. Unequal axes may collide, so every
candidate is confirmed by comparing all `2 * L` packed words.

Use a hash map from fingerprint to the newest event with that fingerprint.
Each event also stores `prev_hash`, producing a chain for collisions and older
events. Walk that chain until the newest live event with an exactly equal axis
is found.

### Rotation events

Store one event for each live non-Clifford rotation:

```text
Event {
    gate:       original gate index,
    offset:     start of its packed axis in the axis arena,
    sign:       +1 or -1,
    angle:      exact multiples of pi/4 plus an f64 residual,
    prev_hash:  previous event in the fingerprint chain,
    live:       whether the event still exists,
}
```

All packed axes live in one flat word arena. This avoids a separate allocation
for every rotation.

An event is live until it is either retained in the output or folded into a
later event. Only live events matter when checking whether two rotations can
commute together.

### Exact angles

Represent an angle as:

```text
Angle {
    eighths: u8,    // multiples of pi/4 (eighths of a turn), modulo 8
    residual: f64,
}
```

T is `(1, 0)`, T† is `(7, 0)`, and exact multiples of $\pi/4$ from Rz are
placed in `eighths`. Other Rz angles are kept in `residual`. This prevents a
long sequence of Clifford+T folds from acquiring floating-point error.

Only classify an Rz as an exact multiple of $\pi/4$ when its stored `f64` value is
exactly equal to `k * (pi/4)`. Do not use an epsilon: rounding a nearby angle
would change the circuit.

## 3. The scan algorithm

The following pseudocode gives the complete control flow. Details such as
allocation caps and edit encoding are omitted here and described below.

```text
if circuit contains an unsupported gate:
    return circuit unchanged

if a cheap fingerprint-only prepass finds no repeated rotation fingerprint:
    return circuit unchanged

initialize exact Clifford frame
initialize fingerprint frame
initialize empty event list, axis arena, and fingerprint map

for each gate at index g:
    if gate is Clifford:
        update both frames
        continue

    convert T/Tdg/Rz to Angle

    if the angle is Clifford:
        update both frames as though its S, Sdg, or Z were here
        keep the original gate in the output
        continue

    hash = current Z fingerprint on the gate's qubit
    row  = current exact signed Z row on the gate's qubit

    candidate = newest live event in fingerprint_map[hash]
    walk prev_hash until candidate.axis == row.unsigned_axis

    if an exact candidate exists:
        scan every live event strictly between candidate and the current event
        if candidate.axis commutes with every intervening axis:
            merge the signed angles
            mark the earlier gate for deletion
            mark the earlier event dead
            place the merged rotation at the current gate index

            if the merged angle is Clifford:
                update both frames with the emitted Clifford
            else:
                append a new live event at the current position
            continue

    append the current rotation as a new live event

rebuild the circuit from the recorded edits
```

### Exact interval check

For packed axes $P=(x_P,z_P)$ and $Q=(x_Q,z_Q)$, compute the symplectic
commutation bit

$$
\omega(P,Q)
=x_P\cdot z_Q+z_P\cdot x_Q\pmod 2.
$$

In packed words:

```text
parity = 0
for word in 0 .. L:
    parity ^= popcount((xP[word] & zQ[word])
                       xor (zP[word] & xQ[word])) & 1
anticommutes = parity != 0
```

Starting immediately after the earlier event, visit every live event up to the
current end of the event list. Reject the fold as soon as one axis
anticommutes. Accept only after the complete interval has been checked.

It is sufficient to try only the newest live event with the same exact axis.
If an intervening event blocks that candidate, the same event also lies after
every older candidate with that axis and blocks those candidates as well.

## 4. Handling a successful fold

When two events merge:

1. Mark the earlier source gate for deletion.
2. Mark its event dead.
3. Record the merged angle as the replacement at the current gate.
4. If the result is still non-Clifford, append a new event at the current
   position. This allows it to fold again later.
5. If the result is S, S†, or Z, do not append an event. Apply that Clifford to
   both frames because it now affects every later rotation axis.
6. If the result is the identity, append nothing and leave the frame unchanged.

Edits are recorded during analysis instead of mutating the circuit in place.
A final linear pass copies unchanged gates, skips deleted gates, and emits each
replacement. This keeps original gate indices stable throughout the analysis.

## 5. Optimizations used by the implementation

### Fingerprint-only prepass

Most circuits or pipeline stages have no repeated Pauli axis. A cheap first
scan maintains only the fingerprint frame and a set of seen rotation hashes.
If no hash repeats, no exact axis can repeat, so the original circuit is
returned without allocating the $O(n^2)$ exact tableau or the axis arena.

The prepass treats exact Clifford-angle Rz gates as Clifford frame updates,
matching the main scan.

### Collision-chain compression

The fingerprint map points into a linked history. When an event is folded
away, later lookups skip it. The lookup rewrites traversed dead links to point
directly at the next live entry, so a long sequence of folds does not repeatedly
walk the same dead prefix.

### Support index for the interval check

Each event's axis has a 64-bit support signature: bit `q / L` for each qubit
`q` it acts on (the same signature the PBC optimizer uses). Per bit, the pass
keeps the ids of the events touching it, in order. Two axes can anticommute
only if they share a qubit, so the interval check visits, for each bit of the
current axis, only the events after the candidate in that bit's list, found by
binary search. Rotations on unrelated qubits cost nothing.

### Flat axis storage

All exact axes are appended to one `Vec<u64>`, `2L` words each, so event `i`'s
axis starts at word `2Li` and events need not store offsets. This reduces
allocator traffic and keeps word comparisons and symplectic products
contiguous in memory.

### Early exits

The pass stops an interval scan at the first anticommuting event. It also
returns the input immediately for malformed circuits, fewer than two
rotations, an impossible exact-tableau allocation, or no repeated candidate
fingerprint.

## 6. Correctness guarantee

Every performed fold satisfies two exact predicates:

1. The full packed unsigned axes of the two rotations are equal.
2. That axis symplectically commutes with every live intervening rotation.

The signed angle formula then applies the Pauli-rotation identity from Section
1. Clifford frame updates are exact stabilizer-tableau operations. Therefore,
assuming those tableau operations and the circuit gate semantics are
implemented correctly, every rewrite preserves the circuit unitary up to
global phase.

The random fingerprints introduce no probability of semantic failure:

- Equal axes always receive equal fingerprints, so they cannot be missed by
  the candidate index.
- A collision between unequal axes is rejected by the packed-axis comparison.
- A collision may consume extra work and can cause the work budget to be
  reached earlier. That changes optimization completeness, not correctness.

Arbitrary Rz angles use `f64` addition and reduction modulo $2\pi$, so their
numeric values have ordinary deterministic floating-point rounding. Exact
Clifford+T angles remain in the integer `eighths` component. The pass never
rounds a merely nearby Rz angle to a Clifford+T value.

## 7. Resource bounds and incomplete optimization

Let $n$ be the number of qubits, $R$ the number of non-Clifford rotation
events, and $L=\lceil n/64\rceil$.

- The exact frame stores `2n` rows of `2L` words: $O(n^2/64)$ words.
- Each retained event axis stores `2L` words: $O(Rn/64)$ words.
- Each Clifford row multiplication and each exact axis comparison costs
  $O(L)$ word operations.
- One interval check costs $O(kL)$ for $k$ intervening events that share a
  support-signature bit with the axis.
- Without a bound, the worst case is quadratic in the number of rotations:
  $O(R^2n/64)$.

Each fold attempt (collision-chain comparisons plus intervening-event checks)
is bounded by $2^{20}$ steps. An attempt that runs out is treated as blocked
and the scan continues, so the pass is $O(G L + R \cdot 2^{20} L)$ for $G$
gates: linear in the circuit, with a large constant in the worst case. Exact
axis storage is capped at 256 MiB; reaching that stops further folding, and
the circuit is rebuilt from the folds already proved valid. These limits may
leave optimizations behind; they never weaken a completed fold.

## 8. Suggested implementation tests

An independent implementation should include at least these tests:

1. Two T gates on one qubit become S.
2. Equal axes fold across a commuting rotation on another qubit.
3. An anticommuting frame-axis rotation blocks a fold.
4. Opposite signed axes subtract their angles.
5. A merge that produces S, S†, or Z updates the later Clifford frame.
6. A merge that produces identity creates no event.
7. Arbitrary Rz angles merge without epsilon classification.
8. Axes spanning qubits 63 and 64 work across packed-word boundaries.
9. Forced fingerprint collisions never authorize a rewrite.
10. CCX and CCZ block the folds they anticommute with and no others, using
    the frame at their position.
11. Measurements and resets are barriers, and circuits with mid-circuit
    measurements keep their exact channel.
12. An attempt that exhausts its budget does not stop later folds.
13. Random small circuits remain unitarily equivalent before and after the
    pass, including circuits with arbitrary finite Rz angles.
