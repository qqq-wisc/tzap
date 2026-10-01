# Pauli-axis propagation as abstract interpretation

These notes develop rotation folding for Clifford+T circuits as an abstract
interpretation. Section 1 states the problem. Section 2 fixes notation for
Pauli operators. Section 3 defines the **concrete domain**: what a concrete
value is, its type, and the exact transformer of every Clifford+T gate,
and of the Toffoli gates CCX and CCZ.
Sections 4 and 5 define the **abstraction** — sets of Pauli axes — with its
Galois connection, abstract transformers, and soundness theorem. Section 6
shows where the abstraction loses precision, Section 7 turns an analysis
result into a sound rewrite, and Section 8 relates the analysis to the domains
used in practice, including the one `PhaseFoldPauli` implements. Section 10
proves that, on Clifford+T circuits, Pauli folding strictly subsumes symbolic
phase folding. Section 11 extends the analysis to measurements, resets,
classically controlled gates, and loops.

Every definition is followed by a small example. All the examples have been
checked numerically.

**Conventions.** A circuit $C=g_1;g_2;\ldots;g_m$ is written in execution
order; its unitary is $U_C=U_{g_m}\cdots U_{g_1}$. Qubits are numbered from
0. Equalities between circuits hold up to a global phase.

## 1. The problem

A *Pauli rotation* about a Pauli operator $P$ by angle $\theta$ is

$$
R_P(\theta)=\exp\!\left(-\tfrac{i\theta}{2}P\right)
=\cos\tfrac{\theta}{2}\,I-i\sin\tfrac{\theta}{2}\,P .
$$

The T gate is a Pauli rotation: $T=R_Z(\pi/4)$ up to global phase. Likewise
$T^\dagger=R_Z(-\pi/4)$, $S=R_Z(\pi/2)$, and $Z=R_Z(\pi)$.

Rotation folding asks when two rotations separated by a circuit segment can be
combined. The answer rests on one identity.

**Lemma 1.1 (transport).** For every unitary $U$, Pauli $P$, and angle
$\theta$,

$$
U\,R_P(\theta)\,U^\dagger=R_{UPU^\dagger}(\theta).
$$

*Proof.* $U\exp(M)U^\dagger=\exp(UMU^\dagger)$ for every matrix $M$. $\square$

Consequently, if a segment $C$ satisfies $U_CPU_C^\dagger=sQ$ with
$s\in\{-1,+1\}$ and $Q$ a Pauli, then $U_C R_P(\alpha)=R_Q(s\alpha)\,U_C$, and

$$
R_P(\alpha);\ C;\ R_Q(\beta)\quad\equiv\quad C;\ R_Q(s\alpha+\beta).
\tag{1.1}
$$

The first rotation moves through $C$ and merges with the second.

**Example 1.2.** In `t q[0]; cx q[0],q[1]; t q[0]`, the segment is the CX and
$P=Q=Z_0$. The CX leaves $Z_0$ unchanged ($s=+1$), so by (1.1) the circuit
equals `cx q[0],q[1]; s q[0]`: two T gates become one S.

**Example 1.3.** In `t q[0]; x q[0]; t q[0]`, conjugation by X gives
$XZX=-Z$, so $s=-1$ and the merged angle is $-\pi/4+\pi/4=0$. The circuit
equals `x q[0]`: the T gates cancel.

**Example 1.4.** In `t q[1]; cx q[0],q[1]; t q[1]`, the CX maps $Z_1$ to
$Z_0Z_1$. The first rotation emerges about $Z_0Z_1$, not about the second
rotation's axis $Z_1$, so (1.1) does not apply.

Folding therefore reduces to one question: **what is $U_CPU_C^\dagger$?** Its
exact computation is the concrete semantics (Section 3). The abstraction
(Sections 4–5) approximates it cheaply and soundly.

*A remark on terminology.* The map $A\mapsto UAU^\dagger$ moves an operator
*forward* through $U$. The Heisenberg picture proper evolves observables by
$A\mapsto U^\dagger AU$, the inverse map. We use the forward map throughout,
since it is the one transport requires.

## 2. Pauli operators

**Definition 2.1 (Pauli strings).** The single-qubit Paulis are

$$
I=\begin{pmatrix}1&0\\0&1\end{pmatrix},\quad
X=\begin{pmatrix}0&1\\1&0\end{pmatrix},\quad
Y=\begin{pmatrix}0&-i\\i&0\end{pmatrix},\quad
Z=\begin{pmatrix}1&0\\0&-1\end{pmatrix}.
$$

An $n$-qubit *Pauli string* is a tensor product $P=P_0\otimes\cdots\otimes
P_{n-1}$ with each $P_q\in\{I,X,Y,Z\}$; $P_q$ is its *letter* at qubit $q$.
The set of Pauli strings is

$$
\mathbb P_n=\{I,X,Y,Z\}^{\otimes n},\qquad |\mathbb P_n|=4^n .
$$

We name a string by its non-identity letters: $X_q$ is $X$ at qubit $q$ and
$I$ elsewhere, and $Z_0Y_1$ is $Z\otimes Y$ on two qubits. We call the
elements of $\mathbb P_n$ *axes*. A *signed Pauli* is $\pm P$ with
$P\in\mathbb P_n$.

**Example 2.2.** On two qubits, $\mathbb P_2$ has 16 axes: $II$, $XI$, …,
$ZZ$. In the naming convention, $XI=X_0$, $IZ=Z_1$, and $ZY=Z_0Y_1$.

Single-qubit Paulis multiply as $XY=iZ$, $YZ=iX$, and $ZX=iY$; each squares
to $I$; distinct non-identity letters anticommute ($XY=-YX$). Strings
multiply letter by letter, so the product of two axes is an axis times a
phase in $\{\pm1,\pm i\}$.

**Definition 2.3 (commutation).** Two axes $P,Q\in\mathbb P_n$ either commute
($PQ=QP$) or anticommute ($PQ=-QP$). They anticommute exactly when the number
of qubits at which both letters are non-identity and different is odd.

**Example 2.4.** $Z_0$ and $X_0Z_1$ differ at qubit 0 only, so they
anticommute. $X_0X_1$ and $Z_0Z_1$ differ at both qubits, so they commute.

**Fact 2.5 (Pauli basis).** The axes form a basis of the $2^n\times2^n$
matrices, orthogonal under the Hilbert–Schmidt inner product:
$\operatorname{tr}(PQ)=2^n$ if $P=Q$ and $0$ otherwise. A Hermitian matrix
$A$ therefore has a unique expansion $A=\sum_P a_PP$ with real coefficients
$a_P=\operatorname{tr}(PA)/2^n$.

## 3. The concrete domain

### 3.1 Concrete values and their type

**Definition 3.1 (concrete values).** A *concrete value* is a Hermitian
operator on $n$ qubits. Their set is

$$
\mathcal H_n=\{\,A\in\mathbb C^{2^n\times2^n} : A=A^\dagger\,\},
$$

a real vector space of dimension $4^n$. The value a rotation axis starts as
is the axis itself, $A=P$.

By Fact 2.5, every $A\in\mathcal H_n$ is determined by its $4^n$ real
*coordinates* in the Pauli basis,

$$
A=\sum_{P\in\mathbb P_n}a_P\,P,\qquad a_P=\tfrac{1}{2^n}\operatorname{tr}(PA).
$$

Throughout, the capital letter is the operator and the lowercase letter with
a subscript is its coordinate: $a_P$ is the coefficient of axis $P$ in $A$,
and $b_P$ the coefficient of $P$ in $B$. The operator is the concrete value;
the coordinates are how we compute with it. (Equivalently, one may take the
coordinate vectors $\mathbb R^{\mathbb P_n}$ as the concrete domain: the map
$A\mapsto(a_P)_P$ is a linear bijection, so nothing below depends on the
choice.)

**Example 3.2.** On one qubit, $\mathcal H_1$ consists of the operators
$a_II+a_XX+a_YY+a_ZZ$ with real coordinates. The operator
$A=\tfrac{1}{\sqrt2}(X+Y)=\tfrac1{\sqrt2}\begin{pmatrix}0&1-i\\1+i&0
\end{pmatrix}$ has coordinates $a_X=a_Y=\tfrac1{\sqrt2}$ and $a_I=a_Z=0$.
The start value for the axis $Z$ is $A=Z$, with $a_Z=1$ and all other
coordinates 0.

**Definition 3.3 (support).** The *support* of a value is the set of axes
with nonzero coordinate:

$$
\operatorname{supp}(A)=\{P\in\mathbb P_n: a_P\neq0\}.
$$

**Example 3.4.** $\operatorname{supp}\bigl(\tfrac1{\sqrt2}(X+Y)\bigr)=\{X,Y\}$
and $\operatorname{supp}(-Z_0Z_1)=\{Z_0Z_1\}$.

### 3.2 Concrete transformers

**Definition 3.5 (concrete transformer).** A gate $g$ with unitary $U_g$ has
the concrete transformer

$$
F_g:\mathcal H_n\to\mathcal H_n,\qquad
F_g(A)=U_gAU_g^\dagger .
$$

A segment $C=g_1;\ldots;g_m$ has
$F_C=F_{g_m}\circ\cdots\circ F_{g_1}$, so $F_C(A)=U_CAU_C^\dagger$.

Each $F_g$ is linear, invertible, maps Hermitian operators to Hermitian
operators, and preserves the norm $\sum_Pa_P^2=\operatorname{tr}(A^2)/2^n$.
Because it is linear, it is determined by its action on the axes:

$$
F_g\Bigl(\sum_Pa_PP\Bigr)=\sum_Pa_P\,F_g(P).
$$

So to transform a value, transform each axis and add up the results with the
value's coordinates as weights. The rest of this section gives $F_g(P)$ for
every gate of the Clifford+T set

$$
\mathcal G=\{X,\ Z,\ H,\ S,\ S^\dagger,\ \mathrm{CX},\ \mathrm{CZ},\ T,\
T^\dagger\},
$$

and for the Rz generalization of T.

#### Clifford gates: signed permutations

**Definition 3.6 (Clifford).** A gate is *Clifford* if it maps every axis to
a signed axis. Then there are a permutation $\pi_g$ of $\mathbb P_n$ and signs
$\sigma_g(P)\in\{-1,+1\}$ with

$$
F_g(P)=\sigma_g(P)\,\pi_g(P).
$$

In coordinates, a Clifford moves each coefficient to a new axis and possibly
flips its sign: if $B=F_g(A)$, then

$$
b_{\pi_g(P)}=\sigma_g(P)\,a_P\qquad\text{for every axis }P.
$$

The single-qubit Cliffords of $\mathcal G$, applied to qubit $q$, change only
the letter at $q$ (and every gate maps $I$ to $I$):

| $g$ | $X\mapsto$ | $Y\mapsto$ | $Z\mapsto$ |
|---|---|---|---|
| $X$ | $X$ | $-Y$ | $-Z$ |
| $Z$ | $-X$ | $-Y$ | $Z$ |
| $H$ | $Z$ | $-Y$ | $X$ |
| $S$ | $Y$ | $-X$ | $Z$ |
| $S^\dagger$ | $-Y$ | $X$ | $Z$ |

The two-qubit Cliffords, with control $c$ and target $t$:

| $g$ | $X_c\mapsto$ | $Y_c\mapsto$ | $Z_c\mapsto$ | $X_t\mapsto$ | $Y_t\mapsto$ | $Z_t\mapsto$ |
|---|---|---|---|---|---|---|
| CX | $X_cX_t$ | $Y_cX_t$ | $Z_c$ | $X_t$ | $Z_cY_t$ | $Z_cZ_t$ |
| CZ | $X_cZ_t$ | $Y_cZ_t$ | $Z_c$ | $Z_cX_t$ | $Z_cY_t$ | $Z_t$ |

Conjugation is multiplicative, $U(PQ)U^\dagger=(UPU^\dagger)(UQU^\dagger)$,
so the image of any string is the product of the images of its letters.

**Example 3.7.** Under CX with control 0 and target 1,

$$
F_{\mathrm{CX}}(Y_0Y_1)
=(Y_0X_1)(Z_0Y_1)=(Y_0Z_0)(X_1Y_1)=(iX_0)(iZ_1)=-X_0Z_1 .
$$

So $\pi_{\mathrm{CX}}(Y_0Y_1)=X_0Z_1$ and $\sigma_{\mathrm{CX}}(Y_0Y_1)=-1$.

#### T and T†: rotations in the X–Y plane

T is not Clifford: it maps X to a combination of two axes.

$$
TXT^\dagger=\tfrac{1}{\sqrt2}(X+Y),\qquad
TYT^\dagger=\tfrac{1}{\sqrt2}(Y-X),\qquad
TZT^\dagger=Z .
$$

For $T^\dagger$, replace $Y$ by $-Y$ on the right-hand sides:
$T^\dagger XT=\tfrac1{\sqrt2}(X-Y)$ and $T^\dagger YT=\tfrac1{\sqrt2}(X+Y)$.

In coordinates, T on qubit $q$ acts on pairs of axes that differ only in the
letter at $q$. For an axis $P$, write $P[q{:=}L]$ for $P$ with its letter at
$q$ replaced by $L$.

**Definition 3.8 (T transformer, in coordinates).** Let $B=F_{T_q}(A)$.

- If the letter of $P$ at $q$ is $I$ or $Z$, then $b_P=a_P$.
- For each pair $P_X=P[q{:=}X]$, $P_Y=P[q{:=}Y]$, with $x=a_{P_X}$ and
  $y=a_{P_Y}$,

$$
\begin{pmatrix}b_{P_X}\\ b_{P_Y}\end{pmatrix}
=\frac{1}{\sqrt2}\begin{pmatrix}1&-1\\1&\phantom{-}1\end{pmatrix}
\begin{pmatrix}x\\ y\end{pmatrix}.
$$

That is, T rotates each (X, Y) coefficient pair by $45^\circ$. $T^\dagger$
rotates it by $-45^\circ$, using the transposed matrix.

**Example 3.9.** Applying T to $A=X+Z$ on one qubit: the $Z$ coefficient is
unchanged, and the pair $(x,y)=(1,0)$ becomes $(\tfrac1{\sqrt2},\tfrac1{\sqrt2})$.
So $F_{T}(X+Z)=\tfrac1{\sqrt2}X+\tfrac1{\sqrt2}Y+Z$.

**Generalization (Rz and Pauli rotations).** $R_Z(\theta)$ on qubit $q$
rotates each pair by $\theta$:
$(x,y)\mapsto(x\cos\theta-y\sin\theta,\ x\sin\theta+y\cos\theta)$. T is
$\theta=\pi/4$. More generally, for a rotation about an axis $Q$,

$$
R_Q(\theta)\,P\,R_Q(\theta)^\dagger=
\begin{cases}
P, & PQ=QP,\\
\cos\theta\,P+\sin\theta\,(iPQ), & PQ=-QP,
\end{cases}
\tag{3.1}
$$

where $iPQ$ is a signed axis whenever $P$ and $Q$ anticommute. For
$R_Z(\theta)$ and $P=X$, (3.1) gives $\cos\theta\,X+\sin\theta\,Y$, since
$iXZ=Y$.

#### Toffoli gates

CCZ on qubits $a,b,c$ is diagonal: $\mathrm{CCZ}\,|z\rangle=(-1)^{z_az_bz_c}|z\rangle$.
It is not Clifford, and unlike T it can spread one axis over four.

**Proposition 3.10 (CCZ transformer).** Let $P$ be an axis, and let
$x\in\{0,1\}^{\{a,b,c\}}$ mark the operands where its letter is X or Y.

- If $x=0$ (every operand letter is I or Z), then $F_{\mathrm{CCZ}}(P)=P$.
- Otherwise $F_{\mathrm{CCZ}}(P)$ is a sum of exactly four signed axes with
  coefficients $\pm\tfrac12$, namely $P\,M$ (up to phase) for the four
  Z-strings $M$ of the table below, up to relabeling the operands.

| X/Y letters on | $M$ ranges over |
|---|---|
| one operand, $c$ | $I,\ Z_a,\ Z_b,\ Z_aZ_b$ |
| two operands, $a$ and $b$ | $I,\ Z_c,\ Z_aZ_b,\ Z_aZ_bZ_c$ |
| all three | $I,\ Z_aZ_b,\ Z_aZ_c,\ Z_bZ_c$ |

*Proof.* Write $D=\mathrm{CCZ}$ and $f(z)=z_az_bz_c$. Z letters commute with
$D$, and X letters flip bits, so $PDP=\operatorname{diag}((-1)^{f(z\oplus x)})$
and $DPD=P\,\Delta$ with $\Delta=\operatorname{diag}((-1)^{f(z\oplus x)+f(z)})$.
For $x\neq0$ the exponent is a product $uv$ of two independent affine
functions: $z_az_b$ for $x=e_c$; $z_c(z_a+z_b+1)$ for $x=e_a+e_b$; and
$(z_a+z_c+1)(z_b+z_c+1)$ for $x=e_a+e_b+e_c$. Then

$$
(-1)^{uv}=\tfrac12\bigl(1+(-1)^u+(-1)^v-(-1)^{u+v}\bigr),
$$

and each $(-1)^{w}$ with $w$ affine is a signed Z-string. Multiplying by $P$
gives the four axes. $\square$

The proposition has also been checked exhaustively over all 64 letter
patterns on the three operands.

**Example 3.11.** For X on the target of a CCZ on qubits 0, 1, 2,

$$
F_{\mathrm{CCZ}}(X_2)=\tfrac12\bigl(X_2+Z_0X_2+Z_1X_2-Z_0Z_1X_2\bigr),
$$

while $F_{\mathrm{CCZ}}(Z_0)=Z_0$ and $F_{\mathrm{CCZ}}(Z_0Z_1Z_2)=Z_0Z_1Z_2$.

CCX with controls $a,b$ and target $t$ is $H_t;\mathrm{CCZ};H_t$, so
$F_{\mathrm{CCX}}=F_{H_t}\circ F_{\mathrm{CCZ}}\circ F_{H_t}$. An axis is
fixed by CCX exactly when its control letters are I or Z and its target
letter is I or X; every other axis spreads over four. For example,
$F_{\mathrm{CCX}}(Z_2)=\tfrac12(Z_2+Z_0Z_2+Z_1Z_2-Z_0Z_1Z_2)$ and
$F_{\mathrm{CCX}}(X_2)=X_2$.

CCZ is also a product of seven commuting Pauli rotations,

$$
\mathrm{CCZ}=
R_{Z_a}(\tfrac\pi4)\,R_{Z_b}(\tfrac\pi4)\,R_{Z_c}(\tfrac\pi4)\,
R_{Z_aZ_b}(-\tfrac\pi4)\,R_{Z_aZ_c}(-\tfrac\pi4)\,R_{Z_bZ_c}(-\tfrac\pi4)\,
R_{Z_aZ_bZ_c}(\tfrac\pi4),
$$

up to global phase, so $F_{\mathrm{CCZ}}$ is also the composition of seven
rotation transformers (3.1). This is how `to_pbc` represents a Toffoli.

### 3.3 A worked example

Propagate the axis $X_1$ through the segment `cx q[0],q[1]; t q[1];
cx q[0],q[1]`:

| after | value | support |
|---|---|---|
| start | $X_1$ | $\{X_1\}$ |
| `cx q[0],q[1]` | $X_1$ | $\{X_1\}$ |
| `t q[1]` | $\tfrac1{\sqrt2}X_1+\tfrac1{\sqrt2}Y_1$ | $\{X_1,Y_1\}$ |
| `cx q[0],q[1]` | $\tfrac1{\sqrt2}X_1+\tfrac1{\sqrt2}Z_0Y_1$ | $\{X_1,Z_0Y_1\}$ |

The first CX fixes $X_t$. The T mixes the letter at qubit 1. The second CX
fixes $X_1$ and maps $Y_1$ to $Z_0Y_1$. The result is not a signed axis, so a
rotation about $X_1$ cannot be moved through this segment.

### 3.4 Two properties of Clifford+T values

**Proposition 3.12 (exact coefficients).** Propagating an axis through a
Clifford+T circuit yields coefficients in the ring
$\mathbb Z[\tfrac1{\sqrt2}]=\{\,u+v\sqrt2 : u,v\text{ dyadic rationals}\,\}$.

*Proof.* The start value has coefficients in $\{0,1\}$. Cliffords only move
coefficients and flip signs. T and $T^\dagger$ replace a pair $(x,y)$ by
$\tfrac1{\sqrt2}(x\mp y,\ x\pm y)$. CCX and CCZ add up coefficients
multiplied by $\pm\tfrac12$. All these operations preserve the ring.
$\square$

This makes exact arithmetic on Clifford+T coefficients possible, and in
particular an exact test for zero (Section 8.3).

**Proposition 3.13 (unit norm).** Propagating an axis keeps
$\sum_Pa_P^2=1$.

*Proof.* Each transformer preserves $\operatorname{tr}(A^2)$, and a single
axis has $\sum_Pa_P^2=1$. $\square$

**Example 3.14.** In Section 3.3, the final value has
$(\tfrac1{\sqrt2})^2+(\tfrac1{\sqrt2})^2=1$.

### 3.5 The collecting domain

An abstraction describes *sets* of concrete values, so we lift the semantics
to sets.

**Definition 3.15 (collecting domain).** The concrete domain is
$\mathcal D=\mathcal P(\mathcal H_n)$, the sets of values, ordered by
$\subseteq$. Transformers lift pointwise:
$F_g(\mathcal A)=\{F_g(A):A\in
\mathcal A\}$.

An analysis of one axis $P$ starts from the singleton $\{P\}$. Sets matter
because an abstract value stands for all the concrete values consistent with
it.

## 4. The abstract domain: Pauli supports

The abstraction keeps only *which* axes may carry a nonzero coefficient, and
forgets the coefficients themselves.

**Definition 4.1 (abstract domain).** The abstract domain is
$\mathcal D^\#=\mathcal P(\mathbb P_n)$, the sets of axes, ordered by
$\subseteq$. It is a complete lattice with

$$
\bot=\varnothing,\qquad\top=\mathbb P_n,\qquad S_1\sqcup S_2=S_1\cup S_2,
\qquad S_1\sqcap S_2=S_1\cap S_2 .
$$

An abstract value $S$ reads: *every axis outside $S$ has coefficient zero.* It
does not claim that the axes in $S$ are present.

**Example 4.2.** $S=\{X,Y\}$ describes every one-qubit operator $xX+yY$:
for instance $X$, $\tfrac1{\sqrt2}(X+Y)$, and $0$. It excludes $X+Z$.

**Definition 4.3 (abstraction and concretization).**

$$
\alpha(\mathcal A)=\bigcup_{A\in\mathcal A}\operatorname{supp}(A),
\qquad
\gamma(S)=\{A\in\mathcal H_n:\operatorname{supp}(A)\subseteq S\}.
$$

**Example 4.4.** $\alpha\bigl(\{\tfrac1{\sqrt2}(X+Y),\ Z\}\bigr)=\{X,Y,Z\}$,
and $\gamma(\{Z\})=\{zZ: z\in\mathbb R\}$, a line through $0$.

**Lemma 4.5 (Galois connection).** For all $\mathcal A\subseteq\mathcal H_n$
and $S\subseteq\mathbb P_n$,

$$
\alpha(\mathcal A)\subseteq S\iff\mathcal A\subseteq\gamma(S).
$$

*Proof.* Both sides say that every $A\in\mathcal A$ has
$\operatorname{supp}(A)\subseteq S$. $\square$

Two consequences are used below: $\alpha$ and $\gamma$ are monotone, and
$\mathcal A\subseteq\gamma(\alpha(\mathcal A))$ (abstracting loses no
possible value).

## 5. Abstract transformers

**Definition 5.1 (soundness).** An abstract transformer
$F^\#_g:\mathcal D^\#\to\mathcal D^\#$ is *sound* if
$F_g(\gamma(S))\subseteq\gamma(F^\#_g(S))$ for every $S$. By Lemma 4.5 this is equivalent to

$$
\operatorname{supp}(F_g A)\subseteq F^\#_g(\operatorname{supp}A)\qquad\text{for every }A\in\mathcal H_n .
$$

The *best* sound transformer is $\alpha\circ F_g\circ
\gamma$. Every transformer below is the best one.

Because each $F_g$ is linear, the abstract transformers act
axis by axis: $F^\#_g(S)=\bigcup_{P\in S}F^\#_g(\{P\})$.

### 5.1 Clifford gates

**Definition 5.2.** For a Clifford $g$, $F^\#_g(S)=
\{\pi_g(P):P\in S\}$. Signs are dropped.

**Example 5.3.** For CX with control 0 and target 1,
$F^\#_{\mathrm{CX}}(\{X_1,Y_1\})=\{X_1,Z_0Y_1\}$, and
$F^\#_{H}(\{X,Y\})=\{Z,Y\}$.

**Lemma 5.4 (Cliffords are exact).** For a Clifford $g$,
$\operatorname{supp}(F_g A)=F^\#_g(\operatorname{supp}A)$.

*Proof.* $F_g A=\sum_Pa_P\,\sigma_g(P)\,\pi_g(P)$. Since
$\pi_g$ is a permutation, distinct input axes reach distinct output axes, so
no two coefficients combine: each $a_P\neq0$ yields a nonzero coefficient at
$\pi_g(P)$, and nothing else. $\square$

### 5.2 T, T†, and rotations

**Definition 5.5.** For $T_q$ and $T_q^\dagger$, per axis:

$$
F^\#_{T_q}(\{P\})=F^\#_{T_q^\dagger}(\{P\})=
\begin{cases}
\{P\}, & \text{the letter of }P\text{ at }q\text{ is }I\text{ or }Z,\\
\{P[q{:=}X],\ P[q{:=}Y]\}, & \text{it is }X\text{ or }Y .
\end{cases}
$$

**Example 5.6.** $F^\#_{T_1}(\{X_1\})=\{X_1,Y_1\}$,
$F^\#_{T_1}(\{Z_0Z_1\})=\{Z_0Z_1\}$, and
$F^\#_{T_0}(\{Y_0Z_1\})=\{X_0Z_1,Y_0Z_1\}$.

For a general rotation $R_Q(\theta)$, (3.1) gives the per-axis rule

$$
F^\#_{R_Q(\theta)}(\{P\})=
\begin{cases}
\{P\}, & PQ=QP,\\
\{P : \cos\theta\neq0\}\cup\{\operatorname{axis}(iPQ):\sin\theta\neq0\},
& PQ=-QP,
\end{cases}
$$

where $\operatorname{axis}(\pm R)=R$ and $\{P:c\}$ means $\{P\}$ if the
condition $c$ holds and $\varnothing$ otherwise. For
$R_Z(\theta)$ this means: if $\theta$ is a multiple of $\pi$, every axis is
kept; if $\theta$ is an odd multiple of $\pi/2$, X and Y are exchanged at
qubit $q$; otherwise both are produced. An implementation that cannot prove a
sine or cosine is exactly zero must include the branch.

**Lemma 5.7 (rotations are sound and best).** Definition 5.5 and the
general rule are sound, and each is the best transformer.

*Proof.* Soundness: by linearity, every output axis with nonzero coefficient
comes from some input axis $P\in\operatorname{supp}(A)$ through (3.1), and the
rule includes every axis that (3.1) can produce from $P$. Coefficients from
different input axes may cancel, which only shrinks the concrete support.
Best: for each $P\in S$, the value $P$ itself lies in $\gamma(S)$, and (3.1)
applied to it produces every axis the rule lists. $\square$

### 5.3 Toffoli gates

**Definition 5.8.** For CCZ, per axis: $F^\#_{\mathrm{CCZ}}(\{P\})=\{P\}$ if
every operand letter of $P$ is I or Z, and otherwise the four axes of
Proposition 3.10. For CCX,
$F^\#_{\mathrm{CCX}}=F^\#_{H_t}\circ F^\#_{\mathrm{CCZ}}\circ F^\#_{H_t}$.

By Proposition 3.10 this is exactly $\operatorname{supp}(F_{\mathrm{CCZ}}(P))$
for each axis, so it is sound and best, by the argument of Lemma 5.7.

**Example 5.9.** $F^\#_{\mathrm{CCZ}}(\{X_2\})=\{X_2,Z_0X_2,Z_1X_2,Z_0Z_1X_2\}$.
Composing the seven rotation transformers instead gives eight axes: these
four, and the same four with $Y_2$ in place of $X_2$. The $Y_2$ terms have
coefficient zero in the concrete value, since the seven rotations' contributions
cancel. As in Section 6, the best transformer of a gate is more precise
than the composition of the best transformers of its parts. An analysis
should use Definition 5.8 directly.

### 5.4 Circuits and soundness

**Definition 5.10.** A segment $C=g_1;\ldots;g_m$ has
$F^\#_C=F^\#_{g_m}\circ\cdots\circ
F^\#_{g_1}$.

**Theorem 5.11 (soundness).** For every segment $C$ and axis $P$,

$$
\operatorname{supp}\bigl(U_CPU_C^\dagger\bigr)\subseteq F^\#_C(\{P\}),
$$

and the same holds at every intermediate position of $C$.

*Proof.* By induction on $m$. For $m=0$ both sides are $\{P\}$. If the claim
holds after $k-1$ gates, the concrete value lies in $\gamma(S_{k-1})$, where
$S_{k-1}$ is the abstract value there. By Definition 5.1, its image under
$g_k$ lies in $\gamma(F^\#_{g_k}(S_{k-1}))=\gamma(S_k)$.
$\square$

**Example 5.12.** The abstract run of Section 3.3 matches the support column
of that table exactly: $\{X_1\}\to\{X_1\}\to\{X_1,Y_1\}\to\{X_1,Z_0Y_1\}$.
Here the abstraction loses nothing.

The theorem bounds one thing: the axes that may carry a nonzero coefficient in
the transported operator $U_CPU_C^\dagger$. It says nothing about the quantum
state, measurement outcomes, or the unitary $U_C$ itself, and the analysis is
deterministic: "may" means "not ruled out once coefficients are forgotten".

## 6. Precision: where the abstraction loses information

Soundness guarantees that the abstract value contains the true support. It
does not guarantee equality.

**Example 6.1.** Propagate $Z$ through `h; t; h; h; tdg; h`: a rotation by
$\pi/4$ about X followed by one by $-\pi/4$ about X. Together they are the
identity, so the exact result is $Z$.

| after | concrete value | concrete support | abstract value |
|---|---|---|---|
| start | $Z$ | $\{Z\}$ | $\{Z\}$ |
| `h` | $X$ | $\{X\}$ | $\{X\}$ |
| `t` | $\tfrac1{\sqrt2}(X+Y)$ | $\{X,Y\}$ | $\{X,Y\}$ |
| `h` | $\tfrac1{\sqrt2}(Z-Y)$ | $\{Y,Z\}$ | $\{Y,Z\}$ |
| `h` | $\tfrac1{\sqrt2}(X+Y)$ | $\{X,Y\}$ | $\{X,Y\}$ |
| `tdg` | $X$ | $\{X\}$ | $\{X,Y\}$ |
| `h` | $Z$ | $\{Z\}$ | $\{Y,Z\}$ |

At `tdg`, the pair $(x,y)=(\tfrac1{\sqrt2},\tfrac1{\sqrt2})$ becomes
$(\tfrac1{\sqrt2}(x+y),\tfrac1{\sqrt2}(y-x))=(1,0)$: the Y coefficient
cancels. The abstract transformer sees only that X and Y are present and
keeps both. The final abstract value $\{Y,Z\}$ is sound, since it contains
$\{Z\}$, but it cannot prove that the rotation about $Z$ may be transported.

The loss is inherent to the domain, not to the transformers. Each transformer
is the best one, but the composition of best transformers need not be the
best transformer of the composition: the best transformer of this whole
segment is the identity. Recovering such cancellations requires keeping
coefficients (Section 8.3).

## 7. From an analysis result to a rewrite

**Theorem 7.1 (fold).** Let $F^\#_C(\{P\})=\{Q\}$. Then
$U_CPU_C^\dagger=sQ$ for some $s\in\{-1,+1\}$, and

$$
R_P(\alpha);\ C;\ R_Q(\beta)\equiv C;\ R_Q(s\alpha+\beta).
$$

*Proof.* By Theorem 5.11, $U_CPU_C^\dagger=cQ$ for a real $c$. By Proposition
3.13, $c^2=1$. Equation (1.1) completes the proof. $\square$

The support domain proves that the result is a single axis but not its sign.
A rewrite needs $s$. An implementation obtains it from an exact signed
tableau while the segment is Clifford, or from the coefficient domain of
Section 8.3.

**Example 7.2.** Take `t q[0]; h q[0]; t q[1]; h q[0]; t q[0]`, with $P=Z_0$
and segment `h q[0]; t q[1]; h q[0]`:
$\{Z_0\}\to\{X_0\}\to\{X_0\}\to\{Z_0\}$. The T on qubit 1 has letter $I$ at
qubit 1, so it keeps $X_0$. The result is the single axis $Z_0$, with sign $+1$
by the signed tableau, so the outer T gates merge into one S at the second
site.

## 8. Refinements and practical domains

The support domain has up to $4^n$ axes per value. Practical analyses use
coarser or bounded domains, and a finer one when cancellations matter.

### 8.1 Singleton or top: the domain of `PhaseFoldPauli`

**Definition 8.1.**
$\mathcal D^\#_{\mathrm{one}}=\{\bot\}\cup\{\mathrm{One}(P):P\in\mathbb
P_n\}\cup\{\top\}$, with $\gamma(\mathrm{One}(P))=\{sP:s\in\mathbb R\}$ and
$\gamma(\top)=\mathcal H_n$. Cliffords map $\mathrm{One}(P)$ to
$\mathrm{One}(\pi_g(P))$. A rotation about $Q$ keeps $\mathrm{One}(P)$ when
$P$ and $Q$ commute and otherwise yields $\top$. Every transformer maps
$\top$ to $\top$.

This domain is the support domain with every non-singleton set widened to
$\top$, so it is sound and loses at least as much. It is what `PhaseFoldPauli`
computes. The pass tracks the Clifford part exactly with a signed tableau,
which also supplies the sign for Theorem 7.1. It declares $\top$ as soon as
the candidate axis anticommutes with one intervening rotation, which is the
rule that the axis must commute with every rotation in between.

**Example 8.2.** In Example 6.1, the first T already anticommutes with the
propagated axis X, so the analysis reaches $\top$ after `t`. For Example 7.2,
it stays at $\mathrm{One}(\cdot)$ throughout and finds the fold.

Toffolis fit the same scheme. A CCZ is the seven commuting rotations
$R_M(\pm\pi/4)$ of Section 3.2, about the Z-strings $M$ on its operands, and a
CCX the same with the target's Z replaced by X; the pass takes their axes from
the frame at the Toffoli's position. The analysis stays at $\mathrm{One}(P)$
exactly when $P$ commutes with all seven, which by Proposition 3.10 is when
$F_{\mathrm{CCZ}}$ (or $F_{\mathrm{CCX}}$) fixes $P$.
Measurements and resets lie outside the unitary semantics, so the pass
treats them as barriers.

### 8.2 Bounded supports

**Definition 8.3.** For a bound $K$,
$\mathcal D^\#_K=\{S\subseteq\mathbb P_n:|S|\le K\}\cup\{\top\}$. Apply the
support transformer, and widen to $\top$ whenever a result would exceed $K$
axes.

Widening only enlarges concretizations, so $\mathcal D^\#_K$ is sound. It
sits between the two previous domains: $K=1$ is
$\mathcal D^\#_{\mathrm{one}}$, and $K=4^n$ is the full support domain.
Propagating one axis through $G$ gates takes $O(KG)$ axis operations.
Propagating every candidate separately is still quadratic in the number of
rotations; sharing work between candidates needs further design.

### 8.3 Coefficients

**Definition 8.4.** A *coefficient* abstract value is a finite map
$A^\#:\mathbb P_n\rightharpoonup\mathcal C^\#$ from axes to elements of a
coefficient domain $\mathcal C^\#$ with concretization
$\gamma_{\mathcal C}$; absent axes have coefficient exactly zero. It
concretizes to $\{\sum_Pa_PP: a_P\in\gamma_{\mathcal C}(A^\#[P])\}$.

The transformers are the concrete ones of Section 3.2, computed in
$\mathcal C^\#$: Cliffords move entries and flip signs, and a rotation adds
its contributions entry by entry, dropping an axis once its coefficient is
proved zero. Soundness follows by the induction of Theorem 5.11, one
coefficient at a time.

For Clifford+T, Proposition 3.12 allows the exact choice
$\mathcal C^\#=\mathbb Z[\tfrac1{\sqrt2}]$, in which zero tests are
decidable. With it, Example 6.1 ends at exactly $\{Z\mapsto1\}$, which gives
both the singleton and the sign. The cost is size: the number of axes can
double at every T.

Other choices trade precision for cost. Outward-rounded intervals are sound
but may fail to prove a coefficient zero. Plain floating point compared
against an epsilon is not a sound zero test.

### 8.4 Summary

| domain | value | a T on an anticommuting axis | cost per axis per gate | proves Example 6.1 |
|---|---|---|---|---|
| $\mathcal D^\#_{\mathrm{one}}$ | one axis or $\top$ | $\top$ | $O(1)$ axis operations | no |
| $\mathcal D^\#_K$ | at most $K$ axes | adds an axis | $O(K)$ | no |
| $\mathcal D^\#$ | any set of axes | adds an axis | up to $4^n$ | no |
| coefficients | axes with coefficients | mixes coefficients | up to $4^n$ | yes |

## 9. Relation to phase folding

Phase folding follows the same pattern: choose an abstract value, give each
gate a transformer, collect candidates keyed by abstract equality, and rewrite
only when the abstract information proves the concrete equality. For circuits
of CNOT, X, and diagonal rotations, the propagated value (a parity of input
bits) never branches, and all the recorded phases commute, so random
fingerprints and a hash map give an expected linear-time analysis.

Across Hadamards, the propagated value is a Pauli axis, and a non-Clifford
rotation can turn one axis into a combination of two. The domains of Section
8 differ in how much of that branching they keep: the singleton domain gives
up at the first branch and stays fast, while coefficients keep everything and
recognize cancellations at exponential worst-case cost.

## 10. Pauli folding subsumes phase folding

This section compares the Pauli abstraction with the *symbolic phase-folding*
abstraction of linear-time T-count optimization (the exact counterpart of
`PhaseFoldRand`, whose random fingerprints only implement its equality test).
The result: on Clifford+T circuits, every fold phase folding finds, Pauli
folding also finds, with the same merged angle, and Pauli folding finds
strictly more.

### 10.1 The phase-folding abstraction

**Definition 10.1 (path variables and wire values).** Start with one Boolean
variable $x_q$ per input qubit. Each H gate introduces one fresh variable.
Every wire carries an *affine* function of the variables seen so far,
$w=\ell\oplus c$ with $\ell$ a linear form (a set of variables) and
$c\in\{0,1\}$. Initially $w_q=x_q$. The gates act as follows:

| gate | effect on wire values |
|---|---|
| $X_q$ | $w_q\mapsto w_q\oplus1$ |
| $\mathrm{CX}_{c\to t}$ | $w_t\mapsto w_t\oplus w_c$ |
| $H_q$ | $w_q\mapsto y$, a fresh variable |
| $Z,\ S,\ S^\dagger,\ \mathrm{CZ},\ T,\ T^\dagger,\ R_Z$ | none (diagonal gates) |

A rotation $R_Z(\theta)$ on qubit $q$ at *site* $k$ is labelled with the
current wire value, $f_k=w_q$. In the path-sum semantics it contributes the
phase $e^{-i\theta(-1)^{f_k}/2}$ to each path, so two rotations whose labels
have the same linear part contribute phases on the same parity.

**Definition 10.2 (phase fold relation).** Sites $k<l$ are *phase-foldable*,
$k\approx_\varphi l$, when $f_k$ and $f_l$ have the same linear part. Their
constants then agree ($s=+1$) or differ ($s=-1$), and the rotations merge into
one of angle $s\theta_k+\theta_l$ at site $l$.

**Example 10.3.** In

```text
cx q[0],q[1]; t q[1]; h q[0]; t q[1]; cx q[0],q[1]; t q[1];
```

the first CX sets $w_1=x_0\oplus x_1$, so the first T is labelled
$x_0\oplus x_1$. The H sets $w_0=y$ but leaves $w_1$ alone, so the second T
has the same label, and the two are phase-foldable *across the H*. The second
CX sets $w_1=x_0\oplus x_1\oplus y$, so the third T is labelled differently
and folds with neither.

### 10.2 The Pauli fold relation

For comparison, take the Pauli analysis as a relation on the same sites. Let
$C_{<k}$ be the product of the *Clifford* gates before site $k$ (the
rotations are omitted), and let $P_k=\pm C_{<k}^\dagger Z_qC_{<k}$ be the
site's axis pulled back to the circuit input, as `PhaseFoldPauli` computes it.

**Definition 10.4 (Pauli fold relation).** Sites $k<l$ are *Pauli-foldable*,
$k\approx_P l$, when $P_k$ and $P_l$ are equal up to sign, and $P_k$ commutes
with $P_m$ for every rotation site $m$ strictly between them.

By Section 8.1 this is exactly when the singleton analysis started from site
$k$ reaches site $l$ as $\mathrm{One}(\pm Z_r)$; Theorem 7.1 then justifies
the merge. (The pass applies the relation greedily, but a relation between
sites is the right level at which to compare abstractions.)

### 10.3 Containment

**Theorem 10.5.** On a circuit over $\{X,Z,S,S^\dagger,H,\mathrm{CX},
\mathrm{CZ},T,T^\dagger,R_Z\}$, if $k\approx_\varphi l$ then
$k\approx_P l$, with the same relative sign $s$, so the two analyses merge
the pair into the same rotation.

Three lemmas prepare the proof. Write $w^{(t)}$ for the wire values just
before gate $t$ of the circuit.

**Lemma 10.6 (independence).** At every point, the linear parts of the $n$
wire values are linearly independent.

*Proof.* Initially they are $x_0,\ldots,x_{n-1}$. X and diagonal gates leave
linear parts unchanged. CX replaces $w_t$ by $w_t\oplus w_c$, an invertible
change. H replaces $w_q$ by a variable that occurs nowhere else. $\square$

Hence every linear form in the span of the wires has a *unique* expression
$\bigoplus_{q\in Q}w_q$ (plus a constant); call $Q$ its *wire set*.

**Lemma 10.7 (old forms stay or leave for good).** Let $g$ be a linear form
over variables that exist at time $t$. If $g$ lies in the span of the wires
at time $t$ and at a later time $t'$, it lies in the span at every time
between.

*Proof.* Consider the span restricted to forms over the variables existing
at time $t$. X, CX, and diagonal gates do not change the span. H on qubit $p$
replaces $w_p$ by a fresh variable, and a fresh variable contributes nothing
to a form over old variables, so the restricted span can only shrink. Once
$g$ leaves it, it cannot return. $\square$

**Lemma 10.8 (Z-string transport).** Let $g$ be in the span at every time in
$[t,t']$, with wire set $Q_s$ at time $s$. Conjugating $Z_{Q_t}=\prod_{q\in
Q_t}Z_q$ forward through the gates from $t$ to $t'$, Clifford or not, gives
$\pm Z_{Q_{t'}}$, and at each time $s$ the transported operator is $\pm
Z_{Q_s}$.

*Proof.* By induction over the gates; in each case the new wire set is the
one the gate's effect on wire values forces.

- *Diagonal gates* commute with every Z-string and leave wire values
  unchanged: $Q_{s+1}=Q_s$.
- *X on $p$* maps $Z_{Q_s}$ to $\pm Z_{Q_s}$ and changes only a constant:
  $Q_{s+1}=Q_s$.
- *CX from $c$ to $t$* maps $Z_t\mapsto Z_cZ_t$ and fixes $Z_c$; since
  $w_t=w'_t\oplus w'_c$ afterwards, the wire set gains or loses $c$ exactly
  when $t\in Q_s$. The two agree.
- *H on $p$*: $g$ lies in the span after the H, which is the span of the other
  wires plus a fresh variable that $g$ does not contain, so by Lemma 10.6
  $p\notin Q_s$. The H does not touch $Z_{Q_s}$, and the other wires are
  unchanged: $Q_{s+1}=Q_s$. $\square$

*Proof of Theorem 10.5.* Let $k\approx_\varphi l$, with rotations on qubits
$q$ and $r$, and let $g$ be the common linear part of $f_k=w_q$ and $f_l=w_r$.
Its wire set is $\{q\}$ at site $k$ and $\{r\}$ at site $l$. By Lemma 10.7,
$g$ is in the span at every time between, so Lemma 10.8 applies to the
segment from $k$ to $l$.

- *Equal axes.* The transported operator is a Z-string at every rotation site
  $m$ in between, and every rotation there has a Z axis, so they commute and
  the rotations can be skipped: transporting through the Clifford gates alone
  gives the same result, $C\,Z_q\,C^\dagger=\pm Z_r$ with $C$ the Clifford
  gates between $k$ and $l$. Since $C_{<l}=C\,C_{<k}$, pulling back gives
  $P_k=\pm P_l$.
- *Commutation.* At each intervening rotation site $m$ on qubit $u$, the
  transported operator $\pm Z_{Q_m}$ and the local axis $Z_u$ are both
  Z-strings, so they commute. Conjugation preserves commutation, so the
  pulled-back axes $P_k$ and $P_m$ commute.
- *Signs.* Label computational basis states by path assignments $v$. At
  site $k$, $Z_q$ multiplies the basis state of path $v$ by
  $(-1)^{f_k(v)}=\varepsilon(-1)^{g(v)}$ with $\varepsilon=(-1)^{c_k}$, where
  $c_k$ is the constant of $f_k$. Every gate between maps each basis state to
  a basis state up to a phase, except H, which acts only off the operator's
  support; so the transported operator keeps multiplying path $v$ by the
  same factor $\varepsilon(-1)^{g(v)}$. At site $l$ it equals $\sigma Z_r$,
  which multiplies path $v$ by $\sigma(-1)^{f_l(v)}=\sigma(-1)^{c_l}
  (-1)^{g(v)}$. Hence $\sigma=(-1)^{c_k\oplus c_l}=s$, the relative sign of
  Definition 10.2. $\square$

### 10.4 Strictness

**Proposition 10.9.** The containment is strict.

*Proof.* In `t q[0]; h q[0]; t q[1]; h q[0]; t q[0]`, the first and last T
are labelled $x_0$ and $y_2$, the variable of the second H, so they are not
phase-foldable. Their axes are $Z_0$ and $H\,H\,Z_0\,H\,H=Z_0$, and the only
rotation between them, on qubit 1, has axis $Z_1$, which commutes; so they
are Pauli-foldable (Example 7.2). $\square$

Phase folding loses here because a fresh variable forgets that $H\,H$ is the
identity; the Pauli frame tracks the Clifford exactly. Running the two passes
on this circuit leaves 3 T gates after `PhaseFoldRand` and 1 after
`PhaseFoldPauli`.

As a check, all 1,560 pairs of T sites in 400 random Clifford+T circuits on
two and three qubits were classified by both relations: all 417
phase-foldable pairs are Pauli-foldable with the same sign, and 11 further
pairs are Pauli-foldable only.

### 10.5 Native Toffolis: the domains are incomparable

Theorem 10.5 needs every gate to be Clifford+T. Symbolic phase folding also
handles a native CCX, with the *nonlinear* update $w_t\mapsto w_t\oplus
w_{c_1}w_{c_2}$, while the Pauli analysis treats a CCX as seven fixed rotations
(Section 8.1). Neither then contains the other. In

```text
t q[2]; ccx q[0],q[1],q[2]; cx q[0],q[3]; ccx q[0],q[1],q[2]; t q[2];
```

both T gates are labelled $x_2$, because the second CCX undoes the first.
Phase folding merges them. The Pauli analysis is blocked: the first CCX's
rotations include the axis $X_2$, which anticommutes with $Z_2$. Here
`PhaseFoldRand` leaves no T gate and `PhaseFoldPauli` leaves 2. Decomposing
the Toffolis into Clifford+T first restores Theorem 10.5, since then both
analyses see only Clifford+T gates.

## 11. Measurements and classical control flow

So far a segment is a unitary circuit. This section extends the analysis to
programs with measurements, resets, classically controlled gates, and loops.
The rewrite stays the same, merging a rotation into a later one. What
changes is the argument that the earlier rotation can be moved.

### 11.1 Measurements are blockers, not barriers

A measurement of Z on qubit $q$ at a point in the circuit measures some
Pauli $M$ in the circuit-input frame, exactly as a rotation's axis is
computed. Its outcome projectors are $\Pi_\pm=\tfrac12(I\pm M)$.

**Proposition 11.1.** If $P$ commutes with $M$, then $R_P(\theta)$ commutes
with both projectors, so moving the rotation across the measurement changes
neither the outcome probabilities nor the post-measurement state. If $P$
anticommutes with $M$, moving it is invalid in general.

*Proof.* $P$ commutes with $M$, hence with $I\pm M$, hence $R_P(\theta)$
does. For the converse, `t q[0]; h q[0]; measure q[0]; h q[0]; t q[0]`
(measuring $X_0$) is not equal to the version with the T gates merged, for
either outcome. $\square$

In the abstract analysis a measurement is therefore one more fixed blocker,
exactly like the rotations of a Toffoli: $\mathrm{One}(P)$ survives it if $P$
commutes with $M$ and becomes $\top$ otherwise. The Clifford frame passes
through a measurement unchanged. `PhaseFoldPauli` currently treats
measurements as full barriers, which is sound but loses folds such as

```text
t q[0]; cx q[0],q[1]; measure q[1] -> c[0]; cx q[0],q[1]; t q[0];
```

where the measured Pauli is $Z_0Z_1$, which commutes with the axis $Z_0$ of
both T gates. The merged circuit equals the original for both outcomes.

**Resets.** A reset is a measurement followed by an X conditioned on its
outcome. The conditional X is a classically controlled gate (Section 11.2),
so a reset contributes two blockers: the measured Pauli, and the X's axis
$C_{<k}^\dagger X_qC_{<k}$. A rotation may move across the reset only if it
commutes with both.

### 11.2 Classical control: joins of branches

A classically controlled gate `if (c) U` runs $U$ on some executions and not
on others. A rotation moved across it must end up about the same signed axis
on both paths, or the merged angle would depend on the run.

The abstract state needs the sign, since it can now differ between paths.

**Definition 11.2 (signed singleton domain).** A candidate's state is
$\bot$, $\mathrm{One}(Q,s)$ with $s\in\{+,-,?\}$, or $\top$. The sign orders
as $+,-\sqsubseteq\ ?$. A path through a Clifford $C$ maps
$\mathrm{One}(Q,s)$ to $\mathrm{One}(\pi_C(Q),\ s\cdot\sigma_C(Q))$, taking
$?\cdot\pm=\ ?$. The join at a control-flow merge is

$$
\mathrm{One}(Q,s)\sqcup\mathrm{One}(Q,s')=\mathrm{One}(Q,s\sqcup s'),\qquad
\mathrm{One}(Q,s)\sqcup\mathrm{One}(Q',s')=\top\ \ (Q\neq Q'),
$$

with $\bot$ as unit and $\top$ absorbing. A fold at a site with local axis
$Z_r$ needs the state $\mathrm{One}(Z_r,s)$ with $s\in\{+,-\}$.

For `if (c) U`, the state after is the join of the state after $U$ and the
state before. Soundness is the usual argument for collecting semantics: each
execution follows one path, and the join overapproximates both.

**Example 11.3.** In `t q[0]; if (c==1) x q[1]; t q[0]`, the X on qubit 1
fixes $Z_0$, so both paths give $\mathrm{One}(Z_0,+)$ and the T gates merge
into an S. In `t q[0]; if (c==1) x q[0]; t q[0]`, the taken path gives
$\mathrm{One}(Z_0,-)$, since the T gates cancel to leave the X, and the other
gives $\mathrm{One}(Z_0,+)$, since they merge into S. The join is
$\mathrm{One}(Z_0,?)$, and no fold is made.

**Control equivalence.** Moving a rotation also changes *when* it runs. The
earlier rotation at site $k$ is deleted and its angle added at site $l$, so
the rewrite is sound only if $l$ runs exactly when $k$ does, once for each
run of $k$. In control-flow-graph terms, $k$ must dominate $l$, $l$ must
post-dominate $k$, and both must lie in the same loop body. A rotation
before an `if` can therefore move *across* the `if`, as in Example 11.3,
but not *into* one of its branches.

### 11.3 Loops: a fixed point

For a loop `while (c) B`, the state at the loop head is the join of the state
on entry and the state after the body. Iterating the body's transformer from
the entry state reaches a fixed point: in the signed singleton domain every ascending chain
has at most four elements ($\bot\sqsubset\mathrm{One}(Q,\pm)\sqsubset
\mathrm{One}(Q,?)\sqsubset\top$), so each candidate's state rises at most
three times, and no widening is needed. The state after the loop is the
fixed point at the head, since the loop may run any number of times,
including zero.

**Example 11.4.** In `t q[0]; for i in 1..k { cx q[0],q[1]; } t q[0]`, the
body maps $Z_0$ to $Z_0$ (a CX fixes its control's Z), so the fixed point is
$\mathrm{One}(Z_0,+)$ and the T gates merge for every $k$. With `x q[0]` as
the body, the sign alternates with the iteration count: the head state goes
from $\mathrm{One}(Z_0,+)$ to $\mathrm{One}(Z_0,?)$, and no fold is made,
unless the trip count's parity is known.

Rotations inside the body fold with each other as before when they are
control-equivalent within one iteration. A rotation at the end of one
iteration cannot fold with one at the start of the next, since neither runs
exactly once per run of the other. Peeling or rotating the loop exposes such
pairs to the same analysis.

### 11.4 The pulled-back frame under control flow

`PhaseFoldPauli` avoids per-candidate propagation by pulling every axis back
to the circuit input through one shared Clifford frame (Section 8.1). With
control flow the frame itself can differ between paths. The shared-frame
method still works on any region where the frame is the same on all paths,
for instance between merge points whose branches apply equal Cliffords.
Across a branch or loop that changes the frame differently on different
paths, the pass would switch to the forward, per-candidate propagation of
Sections 11.2 and 11.3, for the candidates live at that point.
