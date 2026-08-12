# Large-scale circuits with exact expectation-value references

## Purpose

The objective is to benchmark deterministic Pauli truncation (BFS) and
randomized PPS/HT propagation on circuits for which an exact statevector or
exact Pauli vector is too large, while retaining a trustworthy reference value
for one or more observables.

An exact reference need not imply that the complete state is easy to simulate.
It is sufficient that the selected expectation value follows from one of:

1. a planted identity or known net unitary;
2. an exact symmetry;
3. a closed-form expression;
4. a polynomial-size reduced representation, such as a covariance matrix or
   symmetry sector;
5. factorization into independently solvable blocks.

The present simulator assumes the initial state \(|0\rangle^{\otimes n}\) and
propagates a Pauli observable backwards in the Heisenberg picture. For this
initial state,

\[
\langle 0^n|P|0^n\rangle =
\begin{cases}
1, & P\text{ contains only }I\text{ and }Z,\\
0, & P\text{ contains at least one }X\text{ or }Y.
\end{cases}
\]

State-preparation gates can be prepended when another product or entangled
initial state is required. Every proposed family should first be validated at
small \(n\) against `statevector_expectation` before being used beyond the
statevector regime.

---

## 1. Noiseless circuit families

### 1.1 Unitary echo: \(U\) followed by \(U^\dagger\)

#### Construction

Generate an arbitrary deep scrambling circuit

\[
U=G_mG_{m-1}\cdots G_1
\]

and append its inverse in reverse gate order. In forward execution the gate
list is

\[
G_1,G_2,\ldots,G_m,G_m^\dagger,\ldots,G_2^\dagger,G_1^\dagger,
\]

so that the total unitary is

\[
V=U^\dagger U=I.
\]

The exact output state equals the initial state, but at the midpoint the
Heisenberg-evolved observable can have an exponentially large Pauli support.
Exact propagation reconstructs the initial observable in the second half.
Terms discarded or randomly reweighted at the midpoint generally do not
recombine exactly, making this a stringent cancellation and reversibility test.

#### Specific instantiation

Use the existing random Clifford+T brickwork generator for the first half:

- random single-qubit choices from `H`, `S`, and `H` followed by `S`;
- an `RZ(pi/4)` gate independently with density \(p_T\);
- alternating nearest-neighbour CNOT layers.

Append the exact inverse of every generated gate:

| Gate | Inverse to append |
|---|---|
| `H` | `H` |
| `S` | `S S S` (or a new `Sdg` gate) |
| `CNOT` | `CNOT` |
| `RZ(theta)` | `RZ(-theta)` |

The inverse must be generated from the stored gate list, not by drawing a
second random circuit with the same distribution.

#### Repository implementation

Add a helper such as

```cpp
Circuit make_clifford_t_echo(
    int qubits,
    int scrambling_layers,
    double t_density,
    std::uint64_t circuit_seed);
```

It should:

1. generate and store the scrambling-half gates;
2. append them in reverse order using the inverse mapping above;
3. retain all gates explicitly, so no preprocessing cancels the echo;
4. use a circuit-generation seed distinct from every PPS seed;
5. optionally record the exact midpoint event for diagnostics.

Truncation should continue at the normal interval throughout both halves. The
benchmark is uninformative if the inverse gates are simplified away before
Pauli propagation.

#### Exact reference

For any observable \(O\),

\[
\mu_{\mathrm{exact}}
=\langle 0^n|O|0^n\rangle.
\]

Thus a `Z`-only Pauli string has reference \(+1\), whereas the current
`X Z X X` observable has reference zero. It is useful to include both zero and
nonzero observables, because a suite containing only exact zeros cannot test
relative error or multiplicative accuracy.

#### Known-net variant

Append a simple known circuit \(W\) after the echo:

\[
V=WU^\dagger U=W.
\]

Then

\[
\mu_{\mathrm{exact}}
=\langle 0^n|W^\dagger O W|0^n\rangle.
\]

Choose \(W\) to be a product or Clifford circuit, so the reference is obtained
by direct Pauli conjugation or a stabilizer calculation. For example, applying
`H` on every `X` location of the current `X Z X X` observable maps it to a
`Z`-only string and gives reference \(+1\). Product rotations can provide
references strictly between zero and one.

#### Strengths and limitations

- It supports arbitrary width, depth, geometry, and non-Clifford density.
- It directly tests whether truncation preserves delicate cancellations.
- It is intentionally adversarial and time-reversal structured; it should not
  be the only benchmark used to represent generic physical dynamics.

---

### 1.2 Symmetry-protected zero expectation values

#### Construction

Let \(S\) be a unitary symmetry satisfying

\[
S\rho_0S^\dagger=\rho_0,
\qquad [U,S]=0,
\qquad S^\dagger O S=-O.
\]

The expectation is then exactly zero:

\[
\begin{aligned}
\langle O\rangle
&=\operatorname{Tr}(O U\rho_0U^\dagger)\\
&=\operatorname{Tr}(S^\dagger O S\,U\rho_0U^\dagger)\\
&=-\langle O\rangle=0.
\end{aligned}
\]

Unlike an echo, the circuit need not reverse itself. It may describe genuine
long-time, nonintegrable dynamics within one symmetry sector.

#### Specific instantiation: global `Z` parity

Take

\[
S=Z_0Z_1\cdots Z_{n-1},
\qquad \rho_0=|0^n\rangle\langle0^n|.
\]

Construct random layers from parity-preserving rotations such as

\[
e^{-i\alpha_i Z_i/2},\qquad
e^{-i\beta_{ij}Z_iZ_j/2},\qquad
e^{-i\gamma_{ij}X_iX_j/2},\qquad
e^{-i\delta_{ij}Y_iY_j/2}.
\]

An observable containing an odd total number of `X` and `Y` factors is odd
under \(S\), and therefore has exact reference zero. Random `Z`, `ZZ`, `XX`,
and `YY` rotations at generic angles can generate complex operator spreading
while preserving parity.

#### Repository implementation

An arbitrary Pauli rotation

\[
R_P(\theta)=e^{-i\theta P/2}
\]

can be compiled into Clifford basis changes, a CNOT parity ladder, one `RZ`,
and the inverse ladder. Implementing a reusable `append_pauli_rotation` helper
would support this and several later families.

Do not use gates without checking their action on the chosen symmetry. For
example, the current oriented CNOT brickwork does not generally commute with
global `Z` parity.

#### Exact reference

The reference is exactly zero for every circuit realization, width, and depth,
provided the initial state, every gate, and the observable satisfy the symmetry
conditions. No large classical simulation is required.

#### Strengths and limitations

- This gives non-echo dynamics at arbitrary depth.
- The Pauli vector can occupy an exponentially large symmetry sector.
- Only additive error is meaningful for an exact-zero observable. Include
  separate nonzero-reference families in the benchmark suite.
- A mistaken symmetry-breaking gate silently invalidates the reference, so
  small-statevector tests and explicit commutator tests are essential.

---

### 1.3 Commuting diagonal or IQP-type circuits

#### Construction

Prepare \(|+\rangle^{\otimes n}\) and apply commuting `Z` and `ZZ` rotations:

\[
D=
\left(\prod_i e^{-i\phi_i Z_i/2}\right)
\left(\prod_{(i,j)\in E}e^{-i\theta_{ij}Z_iZ_j/2}\right).
\]

Because all generators commute, selected observables have closed-form
expectations even when their Heisenberg Pauli expansions are large.

For one local `X` observable,

\[
\boxed{
\langle X_i\rangle
=\cos\phi_i\prod_{j:(i,j)\in E}\cos\theta_{ij}
}
\]

when each edge has first been combined into one total angle.

The conjugated operator can contain as many as \(2^{\deg(i)}\) Pauli strings.

#### Specific instantiation

Use a high-degree deterministic graph or a seeded random graph. Assign each
edge an angle drawn from a fixed set that avoids Clifford values and avoids
making the reference exponentially indistinguishable from zero. Schedule the
edges into matchings if only disjoint two-qubit gates may be applied
simultaneously.

A bounded-degree one-dimensional graph is usually too easy for a local `X_i`,
because its support is bounded by \(2^{\deg(i)}\). To obtain a difficult
frontier, use one of:

- a graph whose degree grows with \(n\);
- commuting higher-body `Z` rotations;
- a global `X`-string correlator whose exact value is evaluated by a
  low-bond-dimension transfer matrix on a one-dimensional graph.

#### Repository implementation

1. Prepend `H` on every qubit to prepare \(|+^n\rangle\).
2. Implement `RZZ(i,j,theta)` as
   `CNOT(i,j), RZ(j,theta), CNOT(i,j)`.
3. Combine repeated rotations on the same Pauli generator before evaluating
   the analytic reference.
4. Store graph edges and angles in the circuit metadata or result CSV so the
   reference is reproducible.

#### Exact reference

For local `X_i`, evaluate the product above using `long double` or accumulate
logarithms when many cosine factors are present. More general diagonal-circuit
correlators can be mapped to classical Ising partition functions. In one
dimension, a fixed-size transfer matrix gives an exact polynomial-time
reference without constructing a statevector.

#### Strengths and limitations

- The reference is independent of the Pauli propagation implementation.
- The construction directly exposes coefficient branching from Pauli
  rotations.
- Fully commuting repeated layers may collapse into one effective layer, so
  nominal depth alone does not imply increasing physical difficulty.

---

### 1.4 Matchgate or fermionic-Gaussian circuits

#### Construction

Use nearest-neighbour gates generated by Hamiltonians quadratic in fermionic
creation and annihilation operators. In spin language, typical generators are

\[
X_iX_{i+1},\quad Y_iY_{i+1},\quad X_iY_{i+1},\quad
Y_iX_{i+1},\quad Z_i,\quad Z_{i+1}.
\]

At generic angles these circuits are non-Clifford and may have very large spin
Pauli expansions, but their action on Majorana operators is linear.

Define Majoranas

\[
\gamma_{2j}=Z_0\cdots Z_{j-1}X_j,
\qquad
\gamma_{2j+1}=Z_0\cdots Z_{j-1}Y_j.
\]

A Gaussian unitary acts as

\[
U^\dagger\gamma_aU=\sum_b R_{ab}\gamma_b,
\]

where \(R\) is a real \(2n\times2n\) orthogonal matrix.

#### Specific instantiation

Use alternating nearest-neighbour layers of

\[
e^{-i\alpha_\ell(X_iX_{i+1}+Y_iY_{i+1})/2}
\]

with site-dependent `Z` rotations. Generic deterministic or seeded random
angles avoid the Clifford limit. Choose a Gaussian product state such as a
computational-basis occupation state.

Quadratic observables are easy in both representations, so a stronger Pauli
benchmark should include long spin strings or spin correlators whose
Jordan--Wigner representation contains many Majoranas.

#### Repository implementation

- Compile the required two-qubit Pauli rotations using basis changes and the
  existing CNOT--`RZ`--CNOT construction, or add a native matchgate builder.
- Add an independent reference module that builds the Majorana transformation
  matrix for every gate.
- Keep the circuit and reference implementations separate enough that a Pauli
  kernel error cannot be duplicated in the reference solver.

#### Exact reference

Represent the Gaussian state by its covariance matrix

\[
\Gamma_{ab}=\frac{i}{2}\operatorname{Tr}
\left(\rho[\gamma_a,\gamma_b]\right).
\]

After the circuit,

\[
\Gamma\longmapsto R\Gamma R^T.
\]

Quadratic expectations are entries of \(\Gamma\). Expectations of an even
Majorana product follow from Wick's theorem and equal the Pfaffian of the
corresponding covariance submatrix. Thus selected long spin-string
correlations can be computed in polynomial time and to machine precision.

Nearest-neighbour matchgate circuits are a standard polynomial-time
classically simulable family; see
[Valiant's matchgate result](https://epubs.siam.org/doi/10.1137/S0097539700377025).

#### Strengths and limitations

- Both \(n\) and depth can be very large.
- The reference solver is independent and polynomial-time.
- The dynamics are integrable/free, so this should complement rather than
  replace chaotic-circuit benchmarks.
- Care is required with Jordan--Wigner ordering and Pfaffian sign conventions.

---

### 1.5 Dual-unitary brickwork circuits

#### Construction

A two-site gate is dual-unitary when it is unitary both in the usual time
direction and after reshaping indices so that space and time exchange roles.
Brickwork circuits made from such gates can show maximal operator spreading
while selected dynamical correlations reduce to low-dimensional transfer
matrices.

A useful qubit parameterization is

\[
G=(u_1\otimes u_2)
\exp\!\left[-i\left(
\frac{\pi}{4}X\!X+
\frac{\pi}{4}Y\!Y+
JZ\!Z\right)\right]
(v_1\otimes v_2),
\]

with suitable one-qubit unitaries \(u_i,v_i\). Alternate this gate on even and
odd bonds. Generic local dressings and generic \(J\) give interacting,
non-Clifford circuits.

#### Repository implementation

The present gate representation contains only one- and two-qubit elementary
gates, so either:

1. compile `G` into CNOTs and one-qubit Euler rotations; or
2. add a general two-qubit gate type and its Pauli-transfer action.

The first route reuses the current kernel but requires checking that numerical
gate synthesis preserves dual unitarity to the intended tolerance.

The simplest exact references are normalized infinite-temperature
correlations,

\[
C_{AB}(x,t)=2^{-n}\operatorname{Tr}
\left[A_x(t)B_0\right].
\]

The current simulator starts from a pure product state rather than the
maximally mixed state. One implementation is therefore to purify the trace:
prepare a Bell pair between every system qubit and an ancilla, evolve only the
system register, and measure the corresponding system--ancilla Pauli
correlator. This doubles the number of qubits but keeps the reference exact.
Alternatively, implement one of the product or matrix-product initial states
known to be compatible with the selected dual-unitary gate.

#### Exact reference

For appropriate local observables, contract the one-site transfer channel along
the causal-light-cone edge. The cost depends on the local Hilbert-space
dimension and time but not exponentially on the full system width. Boundary
conditions and finite-size wraparound must be included when \(2t\) approaches
\(n\).

Dual-unitary circuits can be interacting and generically nonintegrable while
retaining exact local correlation functions; see
[Bertini, Kos, and Prosen](https://arxiv.org/abs/1904.02140).

#### Strengths and limitations

- This is the strongest non-echo candidate for chaotic operator spreading with
  an exact reference.
- It requires the most additional implementation and validation work.
- Many simple correlations are exactly zero inside the light cone. Select a
  mixture of zero and nonzero transfer-channel eigenmodes.

---

### 1.6 Permutation-symmetric collective circuits

#### Construction

Start from a permutation-symmetric state and apply only collective gates. A
specific kicked-top construction is

\[
U_\ell=
e^{-i\kappa_\ell J_z^2/n}
e^{-i\alpha_\ell J_x},
\qquad
J_a=\frac12\sum_{i=0}^{n-1}\sigma_i^a.
\]

Both the state and evolution remain in the total-spin \(J=n/2\) sector, whose
dimension is only \(n+1\), even though the full Hilbert space has dimension
\(2^n\).

#### Repository implementation

- The collective `Jx` rotation is a product of equal one-qubit `RX` rotations.
- Up to a global phase,

  \[
  e^{-i\kappa J_z^2/n}
  =\prod_{i<j}e^{-i\kappa Z_iZ_j/(2n)},
  \]

  so it can be implemented with all-to-all `RZZ` gates.
- Choose collective observables or fixed-site correlators related to collective
  moments.

This family requires nonlocal connectivity, or a SWAP network if restricted to
nearest neighbours.

#### Exact reference

Build the \((n+1)\times(n+1)\) spin matrices in the Dicke basis
\(|J,m\rangle\). `Jz^2` is diagonal, and collective rotations can be applied
using spin matrices or Wigner rotation matrices. Propagate an \(n+1\)-component
state rather than a \(2^n\)-component state.

Examples of references are

\[
\langle Z_i\rangle=\frac{2}{n}\langle J_z\rangle,
\]

and, for \(i\ne j\),

\[
\langle Z_iZ_j\rangle
=\frac{4\langle J_z^2\rangle-n}{n(n-1)}.
\]

#### Strengths and limitations

- Width and depth can be large with an \(O(n)\)-dimensional reference state.
- The Pauli spectrum contains extensive degeneracies, making this useful for
  testing tie handling and coefficient plateaus.
- It represents collective rather than geometrically local dynamics.

---

### 1.7 Independent hard blocks

#### Construction

Partition the register into \(M=n/b\) blocks of fixed size \(b\). Apply a deep
arbitrary circuit \(U_r\) within every block, with no gates between blocks:

\[
U=\bigotimes_{r=1}^M U_r.
\]

Choose a factorized observable

\[
O=\bigotimes_{r=1}^M O_r
\]

or an additive observable \(O=\sum_r O_r\).

The global Pauli support of a product observable can be the product of the
block supports and therefore exponential in \(M\), even though every block is
small enough for an exact reference calculation.

#### Repository implementation

1. Select a block size such as \(b=8\)--\(16\).
2. Generate one or several seeded deep Clifford+T or Ising-like block circuits.
3. Copy their gates with qubit-index offsets.
4. Construct the global observable from translated block observables.
5. Apply the global support budget, rather than truncating each block
   separately, if the goal is to test the existing simulator unchanged.

#### Exact reference

For a product observable and product input,

\[
\langle O\rangle
=\prod_{r=1}^M
\langle 0^b|U_r^\dagger O_rU_r|0^b\rangle.
\]

For an additive observable, sum the corresponding block expectations. Compute
each block by statevector or exact Pauli propagation. Accumulate logarithms and
signs when a product of many small expectations risks underflow.

#### Strengths and limitations

- It supports arbitrary gates and arbitrarily large total \(n\) and depth.
- It is easy to implement and provides nonzero references.
- There is no entanglement between blocks, so this is primarily a scaling and
  global-budget test rather than a globally chaotic benchmark.

---

### 1.8 Clifford+T circuits with limited T count

#### Construction

Use arbitrarily many Clifford gates but restrict the total number of T gates to
\(t\), independently of the number of qubits and Clifford depth. The circuit
can scramble Pauli locations globally while all non-Clifford complexity remains
controlled by \(t\).

For polynomial scaling, choose \(t=O(\log n)\). For a demanding but still
classically tractable fixed benchmark, choose a larger constant \(t\) supported
by the available reference solver.

#### Repository implementation

- Add a circuit mode that accepts a total T-count budget rather than an
  independent density per site and layer.
- Distribute those T gates at seeded random spacetime locations among deep
  random Clifford layers.
- Record the exact locations and circuit seed.

#### Exact reference

A T or `RZ(pi/4)` gate can be decomposed into a short linear combination of
Clifford operations. Expanding the \(t\) non-Clifford gates gives a sum of
stabilizer branches, evaluated with a stabilizer simulator. A direct expansion
has exponential dependence on \(t\), but no exponential dependence on the
number of Clifford gates or qubits. Stabilizer-rank methods may reduce the
number of required branches.

An exact untruncated Pauli calculation is also possible when \(2^t\) remains
small, but a stabilizer-based reference is preferable because it is
algorithmically independent of the Pauli-streaming implementation under test.

#### Strengths and limitations

- Very large width and Clifford depth are possible.
- The family isolates how error depends on non-Clifford count.
- It cannot reproduce a genuinely extensive non-Clifford density while
  retaining a scalable exact reference.

---

## 2. Noisy circuit families

For noisy evolution the relevant object is a channel, not a unitary. In
general,

\[
\mathcal U^\dagger\circ\mathcal N\circ\mathcal U
\ne\mathcal N,
\]

and the inverse unitary does not undo \(\mathcal N\). Therefore the noiseless
identity reference must not be reused when independent local noise is inserted
inside an arbitrary echo.

### 2.1 Echo blocks with local noise at block boundaries

#### Construction

Use a complete noiseless echo as one block and apply a product channel only
after the block has closed:

\[
\mathcal E
=\mathcal N^{\otimes n}\circ
 \mathcal U^\dagger\circ\mathcal U
=\mathcal N^{\otimes n}.
\]

Repeat this construction \(B\) times. Every block may contain an independently
generated deep scrambling unitary and its exact inverse.

This retains a difficult intermediate Pauli frontier while leaving the
block-to-block physical channel exactly solvable.

#### Repository implementation

Generate the gate list

```text
for block:
    append U_block
    append inverse(U_block)
    append one noise gate on every qubit
```

With the existing depolarizing convention, the adjoint channel multiplies a
nonidentity local Pauli by `1-p`. Do not place these local noise gates between
the gates of `U_block` and its inverse if the closed-form reference below is to
remain valid.

#### Exact reference

Suppose the single-qubit channel is diagonal in the Pauli basis,

\[
\mathcal N^\dagger(P_a)=\lambda_aP_a,
\qquad P_a\in\{X,Y,Z\}.
\]

For a product Pauli \(P=\bigotimes_iP_{a_i}\) and a product input,

\[
\langle P\rangle_B
=\left(\prod_i\lambda_{a_i}^{B}\right)
 \langle P\rangle_0.
\]

For the repository's depolarizing convention, a `Z`-only string of weight
\(w\) on \(|0^n\rangle\) has

\[
\boxed{
\langle Z_S\rangle_B=(1-p)^{Bw}.
}
\]

The formula also extends to arbitrary single-qubit channels by iterating their
small Bloch-vector affine maps. For amplitude damping, use an initial state and
observable that are not already at the channel's fixed point.

#### Strengths and limitations

- It uses the current local depolarizing implementation directly.
- It gives nonzero, tunable noisy references.
- Noise acts only between hard echo blocks, not during the scrambling dynamics.

---

### 2.2 Global depolarizing noise throughout an echo

#### Construction

Define the global \(n\)-qubit depolarizing channel

\[
\mathcal D_\lambda(\rho)
=\lambda\rho+(1-\lambda)\frac{I}{2^n}.
\]

For every traceless observable,

\[
\mathcal D_\lambda^\dagger(O)=\lambda O.
\]

This channel commutes with every unitary channel. It can therefore be inserted
after every layer of an arbitrary echo while preserving an analytic reference.

#### Repository implementation

The existing per-qubit depolarizing gates implement a tensor product of local
channels; that is not the same as global depolarization and does not commute
with a generic entangling unitary. Add a distinct `GlobalDepolarizing` gate or
channel event whose Heisenberg action is:

```text
identity Pauli:      coefficient unchanged
nonidentity Pauli:   coefficient *= lambda
```

#### Exact reference

If \(m\) global depolarizing events occur and the net unitary is the identity,

\[
\boxed{
\langle O\rangle=\lambda^m\langle O\rangle_0
}
\]

for every traceless \(O\). For a known net unitary \(W\), replace the initial
expectation by \(\langle0|W^\dagger O W|0\rangle\).

#### Strengths and limitations

- Noise may be interleaved throughout the difficult circuit.
- Its effect is a uniform rescaling of every nonidentity Pauli, so it does not
  test weight-dependent noise or noise-aware importance sampling.
- It is a collective noise model rather than independent local noise.

---

### 2.3 Symmetry-covariant noisy circuits

#### Construction

Generalize the symmetry argument from a unitary to a channel. Require

\[
\mathcal E(S\rho S^\dagger)
=S\mathcal E(\rho)S^\dagger,
\]

together with

\[
S\rho_0S^\dagger=\rho_0,
\qquad S^\dagger O S=-O.
\]

Then

\[
\operatorname{Tr}[O\mathcal E(\rho_0)]=0
\]

exactly, even with noise after every physical layer.

#### Specific instantiation

Use the global `Z`-parity-preserving random circuit from Section 1.2 and insert
the current single-qubit depolarizing channel after every layer. Local
depolarization is covariant under Pauli conjugation and therefore preserves the
global parity covariance. Local `Z` dephasing also works. An observable with an
odd number of `X/Y` factors retains exact reference zero.

#### Repository implementation

- Add a symmetry-preserving noisy model that reuses the parity-preserving
  unitary builder.
- Apply the chosen noise channel at the desired cadence.
- Add tests that apply the symmetry transformation to every gate/channel and
  verify covariance numerically on one- or two-qubit operator bases.

#### Exact reference

The reference is exactly zero for arbitrary qubit count, circuit depth, noise
strength, and circuit seed, provided every channel event is symmetry covariant.

#### Strengths and limitations

- It supports realistic local noise throughout non-echo dynamics.
- It tests a large noisy Pauli sector.
- It supplies only exact-zero observables unless combined with another
  construction.

---

### 2.4 Diagonal circuits with commuting dephasing

#### Construction

Use the commuting diagonal circuit from Section 1.3 and insert local `Z`
dephasing,

\[
\mathcal Z_\eta^\dagger(I)=I,
\quad
\mathcal Z_\eta^\dagger(Z)=Z,
\quad
\mathcal Z_\eta^\dagger(X)=\eta X,
\quad
\mathcal Z_\eta^\dagger(Y)=\eta Y.
\]

The dephasing channel is a mixture of `I` and `Z` conjugations and commutes with
all diagonal `Z` and `ZZ` unitaries. All dephasing events can therefore be
moved to the end analytically.

#### Repository implementation

Add a `Dephasing` gate kind with the adjoint action above. Prepare
\(|+\rangle^{\otimes n}\), apply the diagonal phase circuit, and insert
dephasing after each scheduled layer.

Independent local depolarization must not be substituted in the analytic
derivation: it does not generally commute through an entangling `ZZ` rotation.

#### Exact reference

If qubit \(i\) experiences \(m_i\) identical dephasing events,

\[
\boxed{
\langle X_i\rangle
=\eta^{m_i}\cos\phi_i
\prod_{j:(i,j)\in E}\cos\theta_{ij}
}
\]

Different dephasing strengths are handled by replacing \(\eta^{m_i}\) with the
product of the individual coherence factors.

#### Strengths and limitations

- It provides a nonzero noisy reference with noise distributed throughout the
  circuit.
- It remains a commuting, structured model.
- Choose angles and dephasing rates so that the reference does not underflow.

---

### 2.5 Noisy fermionic-Gaussian circuits

#### Construction

Combine the matchgate circuit from Section 1.4 with fermionic-Gaussian channels
such as mode loss, gain, or suitable Gaussian thermalization. A Gaussian
channel maps covariance matrices affinely:

\[
\Gamma\longmapsto A\Gamma A^T+B.
\]

The matrices \(A\) and \(B\) must satisfy the complete-positivity constraints
for a fermionic Gaussian channel.

#### Repository implementation

- Add the desired physical channel to the Pauli kernel.
- Independently implement its \((A,B)\) covariance update in the reference
  solver.
- Apply unitary Majorana rotations and noisy affine updates in the same physical
  order as the circuit gates.
- Validate every channel against an exact density-matrix calculation for small
  systems.

Qubit amplitude damping can represent fermionic loss under a consistent
Jordan--Wigner occupation convention, but signs, parity restrictions, and the
chosen mode ordering must be tested explicitly.

#### Exact reference

Iterate

\[
\Gamma_{k+1}=A_kR_k\Gamma_kR_k^TA_k^T+B_k
\]

and evaluate quadratic observables directly or higher even-Majorana
correlations as Pfaffians. The computation remains polynomial in \(n\) and
depth.

#### Strengths and limitations

- It provides large noisy circuits with nonzero, independent references.
- It tests noise interleaved with genuine dynamics.
- Only Gaussian-preserving channels are covered; generic Pauli noise or spin
  dephasing may leave the Gaussian family.

---

### 2.6 Independent noisy blocks

#### Construction

Use the block-factorized circuit from Section 1.7, allowing arbitrary Markovian
noise within each block. Because no gate or channel couples distinct blocks,

\[
\mathcal E=\bigotimes_{r=1}^M\mathcal E_r.
\]

This supports the existing local depolarizing channel at any location in the
block circuit, including after every internal layer.

#### Repository implementation

- Choose a block size small enough for exact density-matrix propagation;
  typically this is smaller than the noiseless statevector block size because a
  dense density matrix has \(4^b\) entries.
- Generate and offset the block gate lists as in the noiseless construction.
- Apply exactly the same noise cadence in the block reference solver and the
  Pauli benchmark.
- Blocks may share parameters or use separately seeded circuits.

#### Exact reference

For product input states and product observables,

\[
\langle O\rangle
=\prod_r\operatorname{Tr}
\left[O_r\mathcal E_r(\rho_{0,r})\right].
\]

Each factor is computed once by exact block density-matrix evolution. Additive
observables give sums instead of products.

#### Strengths and limitations

- It accommodates essentially arbitrary local gates and noise models.
- The reference is simple and independent of global Pauli propagation.
- It does not test noisy entanglement across blocks.

---

## 3. Recommended implementation order

### Phase 1: minimal changes, immediate large-scale benchmarks

1. **Noiseless known-net echo** using the existing Clifford+T gate set.
2. **Echo blocks with boundary local depolarization** using the current noise
   implementation.
3. **Independent noiseless and noisy blocks** using the existing statevector
   and a small exact density-matrix reference.

These require little new gate algebra and immediately provide zero and nonzero
references at widths beyond exact global simulation.

### Phase 2: non-echo analytic families

4. Add a general Pauli-rotation circuit builder.
5. Implement **global-parity symmetry circuits**, with and without local
   depolarization.
6. Implement **commuting diagonal circuits**, then add a `Dephasing` channel.

These distinguish echo-specific behavior from ordinary long-time propagation.

### Phase 3: independent polynomial-time reference solvers

7. Implement a **fermionic covariance/Pfaffian reference solver**.
8. Add **permutation-symmetric kicked-top circuits** and a Dicke-basis solver.
9. Add **dual-unitary circuits** and transfer-channel references if maximally
   spreading nonintegrable dynamics are required.

---

## 4. Validation and reporting requirements

Every new family should record enough information to reconstruct both circuit
and reference:

- model name and version;
- qubit count and layer count;
- circuit-generation seed;
- PPS seed, separately;
- all physical angles and noise strengths;
- observable and state-preparation gates;
- analytic-reference formula or reference-solver name;
- reference precision and numerical tolerance.

For small widths, test all new references against dense statevector or density
matrix evolution. For noisy channels, test both Schrödinger and Heisenberg
implementations on the complete one- or two-qubit Pauli basis.

Use several observables per family:

- an exact zero, to reveal additive leakage;
- an order-one nonzero value, to support relative-error comparisons;
- a moderately attenuated value, to test dynamic range;
- when possible, observables with different Pauli weights and causal cones.

Finally, report the largest exact or untruncated frontier available at small
width, the deterministic retained-support schedule, PPS empirical standard
error, and the distribution of individual PPS estimates. A known reference
alone does not distinguish truncation bias from heavy-tailed sampling variance.
