# Closed-form symmetric quartet probabilities

Use equal sister branch lengths `T_L=T_R=T`, zero terminal lengths, and
shared structured-root duration `T_A`. Define

\[
\Delta=\sqrt{\lambda^2+16m^2},\qquad c_0=(\lambda+4m)/2,
\]

\[
G=e^{-c_0T}\left[\cosh(\Delta T/2)-\frac{\lambda}{\Delta}\sinh(\Delta T/2)\right],
\]

\[
H=e^{-c_0T}\frac{4m}{\Delta}\sinh(\Delta T/2),
\]

\[
K_s=e^{-c_0T}\left[\cosh(\Delta T/2)+\frac{\lambda}{\Delta}\sinh(\Delta T/2)\right].
\]

Also let

\[
\delta=e^{-(\lambda+2m)T},\qquad \varepsilon=e^{-2mT},
\]

\[
a=(G+\delta)/2,\quad b=(G-\delta)/2,\quad c=H/2,
\]

\[
d=(K_s+\varepsilon)/2,\qquad e=(K_s-\varepsilon)/2,
\]

\[
s_S=G+H,\qquad s_D=H+K_s,
\]

and

\[
\alpha=A_m(T_A)=
\frac{\lambda}{\lambda+2m}
\left[1-e^{-2(\lambda+2m)T_A}\right].
\]

Here `s_S` is the survival probability of a same-background sister pair and
`s_D` is the survival probability of a different-background sister pair.

The exact four-class conditional laws are:

## 4:0

\[
q_1^{40}=1-\frac{2s_S^2}{3}+\frac{4\alpha}{3}(ab-c^2),
\]

\[
q_2^{40}=q_3^{40}=\frac{s_S^2}{3}
-\frac{2\alpha}{3}(ab-c^2).
\]

## 3:1

\[
q_1^{31}=1-\frac{2s_Ss_D}{3}
+\frac{2\alpha c}{3}(G-K_s),
\]

\[
q_2^{31}=q_3^{31}=\frac{s_Ss_D}{3}
+\frac{\alpha c}{3}(K_s-G).
\]

## Concordant 2:2, W1

For `0011/1100`, whose arrangement split is `Q1`,

\[
q_1^{22,S}=1-\frac{2s_S^2}{3}
+\frac{2\alpha}{3}(a^2+b^2-2c^2),
\]

\[
q_2^{22,S}=q_3^{22,S}=\frac{s_S^2}{3}
-\frac{\alpha}{3}(a^2+b^2-2c^2).
\]

## Discordant 2:2, W2

For `0101/1010`, whose arrangement split is `Q2`,

\[
q_1^{22,D}=1-\frac{2s_D^2}{3}
+\frac{\alpha}{3}(4c^2-K_s^2),
\]

\[
q_2^{22,D}=\frac{s_D^2}{3}
+\frac{\alpha}{6}(K_s^2+3\varepsilon^2-4c^2),
\]

\[
q_3^{22,D}=\frac{s_D^2}{3}
+\frac{\alpha}{6}(K_s^2-3\varepsilon^2-4c^2).
\]

For W3 (`0110/1001`), exchange `q_2` and `q_3`.

The discordant asymmetry is

\[
q_2^{22,D}-q_3^{22,D}=\alpha e^{-4mT}.
\]

## Topology reversal

For the Q2-associated discordant class,

\[
q_2>q_1
\quad\Longleftrightarrow\quad
\alpha(d^2+e^2-2c^2)>1-s_D^2.
\]

The left side is the structured-coalescent advantage carried by the final
shared population. The right side is the mass already resolved toward the
species topology on the sister branches. We call the parameter set satisfying
this inequality a **conditional topology-reversal region**.

## Unconditional formula

Let `w40` and `w31` be the total 4:0 and 3:1 weights, and let `W1,W2,W3`
be the weights of W1, W2, and W3. Then

\[
q_{\rm MSRC}=w_{40}q_{40}+w_{31}q_{31}+W_1q_{22,S}
+W_2q_{22,D_2}+W_3q_{22,D_3},
\]

with

\[
w_{40}+w_{31}+W_1+W_2+W_3=1.
\]

Therefore

\[
q_2-q_3=(W_2-W_3)\alpha e^{-4mT}.
\]

Directional discordant asymmetry requires both unequal forward sorting weights
and persistent background structure.

## MSC limit

As `m\to\infty`,

\[
\alpha\to0,\qquad s_S,s_D\to e^{-\lambda T/2}.
\]

For equal sister branches,

\[
\tau=\frac{\lambda}{2}(T_L+T_R)=\lambda T,
\]

and every class converges to

\[
\left(1-\frac23e^{-\tau},\frac13e^{-\tau},\frac13e^{-\tau}\right).
\]

