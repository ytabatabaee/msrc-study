# Analytical conditions for anomalous MSRC quartets

This note separates results proved for a fixed-frequency special case from results for the full finite-$M$ Moran-integrated model. Conditioning on terminal arrangement states generally changes the distribution of ancestral arrangement-count histories.

## Definitions and exact criterion

For `(((1,2),3),4)`, let $Q_1=12|34$, $Q_2=13|24$, and $Q_3=14|23$. Define

$$q_i(z;\theta)=P_\theta(Q_i\mid Z=z),\quad
\Delta(z;\theta)=q_2(z;\theta)-q_3(z;\theta),$$

and $\Gamma(z;\theta)=\max\{q_2,q_3\}-q_1$. A conditional quartet is anomalous exactly when $\Gamma>0$.

The conditional recursion returns the unnormalized joint mass

$$J[z,i]=P(Z=z,Q_i)=c_i(z)+
(\pi_{\rm old}\otimes h_i)^T e^{A_4T_A}F_z,$$

where $F_z$ retains the probability of the pattern event. Thus, for $P(Z=z)>0$,

$$q_2>q_1\iff J[z,2]>J[z,1],\qquad
q_3>q_1\iff J[z,3]>J[z,1].$$

In particular,

$$J[z,2]-J[z,1]=c_2(z)-c_1(z)+
\{\pi_{\rm old}\otimes(h_2-h_1)\}^T e^{A_4T_A}F_z,$$

with the analogous expression for $J[z,3]-J[z,1]$. This is a necessary and sufficient criterion for the full finite-state model.

## Correct scope of the crossed-pattern bound

The statement $T_I>M(1-p)\log 2$ is **not** valid for the full stochastic Moran-integrated model from terminal `Z=0101` alone. Terminal conditioning can select ancestral count histories, so the relevant hazard need not be the constant $1/[M(1-p)]$.

The following narrower theorem is valid. Assume: (i) frequency $p\in(0,1)$ is fixed on every branch from sampling through the intermediate branch; (ii) $\rho=0$; (iii) labels are inherited without switching; (iv) same-background rates are $c_0=1/[M(1-p)]$ and $c_1=1/(Mp)$; (v) conditional on `Z=0101`, taxa 1 and 3 enter the intermediate branch in background 0 and taxon 2 in background 1; and (vi) taxa 1 and 2 cannot coalesce on the sister branch.

Let $E_{13}$ be coalescence of taxa 1 and 3 on the intermediate branch. Then

$$P(E_{13}\mid Z=0101)=1-e^{-T_I/[M(1-p)]}.$$

On $E_{13}$ the topology is $Q_2$, while $Q_1$ is impossible. Hence

$$q_2(0101)\ge 1-e^{-T_I/[M(1-p)]},\qquad
q_1(0101)\le e^{-T_I/[M(1-p)]},$$

and

$$\boxed{T_I>M(1-p)\log 2\Longrightarrow q_2(0101)>q_1(0101).}$$

For `Z=0110`, the matching pair is 2 and 3 in background 1, so $T_I>Mp\log 2$ suffices for $Q_3>Q_1$. Complementary patterns use the corresponding matching background. These bounds are sufficient, not necessary, because root-stage events can add discordant mass.

For the full Moran-integrated model, the valid general replacement is $P(E_{13}\mid Z=0101)>1/2$, with that probability calculated from the pattern-specific operator.

## Same-background cases

Under the fixed-frequency, zero-switching assumptions, `0000` and `1111` reduce to a one-background structured coalescent. Their discordant resolutions are exchangeable and the species-tree resolution is at least as probable. Mixed patterns require the exact operator or separate arguments.

## Symmetry

Under neutral dynamics, complementing every arrangement label and replacing $p$ by $1-p$ maps

$$q_i(z;p)=q_i(1-z;1-p)$$

without permuting topology labels. Therefore $\Delta$ is preserved, not sign-reversed. Under the simulator's genic-selection convention $p'=p(1+s)/(1+sp)$, the corresponding label-swap transformation is $s'=-s/(1+s)$ for $s>-1$.

A sign reversal of $\Delta$ requires a separate taxon permutation that exchanges $Q_2$ and $Q_3$ while preserving the relevant tree stages and parameters. The transformation `0101` to `0110` is not a bit complement. Complement symmetry at $p=1/2$ can restore unconditional $q_2=q_3$ when patterns are summed in complementary pairs, but it does not by itself prove unconditional species-tree dominance.

## Short branches and strong switching

For fixed finite state spaces, $e^{AT}=I+AT+R_2(T)$ with $\|R_2(T)\|\le e^{\|A\|T}\|A\|^2T^2/2$. This supplies a local boundary expansion but no global threshold monotonicity.

In the fixed-frequency model, rapid switching has stationary weights $(1-p,p)$ and averaged pair hazard

$$ (1-p)^2c_0+p^2c_1=1/M.$$

This is a leading fast-slow averaging result. Terminal conditioning creates boundary layers, and a complete conditional ordinary-MSC convergence proof for the coupled Moran model is not established here.

## Numerical status and scope

For the base condition $M=10,p=1/2,\rho=.17,T_S=.3,T_I=.2,T_A=.5$, the exact operator gives

$$P(Q\mid Z=0101)=(0.275862,0.449849,0.274289),$$

with $P(Z=0101)=0.044312$, $\Delta=0.175560$, and $\Gamma=0.173987$. Numerical roots in $\rho$ are reported as crossings; uniqueness and global monotonicity are not proved. The tested symmetric settings showed no unconditional anomaly after marginalization, which remains an observation rather than an impossibility theorem.

The exact operator criterion, fixed-frequency crossed-pattern bound, and stated complement mapping are proved under their assumptions. Global threshold uniqueness, unconditional anomaly exclusion, and a full finite-$M$ exponential bound remain open. No result here says that recombination suppression necessarily causes an anomalous quartet or identifies introgression.
