\documentclass[11pt,a4paper]{article}
\usepackage[T1]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage{lmodern,amsmath,amssymb,amsthm,mathtools,mathrsfs}
\usepackage[margin=25mm]{geometry}
\usepackage{booktabs,array,longtable,microtype}
\usepackage[colorlinks=true,linkcolor=blue,citecolor=blue,urlcolor=blue]{hyperref}
\usepackage{enumitem}
\setlist{itemsep=2pt,topsep=4pt}
\setlength{\parskip}{3pt}
\setlength{\emergencystretch}{2em}
\newcommand{\dd}{\mathrm{d}}
\newcommand{\ii}{\mathrm{i}}
\newcommand{\Tr}{\operatorname{Tr}}
\newcommand{\Res}{\operatorname{Res}}
\newcommand{\Wind}{\operatorname{Wind}}
\newcommand{\Lie}{\mathcal{L}}
\newcommand{\eps}{\varepsilon}
\newtheorem{proposition}{Proposition}[section]
\newtheorem{corollary}[proposition]{Corollary}
\theoremstyle{remark}
\newtheorem{remark}[proposition]{Remark}
\title{Thermodynamic Representation and Lorentz Covariance:\\
Entropy Stationarity, Integrability Defects, and Horizon Residues}
\author{Wen-Xiang Chen\\
\small School of Electronic Information, Guangzhou City University of Technology\\
\small \texttt{wxchen4277@qq.com}}
\date{}
\begin{document}
\maketitle
\begin{abstract}
Thermodynamic Legendre transformations and Lorentz transformations both permit changes of description, but act on different structures and do not establish an equivalence between thermodynamics and special relativity. We construct a restricted consistency framework linking three independently defined objects: the first-law one-form on an equilibrium state manifold, the covariant entropy current of a perfect fluid, and a meromorphic radial optical one-form near a static nondegenerate horizon. Constrained entropy maximization identifies the rest-frame grand potential as minus the temperature times the stationary Massieu functional. A first-law integrability defect supplies a calculable loop test, while the optical residue reproduces the Hawking temperature after the normalization of the Killing time is fixed. Explicit counterexamples delimit these results: low-frequency Planck occupation does not cause an infrared energy divergence; an extremal Reissner--Nordstr\"om horizon has a nonzero optical residue although its standard Hawking temperature vanishes; and distinct metric profiles can share the same horizon area, temperature, and asymptotic mass. The resulting construction is a diagnostic synthesis, rather than a derivation of Lorentz symmetry or gravitational field equations from thermodynamic identities. Its contribution is to replace an unrestricted correspondence by explicit assumptions, reproducible calculations, and falsifiable failure criteria.
\end{abstract}
\noindent\textbf{Keywords:} Legendre transformation; Lorentz covariance; constrained entropy; contact geometry; integrability; RVB; horizon residue.

\section{Research question and method}
The question is whether changes of thermodynamic representation, relativistic frame changes, and horizon-temperature calculations can be compared through explicit consistency tests without identifying their distinct mathematical structures. Thermodynamic potentials and their natural variables provide the first anchor \cite{Callen,Mrugala}. Covariant statistical mechanics provides the second \cite{Becattini2016,Becattini2015}. Euclidean horizon regularity provides the third \cite{GibbonsHawking}. These inputs support a limited synthesis; none implies the others.

The working method starts from an established identity, chooses a mathematical language suited to that identity, isolates a discrepancy, constructs the smallest extension that makes the discrepancy calculable, attempts a counterexample, and restricts the conclusion to the surviving domain. Here the languages are contact geometry for the first law, tensor calculus for relativistic fluids, and complex analysis for a horizon one-form. This implements the principle of remaining within a tractable scope while making a small explicit advance.

The distinction between a physical state and its representation is central. A Legendre transformation changes thermodynamic coordinates and generating functions. A passive Lorentz transformation changes tensor components. A change of ensemble additionally changes constraints and need not preserve fluctuations or stability. A change of gravitational theory changes dynamics. Only the first two can be treated as changes of description under the conditions stated below.

\subsection{Conventions and assumptions}
We use $c=\hbar=k_B=1$, keep Newton's constant $G_{\mathrm N}$ explicit, and adopt metric signature $(-,+,+,+)$. Temperature has dimensions of inverse length, and surface gravity is defined relative to a specified Killing time. Greek indices are spacetime indices. The spacetime stress tensor is $\mathcal T^{\mu\nu}$, so it cannot be confused with temperature $T$.

The equilibrium state relation is $C^2$ on a smooth single-phase patch. Whenever natural variables are exchanged, the relevant partial Jacobian must be nonzero. The fluid is simple, extensive, and non-dissipative unless explicitly stated otherwise. The statistical variational calculation assumes a finite partition function and homogeneous equilibrium with constant inverse-temperature vector; rotating equilibrium requires additional angular-momentum constraints. Horizon calculations assume static spherical symmetry, an outer simple zero of the metric function, finite lapse, and a local holomorphic continuation. No modified-gravity field equation is assumed or reconstructed from the residue.

The results concern smooth branches. Phase coexistence, singular response functions, singular coordinate maps, and degenerate horizons are boundary tests, not implicit extensions of those results.

\section{Thermodynamic representation: the certain anchor}
\subsection{Potentials and the local invertibility condition}
For a simple system,
\begin{equation}\label{eq:firstlaw}
\dd U=T\dd S-p\dd V+\mu\dd N,
\quad T=U_S,\quad p=-U_V,\quad\mu=U_N.
\end{equation}
The potential relations and differentials are
\begin{align}
F&=U-TS,&\dd F&=-S\dd T-p\dd V+\mu\dd N,\label{eq:F}\\
H&=U+pV,&\dd H&=T\dd S+V\dd p+\mu\dd N,\\
G_{\mathrm{th}}&=U-TS+pV,&\dd G_{\mathrm{th}}&=-S\dd T+V\dd p+\mu\dd N,\\
\Omega&=U-TS-\mu N,&\dd\Omega&=-S\dd T-p\dd V-N\dd\mu.\label{eq:Omega}
\end{align}
The subscript on $G_{\mathrm{th}}$ distinguishes Gibbs free energy from $G_{\mathrm N}$. The grand-canonical potential is $\Omega$, not $G_{\mathrm{th}}$.

For example, solving $T=U_S(S,V,N)$ for $S$ requires $U_{SS}\ne0$ locally. A simultaneous transformation of several variables requires the corresponding Hessian block to be nonsingular. Convexity or concavity conditions specify stable branches and appropriate extremum principles. At phase transitions a globally invertible smooth chart may fail; a Legendre--Fenchel construction may still exist but need not retain every metastable branch. Consequently, unrestricted information equivalence is too strong.

On a common smooth branch with the same constraints, reconstruction of $U$ from $F$ uses $S=-F_T$ and $U=F+TS$. This preserves the equation of state. It does not imply that the numerical values of $U$ and $F$ coincide, or that canonical and microcanonical fluctuations coincide.

\subsection{Contact structure and equilibrium pullback}
Introduce the seven-dimensional thermodynamic contact space with
\begin{equation}\label{eq:contact}
\eta=\dd U-T\dd S+p\dd V-\mu\dd N.
\end{equation}
An equilibrium embedding $\iota:\mathcal E\to\mathcal C$ is three-dimensional and satisfies $\iota^*\eta=0$. It is the pullback to the equilibrium manifold that vanishes, not the contact form on the full contact space.

Set $q^i=(S,V,N)$ and $y_i=(T,-p,\mu)$. Then
\begin{equation}
\eta=\dd U-y_i\dd q^i,\qquad
\eta\wedge(\dd\eta)^3\ne0.
\end{equation}
In particular, $\eta$ is not an exact one-form on the contact space. For a partial Legendre transformation on an index set $I$,
\begin{equation}
U'=U-\sum_{i\in I}y_iq^i,\qquad
q'^i=y_i,\qquad y'_i=-q^i\quad(i\in I),
\end{equation}
with the remaining pairs unchanged. Substitution gives
\begin{equation}
\dd U'-y'_i\dd q'^i=\dd U-y_i\dd q^i.
\end{equation}
These standard transformations preserve $\eta$ strictly. A general contactomorphism may instead satisfy $L^*\eta=a\eta$ for a nonvanishing function $a$. It preserves the equilibrium contact condition but not necessarily the numerical value of a nonzero loop integral. We will use strict Legendre transformations when comparing such integrals.

\subsection{Maxwell relations}
For fixed $N$, equality of mixed derivatives of $U$ yields
\begin{equation}\label{eq:Maxwell}
\left(\frac{\partial T}{\partial V}\right)_{S,N}
=-\left(\frac{\partial p}{\partial S}\right)_{V,N}.
\end{equation}
For $F$, the corresponding relation is
\begin{equation}
\left(\frac{\partial S}{\partial V}\right)_{T,N}
=\left(\frac{\partial p}{\partial T}\right)_{V,N}.
\end{equation}
These are integrability identities. They do not establish Lorentz symmetry, causality, or any particular gravitational field equation.

\section{A falsifiable difference point: first-law defects}
Fix $N$ and consider candidate fields $T(S,V)$ and $p(S,V)$ on a two-dimensional domain $D$. Define the energy-valued one-form and its defect by
\begin{equation}\label{eq:defect}
\omega=T\dd S-p\dd V,\qquad
K=\partial_V T+\partial_S p.
\end{equation}
Direct differentiation gives
\begin{equation}\label{eq:domega}
\dd\omega=-K\dd S\wedge\dd V.
\end{equation}
Thus the existence of a potential is testable independently of verbal claims of equivalence.

\begin{proposition}[First-law integrability]\label{prop:integrability}
Let $D\subset\mathbb R^2$ be open and simply connected, and let $T,p\in C^1(D)$. A single-valued potential $U$ satisfying $\dd U=\omega$ exists on $D$ if and only if $K=0$. Equivalently, all closed piecewise smooth curves in $D$ have $\oint\omega=0$. The potential is unique up to an additive constant.
\end{proposition}
\begin{proof}
If $\omega=\dd U$, then $\dd\omega=0$, which implies $K=0$. Conversely, a closed $C^1$ one-form is exact on a simply connected planar domain; integrating $\omega$ from a reference point defines $U$ independently of path. Path independence is equivalent to the vanishing of every closed-loop integral.
\end{proof}

For an arbitrary proposed single-valued $U$, define its pulled-back first-law residual
\begin{equation}
\Delta=\dd U-\omega=\dd U-T\dd S+p\dd V.
\end{equation}
For a positively oriented curve $C=\partial A$ with its interior contained in $D$,
\begin{equation}\label{eq:cycle}
\oint_C\Delta=-\oint_C\omega
=\int_A K\dd S\dd V.
\end{equation}
The integral has energy units. It is a thermodynamic loop defect, not automatically a complex residue or an entropy-production rate.

\begin{remark}
$K=0$ ensures the existence of some potential. It does not show that a previously proposed $U$ is the correct potential: $\Delta$ may be a nonzero exact form even though every loop integral vanishes. Testing the actual candidate requires $\Delta=0$ pointwise as well. On a domain with holes, local closedness additionally requires vanishing periods on noncontractible cycles for global exactness.
\end{remark}

\subsection{A controlled inconsistency}
Choose
\begin{align}
T&=T_0+aS,&p&=p_0-bV+\chi S,\label{eq:toy}\\
U_0&=T_0S+\tfrac12aS^2-p_0V+\tfrac12bV^2,
\end{align}
where $a,b>0$, and restrict the state patch to $T>0$. The parameter $\chi$ has dimensions of pressure per entropy. Then
\begin{equation}
K=\chi,\qquad\Delta=\chi S\dd V,
\qquad\oint_C\Delta=\chi\,\Delta S\,\Delta V
\end{equation}
for a counterclockwise rectangular loop. At $\chi=0$ the first law is recovered. At $\chi\ne0$ no energy state function reproduces both proposed fields on an open patch. This is a specific falsification rather than a metaphorical extraction of invariant content.

\subsection{What invariance actually proves}
If $L$ is a strict Legendre transformation and $\widetilde\iota=L\circ\iota$, then
\begin{equation}
\widetilde\iota^*\eta'=\iota^*(L^*\eta')=\iota^*\eta.
\end{equation}
Corresponding curves therefore give the same residual integral. A nonzero mismatch after correctly matching states, contours, and conventions rules out that specified representation change. Equality of one integral is only a necessary check; it is not a sufficient criterion for physical equivalence. Even equality of all integrability periods does not determine the equation of state, since many exact one-forms have zero periods.

\section{Lorentz covariance and constrained entropy stationarity}
\subsection{Fixing the sign at the rest-frame anchor}
Let $u^\mu u_\mu=-1$ and define
\begin{equation}\label{eq:beta}
\beta^\mu=\frac{u^\mu}{T},\qquad
\alpha=\frac{\mu}{T},\qquad
T=(-\beta^\mu\beta_\mu)^{-1/2}.
\end{equation}
For future-directed four-momentum $P^\mu$, the rest-frame energy is
\begin{equation}
E_0=-u_\mu P^\mu>0.
\end{equation}
With the chosen negative-time signature, the homogeneous equilibrium density operator is
\begin{equation}\label{eq:rho}
\rho_* =Z^{-1}\exp\!\left(\beta_\mu\widehat P^\mu+\alpha\widehat N\right).
\end{equation}
At rest, $\beta_\mu\widehat P^\mu=-\widehat E_0/T$, giving the familiar $\exp[-(\widehat E_0-\mu\widehat N)/T]$. Writing $-\beta_\mu P^\mu$ with the same signature and future-directed $\beta^\mu=u^\mu/T$ would reverse the Boltzmann suppression. Formulas quoted from a positive-time convention must be converted consistently \cite{Becattini2016}.

Equation~\eqref{eq:rho} assumes constant $\beta$ and no independent rotational constraints. It is not the most general rotating or accelerating equilibrium density operator. Lorentz transformations act on $u$, $\beta$, $P$, and the operators or states consistently, preserving their contractions.

\subsection{A precise stationary functional}
For normalized density operators, define
\begin{equation}\label{eq:Massieu}
\mathscr M[\rho]=-\Tr(\rho\ln\rho)
+\beta_\mu\Tr(\rho\widehat P^\mu)
+\alpha\Tr(\rho\widehat N),\qquad\Tr\rho=1.
\end{equation}
The multipliers are fixed. Varying $\mathscr M-\lambda(\Tr\rho-1)$ gives
\begin{equation}
\ln\rho=\beta_\mu\widehat P^\mu+\alpha\widehat N-(1+\lambda)\mathbf1,
\end{equation}
so normalization yields \eqref{eq:rho}. Substituting back,
\begin{equation}\label{eq:stationary}
\mathscr M[\rho_*]=\ln Z
=S_*-\frac{U_0-\mu N_*}{T}.
\end{equation}
The statement is stronger than stationarity. For finite $Z$,
\begin{equation}
\mathscr M[\rho]=\ln Z-D(\rho\Vert\rho_*),
\end{equation}
where $D(\rho\Vert\rho_*)=\Tr\rho(\ln\rho-\ln\rho_*)\ge0$. Thus $\rho_*$ maximizes $\mathscr M$. In its rest frame,
\begin{equation}\label{eq:potentialstationary}
\boxed{\Omega=U_0-TS_*-\mu N_*=-T\mathscr M[\rho_*]=-T\ln Z.}
\end{equation}
This realizes the idea of a thermodynamic potential as minus the temperature times a constrained entropy functional evaluated at its extremum. The functional is the entropy \emph{including its multiplier terms}, not the bare entropy alone. At fixed $N$, removal of the number multiplier gives $F=-T\ln Z_N$. The statement does not identify an arbitrary geometric functional with statistical entropy.

\subsection{Two transformations with independent roles}
For a local scalar rest-frame equation of state, one can first express it in another valid thermodynamic chart and then transform tensor components, or reverse those operations. Since the chart change acts on scalar thermodynamic arguments and the Lorentz transformation acts on spacetime indices, the reconstructed fluid tensor agrees. This is a restricted compatibility statement, not an isomorphism between the Lorentz group and thermodynamic contact transformations.

The first law supplies integrability, whereas Lorentz covariance supplies tensorial transformation rules. Neither supplies the microscopic equation of state. Covariance alone also does not guarantee causal matter. For example, at fixed entropy per particle $\sigma=s/n$, choose
\begin{equation}\label{eq:causality}
\eps(n,\sigma)=A(\sigma)n^\Gamma,
\quad p=n\left(\frac{\partial\eps}{\partial n}\right)_\sigma-\eps
=(\Gamma-1)\eps.
\end{equation}
For $A>0$, $A'>0$, and $\Gamma>2$, this is a differentiable scalar constitutive relation with $T=A'(\sigma)n^{\Gamma-1}>0$. It can be placed in a covariant perfect-fluid tensor. Nevertheless,
\begin{equation}
c_s^2=\left(\frac{\partial p}{\partial\eps}\right)_\sigma=\Gamma-1>1.
\end{equation}
The sound characteristics are superluminal. Lorentz covariance of the equations is therefore insufficient for causal admissibility. Constitutive and stability tests remain separate requirements.

\section{A corrected relativistic-fluid example}
\subsection{From extensive variables to densities}
Let $U=\eps V$, $S=sV$, and $N=nV$. Inserting these into \eqref{eq:firstlaw} gives
\begin{equation}
V\dd\eps+\eps\dd V
=TV\dd s+Ts\dd V-p\dd V+\mu V\dd n+\mu n\dd V.
\end{equation}
Extensivity implies Euler's relation $U=TS-pV+\mu N$, hence
\begin{equation}\label{eq:densityfirstlaw}
\eps+p=Ts+\mu n,\qquad\dd\eps=T\dd s+\mu\dd n.
\end{equation}
Pressure work has already entered the Euler identity. Adding $-p\dd(1/n)$ to this density differential double counts the change of variables. For specific energy $e=\eps/n$ and specific entropy $\sigma=s/n$, the correct alternative is
\begin{equation}
\dd e=T\dd\sigma-p\dd(1/n).
\end{equation}
Differentiating the Euler identity and using \eqref{eq:densityfirstlaw} also yields
\begin{equation}
\dd p=s\dd T+n\dd\mu.
\end{equation}

\subsection{Tensor contractions and the entropy current}
The perfect-fluid tensor, particle current, and entropy current are
\begin{equation}\label{eq:fluid}
\mathcal T^{\mu\nu}=(\eps+p)u^\mu u^\nu+pg^{\mu\nu},
\quad N^\mu=nu^\mu,\quad J_s^\mu=su^\mu.
\end{equation}
Contraction with the fluid velocity gives
\begin{equation}
u_\nu\mathcal T^{\mu\nu}=-\eps u^\mu,
\quad u_\mu u_\nu\mathcal T^{\mu\nu}=\eps,
\quad-u_\mu N^\mu=n.
\end{equation}
A useful covariant identity is
\begin{equation}\label{eq:entropycurrent}
\boxed{J_s^\mu=p\beta^\mu-\beta_\nu\mathcal T^{\mu\nu}-\alpha N^\mu
=\frac{\eps+p-\mu n}{T}u^\mu=su^\mu.}
\end{equation}
This follows from the Euler relation and the negative sign in $u_\nu\mathcal T^{\mu\nu}$. The entropy density is $s=-u_\mu J_s^\mu$. For $n>0$,
\begin{equation}
\frac{s}{n}=\frac{-u_\mu J_s^\mu}{-u_\mu N^\mu}.
\end{equation}
This is a scalar ratio of contractions; division of four-vectors is undefined.

An observer with unit velocity $w^\mu$ measures
\begin{equation}
\eps_{(w)}=\mathcal T^{\mu\nu}w_\mu w_\nu
=(\eps+p)(u\cdot w)^2-p.
\end{equation}
For relative speed $v$, $\eps_{(w)}=\gamma^2(\eps+pv^2)$. Proper $\eps,p,s,n,T$ remain scalar fields at the same event, while observed projections change with $w$. A scalar's coordinate independence does not mean it is spatially constant.

\subsection{Dynamical entropy test}
Set $D=u^\mu\nabla_\mu$ and $\theta=\nabla_\mu u^\mu$. Conservation of \eqref{eq:fluid} gives
\begin{equation}
D\eps+(\eps+p)\theta=0,\qquad Dn+n\theta=0.
\end{equation}
Then \eqref{eq:densityfirstlaw} implies
\begin{align}
T Ds+\mu Dn+(\eps+p)\theta&=0,\\
T(Ds+s\theta)&=0.
\end{align}
On the $T>0$ patch,
\begin{equation}\label{eq:entropyconserve}
\nabla_\mu J_s^\mu=0.
\end{equation}
This independently checks the density first law and the tensor contractions. For a dissipative fluid, additional fluxes and constitutive equations are required before an entropy-production inequality is justified. A first-law loop defect cannot simply be renamed entropy production.

Global equilibrium conventionally requires a timelike Killing inverse-temperature field and a constant chemical ratio, with compatible boundaries and no flux \cite{Becattini2016,Becattini2015}:
\begin{equation}\label{eq:killing}
\nabla_{(\mu}\beta_{\nu)}=0,\qquad\nabla_\mu\alpha=0.
\end{equation}
It does not require $\nabla_\mu\beta_\nu=0$; rotating or accelerating equilibrium provides counterexamples. In a static region with Killing vector $\xi=\partial_t$ and $\beta=\xi/T_\infty$,
\begin{equation}\label{eq:Tolman}
T(x)\sqrt{-\xi^2(x)}=T_\infty.
\end{equation}
This is consistent with stationary equilibrium and a nonuniform local temperature.

\section{Radiation: a failed zero-temperature inference}
\subsection{Low frequency is not zero temperature}
At temperature $T_0>0$, the Planck occupation number is
\begin{equation}
n(\omega)=\frac1{e^{\omega/T_0}-1}
=\frac{T_0}{\omega}-\frac12+O(\omega/T_0).
\end{equation}
The occupation of a low-frequency mode diverges, but its mean excitation energy approaches $\omega n(\omega)\to T_0$. More decisively, the two-polarization energy density per angular frequency in three spatial dimensions is
\begin{equation}\label{eq:Planck}
\frac{\dd\rho_\gamma}{\dd\omega}
=\frac{\omega^3}{\pi^2}\frac1{e^{\omega/T_0}-1}
\sim\frac{T_0}{\pi^2}\omega^2.
\end{equation}
Thus $\int_0^\delta(\dd\rho_\gamma/\dd\omega)\dd\omega
\sim T_0\delta^3/(3\pi^2)$ is finite. The photon-number density is also integrable in the infrared. Restoring units multiplies the energy spectrum by $\hbar/c^3$ and replaces $\omega/T_0$ by $\hbar\omega/(k_BT_0)$.

The limit $\omega\to0$ samples the infrared part of a state; it does not set the state's temperature. A nonthermal monochromatic field need not possess any Planck temperature at all. Neither $\dd S=0$ nor absence of entropy flux implies $T=0$: a reversible adiabatic process or a finite-temperature equilibrium state can satisfy those conditions.

\subsection{A directional Doppler limit}
Take isotropic thermal radiation at rest at $T_0$. For an observer boosted at speed $v$ in the $+x$ direction, the radiation-fluid velocity in the observer's coordinates is $u'^\mu=\gamma(1,-v,0,0)$. For a ray with observed angle $\vartheta'$ to $+x$,
\begin{equation}
-\frac{u'_\mu p'^\mu}{T_0}
=\frac{\gamma\omega'(1+v\cos\vartheta')}{T_0}.
\end{equation}
The invariant occupation can therefore be expressed using a directional spectral temperature
\begin{equation}\label{eq:Tdir}
T_{\mathrm{dir}}(\vartheta')
=\frac{T_0}{\gamma(1+v\cos\vartheta')}.
\end{equation}
For the forward ray, $T_{\mathrm{dir}}(0)=\gamma(1-v)T_0\to0$ as $v\to1$. The backward ray instead has $T_{\mathrm{dir}}(\pi)=\gamma(1+v)T_0\to\infty$. There is no isotropic transformed Planck temperature, and the invariant proper $T_0=(-\beta^2)^{-1/2}$ is unchanged.

Since $p_\gamma=\rho_\gamma/3$, the boosted energy density is
\begin{equation}
\rho'_\gamma=\gamma^2\rho_\gamma(1+v^2/3),
\end{equation}
which grows without bound in this limit. A vanishing frequency in one direction cannot establish vanishing total energy or absolute zero. Every finite $v<1$ describes a timelike observer; $v=1$ is not an inertial rest frame for light.

\section{Horizon residues with an explicit domain of validity}
\subsection{A local meromorphic one-form}
The Robson--Villari--Biancalana (RVB) method relates Euclidean curvature integrals and topology to horizon-temperature calculations \cite{RVB}. Chen's residue-based work motivates attention to contour data \cite{ChenResidue}. Here we define a particular optical one-form and derive its relation to temperature. This local construction is not an assertion that every RVB integration constant or every logarithmic derivative is the same observable.

Consider
\begin{equation}\label{eq:metric}
\dd s^2=-e^{2\psi(r)}f(r)\dd t^2+\frac{\dd r^2}{f(r)}+r^2\dd\Omega_2^2,
\end{equation}
with an outer horizon $r_h$ satisfying $f(r_h)=0$, $f'_h>0$, and finite real $\psi_h$. Fix $t$ by a boundary normalization. Surface gravity is
\begin{equation}\label{eq:kappa}
\kappa=\tfrac12 e^{\psi_h}f'_h.
\end{equation}
Assume $f$ and $\psi$ have holomorphic extensions near $z=r_h$, with $e^{\psi_h}\ne0$. Define
\begin{equation}\label{eq:optical}
\vartheta=\frac{e^{-\psi(z)}}{f(z)}\dd z.
\end{equation}
Radial null curves satisfy $\dd t=\pm\vartheta$ on the real static patch. Expanding $x=z-r_h$,
\begin{align}
f(z)&=f'_hx+O(x^2),&e^{-\psi(z)}&=e^{-\psi_h}+O(x),\\
\vartheta&=\left[\frac{e^{-\psi_h}}{f'_h}\frac1x+O(1)\right]\dd z.
\end{align}
The residue is therefore
\begin{equation}\label{eq:residue}
\mathfrak r_h=\Res_{r_h}\vartheta
=\frac{e^{-\psi_h}}{f'_h}=\frac1{2\kappa}.
\end{equation}
It has length units and contains a geometric coefficient, rather than an integer only.

\subsection{Euclidean derivation and contour check}
After $t=-\ii\tau$, set $r-r_h=f'_h\rho^2/4$. To leading order,
\begin{equation}
\dd s_E^2=\dd\rho^2+\kappa^2\rho^2\dd\tau^2+r_h^2\dd\Omega_2^2+\cdots.
\end{equation}
The absence of a conical defect requires period $\beta_H=2\pi/\kappa$ \cite{GibbonsHawking}. Combining this independent regularity condition with \eqref{eq:residue} gives
\begin{equation}\label{eq:TH}
\boxed{T_H=\frac\kappa{2\pi}=\frac1{4\pi\mathfrak r_h}.}
\end{equation}
No tunneling factor is inferred from the radial integral alone. Euclidean regularity fixes the normalization in \eqref{eq:TH}.

For a closed curve in a meromorphic domain, avoiding the poles and enclosing only finitely many of them, the residue theorem states
\begin{equation}\label{eq:winding}
\frac1{2\pi\ii}\oint_\Gamma\vartheta
=\sum_j\Wind(\Gamma,z_j)\Res_{z_j}\vartheta.
\end{equation}
For a contour enclosing only $r_h$ with nonzero winding $n_h$, the single-horizon residue is the integral divided by $n_h$. Winding number changes how often a pole is counted; it does not multiply the physical horizon temperature. If different horizons have unequal surface gravities, their local temperatures remain distinct and a single smooth Euclidean thermal period generally does not exist.

\begin{proposition}[Restricted coordinate invariance]\label{prop:residuecoordinate}
For a local holomorphic coordinate change $z=z(w)$ with $z(w_h)=r_h$ and $z'(w_h)\ne0$, the residue of the pulled-back optical one-form equals $\mathfrak r_h$.
\end{proposition}
\begin{proof}
Write $z-r_h=a_1(w-w_h)+O((w-w_h)^2)$ with $a_1\ne0$. The singular term pulls back as
$\mathfrak r_h\dd z/(z-r_h)
=\mathfrak r_h\dd w/(w-w_h)+$ regular terms. The coefficient of the simple pole is unchanged.
\end{proof}
This statement concerns a differential one-form. The residue of a scalar function written without its differential is not a coordinate-invariant prescription. Singular or branched maps require tracking multiplicity and are outside this proposition.

Time normalization is different. Under $t'=a t$ with $a>0$, $\psi'=\psi-\ln a$ and
\begin{equation}
\vartheta'=a\vartheta,\quad
\kappa'=\kappa/a,\quad T_H'=T_H/a.
\end{equation}
Only temperatures relative to the same normalized generator should be compared. The static local detector temperature is $T_{\mathrm{loc}}=T_H/(e^{\psi}\sqrt f)$, which is unchanged by this constant time rescaling.

\subsection{Logarithmic derivatives and optical residues are different}
If $F(z)=(z-z_0)^m h(z)$ with $h(z_0)\ne0$ and $m\in\mathbb Z$, then
\begin{equation}
\frac{F'}{F}\dd z=\frac{m\dd z}{z-z_0}+\frac{h'}h\dd z,
\qquad\Res_{z_0}\frac{F'}F\dd z=m.
\end{equation}
The argument principle counts zeros minus poles with multiplicities, subject to meromorphicity and a contour free of zeros and poles. It already includes zeros and poles of arbitrary finite order. Winding is part of the residue theorem's contour weighting, not an additional correction to that principle.

For the metric function at a simple horizon, $\Res(f'/f\,\dd z)=1$, irrespective of $f'_h$. For the optical form, $\Res\vartheta=e^{-\psi_h}/f'_h$, which varies continuously. A dimensionless order count cannot equal a temperature without an additional dimensional scale and a derivation of its coupling to the geometry. Multiplying a logarithmic derivative by a regular geometric weight $A(z)$ gives residue $mA(z_0)$, but the weight must be defined physically rather than selected after the desired result is known.

An essential singularity can possess a Laurent residue, but the ordinary zero-minus-pole statement is not generally applicable to it. For example, $F=e^{-a/z}$ has an essential singularity at zero while $F'/F=a/z^2$ has zero residue. The classification of a complex singularity alone does not identify an absolute-zero state or reconstruct gravitational dynamics.

\subsection{Schwarzschild and nonextremal RN checks}
For Schwarzschild, $f=1-2G_{\mathrm N}M/r$ and $\psi=0$. At $r_h=2G_{\mathrm N}M$,
\begin{equation}
\mathfrak r_h=r_h=2G_{\mathrm N}M,\qquad
T_H=\frac1{8\pi G_{\mathrm N}M}.
\end{equation}
For Reissner--Nordstr\"om, use geometric mass $m=G_{\mathrm N}M$ and a geometric charge parameter $q$, with $0\le |q|<m$:
\begin{equation}
f(r)=\frac{(r-r_+)(r-r_-)}{r^2},\qquad
r_\pm=m\pm\sqrt{m^2-q^2}.
\end{equation}
At the outer horizon,
\begin{equation}\label{eq:RN}
f'(r_+)=\frac{r_+-r_-}{r_+^2},\qquad
\mathfrak r_+=\frac{r_+^2}{r_+-r_-},\qquad
T_+=\frac{r_+-r_-}{4\pi r_+^2}.
\end{equation}
As $q\to m^-$, the simple-horizon residue diverges and the temperature tends to zero. A contour isolating $r_+$ must shrink as the horizons merge; holding a contour around both poles is a different calculation.

\subsection{Extremality as a decisive failure test}
At exactly $q=m$, the horizon is a double zero and
\begin{equation}\label{eq:extremal}
\vartheta=\frac{r^2\dd r}{(r-m)^2}
=\left[\frac{m^2}{(r-m)^2}+\frac{2m}{r-m}+1\right]\dd r.
\end{equation}
Its residue is $2m$, although $f'(m)=0$ and the standard surface-gravity temperature is $T_H=0$. Applying \eqref{eq:TH} outside its domain would incorrectly give $1/(8\pi m)$. The degenerate Euclidean throat is not the nondegenerate polar neighborhood used above.

Thus nonzero residue is not a sufficient condition for nonzero temperature. Higher-order poles require their full principal part and a separate physical analysis. The contrast between the divergent simple-pole residue in the near-extremal limit and the finite residue at exact extremality is mathematically consistent: two poles have merged, and the contours and pole structure have changed.

\section{Why horizon thermodynamics does not reconstruct all dynamics}
\subsection{An explicit noninjectivity example}
Let $r_h,L>0$ and $\lambda\ge0$. On the static exterior, compare
\begin{align}
f_0(r)&=1-\frac{r_h}{r},\\
f_\lambda(r)&=1-\frac{r_h}{r}
+\lambda\frac{(r-r_h)^2}{r_h^2}e^{-(r-r_h)/L},\qquad\psi=0.\label{eq:family}
\end{align}
For every member,
\begin{equation}
f_\lambda(r_h)=0,\quad f'_\lambda(r_h)=1/r_h,
\quad A_h=4\pi r_h^2,\quad
\mathfrak r_h=r_h,\quad T_H=1/(4\pi r_h).
\end{equation}
The additional term is nonnegative outside the horizon, and decays faster than $1/r$ at infinity. There is no additional exterior zero. Each metric has the same asymptotic mass $M=r_h/(2G_{\mathrm N})$. In an Einstein--Hilbert area-law setting its horizon entropy is also the same, $S_h=\pi r_h^2/G_{\mathrm N}$.

The metrics nevertheless have different curvature. For \eqref{eq:metric} with $\psi=0$,
\begin{equation}
R=-f''-\frac{4f'}r+\frac{2(1-f)}{r^2}.
\end{equation}
At the horizon,
\begin{equation}\label{eq:curvaturedifference}
f''_\lambda(r_h)-f''_0(r_h)=\frac{2\lambda}{r_h^2},
\qquad R_\lambda(r_h)-R_0(r_h)=-\frac{2\lambda}{r_h^2}.
\end{equation}
Therefore even $(M,A_h,T_H,\mathfrak r_h)$ does not determine the exterior geometry in this metric class.

These are kinematic metric profiles, not a family of vacuum solutions of one specified theory. For $\lambda\ne0$, supporting stress-energy or altered dynamics is required; energy conditions have not been established. In Einstein gravity an effective source can be defined by $8\pi G_{\mathrm N}\mathcal T_{\mu\nu}=G_{\mu\nu}[g_\lambda]$, but this does not provide a microscopic matter model. The example proves noninjectivity before field equations and matter assumptions are supplied. It does not challenge vacuum uniqueness theorems that use those additional assumptions.

\subsection{A one-way map with additional reconstruction data}
The local extraction map is
\begin{equation}
(g,\xi,\text{boundary normalization, entropy law})
\longmapsto(A_h,\kappa,T_H,S_h,\mathfrak r_h).
\end{equation}
Its inverse is underdetermined. Reconstruction can require an action, matter equations, connection data, boundary conditions, constraints at all spacetime points, and a microscopic or effective entropy law. A pole supplies local thermal information under stated assumptions; it does not provide all those data.

Jacobson's field-equation derivation invokes a local horizon relation for all appropriate null directions, area entropy, and additional conservation and geometric assumptions \cite{Jacobson}. That is substantially stronger information than the temperature of one horizon. Entropy-based gravity frameworks also specify the entropy density and the variational constraints \cite{Padmanabhan}. These results cannot be reproduced by an undefined invariant-extraction operator alone.

\subsection{Geometric entropy functionals require a separate definition}
Consider the geometric quantity appearing in the research program,
\begin{equation}\label{eq:Sxi}
\mathcal S_\xi[g,\xi]
=\int_{\mathcal V}\nabla_\mu\xi_\nu\nabla^\mu\xi^\nu\sqrt{-g}\,\dd^4x.
\end{equation}
With the Levi-Civita connection, integration by parts yields
\begin{equation}
\mathcal S_\xi
=\int_{\partial\mathcal V}\xi_\nu\nabla^\mu\xi^\nu\dd\Sigma_\mu
-\int_{\mathcal V}\xi_\nu\nabla_\mu\nabla^\mu\xi^\nu\sqrt{-g}\,\dd^4x.
\end{equation}
For a Killing field and the curvature convention
$[\nabla_\mu,\nabla_\nu]v^\rho=R^\rho{}_{\sigma\mu\nu}v^\sigma$,
the identity $\nabla_\mu\nabla^\mu\xi^\nu=-R^\nu{}_{\lambda}\xi^\lambda$ gives
\begin{equation}\label{eq:SxiR}
\mathcal S_\xi
=\int_{\partial\mathcal V}\xi_\nu\nabla^\mu\xi^\nu\dd\Sigma_\mu
+\int_{\mathcal V}R_{\nu\lambda}\xi^\nu\xi^\lambda\sqrt{-g}\,\dd^4x.
\end{equation}
The relation displays both a boundary contribution and curvature-dependent bulk information. It is not a purely topological identity. Under $\xi\to a\xi$, $\mathcal S_\xi\to a^2\mathcal S_\xi$. The Lorentzian contraction need not be positive. A normalization, entropy units, boundary prescription, matter terms, and precise variation must therefore be supplied before identifying \eqref{eq:Sxi} with physical entropy or deriving field equations from it. No such derivation is assumed here.

Equation~\eqref{eq:potentialstationary} concerns statistical density operators. Equation~\eqref{eq:SxiR} concerns a geometric functional. Their unification is an open construction problem requiring an explicit map, not an established equality.

\subsection{Potential changes and gravitational energy definitions}
In extended AdS thermodynamics, an appropriate mass can be interpreted as enthalpy after introducing pressure associated with the cosmological constant \cite{Kastor}. That setting specifies a first law and its independent variables. It does not imply that ADM, Komar, and Misner--Sharp energies are universally Legendre transforms. Comparisons of gravitational energies must specify asymptotics, generators, boundaries, and ensembles before treating their differences as representational.

\section{Validation ledger and physical scope}
The calculations above define three tests with different dimensions and meanings:
\begin{center}
\begin{tabular}{@{}p{0.22\linewidth}p{0.32\linewidth}p{0.36\linewidth}@{}}
\toprule
Object & Explicit definition & What a mismatch can establish\\
\midrule
First-law defect & $K=\partial_VT+\partial_Sp$ and $\oint\Delta$ & Incompatible candidate thermodynamic derivatives on a smooth patch\\
Entropy-current test & $J_s^\mu$ and $\nabla_\mu J_s^\mu$ & Failure of the stated ideal-fluid conservation assumptions or constitutive identities\\
Horizon residue & $\mathfrak r_h=\Res(e^{-\psi}\dd z/f)$ & Failure of the simple-horizon formula, normalization, or proposed contour comparison\\
\bottomrule
\end{tabular}
\end{center}
They are not one universal residue operator. Their agreement with their own anchors is necessary within their stated domains, but cannot prove complete physical equivalence of two theories.

The accompanying symbolic calculation independently checks the polynomial loop defect, a Lorentz-boosted perfect-fluid tensor, the entropy-current identity, Planck infrared scaling, Schwarzschild and RN residues, exact extremality, and the curvature difference in \eqref{eq:curvaturedifference}. These are equation-level checks. No new observations, simulations of microscopic matter, or solutions of modified-gravity field equations are claimed.

The restricted statements in Propositions~\ref{prop:integrability} and~\ref{prop:residuecoordinate} and the temperature formula \eqref{eq:TH} are derived here from stated inputs. Their constituent mathematical tools are established. The contribution is their organization into an explicit discrepancy-testing framework and the associated boundary counterexamples. Priority or a comprehensive absence of similar syntheses in the literature is not established by these calculations.

For the scope, definitions, algebraic consistency, independent reproduction, and accessible counterexample tests, the present manuscript supplies explicit evidence. Broader novelty comparisons and a dynamical completion of \eqref{eq:Sxi} remain open. Journal-specific submission requirements are outside the present scope.

\section{Conclusion}
Thermodynamic representation changes and Lorentz frame changes are compatible operations in a restricted covariant equilibrium model, but they are not the same symmetry. Constrained entropy maximization yields a precise potential relation, $\Omega=-T\mathscr M[\rho_*]$, once the multiplier terms and rest-frame convention are fixed. The first-law loop defect provides an independent integrability test, and the optical residue yields $T_H=(4\pi\mathfrak r_h)^{-1}$ at a static simple horizon with fixed Killing normalization.

The strongest conclusions come from the failure boundaries. Low-frequency occupation does not require a zero-temperature radiation state. The exact extremal RN example invalidates an unrestricted identification of nonzero residue with nonzero temperature. The equal-horizon-data metric family demonstrates that thermal boundary data alone do not reconstruct exterior curvature or dynamics. These results preserve a tractable route from certain identities to testable extensions while keeping speculative geometric entropy constructions distinct from established statistical entropy.

\section*{Declarations}
\noindent\textbf{Data and code availability.} No new empirical data were generated. 
\noindent\textbf{Conflict of interest.} The author declares no conflict of interest.

\appendix
\section{Representation dictionary with explicit limits}
\begin{longtable}{@{}p{0.25\linewidth}p{0.31\linewidth}p{0.34\linewidth}@{}}
\toprule
Thermodynamic object & Relativistic object & Scope of the comparison\\
\midrule\endhead
Potential chart & Inertial coordinate frame & Both change a description; their spaces and transformation groups differ.\\
Strict Legendre map & Passive Lorentz map & Each preserves its own defining structure, not the other's.\\
Contact one-form $\eta$ & Metric $g_{\mu\nu}$ & A one-form and a symmetric rank-two tensor are not interchangeable invariants.\\
Maxwell integrability & Lorentz covariance & Independent checks; integrability does not derive spacetime symmetry.\\
Scalar equation of state & Covariant constitutive tensor & Scalar thermodynamic arguments can define a tensor; dynamics and causality need additional input.\\
$U=F+TS$ & $E^2-\boldsymbol P^2=m^2$ & The first is a potential relation; the second is a mass-shell relation. Neither derives the other.\\
Loop defect & Horizon contour residue & Distinct differential forms and dimensions; both support specified consistency tests.\\
\bottomrule
\end{longtable}



\begin{thebibliography}{99}
\bibitem{Callen} H.~B.~Callen, \emph{Thermodynamics and an Introduction to Thermostatistics}, 2nd ed. (Wiley, 1985).
\bibitem{Mrugala} R.~Mruga\l a, ``Geometrical formulation of equilibrium phenomenological thermodynamics,'' \emph{Rep. Math. Phys.} \textbf{14}, 419--427 (1978), \href{https://doi.org/10.1016/0034-4877(78)90010-1}{doi:10.1016/0034-4877(78)90010-1}.
\bibitem{Becattini2016} F.~Becattini, ``Thermodynamic equilibrium in relativity: four-temperature, Killing vectors and Lie derivatives,'' \emph{Acta Phys. Pol. B} \textbf{47}, 1819 (2016), \href{https://doi.org/10.5506/APhysPolB.47.1819}{doi:10.5506/APhysPolB.47.1819}, \href{https://arxiv.org/abs/1606.06605}{arXiv:1606.06605}.
\bibitem{Becattini2015} F.~Becattini, L.~Bucciantini, E.~Grossi, and L.~Tinti, ``Local thermodynamical equilibrium and the $\beta$ frame for a quantum relativistic fluid,'' \emph{Eur. Phys. J. C} \textbf{75}, 191 (2015), \href{https://doi.org/10.1140/epjc/s10052-015-3384-y}{doi:10.1140/epjc/s10052-015-3384-y}.
\bibitem{GibbonsHawking} G.~W.~Gibbons and S.~W.~Hawking, ``Action integrals and partition functions in quantum gravity,'' \emph{Phys. Rev. D} \textbf{15}, 2752 (1977), \href{https://doi.org/10.1103/PhysRevD.15.2752}{doi:10.1103/PhysRevD.15.2752}.
\bibitem{RVB} C.~W.~Robson, L.~Di Mauro Villari, and F.~Biancalana, ``Topological nature of the Hawking temperature of black holes,'' \emph{Phys. Rev. D} \textbf{99}, 044042 (2019), \href{https://doi.org/10.1103/PhysRevD.99.044042}{doi:10.1103/PhysRevD.99.044042}, \href{https://arxiv.org/abs/1810.09322}{arXiv:1810.09322}.
\bibitem{ChenResidue} W.-X.~Chen, ``Calculating the Hawking temperature of black holes in $f(Q)$ gravity using the RVB method: a residue-based approach,'' \href{https://doi.org/10.1139/cjp-2024-0214}{doi:10.1139/cjp-2024-0214}, \href{https://arxiv.org/abs/2501.17485}{arXiv:2501.17485} (2025).
\bibitem{Jacobson} T.~Jacobson, ``Thermodynamics of spacetime: The Einstein equation of state,'' \emph{Phys. Rev. Lett.} \textbf{75}, 1260--1263 (1995), \href{https://doi.org/10.1103/PhysRevLett.75.1260}{doi:10.1103/PhysRevLett.75.1260}, \href{https://arxiv.org/abs/gr-qc/9504004}{arXiv:gr-qc/9504004}.
\bibitem{Padmanabhan} T.~Padmanabhan, ``Thermodynamical aspects of gravity: New insights,'' \emph{Rep. Prog. Phys.} \textbf{73}, 046901 (2010), \href{https://doi.org/10.1088/0034-4885/73/4/046901}{doi:10.1088/0034-4885/73/4/046901}, \href{https://arxiv.org/abs/0911.5004}{arXiv:0911.5004}.
\bibitem{Kastor} D.~Kastor, S.~Ray, and J.~Traschen, ``Enthalpy and the mechanics of AdS black holes,'' \emph{Class. Quantum Grav.} \textbf{26}, 195011 (2009), \href{https://doi.org/10.1088/0264-9381/26/19/195011}{doi:10.1088/0264-9381/26/19/195011}, \href{https://arxiv.org/abs/0904.2765}{arXiv:0904.2765}.
\end{thebibliography}
\end{document}
