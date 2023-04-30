E XPLORING P OSITIVE N OISE IN E STIMATION T HEORY

arXiv:1910.01569v2 [eess.SP] 24 Nov 2019

A P REPRINT
Kamiar Radnosrati
Department of Electrical Engineering
Linköping University
Linköping, Sweden
kamiar.radnosrati@liu.se

Gustaf Hendeby
Department of Electrical Engineering
Linköping University
Linköping, Sweden
gustaf.hendeby@liu.se

Fredrik Gustafsson
Department of Electrical Engineering
Linköping University
Linköping, Sweden
fredrik.gustafsson@liu.se

November 26, 2019

A BSTRACT
Estimation of a deterministic quantity observed in non-Gaussian additive noise is explored via order
statistics approach. More specifically, we study the estimation problem when measurement noises
either have positive supports or follow a mixture of normal and uniform distribution. This is a problem of great interest specially in cellular positioning systems where the wireless signal is prone to
multiple sources of noises which generally have a positive support. Multiple noise distributions are
investigated and, if possible, minimum variance unbiased (MVU) estimators are derived. In case of
uniform, exponential and Rayleigh noise distributions, unbiased estimators without any knowledge
of the hyper parameters of the noise distributions are also given. For each noise distribution, the
proposed order statistic-based estimator’s performance, in terms of mean squared error, is compared
to the best linear unbiased estimator (BLUE), as a function of sample size, in a simulation study.
Keywords Order statistics · Estimation · Non-Gaussian noise · Mean squared error

1

Introduction

We consider the problem of estimating the mean x observed in noise as yk = x + ek , for k = 1, 2, . . . , N , also known
as “estimation of location" [Kassam and Poor, 1985], where the noise ek has positive support. We will refer to such
distributions as positive noise. Examples of distributions we will study include uniform, exponential, Rayleigh, Pareto.
A bias compensated linear estimator as the sample mean has a variance that decays as 1/N , while it is well-known
from the statistical literature, see for example [Kay, 1993, Lehmann and Casella, 1998], that the minimum has a
variance that decays as 1/N 2 . The minimum is the simplest example of order statistics. Certain care has to be taken
for the cases where the parameters in the distributions are unknown, in which case bias compensation becomes tricky.
This paper derives all combinations of known/unknown parameters for order statistics/BLUE (best linear unbiased
estimator) for some selected and common distribution that allow for analytical solutions.
Problems involving positive noise can be motivated from applications where the arrival times of radio or sound waves
are used. Such waves travel with the speed of the medium, and non line of sight conditions give rise to delayed arrival
times. Physics does simply not allow for negative noise, only positive one. This case occur in a variety of applications
such as target tracking using radar or lidar, and localisation using radio waves such as is done in for instance global
satellite navigation systems [Kok et al., 2015, Chen et al., 2009, Gustafsson and Gunnarsson, 2005, Eling, 2012].

A PREPRINT - N OVEMBER 26, 2019

Figure 1: The error histograms of time-of-arrival measurements collected from three separate cellular antennas described in [Medbo et al., 2009].

For example, the error histograms of time-of-arrival measurements collected from three separate cellular antennas are
given in Figure 1. For detailed description of hardware and the measurement campaign see [Medbo et al., 2009]
To deal with the estimation’s performance degradation in non-Gaussian error conditions, conventional estimation techniques which are developed based on Gaussian assumptions need to be adjusted properly. As discussed in [Yin et al.,
2013], “identify and discard”, “mathematical programming”, and “robust estimation” are the three broad categories of
estimation methods which are robust against non-Gaussian errors. Robustness of the estimator has been a concern for
many years in both research [Stigler, 1973] and different engineering topics [Kassam and Poor, 1985, Kassam, 1988,
Stewart, 1999, Arce, 2004] for a long time now. A more recent survey on this topic containing more references can be
found in [Zoubir et al., 2012].
The maximum likelihood estimator (MLE), developed under Gaussian assumptions, can be modified to become robust
in presence of non-Gaussian noises. The authors in [Eskin, 2000] first detect and then reject the outliers by learning
the probability density function (PDF) of the measurements and develop a mixture model of outliers and clean data. A
similar idea to k-nearest-neighbor approach is used in [Chawla et al., 2010] to classify outliers as the data points that
are less similar to a number of neighboring measurements. Surveys of advances in clustering the data into outliers and
clean data can be found in [Hodge and Austin, 2004, Yin et al., 2013, Fritsche et al., 2009]. While these approaches
might result in high estimation accuracy, they typically require large datasets [Zoubir et al., 2012].
M-estimators [Huber and Ronchetti, 2009], in the contrary to identification-based methods, do not require preprocessing and can be used in non-Gaussian noise conditions. In principle, M-estimators can be seen as generalization
of MLE and rely on solving a minimization problem of some loss function. For a detailed discussion on different loss
functions, see [Huber and Ronchetti, 2009]. Since minimization problems are typically solved numerically based on
the derivative of the loss function [Maronna et al., 2006], they might converge to local minima.
In this work, we strive to find minimum variance unbiased (MVU) estimators for the location of estimation problems
for non-Gaussian noise distributions where multiple distributions with positive support are considered. In case where
MVU is not found, we introduce unbiased order-statistic-based estimators and compare their variances against the
BLUE. The MVU estimators without any knowledge of the hyper parameters of the noise distributions are also derived,
if possible. Finally, we derive an estimator for the case in which the noise follows a mixture of normal and uniform
distribution.
The rest of this paper is structured as follows. In Section 2 the marginal distribution of order statistics is introduced.
In Section 3 the location estimation problem is formulated. The problem is then investigated for different noise
distributions and estimators for each distribution are derived in Sections 4–6. The proposed estimators are evaluated
in a simulation study in Section 8 followed by the concluding remarks given in Section 9.

2

Marginal Distribution of Order Statistics

The marginal distribution of order statistics, in this work is computed by differentiating the corresponding cumulative
distribution function (CDF). In this section, we first introduce the minimum, also know as first or extreme, order
statistic and then give the generalization to any statistics of order k. Let F denote the common CDF of N independent
and identically distributed sample of random variables y1 , . . . , yN . We let y(k) denote the k:th order statistic of the
sample, defined as the k:th smallest value of the set {yi }N
i=1 . We define f(k,N ) (y) as the marginal PDF of the k:th
2

A PREPRINT - N OVEMBER 26, 2019

order statistics corresponding to a sample of size N . The PDF f(k,N ) (y) is then calculated by differentiating F(k,N ) (y)
with respect to y.
2.1

Marginal distribution of minimum order statistic

To further illustrate the problem, consider fist an example in which we have drawn N = 5 independent random
variables {yi }5i=1 each from a common distribution with PDF f (y). Assume that we are interested in the PDF of the
first order statistic, f(1,5) (y). The CDF F(1,5) (y) is defined as P (y(1) < y). We note that the minimum order statistic
Y(1) would be less than y if at least 1 of the random variables y1 , y2 , y3 , y4 , y5 are less than y. In other words, we
need to count the number of ways that can happen such that at least one random variable is less than y. This leads to a
binomial probability calculation. The ’success’ is considered to be the event {yi < y}, i = 1 and we let ζ denote the
number of successes in five trials, then
F(1,5) (y) = P (y(1) < y) = P (ζ = 1) + . . . + P (ζ = 5),
d
f(1,5) (y) =
F(1,5) (y).
dy
To generalize the example, let y(1) < y(2) < . . . < y(N ) be the order statistics of N independent observations from
a continuous distribution with cumulative distribution function F (y) and probability density function f (y) = F 0 (y).
The marginal PDF f(1,N ) (y) of the minimum order statistic can be obtained by considering the event {Yi ≤ y}, i = 1
as a "success," and letting ζ = the number of such successes in N mutually independent trials. ζ is a binomial random
variable with N trials and probability of success P (yi ≤ y). Hence, the CDF of the minimum order statistic is given
by,
F(1,N ) (y) =

N
X

P (ζ = n).

(1a)

n=1

Noting that the probability mass function of this binomial distribution is given by,
 
N
P (ζ = n) =
[F (y)]n [1 − F (y)]N −n .
n

(1b)

Substituting (1b) into (1a) and taking the last term out of the sum, we get
N
−1  
X
N
F(1,N ) (y) =
[F (y)]n [1 − F (y)]N −n + [F (y)]N .
n

(1c)

n=1

Differentiating (1c) with respect to y gives a telescoping sum of the form,
f(1,N ) (y) =

N
−1
X
n=1

+

N
−1
X
n=1

N!
[F (y)]n−1 f (y)[1 − F (y)]N −n
(n − 1)!(N − n)!
N!
[F (y)]n [1 − F (y)]N −n−1 (−f (y))
n!(N − n − 1)!

+ N [F (y)]N −1 f (y),

(2)

in which, except the first term, all other terms cancel each other out. Hence, the marginal probability density function
of the minimum order statistic of a set of N independent and identically random variables with common CDF F (y)
and PDF f (y) is given by,
f(1,N ) (y) = N f (y) (1 − F (y))
2.2

N −1

.

(3)

Marginal distribution of general order statistic

The marginal PDF f(k,N ) (y) of the general order statistic k can be obtained by generalizing the results of the previous
section, and considering the event {yi ≤ y}, i = 1, 2, . . . , k as a "success," and letting ζ = the number of such
successes in N mutually independent trials,
N
−1  
X
N
F(k,N ) (y) =
[F (y)]n [1 − F (y)]N −n + [F (y)]N .
(4)
n
n=k

3

A PREPRINT - N OVEMBER 26, 2019

Table 1: Notation.
noisy measurements of the unknown parameter x

{yk }N
k=1
N
y(m) m=1

ordered measurement sequence

θ

parameters of the noise distribution

δ(θ)

bias compensation term

e
x̂pBLUE
(y1:N , θ)
pe
x̂MVU (y1:N , θ)
e
x̂pMVU
(y1:N )
pe
x̂ (y1:N , θ)

BLUE when ek ∼ pe for known θ
MVU estimator when ek ∼ pe for known θ
MVU estimator when ek ∼ pe for unknown θ
unbiased estimator when ek ∼ pe for known θ

pe

unbiased estimator when ek ∼ pe for unknown θ

x̂ (y1:N )

Differentiating (4) with respect to y gives a telescoping sum of the form,
f(k,N ) (y) =

N
−1
X

N!
[F (y)]n−1 f (y)[1 − F (y)]N −n
(n − 1)!(N − n)!

n=k

+

N
−1
X
n=k

N!
[F (y)]n [1 − F (y)]N −n−1 (−f (y))
n!(N − n − 1)!

+ N [F (y)]N −1 f (y),
(5)
in which, except the first term, all other terms cancel each other. Hence, the marginal probability density function of
the k:th order statistic of a set of N independent and identically random variables with common CDF F (y) and PDF
f (y) is given by,


N −1
N −k
f(k,N ) (y) = N f (y)
F (y)k−1 (1 − F (y))
.
(6)
k−1

3

Location Estimation Problem

Consider the location estimation problem in which we have measurements yk , k = 1, . . . , N of the unknown parameter
x. Assuming that the measurements are corrupted with additive noise ek ∼ pe (θ), where θ denotes the parameter(s)
of the noise distribution, the measurement model is given by
yk = x + ek , k = 1, . . . , N.
(7)
The BLUE for the estimation problem (7) is given by
p (θ)

e
x̂BLUE
(y1:N , θ) =

N
1 X
yk − δ(θ),
N

(8)

k=1

where y1:N = {yk }N
k=1 and δ(θ) = E(ek ) is the bias compensation term.
In the following sections, closed-form expressions for the mean squared error (MSE) of the BLUE estimator for
multiple noise distributions with positive support are provided. Given hyperparameter θ, the MVU estimator for
e
each noise distribution pe is denoted by x̂pMVU
(y1:N , θ). MVU estimators with unknown hyperparameter are denoted
pe
by x̂MVU (y1:N ). If the MVU cannot be found, an unbiased order-statistics-based estimator is derived and denoted by
x̂pe (y1:N , θ) and x̂pe (y1:N ) for known and unknown hyperparameter cases, respectively. For example, x̂U
MVU (y1:N , β)
denotes the MVU estimator when ek ∼ U[0, β] and β is known. x̂U
(y
),
on
the
other
hand,
corresponds
to the
1:N
MVU
MVU estimator of uniform noise with unknown hyper parameters of the distribution. Table 1 summarizes the notation
used throughout this work.
e
For each noise distribution, we also consider the minimum order statistic estimator, denoted by x̂pmin
(y1:N ). Let
N
pe
y(m) m=1 denote the ordered sequence obtained from sorting y1:N in an ascending order, x̂min (y1:N ) is defined as

pe
x̂min
(y1:N ) = y(1) , min yk .
k

4

(9)

A PREPRINT - N OVEMBER 26, 2019

Noting that for any generic estimator x̂ the MSE is given by
MSE(x̂) = Var(x̂) + b2 (x̂).

(10)

The MSE for BLUE and MVU or any other bias compensated estimator coincides with the estimator’s variance. In
e
case of x̂pmin
, the existing bias enters the MSE.
In order to find the MVU estimator, the first step is to find the PDF f (y1:N ; θ) with θ denoting the parameters of the
distribution. If the PDF satisfies regularity conditions, the CRLB can be determined. Any unbiased estimator that
satisfies CRLB is thus the MVU estimator. However, the considered PDFs do not satisfy the regularity conditions,


∂ ln f (yk ; θ)
E
6= 0.
(11)
∂θ
Hence, the CRLB approach is not applicable. Instead, we rely on the RBLS theorem [Lehmann and Scheffé, 1950,
1955, Kay, 1993], to find the MVU estimator. The RBLS theorem [Lehmann and Scheffé, 1950] states that for any
unbiased estimator x̃ and sufficient statistics T (y1:N ), x̂ = E(x̃ | T (y1:N )) is unbiased and Var(x̂) ≤ Var(x̃).
Additionally, if T (y1:N ) is complete, then x̂ is MVU.
As shown in [Kay, 1993], if the dimension of the sufficient statistics is equal to the dimension of the parameter, then
the MVU estimator is given by x̂ = g(T (y1:N )) for any function g(·) that satisfies
E(g(T )) = x.

(12)

Hence, the problem of MVU estimator turns into the problem of finding a complete sufficient statistic. The NeymanFisher theorem [Fisher and Phil, 1922, Halmos and Savage, 1949] gives the sufficient statistic T (y1:N ), if the PDF can
be factorized as follows
f (y1:N ; Ψ) = g(T (y1:N ), Ψ)h(y1:N ),

(13)

where ψ is the union of the noise hyper parameter(s) θ and x. The estimators in this work are derived in the order
statistics framework.

4

Uniform Distribution

As the first scenario, consider the case in which the additive noise ek has a uniform distribution with a positive support,
pe (θ) = U[0, β], β > 0 and θ = β. The BLUE is given by
x̂U
BLUE (y1:N , β)

N
1 X
β
=
yk − .
N
2

(14a)

k=1

The MSE of BLUE for this case is given by


N

β
1 X
β2
MSE x̂U
(y
,
β)
=
Var
y
−
=
.
1:N
k
BLUE
N2
2
12N

(14b)

k=1

In order to find the MSE of the minimum order statistics estimator, x̂U
min (y1:N ), we need to find the first two moments
of the estimator. Let ỹk = β1 yk . Since yk ∼ U[x, x + b], then for any constant β > 0, y˜k ∼ U[ βx , βx + 1]. Hence,
f (ỹk ) = 1 and F (ỹk ) = β1 (yk − x). From (6) we get,


k−1 
N −k
ỹ − x
β − (ỹ − x)
N −1
U [0,β]
f (k,N ) (ỹ) = N
k−1
β
β

k−1 
N −k
ỹ − x
β − (ỹ − x)
(N )!
.
(15a)
=
(k − 1)!(N − k)!
β
β
since N ∈ N+ , k ∈ N+ , and k ∈ [1, N ] we can the change the factorials to gamma functions,

k−1 
N −k
Γ(N + 1)
ỹ − x
β − (ỹ − x)
U [0,β]
f(k,N ) (ỹ) =
.
Γ(k)Γ(N − k + 1)
β
β

(15b)

The marginal distribution (15b) is a generalized beta distribution, also known as four parameters beta distribution [McU [0,β]
U [0,1]
Donald and Xu, 1995]. The support of this distribution is from 0 to β > 0 and f(k,N ) (·) = β1 f(k,N ) (·). The bias
5

A PREPRINT - N OVEMBER 26, 2019

and variance of the general k:th order statistic estimator x̂U
(k) (y1:N ) in case of uniform noise with support on [0, β] are
given by
βk
,
N +1
k(N − k + 1)β 2
Var(x̂U
.
(k) (y1:N )) =
(N + 1)2 (N + 2)
E(x̂U
(k) (y1:N )) =

(16a)
(16b)

The first two moments of the minimum order statistic estimator are obtained by letting k = 1 in (16)

β
b x̂U
min (y1:N ) =
N +1

N β2
.
Var x̂U
min (y1:N ) =
(N + 1)2 (N + 2)

(17a)
(17b)

The MSE of x̂U
min (y1:N ) is then given by

MSE x̂U
min (y1:N ) =
4.1

2β 2
.
(N + 1)(N + 2)

(18)

MVU estimator

In order to find the MVU estimator, we note that the PDF can be written in a compact form using the step function
σ(·) as
1
f (yk ; x, β) = [σ(yk − x) − σ(yk − x − β)] .
(19a)
β
which gives
f (y1:N ;x, β) =

N
1 Y
[σ(yk − x) − σ(yk − x − β)]
βN
k=1


1 
= N σ(y(1) − x) − σ(y(N ) − x − β) ,
β

(19b)

where y(N ) , maxk yk , k = 1, . . . , N . The expressions for the MVU estimator is derived for two different
scenarios. We first assume that the hyper parameter β of the noise distribution is known and then further discuss
the unknown hyper parameter case. In the general case, let Ψ = [x, β]> denote the unknown parameter vector, the
Neyman-Fisher factorization gives h(y1:N ) = 1 and
 


y(1)
T1 (y1:N )
=
.
(20)
T (y1:N ) =
T2 (y1:N )
y(N )
4.1.1

Known hyper parameter β

When the maximum support of the uniform noise β is known, the dimensionality of the sufficient statistic is larger
than that of the unknown parameter x. As discussed in [Kay, 1993], the RBLS theorem can be extended to address
this case if the form of a function g(T1 (y1:N ), T2 (y1:N )) can be found that combines T1 and T2 into a single unbiased
estimator of x.
Let Z = T1 (y1:N ) + T2 (y1:N ) = u + v. Since T1 and T2 are dependent,
Z ∞
fZ (z) =
fy(1) ,y(N ) (u, z − u) du ,

(21a)

−∞

where fy(1) ,y(N ) (u, z − u) is the joint density of minimum and maximum order statistics. As shown in [David and
Nagaraja, 2004], for −∞ < u < v < ∞, the joint density of two order statistics y(i) and y(j) is given by
fy(i) ,y(j) (u, v) =

N!
(i − 1)!(j − 1 − i)!(N − j)!
i−1

× fY (u)fY (v) [FY (u)]

j−1−i

× [FY (v) − FY (u)]
6

N −j

[1 − FY (v)]

,

(21b)

A PREPRINT - N OVEMBER 26, 2019

that for the extreme orders, i = 1 and j = N can be simplified such that for u < v
N −2

fy(1) ,y(N ) (u, v) = N (N − 1) [FY (v) − FY (u)]

fY (u)fY (v).

(21c)

and zero otherwise. Substituting (21c) into (21a), we get
fZ (z) =

1
N β −N (2x + 2β − z)N −1 ,
2

(21d)

for 2x + β < z < 2(x + β) and
(z − 2x)N
1
,
fZ (z) = − N β −N
2
2x − z

(21e)

for 2x < z ≤ 2x + β and zero otherwise. It can be shown that
E(fZ (z)) = 2x + β.

(21f)

Hence, noting that β is known, the function g(T1 (y1:N ), T2 (y1:N )) that gives an unbiased estimator should be of the
form of
x̂U
MVU (y1:N , β) = g(T1 (y1:N ), T2 (y1:N ))
β
1
= (y(1) + y(N ) ) − .
2
2

(22a)

The MSE of the MVU estimator is given by

MSE x̂U
MVU (y1:N , β) =

β2
.
2N (N + 3) + 4

(22b)

Comparing to (14b), the order statistics based MVU estimator outperforms the BLUE one order of magnitude.
4.1.2

Unknown hyper parameter β

In this case, the MVU estimators for the parameter vector Ψ = [x, β]> can be derived from sufficient statistics (20),
Ψ̂ = g(T (y1:N )),

s.t. E (g (T (y1:N ))) = Ψ.

(23)

In this case, we have

x+
E(T (y1:N )) = 
x+

β 
N +1
Nβ
N +1



To find the transformation that makes (24) unbiased, we define
 1

N −1 (N T1 (y1:N ) − T2 (y1:N ))

g(T (y1:N )) = 
N +1
N −1 (T2 (y1:N ) − T1 (y1:N ))

(24)

(25a)

that gives
 
x
E (g(T (y1:N ))) =
.
β

(25b)

Finally, the MVU estimator of x when the hyper parameter β is unknown is given by
x̂U
MVU (y1:N ) =

N
1
y(1) −
y(N ) .
N −1
N −1

(26a)

N β2
.
(N + 2)(N 2 − 1)

(26b)

and its MSE is

MSE x̂U
MVU (y1:N ) =
This is naturally slightly larger than (22b) for finite N .
7

A PREPRINT - N OVEMBER 26, 2019

5

Distributions in the exponential family

The exponential family of probability distributions, in their most general form, is defined by
f (y; θ) = h(y)g(θ) exp {A(θ) · T (y)} ,
(27)
where θ is the parameter of the distribution, and h(y), g(θ), A(θ), and T (y) are all known functions. In this section,
we only consider some example distributions of this family and show that the minimum order statistic estimator gets
the same form of distribution as the noise distribution but with modified parameters. For the selected distributions, if
possible, MVU estimators for both cases of known and unknown hyperparameter are derived. Otherwise, unbiased
estimators with less variance than BLUE are proposed.
5.1

Exponential distribution

Exponential distributions are members of the gamma family with shape parameter 1; strongly skewed with no left
sided tail (yk ∈ [x, ∞]). Let β > 0 denote the scale parameter, the PDF of an exponential distribution is then given by
f Exp (yk ; x, β) =

1
β

exp(− (ykβ−x) ) yk ≥ x,
0
yk < x.

(28a)

and the CDF, for y ≥ x, is given by
F Exp (yk ; x, β) = 1 − exp(−

(yk − x)
).
β

(28b)

For the BLUE estimator (8), from the properties of exponential distribution, we have
x̂Exp
BLUE (y1:N , β)

N
1 X
=
yk − β,
N

MSE(x̂Exp
BLUE ) =

k=1

β2
.
N

Substituting (28) into (6), the marginal density of the k:th order statistic is given by


k−1


(y − x)
(N − k + 1)(y − x)
N N −1
Exp
1 − exp(−
)
exp −
.
f(k,N ) (y; x, β) =
β k−1
β
β

(29)

(30)

The first order statistic density is then given by letting k = 1 in (30) that results in another exponential distribution,
Exp
Exp
f(1,N
(y; x, β̄),
) (y; x, β) = f

where β̄ =

(31)

β
N.

Hence, the MSE of the minimum order statistics estimator is given by

 2β 2
MSE x̂Exp
(y
)
= 2.
min 1:N
N
In order to find the MVU estimator, we re-write the PDF as
#
"


N
1
1X
N
f (y1:N ; x, β) = N exp −
yk exp − x × σ(y(1) − x).
β
β
β

(32)

(33)

k=1

5.1.1

Known hyper parameter β

In case of the known hyper parameter β, the Neyman-Fisher factorization of PDF (33) gives
T (y1:N ) = y(1)
"
#
N
1X
1
h(y1:N ) = N exp −
yk .
β
β

(34a)
(34b)

k=1

The MVU estimator can then be obtained from a transformation of the minimum order statistic that makes it an
unbiased estimator. Finally, in case of exponential noise with known hyper parameter of the distribution, the MVU
estimator and its MSE are given by
β
x̂Exp
(35a)
M V U (y1:N , β) = y(1) −
N


2
β
MSE x̂Exp
.
(35b)
M V U (y1:N , β) =
N2
8

A PREPRINT - N OVEMBER 26, 2019

5.1.2

Unknown hyper parameter β

If the hyper parameter β is unknown, the factorization gives

 

y(1)
T (y )
T (y1:N ) = PN
= 1 1:N .
T2 (y1:N )
k=1 yk

(36a)

Noting that sum of exponential random variables results in a Gamma distribution, we have T2 (y1:N ) ∼ Γ(N, β).
Hence,


β
x+ N
.
E(T (y1:N )) = 
(36b)
N (x + β)
Following the same line of reasoning as in Section 4.1.2, the unbiased estimator is given by the transformation

 1
1
N −1 N T1 (y1:N ) − N T2 (y1:N )
,
g(T (y1:N )) = 
(37a)
1
(T
(y
)
−
N
T
(y
))
2
1:N
1
1:N
N −1
that gives
E (g(T (y1:N ))) =

 
x
.
β

(37b)

Finally, the MVU estimator when the hyper parameter β is unknown, is given by
N

x̂Exp
MVU (y1:N )

X
1
N
y(1) −
yk
=
N −1
N (N − 1)
k=1

N
1
=
y(1) −
ȳ,
(38a)
N −1
N −1
where ȳ is the sample mean. Assuming that N is large mink yk and ȳ are independent and the MSE of the estimator,
asymptotically, is given by


β 2 (N + 1)
.
(38b)
MSE x̂Exp
(y
)
=
1:N
MVU
N (N − 1)2
5.2

Rayleigh distribution

One generalization of the exponential distribution is obtained by parameterizing in terms of both a scale parameter β
and a shape parameter α. Rayleigh distribution is a special case obtained by setting α = 2
(
(yk −x)2
yk −x
Rayleigh
β 2 exp(− 2β 2 ) yk > x,
f
(yk ; x, β) =
(39a)
0
yk ≤ x.
and the CDF, for yk > x is given by
F Rayleigh (yk ; x, β) = 1 − exp(−

(yk − x)2
).
2β 2

(39b)

Hence, the BLUE estimator (8), becomes
x̂Rayleigh
BLUE (y1:N , β) =

N
1 X
yk −
N
k=1

r

π
β,
2


 (4 − π)β 2
MSE x̂Rayleigh
.
BLUE (y1:N , β) =
2N
The marginal density of the k:th order statistic is given by



k−1


2
2
 Ny N − 1
)
exp − (N −k+1)(y−x)
y > x,
1 − exp(− (y−x)
Rayleigh
2
2
2
β
2β
2β
k−1
f(k,N ) (y; x, β) =

0
y ≤ x.

9

(40a)
(40b)

(41)

A PREPRINT - N OVEMBER 26, 2019

Hence, the minimum order statistics density also is Rayleigh distributed
Rayleigh
(y; x, β) = f Rayleigh (y; x, β̄),
f(1,N
)

where β̄ =

√β .
N

(42)

The MSE of the minimum order statistics is given by
 2β 2

(y1:N ) =
MSE x̂Rayleigh
.
min
N

The joint PDF of N independent observations y1:N is given by
"N
#
QN
X (yk − x)2
(y
−
x)
k
exp
−
σ(y(1) − x).
f (y1:N ; x, β) = k=1 2N
β
2β 2

(43)

(44a)

k=1

Noting that
N
X

(yk − x)2 =

k=1

N
X

(yk )2 − 2x

k+1

N
X

yk + N x2 ,

(44b)

k=1

the PDF becomes
"

#
N
−1 X 2
yk
f (y1:N ; x, β) = β
(yk − x) exp
2β 2
k=1
k=1
#
"


N
N x2
x X
× exp − 2 exp 2
yk σ(y(1) − x).
2β
β
−2N

N
Y

(44c)

k=1

5.2.1

Known hyper parameter β

Since (44c) cannot be factorized in the form of f (y1:N ; x, β) = g(T (y1:N ), x)h(y1:N ), the RBLS theorem cannot
be used. Hence, even if an MVU estimator exists for this problem, we may not be able to find it. Thus, in case of
Rayleigh-distributed measurement noise, we propose unbiased estimators based on order statistics.
If the hyper parameter of the distribution is known, the unbiased order statistic based estimator x̂Rayleigh (y1:N , β) is
then given by,
√
πβ
x̂Rayleigh (y1:N , β) = y(1) − √
,
(45a)
2N
 (4 − π)β 2
MSE x̂Rayleigh (y1:N , β) =
.
(45b)
2N
which has the same variance as the BLUE estimator.
5.2.2

Unknown hyper parameter β

In case of unknown hyper parameters, as for the known case, no factorization that enables us to use the RBLS theorem
can be found. In this case, we propose the following unbiased estimator
√
N
X
N
1
Rayleigh
√
x̂
(y1:N ) = √
y(1) −
yk
N −1
N ( N − 1) k=1
√
1
=√
( N y(1) − ȳ).
(46)
N −1
Asymptotically, for large N , the sample mean and minimum order statistic are independent and the estimator MSE is
given by
 (1 + N )(4 − π)β 2
√
MSE x̂Rayleigh (y1:N ) =
.
2N ( N − 1)2
10

(47)

A PREPRINT - N OVEMBER 26, 2019

5.3

Weibull distribution

Weibull distribution is a generalization of the Rayleigh, distribution that is parameterized by two parameters–scale
parameter β and shape parameter α > 0. In fact Weibull distribution is obtained by relaxing the assumption α = 2 in
the Rayleigh distribution and its density function is given by
( 
α−1
α yk −x
exp(−( ykβ−x )α ) yk > x,
Weibull
β
f
(yk ; x, β, α) = β
(48a)
0
yk ≤ x.
and the CDF, for yk ≥ x is given by
F Weibull (yk ; x, β, α) = 1 − exp(−(

yk − x α
) ).
β

(48b)

The BLUE estimator, in case of Weibull-distributed measurement noises is given by
x̂Weibull
BLUE (y1:N , β, α)

N
1 X
1
=
yk − βΓ(1 + ),
N
α

MSE(x̂Weibull
BLUE (y1:N , β, α))

k=1

"

2 #
β2
α+2
α+1
=
Γ(
) − Γ(
)
.
N
α
α
(48c)

The marginal density of the k:th order statistic is given by





α k−1
y−x α
N α N − 1 y − x α−1 
−( y−x
Weibull
β )
(
)
exp
−(N
−
k
+
1)(
)
. (49)
f(k,N
(y;
x,
β,
α)
=
1
−
exp
)
k−1
β
β
β
Hence, the first order statistic density in case of ek ∼ Weibull(β, α), is another Weibull distribution,
Weibull
Weibull
f(1,N
(y; x, β̄, α),
) (y; x, β, α) = f

where β̄ =

(50)

√

−α

N β. This gives the MSE of the minimum order statistic estimator as
−2
α+2
MSE(x̂Weibull
(y1:N )) = β 2 N α Γ(
)
min
α
Given N independent observations, the joint density is given by
α−1
N 
N
X
α N Y yk − x
yk − x α
Weibull
f
(y1:N ; x, β, α) = ( )
(
) )σ(y(1) − x)
exp(−
β
β
β
k=1

(51)

(52)

k=1

Since (52) cannot be factorized using Neyman-Fisher factorization, RBLS is not applicable. Additionally, in this case,
it is not possible to find an unbiased estimator when the hyper parameters α and β are unknown. In case of known
hyper parameters, the unbiased minimum order statistic estimator, however, can be computed. The unbiased estimator
based on minimum order statistic is given by,
1
1
x̂Weibull (y1:N , β, α) = y(1) − βN − α Γ(1 + ),
α
"

2 #
−2
α+2
α+1
Weibull
2
α
MSE(x̂
(y1:N , β, α)) = β N
Γ(
) − Γ(
)
.
(53)
α
α
An order-statistics-based unbiased estimator with unknown hyper parameters of the distribution could not be obtained.

6

Other Distributions

In this section, we further study the location estimation problem for two other noise distributions. In the rest, the
Pareto distribution with positive support is first studied followed by the mixture of uniform and normal distribution.
6.1

Pareto distribution

Let the scale parameter β (necessarily positive) denote the minimum possible value of yk , and α > 0 denote the shape
parameter. The Pareto Type I distribution is characterized by β and α

αβ α (yk − x)−(α+1) yk ≥ x + β,
f Pareto (yk ; x, β, α) =
(54a)
0
yk < x + β.
11

A PREPRINT - N OVEMBER 26, 2019

and the CDF is given by
F

Pareto


(yk ; x, β, α) = 1 −

β
y−x

α
.

(54b)

For the BLUE we get,
x̂Pareto
BLUE (y1:N , β, α) =

N
1 X
αβ
yk −
,
N
α−1

α>1

(55a)

αβ 2
,
N (α − 1)2 (α − 2)

α>2

(55b)

k=1

MSE(x̂Pareto
BLUE (y1:N , β, α)) =

The RBLS theorem cannot be used in case Pareto-distributed noises. We provide an unbiased estimator using minimum
order statistics and its variance. The marginal density of the k:th order statistic, for y ≥ x + β is given by


k−1 

β α
β
Pareto
α
−(α+1) N − 1
f(k,N
(y;
x,
β,
α)
=
N
αβ
(y
−
x)
1
−
(
)
α(N − k).
(56)
)
k−1
y−x
y−x
The minimum order statistic has the same form of distribution
Pareto
Pareto
f(1,N
(y; β, ᾱ),
) (y; β, α) = f

(57)

where ᾱ = N α. The MSE of the minimum order statistic estimator is
MSE(x̂Pareto
min (y1:N )) =

N αβ 2
Nα − 2

(58)

The unbiased estimator is thus given by,
N αβ
,
Nα − 1
N αβ 2
MSE(x̂Pareto (y1:N , β, α)) =
,
(N α − 1)2 (N α − 2)
x̂Pareto (y1:N , β, α) = y(1) −

Nα > 1

(59a)

Nα > 2

(59b)

No unbiased estimator for unknown hyper parameter case could be found for Pareto distribution.

7

Mixture of Normal and Uniform Noise Distribution

Suppose the error is distributed as
ek ∼ αN (0, σ 2 ) + (1 − α)U[0, β],
where α is the mixing probability of the mixture distribution. Define f U ,N (yk ) as the probability density function of
the considered mixture distribution given by
f U ,N (yk ; x, α, σ 2 , β) =



(yk − x)2
1−a

√ a
exp −
+
2σi2
β
2πσ 2 h
2

 √ a exp − (yk −x)
2σ 2
2πσ 2

0 ≤ yk − x ≤ β

(60)

Otherwise.

The BLUE, in case of the mixture of normal and uniform measurement noises is given by
N
1 X
β(1 − α)
yk −
,
N
2

(61a)

β 2 (1 + (2 − 3α)α) + 12ασ 2
.
12N

(61b)

,N
2
x̂U
BLUE (y1:N , α, β, σ ) =

k=1



,N
2
MSE x̂U
BLUE (y1:N , α, β, σ ) =

12

A PREPRINT - N OVEMBER 26, 2019

Noting that at yk − x = 0 contributions of the uniform distribution and the mean (mode) of the normal distribution are
added together, (60) is maximized at this point. The order statistics PDF for 0 ≤ y − x ≤ β is given by
U ,N
2
f(k,N
) (y; α, β, σ , x, k) =
!
2


α exp(− (y−x)
)
1
−
α
2
N −1
2σ
√
N
+
k−1
β
2πσ 2
 k−1


(1 − α)(y − x) α
y−x
)
×
+ (1 + Erf √
β
2
2σ 2
 N −k


(α − 1)(y − x) α
y−x
,
)
× 1+
− (1 + Erf √
β
2
2σ

(62)

R·
2
where Erf(·) = √2π 0 e−t dt is the error function. In order to find the best order statistic estimator, we maximize the
likelihood function `(k | y = x, a, β, σ 2 )


α
N −1
`(k | y = x, α, β, σ 2 ) = N
2(2 − α)−k (1 − )N αk−1
k−1
2


1−α
α
×
+√
.
(63a)
β
2πσ


√α
Noting that 1−α
β + 2πσ is always positive and independent of k, we extract it from the likelihood function. Simplifying (63a) by means of manipulating the terms, we get
2(2 − α)−k = 21−k (1 −

α −k
) ,
2

(63b)

α
α
αk−1 = (2 )k−1 = 2k−1 ( )k−1 .
2
2
the likelihood function to be maximized can be re-written as


α
N − 1 α k−1
2
`(k | y = x, α, , σ ) ∝
( )
(1 − )N −k .
k−1
2
2

(63c)

(63d)

In order to find the maximum likelihood estimate k̂ = arg maxk `(k | y − x = 0), we note that (63d) is a binomial
distribution with probability of success α2 . Hence, the maximum of the function is given at the mode of the distribution,
Nα
Nα
c + 1 or d
e.
(64)
2
2
This gives the best order statistic estimator for the case when noise is a mixture of normal and uniform distribution as
k̂ = b

x̂U ,N (y1:N , α) = y(k̂) .

(65)

e
Table 2: Bias and MSE of minimum order statistics estimators x̂pmin
.
distribution
bias
MSE
β
2β 2
U[0, β]
N +1
(N +1)(N +2)

Exp(β)
Rayleigh(β)

β
N
√
√ πβ
2N
1

2β 2
N2
2β 2
N
2

Weibull(β, α)

βN − α Γ(1 + α1 )

β 2 N − α Γ(1 + α2 )

Pareto(β, α)

N αβ
N α−1

N αβ 2
N α−2

13

A PREPRINT - N OVEMBER 26, 2019

Table 3: Estimators and their MSEs derived for multiple noise distributions.
noise
estimator
MSE
distribution
P
N
β
β2
1
x̂U
BLUE (y1:N , an ) = N
k=1 yk − 2
12N
U[0, β]

1
x̂U
MVU (y1:N , an ) = 2 (y(1) + y(N ) ) −

x̂U
MVU (y1:N ) =

N
N −1 y(1)

x̂Exp
BLUE (y1:N , β) =
Exp(β)

1
N

x̂Rayleigh
BLUE (y1:N , β) =

1
N

αN (0, σ 2 ) + (1 − α)U[0, β]

β2
N
β2
N2

β
N
PN

1
N

P

√
√ πβ
2N

(N +1)β 2
N (N −1)2
(4−π)β 2
2N
(4−π)β 2
2N

PN
k=1 yk
√
N ( N −1)

(1+N )(4−π)β 2
√
2N ( N −1)2

yk − βΓ(1 + α1 )

h
2 i
β 2 N −1 Γ(1 + α2 ) − Γ(1 + α1 )
h
2 i
2
β 2 N − α Γ(1 + α2 ) − Γ(1 + α1 )

−

1

x̂Weibull (y1:N , β, α) = y(1) − βN − α Γ(1 + α1 )
P
αβ
1
x̂Pareto
yk − α−1
,α>1
BLUE (y1:N , β, α) = N
N αβ
N α−1 ,

x̂Pareto (y1:N , β, α) = y(1) −
Nα > 1
PN
β(1−α)
,N
1
x̂U
BLUE (y1:N , α, β) = N
k=1 yk −
2
x̂U ,N (y1:N , α, β) = y(b N α c+1)
2

8

y

k=1 k
− N (N
−1)
pπ
PN
k=1 yk −
2β

√
√ N y(1)
N −1

x̂Weibull
BLUE (y1:N , β, α) =

Pareto(β, α)

yk − β

x̂Rayleigh (y1:N , β) = y(1) −
x̂Rayleigh (y1:N ) =

Weibull(β, α)

N β2
(N +2)(N 2 −1)

k=1

N
N −1 y(1)

β2
2N (N +3)+4

1
N −1 y(N )

PN

x̂Exp
MVU (y1:N , β) = y(1) −
x̂Exp
MVU (y1:N ) =

Rayleigh(β)

−

β
2

αβ 2
N (α−1)2 (α−2) ,

α>2

N αβ 2
(N α−1)2 (N α−2)
β 2 (1+(2−3α)α)+12ασ 2
12N

Unknown

Performance Evaluation

The estimators (both unbiased and the ones without bias compensation) derived in sections 5–6 for different noise
distributions together with their MSE are summarized in Tables 3 and 2. The biased minimum order statistics based
estimators and their MSE are also The estimators derived for each noise distribution are compared against each other
as a function of the sample size N ∈ [2, . . . , 2000]. Additionally, in order to verify the analytical derivations of the
estimator variances, they are compared against the numerical variances obtained from M = 5000 Monte Carlo runs.
8.1

Simulation Setup

For each sample size, N noisy measurements of the unknown parameter x are generated. The hyper parameters of the
noise distributions are randomly selected in each repetition. In order to have a fair comparison, the hyper parameters
are randomly drawn such that the error densities are mostly in the same range for all scenarios. The noise realizations
are generated from the six considered distributions with the following hyper parameters
•
•
•
•
•
•

Uniform noise: β ∼ U[6, 50]
Exponential noise: β ∼ U[5, 14]
Rayleigh noise: β ∼ U[5, 12]
Weibull noise:β = 1, α ∼ U[5, 10]
Pareto noise: β = 6, α ∼ U[2.1, 2.5]
Mixture noise: σ ∼ U[1, 9], β ∼ U[1, 50]

The empirical CDF of the error values used in the simulations are presented in Figure 2. The support of the noise
values, as can be read from the figure, is em ∈ [0, 60] unit.
(m)

Let x̂N denote the estimated value of the unknown parameter x in the m:th repetition obtained from a sample of
size N . For each noise distribution, the estimators’ performances are evaluated in terms of the obtained MSEs. The
14

A PREPRINT - N OVEMBER 26, 2019

BLUE, known hyperparameter
proposed estimator, known hyperparameter
proposed estimator, unknown hyperparameter
minimum order statistics estimator

100

102

80

10

1

100

CDF

60
10-1

40

20

0

10-2

Uniform
Exponential
Rayleigh
Pareto
Weibull

0

10

20

30

40

50

10-3
10-4

60

0

Measurement Error

Figure 2: Empirical CDF of measurement erros computed from noise realizations used in the simulations.

500

1000
sample size N

1500

2000

Figure 3: Analytical (marked with solid lines) and numerical (marked with crosses) MSE for uniform noise
distribution as a function of the sample size N .

theoretical MSE of each estimator, as defined in Table 3 and Table 2, is compared against the numerical MSE obtained
in simulations.
PM
(m)
1
We let E[x̂N ] = M
m=1 x̂N and define
b̂N = E[x̂N ] − x
2
σ̂N
=

1
M

M
X

(m)

(x̂N − E[x̂N ])2 .

(66a)
(66b)

m=1

The numerical MSE for each sample size N is then computed by
2
2
ˆ
MSE(x̂
N ) = σ̂N + b̂N .

8.2

(66c)

Simulation Results

Figure 3 presents the performance of the four estimators when the noise is uniformly distributed. The solid lines
correspond to the theoretical MSEs and the crosses are the numerical MSEs obtained from M = 5000 repetitions.
Both MVU estimators, with and without any knowledge of the hyper parameters of the underlying noise, result in
noticeably less MSE compared to the BLUE estimator. The minimum order statistics estimator also outperforms
BLUE when measurements are corrupted with additive, uniformly distributed, noise. It can be further observed that if
the hyper parameter β is unknown, the MSE of the proposed estimator is negligibly larger than the case with known
β.
For the exponential noise distribution, as shown in Figure 4a, there is still a non-negligible difference between BLUE
and the other three estimators in terms of estimators’ MSE. However, the two MVU estimators, specially for large
values of N , behave similarly. In order to verify their performance for smaller sample sizes, Figure 4b illustrates the
variances of all estimators for N ≤ 20. At the beginning, N ∈ [2, 4] the estimator with unknown hyperparameter
has the largest MSE. However, for larger sample sizes, the two MVU estimators are almost equal and both have less
MSE than the BLUE estimator. As in case of uniformly distributed measurement noise, the minimum order statistics
estimator outperforms BLUE specially for large sample sizes.
In case of Rayleigh noise distribution, as given in Table 3, the minimum order statistics estimator has the largest
MSE while the BLUE and the proposed unbiased estimator with known hyper parameter, result in similar estimation
15

A PREPRINT - N OVEMBER 26, 2019

BLUE, known hyperparameter
proposed estimator, known hyperparameter
proposed estimator, unknown hyperparameter
minimum order statistics estimator
104

103

102

102

100

101

10-2

100

10-4

10-1

10-6

0

500

1000
sample size N

1500

10-2

2000

2

4

6

8

10
12
sample size N

14

16

18

20

(b) Numerical and analytical MSE for N ≤ 20.

(a) Numerical and analytical MSE.

Figure 4: Analytical (marked with solid lines) and numerical (marked with crosses) MSE for exponential noise distribution as a function of the sample size N . The blue and the proposed estimator with known hyperparameters have
equal variances, hence the blue line is invisible in these plots. Since the MVU estimators have similar results for large
sample sizes, the MSE of the four estimators for smaller sample sizes are presented separately.
BLUE, known hyperparameter
proposed estimator, known hyperparameter
proposed estimator, unknown hyperparameter
minimum order statistics estimator
103

103

102

102

101

101

100

100

10-1

10-1

10-2

0

500

1000
sample size N

1500

10-2

2000

2

4

6

8

10
12
sample size N

14

16

18

20

(b) Numerical and analytical MSE for N ≤ 20.

(a) Numerical and analytical MSE.

Figure 5: Analytical (marked with solid lines) and numerical (marked with crosses) MSE for Rayleigh noise distribution as a function of the sample size N . Since the BLUE and the proposed estimators have similar results for large
sample sizes, the MSE of the four estimators for smaller sample sizes are presented separately.

variance. This can be verified also in the simulation results presented in Figure 5a. For large sample sizes, N > 20,
these two estimators and the proposed estimator with unknown hyperparameter have similar values. However, for the
smaller sample sizes, as illustrated in Figure 5b, the BLUE (and order statistic with known hyper parameter) estimator
has smaller variance compared to the case with unknown hyper parameter. The minimum order statistics estimator
results in larger MSE compared to the other three estimators in case of Rayleigh noise distribution.
16

A PREPRINT - N OVEMBER 26, 2019

BLUE, known hyperparameter
proposed estimator, known hyperparameter
minimum order statistics estimator

BLUE, known hyperparameter
proposed estimator, known hyperparameter
minimum order statistics estimator

102

102

100

100

10

-2

10

10-4

10-6

-2

10-4

0

500

1000
sample size N

1500

10-6

2000

0

500

(a) Weibull noise distribution.

1000
sample size N

1500

2000

(b) Pareto noise distribution.

Figure 6: Analytical (marked with solid lines) and numerical (marked with crosses) MSE for Weibull and Pareto noise
distributions as a function of the sample size N .

As Table 3 suggests, for Pareto and Weibull noise distributions, we only derived BLUE and an unbiased order statistics
based estimators when the two hyperparameters of the distributions are known. For both noise distributions, the MSE
of the two unbiased estimators as well as the MSE of the minimum order statistics estimator are compared and the results are presented in Figure. 6. In both cases, the proposed estimators outperform the BLUE in terms of variance. The
minimum order statistics estimator results in a lower MSE than the BLUE for Weibull noise distributions. However,
in case of Pareto noise, the BLUE has a better performance compared to the minimum order statistics estimator.
In case of mixture noise distribution, we consider three different scenarios based on the mixing probabilities; two
extreme cases with dominant contribution from uniform noise, α = 0.01, and dominant contribution from normal
noise, α = 0.99, and the case with α = 0.5. Fig. 7a illustrates the histogram of the noise realizations of the considered

100

80

CDF

60

40

20

0
-20

(a) Histogram of noise realizations and the fitted PDF.

0

20

40

60

(b) Empirical CDF.

Figure 7: Noise realizations, PDF, and empirical CDF of the mixture noise distribution ek ∼ αN (0, 82 ) + (1 −
α)U(0, 60) for three different values of α.
17

A PREPRINT - N OVEMBER 26, 2019

BLUE, known hyperparameter
proposed estimator, known hyperparameter
103

103
101
102

102
100
10

101

1

100
10-1
100

10-1

10-1

0

500

1000
sample size N

(a) α = 0.01.

1500

2000

10-2

0

500

1000
sample size N

(b) α = 0.99.

1500

2000

10-2

0

500

1000
sample size N

1500

2000

(c) α = 0.5.

Figure 8: Analytical (marked with solid lines) and numerical (marked with crosses) MSE for three different values of
mixing probability α, when ek ∼ αN (0, 82 ) + (1 − α)U(0, 60).
mixture noise distributions ek ∼ αN (0, 82 ) + (1 − α)U(0, 60) and the fitted densities. The empirical CDFs of the
errors for the three cases are presented in Figure 7b.
In order to estimate the unknown parameter x, in each Monte Carlo run, we sort the measurements and then find the
(b N2α c+1):th component. Figure 8 presents the estimation MSE for the three different scenarios with different mixing
probabilities. As the results indicate, when the main contribution of the noise is from uniform distribution, α = 0.01,
BLUE outperforms the proposed estimator. In this case, a periodic behavior for the MSE can be observed. The jumps
in the MSEe occur exactly at pints where b N2α c + 1 switches from the k:th measurement to the k + 1:th measurement.
For instance, for N ∈ [1, 199], b N2α c = 0, hence x̂ = y(1) . However, at N = 200, b N2α c = 1, resulting in x̂ = y(2) .
The proposed estimator and the BLUE result in similar estimation MSE for α = 0.99, as shown in Figure 8b, in which
the normal component is the dominant source of error. However, the most interesting results are obtained when both
distributions have equal contributions in the measurement noise, i.e α = 0.5. In this case, as Figure 8c suggests, the
proposed estimator outperforms the BLUE.

9

Conclusions

In this work, the location estimation problem was studied in which an unknown parameter was estimated from observations under additive noise. Multiple noise distributions were considered and, in some cases, MVU estimators were
proposed. In other cases an unbiased estimator based on minimum order statistic was derived. Furthermore, if applicable, MVU and minimum order statistic estimators without any knowledge of the hyper parameters of the underlying
noise distributions were provided. The results of all the estimators were compared with BLUE in terms of variance
for various measurement sample sizes. The results indicate better performance of the proposed estimators compared
to BLUE. Additionally, the location estimation problem under mixture of normal and uniform noise distribution was
studied and the numerical MSE of the proposed estimator were evaluated. The simulation results indicate that for
the extreme cases where either of the two components, Gaussian or uniform, are dominant, the proposed estimator
cannot beat the BLUE. However, when the mixing probability is not in the extreme region, e.g larger than 1 percent,
the proposed estimator has a noticeably less MSE compared to the BLUE.

References
S. A. Kassam and H. V. Poor. Robust techniques for signal processing: A survey. Proceedings of the IEEE, 73(3):
433–481, March 1985.
S. M. Kay. Fundamentals of Statistical Signal Processing: Estimation Theory. Prentice-Hall, Inc., Upper Saddle
River, NJ, USA, 1993.
E. L. Lehmann and G. Casella. Theory of Point Estimation. Springer-Verlag New York, 1998.
M. Kok, J. D. Hol, and T. B. Schön. Indoor positioning using ultrawideband and inertial measurements. IEEE
Transactions on Vehicular Technology, 64(4):1293–1303, April 2015.
18

A PREPRINT - N OVEMBER 26, 2019

B. Chen, C. Yang, F. Liao, and J. Liao. Mobile location estimator in a rough wireless environment using extended
kalman-based IMM and data fusion. IEEE Transactions on Vehicular Technology, 58(3):1157–1169, March 2009.
F. Gustafsson and F. Gunnarsson. Mobile positioning using wireless networks: possibilities and fundamental limitations based on available wireless network measurements. IEEE Signal Processing Magazine, 22(4):41–53, July
2005.
M. Eling. Fitting insurance claims to skewed distributions: Are the skew-normal and skew-student good models?
Insurance: Mathematics and Economics, 51(2):239–248, 2012.
J. Medbo, I. Siomina, A. Kangas, and J. Furuskog. Propagation channel impact on LTE positioning accuracy: A study
based on real measurements of observed time difference of arrival. In Proc. of 20th IEEE International Symposium
on Personal, Indoor and Mobile Radio Communications, pages 2213–2217, Westin Toyko, Toyko, Japan, September
2009.
F. Yin, C. Fritsche, F. Gustafsson, and A. M. Zoubir. TOA-based robust wireless geolocation and cramér-rao lower
bound analysis in harsh LOS/NLOS environments. IEEE Transactions on Signal Processing, 61(9):2243–2255,
May 2013.
S. M. Stigler. Simon Newcomb, Percy Daniell, and the history of robust estimation 1885-1920. Journal of the
American Statistical Association, 68(344):872–879, 1973.
S. A. Kassam. Signal Detection in Non-Gaussian Noise. Springer-Verlag New York, 1988.
C. Stewart. Robust parameter estimation in computer vision. SIAM Review, 41(3):513–537, 1999.
G. R. Arce. Nonlinear Signal Processing: A Statistical Approach. Hoboken, NJ: Wiley, 2004.
A. M. Zoubir, V. Koivunen, Y. Chakhchoukh, and M. Muma. Robust estimation in signal processing: A tutorial-style
treatment of fundamental concepts. IEEE Signal Processing Magazine, 29(4):61–80, July 2012.
E. Eskin. Anomaly detection over noisy data using learned probability distributions. In In Proc. of the International
Conference on Machine Learning, pages 255–262, Stanford, CA, USA, June 2000.
S. Chawla, D. Hand, and V. Dhar. Outlier detection special issue. Data Mining and Knowledge Discovery, 20(2):
189–190, March 2010.
V. J. Hodge and J. Austin. A survey of outlier detection methodologies. Artificial Intelligence Review, 22(2):85–126,
October 2004.
C. Fritsche, U. Hammes, A. Klein, and A. M. Zoubir. Robust mobile terminal tracking in NLOS environments
using interacting multiple model algorithm. In Proc. of International Conference on Acoustics, Speech and Signal
Processing (ICASSP), pages 3049–3052, Taipei, Taiwan, April 2009.
P.J. Huber and E.M. Ronchetti. Robust Statistics. Hoboken, NJ: Wiley„ 2009.
R. A. Maronna, R. D. Martin, and V. J. Yohai. Robust Statistics: Theory and Methods. Hoboken, NJ: Wiley„ 2006.
E. L. Lehmann and H. Scheffé. Completeness, similar regions, and unbiased estimation: Part I. The Indian Journal of
Statistics, 10(4):305–340, 1950.
E. L. Lehmann and H. Scheffé. Completeness, similar regions, and unbiased estimation: Part II. The Indian Journal
of Statistics, 15(3):219–236, July 1955.
R. A. Fisher and M. A. Phil. On the mathematical foundations of theoretical statistics. Philosophical Transactions of
the Royal Society of London A: Mathematical, Physical and Engineering Sciences, 222(594-604):309–368, January
1922.
P. R. Halmos and L. J. Savage. Application of the Radon-Nikodym theorem to the theory of sufficient statistics. The
Annals of Mathematical Statistics, 20(2):225–241, June 1949.
J. B. McDonald and Y. J. Xu. A generalization of the beta distribution with applications. Journal of Econometrics, 66
(1):133–152, March 1995.
H. A. David and H. N. Nagaraja. Order Statistics. John Wiley & Sons, 2004.
19


