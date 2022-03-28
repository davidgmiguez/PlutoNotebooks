### A Pluto.jl notebook ###
# v0.17.7

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ d096a6be-65a1-428d-9bfb-da7fe89f4c19
using Plots 

# ╔═╡ f0413f68-1eb4-4a3d-bdcf-61e8aa96c0e7
using DifferentialEquations

# ╔═╡ 6d9759c3-63d9-4166-8884-cfaa99ee33c4
using HypertextLiteral

# ╔═╡ 54dd5337-c76f-4a23-8bec-a0cd57b433dd
using PlutoUI

# ╔═╡ 7ba042d6-107d-4188-a877-c7f7d7b35c46
md" # Continuous Models 

## 4 Unconstrained growth Continuous

Stem Cell Proliferation Growth In general stem cells undergo a phase of proliferation that expands the population, then they start to differentiate to produce the different types of progeny. Let’s calculate how a population of stem cells grow in time. Lets assume we have an initial number of cells `N` in a population 

After a given time we have a number of new stem cells `∆N` (number of newly produced cells) during a given time interval, `∆t`, is proportional to the initial  number of cells `N`. (If a population of 20,000 cells produces 1200 new cells in 1 h, a 4-fold bigger population of 80,000 cells of the same type of microorganism will produce 4 times as many, viz. 4800, new cells in 1 h): 

```math
\frac{\Delta N}{\Delta t} \sim N\tag{1}
```

Because the newly produced cells always add to the population, i.e. N increases steadily, we have to regard time intervals (and accordingly numbers of newly formed cells) that are as small as possible. Mathematically we thus have to deal with infinitesimal increases or differentials: 

```math
\frac{\mathrm{d} N}{\mathrm{d} t} \sim N\tag{2}
```

To get from proportionality  to an equation, a proportionality factor, μ, is introduced. This is the specific growth rate or often termed simply growth rate. Sorting of variables yields: 

```math
\frac{\mathrm{d} N}{N}  = \mu \mathrm{d} t\tag{3}
```

If we integrate this differential equation:

```math
\int \frac{\mathrm{d} N}{N}  = \int \mu \mathrm{d} t \tag{4}
```

 we obtain:
 
 ```math
\ln N = \mu t + C \tag{5}
```
 
 The undefined integration constant `C` can be fixed if we define the initial conditions: $N(t = 0) = N_0$

 ```math
\ln N_0 = \mu 0 + C =C \tag{6}
```
 
 so the equation becomes:
 
  ```math
\ln N = \mu t + \ln N_0 \tag{7}
```
  
  and after applying the properties of logatithms:
  
```math
\begin{align*}
\ln N - \ln N_0 &= \mu t  \tag{8}\\
\ln \frac{N}{N_0}&= \mu t  \tag{9}\\
\frac{N}{N_0} &= e^{\mu t}  \tag{10}\\
N &= N_0 e^{\mu t}  \tag{11}\\
\end{align*}
```
  
   If the time that has passed is exactly the average length of the cell cycle of the stem cell population, we can write that in $t = T$ the population has doubled, so $N_0$ has increased to $2 N_0$:

  
```math
\begin{align*}
2 N &= N_0 e^{\mu T}  \\
\frac{2 N_0}{N_0 }& =2= e^{\mu T}  \tag{13}\\
 \ln 2&= \mu T  \tag{14}\\
\end{align*}
```
   
   So, the proportionality constant is related to the cell cycle length as:
```math
\begin{align*}
 \mu = \frac{\ln 2}{T}  \tag{15}
\end{align*}
```


"

# ╔═╡ 4b516468-55ff-492b-a026-61154c8062a0
function unconstrained_growth_continuous(N₀,t,T)
	μ=log(2)/T
    N₀ .* exp.(μ .* t)
end

# ╔═╡ 51945b7a-f1ac-423b-affd-f1676987de7b
begin
	T_slide = @bind TT html"<input type=range min=5 max=15 step=1>"
	
	md"""
	**Set the Cell cycle length**
	
	value of T: $(T_slide)
	
	"""
end

# ╔═╡ 1c24a2f1-2b15-48e6-b0d8-c3278103036d

	plot(collect(0:1:100),unconstrained_growth_continuous(20,collect(0:1:100),TT),label="T= $TT",seriestype=:line,xlabel=("Time"),ylabel=("Number of cells"),ylims = (0,400))


# ╔═╡ 85bc871e-6067-4a1d-bbb9-e3a9b0dfbd85
md"
# Logistic Growth 

We can introduce in the original equation a term that accounts for the limitation of resources 

```math
\begin{align*}
1 -\frac{N}{K} \tag{16}
\end{align*}
```

This factor is close to 1 (i.e., has no effect) when N is much smaller than K, and which is close to 0 when N is close to K. The resulting differential equation is called is the logistic growth model.


```math
\begin{align*}
\frac{\mathrm{d} N}{\mathrm{d} t}=\mu N(1 -\frac{N}{K}) \tag{17}
\end{align*}
```

separate variables

```math
\begin{align*}
\frac{\mathrm{d} N}{N(1 -\frac{N}{K})}=\mu \mathrm{d} t   \tag{18}
\end{align*}
```

Our next goal would be to integrate both sides of this equation, but the form of the right hand side doesn't look elementary and will require a partial fractions expansion. That is, we wish to write 


```math
\begin{align*}
\frac{1}{N(1 -\frac{N}{K})} = \frac{A}{N}+\frac{B}{1 -\frac{N}{K}} \tag{19}
\end{align*}
```

where $A$ and $B$ are unknown constants. If we multiply on the left and right hand sides by  $N \left( 1- \frac{N}{K} \right)$ (which is equivalent to putting the right hand side over a common denominator) we arrive at the equation 

```math
\begin{align*}
1 = A \ \left( 1 -\frac{N}{K}\right) + B \cdot N  = A + N (B - \frac{A}{K}) \tag{20}
\end{align*}
```

Since there is no term with $N$ on the left hand side, we see that 

```math
\begin{align*}
B - \frac{A}{K} = 0 \quad \mbox{ or } \quad B = \frac{A}{K} \tag{21}
\end{align*}
```

If we set $B = \frac{A}{K}$ then we are left with $A=1$, and thus the partial fraction decomposition is 

```math
\begin{align*}
\frac{1}{N(1 -\frac{N}{K})} = \frac{1}{N}+\frac{\frac{1}{K}}{1 -\frac{N}{K}} \tag{22}
\end{align*}
```

so the integral becomes

```math
\begin{align*}
\frac{\mathrm{d} N}{N}+\frac{\frac{\mathrm{d} N}{K}}{1 -\frac{N}{K}}=\mu \mathrm{d} t   \tag{23}
\end{align*}
```

the first part is simply:

```math
\begin{align*}
\int\frac{dN}{N} = \ln (N) \tag{24},
\end{align*}
```

For the second term, we must use a substitution $u=1-\frac{N}{K}$, which gives a differential  $du = \frac{-1}{K} \ dN$. Thus we may write the second term on the right hand side as:

```math
\begin{align*}
\int \frac{dN/K}{\left( 1-\frac{N}{K} \right)} = \int \frac{-du}{u} = -\ln (u) = -\ln (1-N/K) \tag{25}
\end{align*}
```

Putting all these terms together gives us:

```math
\begin{align*}
\mu t + c = \ln (N)-\ln (1-\frac{N}{K}) = \ln \left[\frac{N}{1-N/K} \right] \tag{26}
\end{align*}
```

Here we have used the property of logarithms to equate the difference of the logs with the log of the quotient. The additional term, $c$, on the left hand side is the free constant of integration, which will be determined by considering initial conditions to the differential equation. Exponentiating both sides of the equation gives 

```math
\begin{align*}
e^{\mu t + c} = \frac{N}{1-\frac{N}{K}} \tag{27},
\end{align*}
```

so

```math
\begin{align*}
 e^{\mu t} e^c = \frac{N}{1-\frac{N}{K}} \tag{28}
\end{align*}
```

```math
\begin{align*}
e^{\mu t} C = \frac{N}{1-\frac{N}{K}} \tag{29}
\end{align*}
```

to find $C$ we use teh inital condition that $N(t=0)=N_0$, and substituting gives 

```math
\begin{align*}
C = C e^0 = \frac{N_0}{1-\frac{N_0}{K}} = \frac{N_0}{1-\frac{N_0}{K}} \frac{K}{K} =\frac{K N_0}{K-N_0} \tag{30}
\end{align*}
```

Solving now for $P$, we first cross-multiply to arrive at 

```math
\begin{align*}
\left(1-\frac{N}{K} \right) C e^{\mu t} = N \tag{31}
\end{align*}
```

and putting all terms including $N$ on one side of the equation, 

```math
\begin{align*}
C e^{ \mu t} = N \left[1 + \frac{C e^{ \mu t}}{K} \right]  \tag{32}
\end{align*}
```

Solving now for $N$, 

 
```math
\begin{align*}
N = \frac{C e^{\mu t}}{1 + \frac{C e^{\mu t}}{K} } =
\frac{\frac{K \cdot N_0}{K-N_0} e^{\mu t}}{1 + \frac{\frac{K\cdot N_0}{K-N_0} e^{\mu t}}{K} }\tag{33}
\end{align*}
```

Simplifying this expression by multiplying numerator and denominator by  $(K-N_0) e^{-\mu t}$ gives 

```math
\begin{align*}
N = \frac{K N_0}{N_0 +(K-N_0) e^{- \mu t}} \tag{34}
\end{align*}
```

if we monitor hwo the value of $\mu$ changes 

```math
\begin{align*}
(K-N_0) e^{- \mu t} &=\frac{K N_0}{N} - N_0\\
(K-N_0) e^{- \mu t} &=\frac{N_0 (K - N)}{N}\\
 e^{- \mu t} &=\frac{N_0 (K - N)}{N (K-N_0)}\\
 - \mu t &= Log(\frac{N_0 (K - N)}{N (K-N_0)})\\
  \mu  &= \frac{1}{t}Log(\frac{N (K-N_0)}{N_0 (K - N)})\\
    \frac{Log 2}{T}  &= \frac{1}{t}Log(\frac{N (K-N_0)}{N_0 (K - N)})\\
    T  &=  \frac{t \cdot Log 2}{Log\frac{N (K-N_0)}{N_0 (K - N)}} \\
\end{align*}
```
"

# ╔═╡ e1601070-d88e-4515-8ad7-af810d98e46d
function logistic_growth_continuous(N₀,K,t,T)
	μ=log(2)/T
    #N₀ .* exp.(μ .* t)
	N₀.* K  ./ (N₀ .+(K .- N₀) .* exp.(-μ .* t))
end

# ╔═╡ c871a80d-113c-47ee-8292-3a857752b094
begin
	TTT_slide = @bind TTT html"<input type=range min=5 max=15 step=1>"
	K_slide = @bind K html"<input type=range min=100 max=400 step=1>"
	
	md"""
	**Set the Cell cycle length and the Carrying Capacity **
	
	value of T: $(TTT_slide) 
	
	value of K: $(K_slide)
	
	"""
end

# ╔═╡ 25caf1e9-788d-41a2-add3-f0cb451420aa
	plot(collect(0:1:100),logistic_growth_continuous(20,K,collect(0:1:100),TTT),label="T= $TTT",seriestype=:line,xlabel=("Time"),ylabel=("Number of cells"),ylims = (0,400))

# ╔═╡ 5a387568-2284-4d4e-85fc-f106db7377c9

question: how long until half of the population  N=K/2
how long until 10% of the population  N=K/2

# ╔═╡ 61fe29b8-065c-430f-a28d-b463780f8b00
md" ## Conclusions

Fitness is a term that describes the survival and reproductive output of an individual in a population relative to other members of the population. In other words, it represnets how well an organism is adapted to its environment. The fitness of an individual animal is a measure of its ability, relative to others, to leave viable offspring.

Fitness can fundamentally be achieved by two different strategies: long life (stability) or fast reproduction (multiplication, replication). These strategies are to some degree dependent: since no organism is immortal, a minimum amount of reproduction is needed to replace the organisms that have died; yet, in order to reproduce, the system must live long enough to reach the degree of development where it is able to reproduce. On the other hand, the two strategies cannot both be maximally pursued: the resources used for fast reproduction cannot be used for developing a system that will live long, and vice-versa. This means that all evolutionary systems are confronted with a development-reproduction trade-off: they must choose whether they invest more resources in the one or in the other.
How much a given system will invest in one strategy at the expense of the other one depends on the selective environment. In biology, this is called r-K selection: in an r-situation, organisms will invest in quick reproduction, in a K-situation they will rather invest in prolonged development and long life. Typical examples of r-species are mice, rabbits, weeds and bacteria, which have a lot of offspring, but a short life expectancy. Examples of organisms undergoing K-selection are tortoises, elephants, people, and sequoia trees: their offspring are few but long-lived. In summary, r-selection is selection for quantity, K-selection for quality of offspring.

| r-organisms          | K-organisms | 
|--------------|:-----:|
| short-lived |  long-lived | 
| small      |  Large |  
| weak      |  strong or well-protected |  
| waste a lot of energy      |  energy efficient |  
| less intelligent, experienced...      |  intelligent, experienced... |
|have large litters| have small litters|
|reproduce at an early age|reproduce at a late age|
|fast maturation|slow maturation|
|little care for offspring|much care for offspring|
|strong sex drive|weak sex drive|
|small size at birth|large size at birth|

"

# ╔═╡ 2ecd106c-faee-469e-af4c-0f3b6c28bc02
md"
| r-organisms          | K-organisms | 
|--------------|:-----:|
| short-lived |  long-lived | 
| small      |  Large |  
| weak      |  strong or well-protected |  
| waste a lot of energy      |  energy efficient |  
| less intelligent, experienced...      |  intelligent, experienced... |
|have large litters| have small litters|
|reproduce at an early age|reproduce at a late age|
|fast maturation|slow maturation|
|little care for offspring|much care for offspring|
|strong sex drive|weak sex drive|
|small size at birth|large size at birth|
"

# ╔═╡ 5220090a-c8ef-4b9b-ae17-6692fc4b19e3


# ╔═╡ Cell order:
# ╟─7ba042d6-107d-4188-a877-c7f7d7b35c46
# ╠═4b516468-55ff-492b-a026-61154c8062a0
# ╠═51945b7a-f1ac-423b-affd-f1676987de7b
# ╠═1c24a2f1-2b15-48e6-b0d8-c3278103036d
# ╟─85bc871e-6067-4a1d-bbb9-e3a9b0dfbd85
# ╠═e1601070-d88e-4515-8ad7-af810d98e46d
# ╟─c871a80d-113c-47ee-8292-3a857752b094
# ╠═25caf1e9-788d-41a2-add3-f0cb451420aa
# ╠═5a387568-2284-4d4e-85fc-f106db7377c9
# ╟─61fe29b8-065c-430f-a28d-b463780f8b00
# ╠═2ecd106c-faee-469e-af4c-0f3b6c28bc02
# ╠═d096a6be-65a1-428d-9bfb-da7fe89f4c19
# ╠═f0413f68-1eb4-4a3d-bdcf-61e8aa96c0e7
# ╠═6d9759c3-63d9-4166-8884-cfaa99ee33c4
# ╠═54dd5337-c76f-4a23-8bec-a0cd57b433dd
# ╠═5220090a-c8ef-4b9b-ae17-6692fc4b19e3
