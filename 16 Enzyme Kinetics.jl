### A Pluto.jl notebook ###
# v0.19.22

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

# ╔═╡ 507bfc5b-2f6b-43ea-be67-f224aa842afa
using DifferentialEquations

# ╔═╡ 8e41bbfa-7241-44c3-bc2d-1a018549693e
using ParameterizedFunctions

# ╔═╡ 9914fc0f-4fd3-416d-8c3f-9b1c26999b20
using Plots

# ╔═╡ 057ba01c-7f78-456d-a92e-9c57bf183add
using PlutoUI

# ╔═╡ 5ed95079-b774-4cbe-8c00-1a3c1388539c
html"<button onclick='present()'>present</button>"

# ╔═╡ e3979b13-75a1-4a6c-9e77-40bd01e89404
md" # 8. Enzyme-substrate kinetics

## 8.1 Reminder of concept of Gibbs Free Enery 

A fundamental task of proteins is to act as catalysts of chemical reactions in the cell. A catalyst is actually a substance or factor to speed up our (chemical or biological) processes/reactions without taking part in them. 

In biochemical systems, these catalysts are called enzymes, and they increase the rate of virtually all the chemical reactions within cells.
##
To understand how enzymes speed up biochemical reactions we will introduce the concept of __Gibbs Free Energy__ (originally called _available energy_): a thermodynamic potential proportional to the maximum of reversible work that may be performed by a thermodynamic system at a constant temperature and pressure 

```math
\begin{align}
 ΔG^0 = ΔH^0 -T \cdot ΔS^0
\end{align}
```
$ΔG^0$ is the change in Gibbs free energy, $T$ is the temperature (SI unit: kelvin),  $ΔS^0$ is the change in entropy (SI unit: joule per kelvin), $ΔH^0$ is the change in enthalpy (SI unit: joule).
##



"

# ╔═╡ fee3f3ef-a9a2-4bb7-965d-46c0ad7d289c
Resource("https://schoolbag.info/chemistry/mcat_2/mcat_2.files/image357.jpg",:width => 400)

# ╔═╡ 96208e5e-2937-41b1-a270-f3f06e1e7ad9
md"
##
A consequence of the second law of thermodynamics, is that every reaction tends to achieve a minimum in its Gibbs free energy. A quantitative measure of the favorability of a given reaction is the change ΔG caused by the reaction. ∆G can be understood as the amount of _free_ or _useful_ energy available to do work.

##
when the value of ∆G>0 means that the reaction requires free energy, so it cannot take place spontaneously (endergonic).

On the contrary, when the value of ∆G<0 means that the reaction liberates free energy as its proceeds, so its can take place spontaneously (exergonic).

"

# ╔═╡ a40571e8-2112-420f-94a0-dc31706e3e08
md"
## 8.2 The Transition State and the role of enzymes

Between the state where molecules are reactants and the state where molecules are products in biochemical reactions, there is often a _bump_ in the value of ∆G, which is often refereed as the __transition state__. In the transition state, the reactants are combined in a intermediate species called the activated complex. The stability of this intermediate activated complex affects the rate of the reaction (if it goes towars the producs or it goes towards the reactants). The higher the bump, the more difficult for the reaction to take place, even in conditions of strongly exergonic biochemical reactions (∆G<<0). 

"

# ╔═╡ e02fd130-545e-4983-b12f-14eddb6f4bf3
function gauss(x,E_a2)
E_a2*exp(-(x - 5)^2/(2 *  2^2))/(sqrt(2*  π) * 2);
end;

# ╔═╡ 585ed61f-8c22-4deb-91b3-7dac9789fe25
function sigm(x)
    x^2 / (x^2+5^2) /10
end;

# ╔═╡ 0984844f-868a-4c58-9650-cddb4bd75bba
begin
	E2_slide = @bind E_a2 html"<input type=range min=1.0 max=2.0 step=0.1>"

	
	md"""
	##
	**Move the silder to change the Activation Energy**
	
	Ea: $(E2_slide)
	
	"""
end

# ╔═╡ a28457ea-be99-4cb8-a87e-797f5fee7421
begin
	P1=plot(x -> sigm(x)+gauss(x,E_a2),collect(1:0.1:12),title=("endergonic"),seriestype=:line,ylabel=("Gibbs free energy(G)"),xlabel=("reaction coordinate"),label="DeltaG > 0",yaxis=[0,0.5])
	P2=plot(x -> 0.1-sigm(x)+gauss(x,E_a2),collect(1:0.1:12),title=("exergonic"),seriestype=:line,ylabel=("Gibbs free energy(G)"),xlabel=("reaction coordinate"),yaxis=[0,0.5],label="DeltaG < 0")
	plot(P1,P2,layout=(1,2),legend=true,size = (800, 300))
end

# ╔═╡ cf486828-89ef-4698-b658-4205e847ede1
md"
##
A way to speed up (i.e., catalyze) a chemical reaction is by stabilization of this intermediate activated complex. A fundamental task of proteins in the cell is to act as these catalysts that increase the rate of virtually all the chemical reactions within cells. 
##
These type of proteins are called enzymes. In the absence of enzymatic catalysis, most biochemical reactions are so slow that they would not occur under the mild conditions of temperature and pressure that are compatible with life. Enzymes accelerate the rates of such reactions by well over a million-fold, so reactions that would take years in the absence of catalysis can occur in fractions of seconds if catalyzed by the appropriate enzyme.
##
Enzymes work basically as a key that brings molecules together in the right orientation, 
"

# ╔═╡ 2f66f76c-01fb-4176-93bc-3de8d5831e8f
Resource("https://1.bp.blogspot.com/-ttv1gwKAvsc/X8Y4zo6S6lI/AAAAAAAAAQg/LkU50ozS8lsKPsnZ6keV_FU1binSjIAxwCLcBGAsYHQ/w1200-h630-p-k-no-nu/what%2Bis%2Benzyme%2Bin%2Bbiology.PNG")

# ╔═╡ dac8f83a-bab8-4ead-9a64-066218e77d5b
md"
##
In brief, the main features of an enzyme are:

- Increase the rate of reactions without being consumed or changed by the reaction. 
- Increase reaction rates without altering the chemical equilibrium between reactants and products (they do not change ∆G).
- Enzymes work by lowering the energy of the transition state. They stabilise the intermediate activated complex (they reduce the entropy)"

# ╔═╡ a1dab4c7-2629-4427-8bff-e380dfae5bc8
md" ## 
Another very important effect of enzymes is that they reduce the dependence on temperature, compared to uncatalyzed reactions. 

We have seen that the speed of a reaction is proportional to the concentration of reactants, and to the kinetic rate constant. The Arrhenius Equation 
The dependence of the rate constant with the activation energy is set by the Arrhenius Equation:
```math
k = A e^{-E_a/RT} 
```
The parameter of $E_a$ is the activation energy, and $A$ is a parameter that can be identified as the rate of frequency of collisions and the orientation of a favorable collision probability. $A$ is rate of the reaction if either the activation energy were zero or the temperature was infinite.

##
The activation energy is the energy required for two molecules to interact, it is like a barrier of energy.  More temperature means more energy to cross the barrier.

Enzyme have a key effect in how the rate constant of biochemical reactions becomes more indeependent on the temperature, since the exponential term includes the activation energy as the numerator and the temperature as the denominator, a smaller activation energy will have less of an impact on the rate constant compared to a larger activation energy. 

Hence, the rate of an uncatalyzed reaction is more affected by temperature changes than a catalyzed reaction.

##
the Arrienius equation in linear form was:

```math
\begin{align}  
\ln k &= \ln \left(Ae^{-E_a/RT} \right) \\
&= \ln A + \ln \left(e^{-E_a/RT}\right) \\
&= \left(\dfrac{-E_a}{R}\right) \left(\dfrac{1}{T}\right) + \ln A \\
&=\ln A - \dfrac{E_{a}}{RT}
\end{align}
```

We can see that by lowering the $E_a$, enzymes lower the dependence on temperature of the rates of a reaction"

# ╔═╡ 3812d970-60fb-4fa7-810f-0f92dc4eb28b
begin
	E_slide = @bind E_a html"<input type=range min=1.0 max=20 step=1>"

	
	md"""
	##
	**Move the silder to change the Activation Energy**
	
	Ea: $(E_slide)
	
	"""
end

# ╔═╡ 55434867-6774-4863-926c-990c9deac192
begin
	A=1
	T= LinRange(100,300,100)
	p1=plot(T,A.*exp.(-E_a./(0.082.*T)),ylims = (0,1),title=("Without Enzyme"))
	p2=plot(T,A ./2 .*exp.(-E_a./10 ./ (0.082.*T)),ylims = (0,1),title=("With Enzyme"))
	plot(p1,p2,layout=(1,2),legend=true,size = (800, 500))
	#title!("Rate constant change with temperature")
	ylabel!("k")
	xlabel!("Temperature [K]")
end

# ╔═╡ 20bbb617-1d2f-4210-8a3b-19a1631ba70d
md"##
This is one of the reasons underlying the observed robustness of biological systems to changes in the temperature, since most of the reactions that take place in biologoical processes are calalyzed by enzymes."

# ╔═╡ b87efce9-71a9-4407-9d77-bd029e6aa557
md"## 8.3 The enzyme-substrate interaction scheme 

We will now solve the differential equations that result from a generic enzimatic reaction. The scheme of the reaction is as follows:

```math
s + e \overset{k_1}{\underset{k_2}{\longleftrightarrow}}  c \overset{k_3}{\longrightarrow} p + e \tag{1}
```
##
where $s$ represents the substrate, $e$ represents the enzyme, $c$ represents the complex, and $p$ the product of the reaction. 

"

# ╔═╡ ddcc3158-a0f7-4d56-80d2-1b6d53a68d96
md" 
##
> __Task 1__: Generate the Differential equations for the Enzyme-substrate system.
"

# ╔═╡ dc68a2ab-773b-40bd-a5d3-dfdd2c987489
md"
##
__Solution__: The reaction can be written as three irreversible reactions: 

```math
\begin{align*}
s + e \rightarrow  c \tag{2}\\
c \rightarrow s + e \tag{3}\\
c \rightarrow p + e  \tag{4}
 \end{align*}
```

##
with kinetic constants $k_1$, $k_2$ and $k_3$, respectively. If we rename variables as $X_1$=s,$X_2$=c, $X_3$=e, $X_4$=p to use the state vector approach.

```math
\begin{align*}
X_1 + X_3 \rightarrow  X_2 \tag{2}\\
X_2 \rightarrow X_1 + X_3 \tag{3}\\
X_2 \rightarrow X_3 +  X_4 \tag{4}
\end{align*}
```

"

# ╔═╡ 30afb560-a44a-11ec-1ebf-6fbfadad9353
md"
##

The kinetics of these chemical equations is described by the __Mass action Law__ 

```math
A=\begin{bmatrix}
 1 & 0 & 1 & 0\\
 0 & 1 & 0 & 0 \\
 0 & 1 & 0 & 0\end{bmatrix} ;
B=\begin{bmatrix}
 0 & 1 & 0 & 0\\
 1 & 0 & 1 & 0 \\
 0 & 0 & 1 & 1\end{bmatrix} ; (B-A)^T= \begin{bmatrix}
-1 & 1 & 0\\
 1 & -1  & -1\\
 -1 & 1 & 1 \\
 0 & 0 & 1\end{bmatrix} 
```

##

 in this particular case

```math
K=\begin{pmatrix}
 k_1 & 0 & 0 \\
 0 & k_2 & 0 \\
 0 & 0 & k_3 
\end{pmatrix}
```

##
```math
X^A=\begin{pmatrix}
X_1^1\cdot X_2^0 \cdot X_3^1 \cdot X_4^0 \\
X_1^0\cdot X_2^1 \cdot X_3^0 \cdot X_4^0 \\
X_1^0\cdot X_2^1 \cdot X_3^0 \cdot X_4^0 \\
\end{pmatrix} = \begin{pmatrix}
X_1 \cdot X_3\\
 X_2\\
  X_2
\end{pmatrix}
```
##
so 

```math
\begin{align}
 \begin{bmatrix}
\frac{\mathrm{d} X_1}{\mathrm{d} t}\\ \frac{\mathrm{d} X_2}{\mathrm{d} t} \\ \frac{\mathrm{d} X_3}{\mathrm{d} t} \\ \frac{\mathrm{d} X_4}{\mathrm{d} t}\end{bmatrix}= \begin{bmatrix}
-1 & 1 & 0\\
 1 & -1  & -1\\
 -1 & 1 & 1 \\
 0 & 0 & 1\end{bmatrix}\begin{pmatrix}
 k_1 & 0 & 0 \\
 0 & k_2 & 0 \\
 0 & 0 & k_3 
\end{pmatrix} \begin{pmatrix}
X_1 \cdot X_3\\
 X_2\\
  X_2
\end{pmatrix} 
\end{align}
```
##
that when we operate the matrices becomes

```math
\begin{align}
 \begin{bmatrix}
\frac{\mathrm{d} X_1}{\mathrm{d} t}\\ \frac{\mathrm{d} X_2}{\mathrm{d} t} \\ \frac{\mathrm{d} X_3}{\mathrm{d} t} \\ \frac{\mathrm{d} X_4}{\mathrm{d} t}\end{bmatrix}= \begin{bmatrix}
-1 & 1 & 0\\
 1 & -1  & -1\\
 -1 & 1 & 1 \\
 0 & 0 & 1\end{bmatrix} \begin{pmatrix}
k_1 \cdot X_1 \cdot X_3\\
k_2 \cdot X_2\\
k_3 \cdot  X_2
\end{pmatrix} 
\end{align}
```
##
and operating further,

```math
\begin{align}
 \begin{bmatrix}
\frac{\mathrm{d} X_1}{\mathrm{d} t}\\ \frac{\mathrm{d} X_2}{\mathrm{d} t} \\ \frac{\mathrm{d} X_3}{\mathrm{d} t} \\ \frac{\mathrm{d} X_4}{\mathrm{d} t}\end{bmatrix}= \begin{bmatrix}
-k_1 \cdot X_1 \cdot X_3 + k_2 \cdot X_2\\
 k_1 \cdot X_1 \cdot X_3 - k_2 \cdot X_2  - k_3 \cdot  X_2\\
-k_1 \cdot X_1 \cdot X_3 + k_2 \cdot X_2 + k_3 \cdot  X_2\\
 k_3 \cdot  X_2\end{bmatrix} 
\end{align}
```
##
which rearranging terms provides the following set of coupled differential equations 

```math
\begin{align*}
  \dot{X_1}&= k_2 \cdot X_2 - k_1 X_3 \cdot X_1 \tag{8}\\ 
  \dot{X_2}&= k_1 \cdot X_3 \cdot X_1 - (k_2+k_3) \cdot X_2 \tag{9}\\
  \dot{X_3}&= (k_2+k_3) \cdot X_2 - k_1 \cdot X_3(t) \cdot X_1 \tag{10}\\ 
  \dot{X_4}&= k_3 \cdot X_2 \tag{11}
\end{align*}
```
##
To solve it numericaly, we define the ODE problem as follows:
"

# ╔═╡ 8fca7a96-b7d5-4bac-86f4-400f260db607
function simpleODEEnzyme!(du,u,p,t)
    k1,k2,k3 = p
    du[1] = k2*u[2] - k1*u[1]*u[3]
    du[2] = k1*u[1]*u[3] - (k2 + k3) * u[2]
    du[3] = -k1*u[1]*u[3]+(k2 + k3) * u[2]
    du[4] = k3 * u[2] 
end

# ╔═╡ 3fe78776-a90a-457a-9176-7a495d0d3e30
md"
##
We can also use the simpler DSL direct notation, this is useful because we can use directly the differential equations with the original names of the variables:

```math
\begin{align*}
\dot{s}  &= k_2 \cdot c - k_1 \cdot e \cdot s \tag{5}\\
\dot{c} &= k_1 \cdot e\cdot s - (k_2+k_3) \cdot c \tag{6}\\
\dot{e} &= (k_2+k_3) \cdot c -k_1 \cdot e \cdot s \tag{7}\\
\dot{p} &= k_3 \cdot c
\end{align*}
```


  
(in Julia, we cannot use $e$ as the name of a variable, so we use $en$)
	"

# ╔═╡ 3f4053b5-f0ed-4b6b-8a13-eb7b039a8fad
enzyme_kinetics! = @ode_def ab1 begin
  ds = -k1*en*s+k2*c
  dc = k1*en*s-(k2+k3)*c
  den = -k1*en*s+(k2+k3)*c
  dp = k3*c
end k1 k2 k3;

# ╔═╡ 8801582a-dabd-4444-832c-8cc137225396
begin
	k1=1e3
	k2=0.1
	k3=0.05
	e0=0.002
	s0=0.002
	prob1 = ODEProblem(enzyme_kinetics!,[s0,0,e0,0],(0.0,100),[k1,k2,k3])
	sol1 = solve(prob1)
	plot(sol1)
	title!("Enzyme-substrate kinetics using DSL notation")
	xlabel!("Time [s]")
	ylabel!("Concentration [M]")
end

# ╔═╡ 435a8401-1c27-4a85-bb78-ae8a9c39eeb4
md"
##
Next, we will take advantage of the __Mass Conservation Law__ to reduce the number of independent variables (i.e, the number of differential equations to solve). 

```math
\begin{align*}
C_1 \cdot X_1 + C_2 \cdot X_2 + C_3 \cdot X_3 + C_4 \cdot X_4 & = cte \\
\end{align*}
```
##
so
```math
\begin{align*}
C \cdot (B-A)^T &=0\\ 
\begin{pmatrix}C_1 & C_2 & C_3 & C_4 \end{pmatrix} \begin{bmatrix}
-1 & 1 & 0\\
 1 & -1  & -1\\
 -1 & 1 & 1 \\
 0 & 0 & 1\end{bmatrix} &=0
 \end{align*}
```
##
the solution of the system of equations is: 

```math
\begin{align*}
 -C_1 + C_2 - C_3= 0 &\implies  C_2 = C_3 + C_1\\
 -C_2 + C_3 + C_4= 0 &\implies  C_2 = C_3 + C_4 
\end{align*}
```
##
Therefore, from the two equations above, $C_1 = C_4$, and the conservation is:

```math
\begin{align*}
C_1 \cdot X_1(t) + (C_3 + C_1) \cdot X_2(t) + C_3 \cdot X_3(t) + C_1 \cdot X_4(t) &= cte\\
C_1 \cdot X_1(t) + C_3 \cdot X_2(t) + C_1 \cdot X_2(t)  + C_3 \cdot X_3(t) + C_1 \cdot X_4(t) &= cte\\
\end{align*}
```
##
which is valid for any value of $C_1$ and $C_3$, so, separating these two dependencies

```math
\begin{align*}
C_1 \cdot X_1(t) + C_1 \cdot X_2(t)  + C_1 \cdot X_4(t) &= cte\\
\end{align*}
```

which makes sense, since the $s$, $c$ and $p$ are three forms of the same protein, and their sum has to be constant,a nd equal to the initial value of substrate used. 
##
Using the value $C_1=1$ we obtain:

```math
\begin{align*}
X_1(t) + X_2(t) + X_4(t) &= cte = X_1(0) + X_2(0) + X_4(0) =  X_1(0) \\
\end{align*}
```

is our first conservation law (typically we start a biochemical reaction adding only $e$ and $s$). 

##
For the second, we use the parameters that depend on the other constant $C_2$ 

```math
\begin{align*}
C_3 \cdot X_2(t) + C_3 \cdot X_3(t) &= cte
\end{align*}
```

which also makes sense, because the amoun of enzyme is not consumed and is only in teh form of the reactans $e$ and $c$. 
##
So, giving the value $C_3=1$, we obtain the conservation law, 

```math
\begin{align*}
X_2(t) + X_3(t) & = cte &= X_2(0) + X_3(0) = X_3 (0) 
\end{align*}
```

which is true because the enzyme is a catalyst that facilitates the reaction but does not react itself, and if we assume that there is no complex $c(0)=0$ before the reacton starts, we can assume ($e(0)+c(0)=e_0$). 
##
This conservation law allows us to reduce the four differential equations into the following three coupled ordinary differential equations:

```math
\begin{align*}
  \dot{s}  &= k_2 \cdot c - k_1 \cdot (e_0-c) \cdot s \tag{12}\\
  \dot{c} &= k_1 \cdot (e_0-c)\cdot s - (k_2+k_3) \cdot c \tag{13}\\
  \dot{p} &= k_3 \cdot c \tag{14}
  \end{align*}
```
##
which rearranging terms become:

```math
\begin{align*}
  \dot{s}  &= (k_2+k_1 \cdot s) \cdot c -k_1 \cdot e_0 \cdot s \tag{15}\\ 
  \dot{c} &= k_1\cdot e_0\cdot s-(k_3+k_2+k_1 \cdot s) \cdot c \tag{16}\\
  \dot{p} &= k_3 \cdot c \tag{17}
  \end{align*}
```
  

"

# ╔═╡ a4034df8-7185-41dd-8c4e-a46e4ac6028f
md" 
##
> __Task 2__: Sove numerically the full system and comprate it with the solution obtained using mass conservation. Use the following values:
>	k1=1e3
>	k2=0.1
>	k3=0.05
>	e0=0.002
>	s0=0.002
>    time= 100 s
"

# ╔═╡ 75703719-6fc1-487b-b0c9-14ba9e339239
enzyme_kinetics2! = @ode_def ab begin
  ds = (k2+s*k1)*c-k1*e0*s
  dc = k1*e0*s-(k2+k3+k1*s)*c
  dp = k3*c
end k1 k2 k3 e0;

# ╔═╡ 3021edfe-198c-49d0-8006-667b83c568c5
begin
	prob2 = ODEProblem(enzyme_kinetics2!,[s0,0,0],(0.0,100.0),[k1,k2,k3,e0])
	sol2 = solve(prob2)
	plot(sol2)
	e=[e0-u[2] for (u,t) in tuples(sol2)]
	plot!(sol2.t,e,
	        label="e",
	        linealpha = 1,
	        linewidth = 3,
	        linecolor = :purple)
	title!("Enzyme-substrate kinetics simplified using mass action")
	xlabel!("Time [s]")
	ylabel!("Concentration [M]")
end

# ╔═╡ a04787f9-8262-43dc-921f-be23863a019c
md"## 8.4 Quasi steady state approximation:

The other conservation law could be used to reduce the system to just two equations, but instead, we will use another approiximation that allows us to go further and reduce the system to a single differential equation. This aproximation si swidely used in all biochemistry books and has been famously proposed by Miachelis and Menten. 
##
In brief, it assumes that the formation and degradation of the intermediate complex is often very fast and reaches equilibrium early in the reaction. 

```math
k_1\cdot x_1(0) ≈ k_2 >> k_3 \tag{18}
```


The substrate-enzyme binding occurs at much faster time scales than the turnover into product (often the case in biologically relevant biochemical reactions)
##
Under these circumstances one expects that after an initial short transient period there will be a balance between the formation of the enzyme-substrate complex and the breaking apart of complex (either to enzyme and substrate, or to enzyme and product) 

```math
\frac{d c}{dt}= \frac{d e}{dt}=0 \tag{19}
```
##
In these conditions, we can calculate the equilibrium for $c$ from eq. 2.16 as 

```math
c_{eq}=\frac{k_1\cdot e_0\cdot s}{k_3+k_2+k_1 \cdot s} \tag{20}
```

This is called the __pseudo-steady state aproximation__ for biochemical systems"

# ╔═╡ 7eb42169-410b-4f1a-bc4c-f9230cd3a904
begin
	c_eq=[(k1 * e0 * u[1])/(k3+k2+k1 * u[1]) for (u,t) in tuples(sol2)]
	plot!(sol2.t,c_eq,
	        label="c_eq",        
	        linealpha = 0.5,
	        linewidth = 3,
	        linestyle= :dash,
	        linecolor = :red)
end

# ╔═╡ c4e8b093-aed0-4d49-8cbd-b57f950a29fa
md"
##
As we can see, the equilibrium solution $c_{eq}$ is quite close to the concentration of $c$, in this conditions we can simplify all equations as:
  
```math
\begin{align*}
  \dot{s}  &= (k_2+k_1 \cdot s) \cdot \frac{k_1\cdot e_0\cdot s}{k_3+k_2+k_1 \cdot s} -k_1 \cdot e_0 \cdot s \tag{21}\\ 
  \dot{p}  &= k_3 \frac{k_1\cdot e_0\cdot s}{k_3+k_2+k_1 \cdot s} \tag{22}
\end{align*}
```
##
Eq.  can be rewritten in a more simplified form, after some some algebraic manipulation:

```math
\begin{align*}
  \dot{s}  &= (k_2+k_1 \cdot s) \cdot \frac{k_1\cdot e_0\cdot s}{k_3+k_2+k_1 \cdot s} -k_1 \cdot e_0 \cdot s \tag{23}\\ 
   \dot{s}  &= k_1 \cdot e_0 \cdot s  \left [ \frac{k_2+k_1 \cdot s}{k_3+k_2+k_1 \cdot s} -1 \right ]  \tag{24}\\ 
\dot{s}  &= k_1 \cdot e_0 \cdot s  \left [ \frac{k_2+k_1 \cdot s - k_3-k_2-k_1 \cdot s}{k_3+k_2+k_1 \cdot s} \right ]  \tag{25}\\ 
\dot{s}  &=   \left [ \frac{- k_3\cdot k_1 \cdot e_0  \cdot s}{k_3+k_2+k_1 \cdot s} \right ]  \tag{26}\\ 
\end{align*}
```
##
So the final set of equations takes the form:

```math
\begin{align*}
  \dot{s}  &=- \frac{k_1 \cdot k_3 \cdot e_0 \cdot s}{k_3+k_2+k_1 \cdot s} \tag{27}\\ 
  \dot{p}  &=  \frac{ k_3 \cdot k_1\cdot e_0\cdot s}{k_3+k_2+k_1 \cdot s} \tag{28}
\end{align*}
```


"

# ╔═╡ 75e15198-e559-4201-a292-0d5e85898987
enzyme_kinetics10! = @ode_def ab10 begin
     ds = - (k3 * e0 * s)/((k3+k2)/k1 +  s)
     dp = (k3 * e0 * s)/((k3+k2)/k1 +  s)
end k1 k2 k3 e0;



# ╔═╡ feeb275b-f7f2-4d10-9449-48657dc3599a
begin
	prob10 = ODEProblem(enzyme_kinetics10!,[s0,0.00001],(0.0,100.0),[k1,k2,k3,e0])
	sol10 = solve(prob10)
	plot(sol10)
	title!("Enzyme-substrate kinetics simplified using mass action")
	xlabel!("Time [s]")
	ylabel!("Concentration [M]")
end

# ╔═╡ 25aa3755-d698-4e01-beff-90567b54fa62
md"
Now, the conservation of mass $s+p=s_0$ allows us to reduce everything to a single ODE:
```math
  \dot{p}  =  \frac{ k_3 \cdot k_1\cdot e_0\cdot (s_0-p)}{k_3+k_2+k_1 \cdot (s_0-p)} \tag{29}
```

Or, if the divide every term by $k_1$

```math
  \dot{p}  =  \frac{ k_3 \cdot e_0\cdot (s_0-p)}{\frac{k_3+k_2}{k_1}+(s_0-p)} \tag{30}
```
"

# ╔═╡ 3ae1813b-4741-44b9-9e96-0ddb8246d57e
begin
	enzyme_kinetics3! = @ode_def ab3 begin
	  dp = (k3 * e0 * (s0 - p))/((k3+k2)/k1 +  (s0 - p))
	    end k1 k2 k3 e0 s0
	
	
	prob3 = ODEProblem(enzyme_kinetics3!,[0.0],(0.0,100.0),[k1,k2,k3,e0,s0])
	sol3 = solve(prob3)
	plot(sol3)

	prob5 = ODEProblem(enzyme_kinetics2!,[s0,0,0],(0.0,100.0),[k1,k2,k3,e0])
	sol5 = solve(prob5)
	pp=[u[3] for (u,t) in tuples(sol5)]
	
	plot!(sol5.t,pp,
	        label="p full mass action",
	        linealpha = 1,
	        linewidth = 3,
	        linecolor = :purple)
	
	title!("Enzyme-substrate kinetics simplified using mass action")
	xlabel!("Time [s]")
	ylabel!("Concentration [M]")
end

# ╔═╡ 5847adf5-8444-4563-9562-81de6ee87d7b
md"
##
The steady state solution aproximates quite well to the rate of the reaction at the starting point of the reaction, but not at the intermediate time points. This, of course, will depend on the reaction parameters and the initial conditions. We can set up a program to test when the steady state aproximation is more correct."

# ╔═╡ 708b2261-19ca-4d90-839a-65eedf9458e3
md" we can see that increasing the $k_2$ and reducing the amount of enzyme $e_0$ results in a much better approximation

## 8.6 Michaelis-Menten equation: 

The previous equation 

```math
  \dot{p}  =  \frac{ k_3 \cdot e_0\cdot (s_0-p)}{\frac{k_3+k_2}{k_1}+(s_0-p)} \tag{31}
```

leads to the traditional Michaelis-Menten equation, which predicts the initial turnover rate of the enzymatic reaction $V_0$ as a function of initial substrate concentration $s_0$. 
##
So at the initial state of the reaction ($p=0$):

```math
  V_0  =  \frac{ k_3 \cdot e_0\cdot s_0}{\frac{k_3+k_2}{k_1}+s_0}=\frac{ v_{max} \cdot s_0}{K_M+s_0} \tag{32}
```

where the constant $K_M = \frac{k_3+k_2}{k_1}$ is called the Michaelis-Menten constant and $v_{max}= k_3 \cdot e_0$ is the maximum turn-over rate. 
The $K_M$  reflects the affinity of the reaction. Strong affinity means small $K_M$. At $s_0=K_M$ the turn-over rate is half maximal, i.e., $\dot{p}=\frac{v_{max}}{2}$

### Michaelis-Menten curve

The famous Michaelis-Menten plot is a culve that shows the dependence of the velocity of the reaction in terms of production of product (i.e., $\dot p$) on teh concetration of substrate. 
##
It showns three regimes:
- __linear__: there is a lot of enzyme available to bind the substrate, so the velocity of teh reaction increases linearly with the concetration of subtrate, as it ocurs in noncatalized reactions. 
- __Constrained__: the intermediate regime 
- __Saturated__: there is a lot of substrate available, and the enzyme is the limitant factor for the catalisis, so the speed of teh reactioon does not increase if the increase further the concetration of substrate.


"

# ╔═╡ 8f1a6377-1894-4c0b-8dfe-727c3d3c1c7e
md" ## 8.6 Validity of the Michaelis-Menten Approach

We can compare this analytical equation with the result of the numerical simulation of the full ODE system. To do that we solve the numerical system, and compare the rates of the reactions predicted by the Mass action and the Michaelis-Menten approach."

# ╔═╡ 1215f027-d12a-4fc1-ab1c-82c27fd75d95
begin
	s0_slide = @bind s0_ html"<input type=range min=0.001 max=0.005 step=0.001>"
	e0_slide = @bind e0_ html"<input type=range min=0.0005 max=0.005 step=0.001>"
	
	md"""
	**How many molecules do you have?**
	
	Initial concentration of substrate: $(s0_slide)
	
	Initial concentration of enzyme: $(e0_slide)
	"""
end

# ╔═╡ fac136db-3cfa-47d3-8329-b21c2cca74bb
begin
	s_vector= s0_*LinRange(0,1,100)
	plot(s_vector,k3 * e0_ * s_vector ./(s_vector.+(k3+k2)/k1),label="Michaelis Menten",)
	println("K_m = ",(k3+k2)/k1)
	title!("Michaelis-Menten Plot")
	xlabel!("concentration of substrate")
	ylabel!("speed of the reaction")
end

# ╔═╡ 78f83d23-863f-4113-a2ba-07282397e9f2
begin
	speed_vector = similar(s_vector)
	plot(s_vector,k3 * e0_ * s_vector ./(s_vector.+(k3+k2)/k1),label="Michaelis Menten",)
	for i in 1:100 
	    prob6 = ODEProblem(enzyme_kinetics2!,[s_vector[i],0,0],(0.0,10),[k1,k2,k3,e0_])
	    sol6 = solve(prob6)
	    # we shoudl compute the speed at the maximum formation of Complex. 
	    speed=[k3*u[2] for (u,t) in tuples(sol6)]
		speed_vector[i]=maximum(speed)
	end
	plot!(s_vector,speed_vector,label="rate full mass action")
	title!("Michaelis-Menten vs. Mass action")
	xlabel!("Concentration of substrate")
	ylabel!("Speed of the reaction")
end

# ╔═╡ f46dd83a-acb9-4e3f-88d9-5ab068478df7
md" 
##
> __Task 1__: (a) Generate the Michaelis-Menten curve ($v_0$ versus $e_0$) for the full system of equations of interactions enzyme-substrate obtained using the Mass Action Law.  Use the following values for the parameters: 
>	k1=1e3
>	k2=0.1
>	k3=0.05
>	e0=0.002
>	s0=0.002
>    time= 100 s
> (b) Use the same values to generate the Michaelis-menten curve ($v_0$ versus $e_0$) for the equation obtained after applying the pseudo-steady state aproximation.
> (c) Plot the analytical solution obtained by Michaelis-Menten and find parameters where the pseudo-steady state aproximation is a better approximation to the full system.
"

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DifferentialEquations = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
ParameterizedFunctions = "65888b18-ceab-5e60-b2b9-181511a3b968"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
DifferentialEquations = "~7.1.0"
ParameterizedFunctions = "~5.13.1"
Plots = "~1.27.1"
PlutoUI = "~0.7.37"

[extras]
CPUSummary = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.5"
manifest_format = "2.0"
project_hash = "b3c55e4cf515908a4dd62657799508bba2d9b552"

[[deps.AbstractAlgebra]]
deps = ["GroupsCore", "InteractiveUtils", "LinearAlgebra", "MacroTools", "Markdown", "Random", "RandomExtensions", "SparseArrays", "Test"]
git-tree-sha1 = "c0750f99036c12a14c93a7b10609cb892ffbd092"
uuid = "c3fe647b-3220-5bb0-a1ea-a7954cac585d"
version = "0.24.1"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.AbstractTrees]]
git-tree-sha1 = "03e0550477d86222521d254b741d470ba17ea0b5"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.3.4"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

[[deps.ArgCheck]]
git-tree-sha1 = "a3a402a35a2f7e0b87828ccabbd5ebfbebe356b4"
uuid = "dce04be8-c92d-5529-be00-80e4d2c0e197"
version = "2.3.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "62e51b39331de8911e4a7ff6f5aaf38a5f4cc0ae"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.2.0"

[[deps.ArrayInterface]]
deps = ["Compat", "IfElse", "LinearAlgebra", "Requires", "SparseArrays", "Static"]
git-tree-sha1 = "1ee88c4c76caa995a885dc2f22a5d548dfbbc0ba"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "3.2.2"

[[deps.ArrayLayouts]]
deps = ["FillArrays", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "56c347caf09ad8acb3e261fe75f8e09652b7b05b"
uuid = "4c555306-a7a7-4459-81d9-ec55ddd5c99a"
version = "0.7.10"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AutoHashEquals]]
git-tree-sha1 = "45bb6705d93be619b81451bb2006b7ee5d4e4453"
uuid = "15f4f7f2-30c1-5605-9d31-71845cf9641f"
version = "0.2.0"

[[deps.BandedMatrices]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "Random", "SparseArrays"]
git-tree-sha1 = "0a48e1230cbbf1dce4a6961e0c92615c8be20ddd"
uuid = "aae01518-5342-5314-be14-df237901396f"
version = "0.16.12"

[[deps.BangBang]]
deps = ["Compat", "ConstructionBase", "Future", "InitialValues", "LinearAlgebra", "Requires", "Setfield", "Tables", "ZygoteRules"]
git-tree-sha1 = "b15a6bc52594f5e4a3b825858d1089618871bf9d"
uuid = "198e06fe-97b7-11e9-32a5-e1d131e6ad66"
version = "0.3.36"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Baselet]]
git-tree-sha1 = "aebf55e6d7795e02ca500a689d326ac979aaf89e"
uuid = "9718e550-a3fa-408a-8086-8db961cd8217"
version = "0.1.1"

[[deps.Bijections]]
git-tree-sha1 = "705e7822597b432ebe152baa844b49f8026df090"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.1.3"

[[deps.BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "28bbdbf0354959db89358d1d79d421ff31ef0b5e"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.3"

[[deps.BoundaryValueDiffEq]]
deps = ["BandedMatrices", "DiffEqBase", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "NLsolve", "Reexport", "SparseArrays"]
git-tree-sha1 = "fe34902ac0c3a35d016617ab7032742865756d7d"
uuid = "764a87c0-6b3e-53db-9096-fe964310641d"
version = "2.7.1"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CEnum]]
git-tree-sha1 = "215a9aa4a1f23fbd05b92769fdd62559488d70e9"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.1"

[[deps.CPUSummary]]
deps = ["Hwloc", "IfElse", "Preferences", "Static"]
git-tree-sha1 = "68150205edbf60f0410ba2463b5b38eae44cad1f"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.1.15"

[[deps.CSTParser]]
deps = ["Tokenize"]
git-tree-sha1 = "b66abc140f8b90a1d6bc7bfad5c80070f8c1ddc6"
uuid = "00ebfdb7-1f24-5e51-bd34-a7502290713f"
version = "3.3.3"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c9a6160317d1abe9c44b3beb367fd448117679ca"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.13.0"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "bf98fa45a0a4cee295de98d4c1462be26345b9a1"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.2"

[[deps.CloseOpenIntervals]]
deps = ["ArrayInterface", "Static"]
git-tree-sha1 = "f576084239e6bdf801007c80e27e2cc2cd963fe0"
uuid = "fb6a15b2-703c-40df-9091-08a04967cfa9"
version = "0.1.6"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "12fc73e5e0af68ad3137b886e3f7c1eacfca2640"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.17.1"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonMark]]
deps = ["Crayons", "JSON", "URIs"]
git-tree-sha1 = "4cd7063c9bdebdbd55ede1af70f3c2f48fab4215"
uuid = "a80b9123-70ca-4bc0-993e-6e3bcb318db6"
version = "0.8.6"

[[deps.CommonSolve]]
git-tree-sha1 = "68a0743f578349ada8bc911a5cbd5a2ef6ed6d1f"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.0"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "96b0bc6c52df76506efc8a441c6cf1adcb1babc4"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.42.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.1+0"

[[deps.CompositeTypes]]
git-tree-sha1 = "d5b014b216dc891e81fea299638e4c10c657b582"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.2"

[[deps.CompositionsBase]]
git-tree-sha1 = "455419f7e328a1a2493cabc6428d79e951349769"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.1"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f74e9d5388b8620b4cee35d4c5a618dd4dc547f4"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.3.0"

[[deps.Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DEDataArrays]]
deps = ["ArrayInterface", "DocStringExtensions", "LinearAlgebra", "RecursiveArrayTools", "SciMLBase", "StaticArrays"]
git-tree-sha1 = "5e5f8f363c8c9a2415ef9185c4e0ff6966c87d52"
uuid = "754358af-613d-5f8d-9788-280bf1605d4c"
version = "0.2.2"

[[deps.DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3daef5523dd2e769dad2365274f760ff5f282c7d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.11"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DefineSingletons]]
git-tree-sha1 = "0fba8b706d0178b4dc7fd44a96a92382c9065c2c"
uuid = "244e2a9f-e319-4986-a169-4d1fe445cd52"
version = "0.1.2"

[[deps.DelayDiffEq]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "LinearAlgebra", "Logging", "NonlinearSolve", "OrdinaryDiffEq", "Printf", "RecursiveArrayTools", "Reexport", "UnPack"]
git-tree-sha1 = "52f54bd7f7bc1ce794add0ccf08f8fa21acfaed9"
uuid = "bcd4f6db-9728-5f36-b5f7-82caef46ccdb"
version = "5.35.1"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[deps.DiffEqBase]]
deps = ["ArrayInterface", "ChainRulesCore", "DEDataArrays", "DataStructures", "Distributions", "DocStringExtensions", "FastBroadcast", "ForwardDiff", "FunctionWrappers", "IterativeSolvers", "LabelledArrays", "LinearAlgebra", "Logging", "MuladdMacro", "NonlinearSolve", "Parameters", "PreallocationTools", "Printf", "RecursiveArrayTools", "RecursiveFactorization", "Reexport", "Requires", "SciMLBase", "Setfield", "SparseArrays", "StaticArrays", "Statistics", "SuiteSparse", "ZygoteRules"]
git-tree-sha1 = "df03eb34293066d699f8a535d1ccdcff94cb9765"
uuid = "2b5f629d-d688-5b77-993f-72d75c75574e"
version = "6.82.1"

[[deps.DiffEqCallbacks]]
deps = ["DataStructures", "DiffEqBase", "ForwardDiff", "LinearAlgebra", "NLsolve", "OrdinaryDiffEq", "Parameters", "RecipesBase", "RecursiveArrayTools", "SciMLBase", "StaticArrays"]
git-tree-sha1 = "c4b99e3a199e293e7290eea94ba89364d47ee557"
uuid = "459566f4-90b8-5000-8ac3-15dfb0a30def"
version = "2.22.0"

[[deps.DiffEqJump]]
deps = ["ArrayInterface", "Compat", "DataStructures", "DiffEqBase", "FunctionWrappers", "Graphs", "LinearAlgebra", "PoissonRandom", "Random", "RandomNumbers", "RecursiveArrayTools", "Reexport", "StaticArrays", "TreeViews", "UnPack"]
git-tree-sha1 = "eec5fd03c26dadc6b20f84d815309d060358e95b"
uuid = "c894b116-72e5-5b58-be3c-e6d8d4ac2b12"
version = "8.3.0"

[[deps.DiffEqNoiseProcess]]
deps = ["DiffEqBase", "Distributions", "LinearAlgebra", "Optim", "PoissonRandom", "QuadGK", "Random", "Random123", "RandomNumbers", "RecipesBase", "RecursiveArrayTools", "Requires", "ResettableStacks", "SciMLBase", "StaticArrays", "Statistics"]
git-tree-sha1 = "d6839a44a268c69ef0ed927b22a6f43c8a4c2e73"
uuid = "77a26b50-5914-5dd7-bc55-306e6241c503"
version = "5.9.0"

[[deps.DiffResults]]
deps = ["StaticArrays"]
git-tree-sha1 = "c18e98cba888c6c25d1c3b048e4b3380ca956805"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.0.3"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "dd933c4ef7b4c270aacd4eb88fa64c147492acf0"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.10.0"

[[deps.DifferentialEquations]]
deps = ["BoundaryValueDiffEq", "DelayDiffEq", "DiffEqBase", "DiffEqCallbacks", "DiffEqJump", "DiffEqNoiseProcess", "LinearAlgebra", "LinearSolve", "OrdinaryDiffEq", "Random", "RecursiveArrayTools", "Reexport", "SteadyStateDiffEq", "StochasticDiffEq", "Sundials"]
git-tree-sha1 = "3f3db9365fedd5fdbecebc3cce86dfdfe5c43c50"
uuid = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
version = "7.1.0"

[[deps.Distances]]
deps = ["LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "3258d0659f812acde79e8a74b11f17ac06d0ca04"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.7"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "c43e992f186abaf9965cc45e372f4693b7754b22"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.52"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.DomainSets]]
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "StaticArrays", "Statistics"]
git-tree-sha1 = "5f5f0b750ac576bcf2ab1d7782959894b304923e"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.5.9"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "90b158083179a6ccbce2c7eb1446d5bf9d7ae571"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.7"

[[deps.DynamicPolynomials]]
deps = ["DataStructures", "Future", "LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Pkg", "Reexport", "Test"]
git-tree-sha1 = "d0fa82f39c2a5cdb3ee385ad52bc05c42cb4b9f0"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.4.5"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[deps.EllipsisNotation]]
deps = ["ArrayInterface"]
git-tree-sha1 = "d7ab55febfd0907b285fbf8dc0c73c0825d9d6aa"
uuid = "da5c29d0-fa7d-589e-88eb-ea29b0a81949"
version = "1.3.0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ae13fcbc7ab8f16b0856729b050ef0c446aa3492"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.4+0"

[[deps.ExponentialUtilities]]
deps = ["ArrayInterface", "LinearAlgebra", "Printf", "Requires", "SparseArrays", "libblastrampoline_jll"]
git-tree-sha1 = "b026981973ccbe38682fbb4ccb0732fd6b1e1207"
uuid = "d4d017d3-3776-5f7e-afef-a10c40355c18"
version = "1.13.0"

[[deps.ExprTools]]
git-tree-sha1 = "56559bbef6ca5ea0c0818fa5c90320398a6fbf8d"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.8"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "d8a578692e3077ac998b50c0217dfd67f21d1e5f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.0+0"

[[deps.FastBroadcast]]
deps = ["LinearAlgebra", "Polyester", "Static"]
git-tree-sha1 = "f39bcc05eb0dcbd2c0195762df7a5737041289b9"
uuid = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
version = "0.1.14"

[[deps.FastClosures]]
git-tree-sha1 = "acebe244d53ee1b461970f8910c235b259e772ef"
uuid = "9aa1b823-49e4-5ca5-8b0f-3971ec8bab6a"
version = "0.3.2"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "deed294cde3de20ae0b2e0355a6c4e1c6a5ceffc"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.12.8"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "56956d1e4c1221000b7781104c58c34019792951"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.11.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "1bd6fc0c344fc0cbee1f42f8d2e7ec8253dda2d2"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.25"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.FunctionWrappers]]
git-tree-sha1 = "241552bc2209f0fa068b6415b1942cc0aa486bcc"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.2"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "51d2dfe8e590fbd74e7a842cf6d13d8a2f45dc01"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.6+0"

[[deps.GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "RelocatableFolders", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "9f836fb62492f4b0f0d3b06f55983f2704ed0883"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.64.0"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a6c850d77ad5118ad3be4bd188919ce97fffac47"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.64.0+0"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "83ea630384a13fc4f002b77690bc0afeb4255ac9"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.2"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "a32d672ac2c967f3deb8a81d828afc739c838a06"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "57c021de207e234108a6f1454003120a1bf350c4"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.6.0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.Groebner]]
deps = ["AbstractAlgebra", "Combinatorics", "Logging", "MultivariatePolynomials", "Primes", "Random"]
git-tree-sha1 = "1a13fd88bfaf0c5181da639edf4e773f756236c7"
uuid = "0b43b601-686d-58a3-8a1c-6623616c7cd4"
version = "0.2.0"

[[deps.GroupsCore]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "9e1a5e9f3b81ad6a5c613d181664a0efc6fe6dd7"
uuid = "d5909c97-4eac-4ecc-a3dc-fdd0858a4120"
version = "0.4.0"

[[deps.HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "0fa77022fe4b511826b39c894c90daf5fce3334a"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.17"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.HostCPUFeatures]]
deps = ["BitTwiddlingConvenienceFunctions", "IfElse", "Libdl", "Static"]
git-tree-sha1 = "18be5268cf415b5e27f34980ed25a7d34261aa83"
uuid = "3e5b6fbb-0976-4d2c-9146-d79de83f2fb0"
version = "0.1.7"

[[deps.Hwloc]]
deps = ["Hwloc_jll"]
git-tree-sha1 = "92d99146066c5c6888d5a3abc871e6a214388b91"
uuid = "0e44f5e4-bd66-52a0-8798-143a42290a1d"
version = "2.0.0"

[[deps.Hwloc_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d8bccde6fc8300703673ef9e1383b11403ac1313"
uuid = "e33a78d0-f292-5ffc-b300-72abe9b543c8"
version = "2.7.0+0"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "SpecialFunctions", "Test"]
git-tree-sha1 = "65e4589030ef3c44d3b90bdc5aac462b4bb05567"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.8"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
git-tree-sha1 = "2b078b5a615c6c0396c77810d92ee8c6f470d238"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.3"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.Inflate]]
git-tree-sha1 = "f5fc07d4e706b84f72d54eedcc1c13d92fb0871c"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.2"

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.InitialValues]]
git-tree-sha1 = "4da0f88e9a39111c2fa3add390ab15f3a44f3ca3"
uuid = "22cec73e-a1b8-11e9-2c92-598750a2cf9c"
version = "0.3.1"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IntervalSets]]
deps = ["Dates", "EllipsisNotation", "Statistics"]
git-tree-sha1 = "bcf640979ee55b652f3b01650444eb7bbe3ea837"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.5.4"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "91b5dcf362c5add98049e6c29ee756910b03051d"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.3"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.IterativeSolvers]]
deps = ["LinearAlgebra", "Printf", "Random", "RecipesBase", "SparseArrays"]
git-tree-sha1 = "1169632f425f79429f245113b775a0e3d121457c"
uuid = "42fd0dbc-a981-5370-80f2-aaf504508153"
version = "0.9.2"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[deps.JuliaFormatter]]
deps = ["CSTParser", "CommonMark", "DataStructures", "Pkg", "Tokenize"]
git-tree-sha1 = "7a6717e8055dc3ee240c8d44668e5908d40efd2a"
uuid = "98e50ef6-434e-11e9-1051-2b60c6c9e899"
version = "0.22.7"

[[deps.KLU]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse_jll"]
git-tree-sha1 = "cae5e3dfd89b209e01bcd65b3a25e74462c67ee0"
uuid = "ef3ab10e-7fda-4108-b977-705223b18434"
version = "0.3.0"

[[deps.Krylov]]
deps = ["LinearAlgebra", "Printf", "SparseArrays"]
git-tree-sha1 = "a024280a69c49f51ba29d2deb66f07508f0b9b49"
uuid = "ba0b0d4f-ebba-5204-a429-3ac8c609bfb7"
version = "0.7.13"

[[deps.KrylovKit]]
deps = ["LinearAlgebra", "Printf"]
git-tree-sha1 = "0328ad9966ae29ccefb4e1b9bfd8c8867e4360df"
uuid = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"
version = "0.5.3"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.LabelledArrays]]
deps = ["ArrayInterface", "ChainRulesCore", "LinearAlgebra", "MacroTools", "StaticArrays"]
git-tree-sha1 = "fbd884a02f8bf98fd90c53c1c9d2b21f9f30f42a"
uuid = "2ee39098-c373-598a-b85f-a56591580800"
version = "1.8.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "4f00cc36fede3c04b8acf9b2e2763decfdcecfa6"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.13"

[[deps.LayoutPointers]]
deps = ["ArrayInterface", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static"]
git-tree-sha1 = "b651f573812d6c36c22c944dd66ef3ab2283dfa1"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.6"

[[deps.LevyArea]]
deps = ["LinearAlgebra", "Random", "SpecialFunctions"]
git-tree-sha1 = "56513a09b8e0ae6485f34401ea9e2f31357958ec"
uuid = "2d8b4e74-eb68-11e8-0fb9-d5eb67b50637"
version = "1.0.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c7cb1f5d892775ba13767a87c7ada0b980ea0a71"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+2"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "c9551dd26e31ab17b86cbd00c2ede019c08758eb"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+1"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "f27132e551e959b3667d8c93eae90973225032dd"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.1.1"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LinearSolve]]
deps = ["ArrayInterface", "DocStringExtensions", "IterativeSolvers", "KLU", "Krylov", "KrylovKit", "LinearAlgebra", "RecursiveFactorization", "Reexport", "Requires", "SciMLBase", "Setfield", "SparseArrays", "SuiteSparse", "UnPack"]
git-tree-sha1 = "a25bc80647e44d0e1e1694b47000603497700b18"
uuid = "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"
version = "1.13.0"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "58f25e56b706f95125dcb796f39e1fb01d913a71"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.10"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoopVectorization]]
deps = ["ArrayInterface", "CPUSummary", "ChainRulesCore", "CloseOpenIntervals", "DocStringExtensions", "ForwardDiff", "HostCPUFeatures", "IfElse", "LayoutPointers", "LinearAlgebra", "OffsetArrays", "PolyesterWeave", "SIMDDualNumbers", "SLEEFPirates", "SpecialFunctions", "Static", "ThreadingUtilities", "UnPack", "VectorizationBase"]
git-tree-sha1 = "077c7c9d746cbe30ac5f001ea4c1277f64cc5dad"
uuid = "bdcacae8-1622-11e9-2a5c-532679323890"
version = "0.12.103"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.ManualMemory]]
git-tree-sha1 = "bcaef4fc7a0cfe2cba636d84cda54b5e4e4ca3cd"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.8"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[deps.Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[deps.Metatheory]]
deps = ["AutoHashEquals", "DataStructures", "Dates", "DocStringExtensions", "Parameters", "Reexport", "TermInterface", "ThreadsX", "TimerOutputs"]
git-tree-sha1 = "0886d229caaa09e9f56bcf1991470bd49758a69f"
uuid = "e9d8d322-4543-424a-9be4-0cc815abe26c"
version = "1.3.3"

[[deps.MicroCollections]]
deps = ["BangBang", "InitialValues", "Setfield"]
git-tree-sha1 = "6bb7786e4f24d44b4e29df03c69add1b63d88f01"
uuid = "128add7d-3638-4c79-886c-908ea0c25c34"
version = "0.1.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.ModelingToolkit]]
deps = ["AbstractTrees", "ArrayInterface", "ConstructionBase", "DataStructures", "DiffEqBase", "DiffEqCallbacks", "DiffEqJump", "DiffRules", "Distributed", "Distributions", "DocStringExtensions", "DomainSets", "Graphs", "IfElse", "InteractiveUtils", "JuliaFormatter", "LabelledArrays", "Latexify", "Libdl", "LinearAlgebra", "MacroTools", "NaNMath", "NonlinearSolve", "RecursiveArrayTools", "Reexport", "Requires", "RuntimeGeneratedFunctions", "SafeTestsets", "SciMLBase", "Serialization", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicUtils", "Symbolics", "UnPack", "Unitful"]
git-tree-sha1 = "cdadda283d438908ca4840c0a46ea61c4959e8a0"
uuid = "961ee093-0014-501f-94e3-6117800e7a78"
version = "8.5.3"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.MuladdMacro]]
git-tree-sha1 = "c6190f9a7fc5d9d5915ab29f2134421b12d24a68"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.2"

[[deps.MultivariatePolynomials]]
deps = ["DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "13ecec3d52a409be2e8653516955ec58f1c5a847"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.4.4"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "ba8c0f8732a24facba709388c74ba99dcbfdda1e"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.0.0"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "50310f934e55e5ca3912fb941dec199b49ca9b68"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.2"

[[deps.NLsolve]]
deps = ["Distances", "LineSearches", "LinearAlgebra", "NLSolversBase", "Printf", "Reexport"]
git-tree-sha1 = "019f12e9a1a7880459d0173c182e6a99365d7ac1"
uuid = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
version = "4.5.1"

[[deps.NaNMath]]
git-tree-sha1 = "b086b7ea07f8e38cf122f5016af580881ac914fe"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.7"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.NonlinearSolve]]
deps = ["ArrayInterface", "FiniteDiff", "ForwardDiff", "IterativeSolvers", "LinearAlgebra", "RecursiveArrayTools", "RecursiveFactorization", "Reexport", "SciMLBase", "Setfield", "StaticArrays", "UnPack"]
git-tree-sha1 = "aeebff6a2a23506e5029fd2248a26aca98e477b3"
uuid = "8913a72c-1f9b-4ce2-8d82-65094dcecaec"
version = "0.3.16"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "043017e0bdeff61cfbb7afeb558ab29536bbb5ed"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.10.8"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ab05aa4cc89736e95915b01e7279e61b1bfe33b8"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.14+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Optim]]
deps = ["Compat", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "bc0a748740e8bc5eeb9ea6031e6f050de1fc0ba2"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.6.2"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.OrdinaryDiffEq]]
deps = ["Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "ExponentialUtilities", "FastClosures", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "LinearSolve", "Logging", "LoopVectorization", "MacroTools", "MuladdMacro", "NLsolve", "NonlinearSolve", "Polyester", "PreallocationTools", "RecursiveArrayTools", "Reexport", "SciMLBase", "SparseArrays", "SparseDiffTools", "StaticArrays", "UnPack"]
git-tree-sha1 = "509aa6d3b2773e5109d4a4dd9a300259ac727961"
uuid = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"
version = "6.7.1"

[[deps.PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "e8185b83b9fc56eb6456200e873ce598ebc7f262"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.7"

[[deps.ParameterizedFunctions]]
deps = ["DataStructures", "DiffEqBase", "DocStringExtensions", "Latexify", "LinearAlgebra", "ModelingToolkit", "Reexport", "SciMLBase"]
git-tree-sha1 = "2f48f745e976dc5575bbc301e6c63b8fb5f12155"
uuid = "65888b18-ceab-5e60-b2b9-181511a3b968"
version = "5.13.1"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "85b5da0fa43588c75bb1ff986493443f821c70b7"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.2.3"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Requires", "Statistics"]
git-tree-sha1 = "a3a964ce9dc7898193536002a6dd892b1b5a6f1d"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "2.0.1"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "bb16469fd5224100e422f0b027d26c5a25de1200"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.2.0"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "1690b713c3b460c955a2957cd7487b1b725878a7"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.27.1"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "bf0a1121af131d9974241ba53f601211e9303a9e"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.37"

[[deps.PoissonRandom]]
deps = ["Random", "Statistics", "Test"]
git-tree-sha1 = "44d018211a56626288b5d3f8c6497d28c26dc850"
uuid = "e409e4f3-bfea-5376-8464-e040bb5c01ab"
version = "0.4.0"

[[deps.Polyester]]
deps = ["ArrayInterface", "BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "ManualMemory", "PolyesterWeave", "Requires", "Static", "StrideArraysCore", "ThreadingUtilities"]
git-tree-sha1 = "ad769d3f29cffb33380ab28318a10c1ccb19c827"
uuid = "f517fe37-dbe3-4b94-8317-1923a5111588"
version = "0.6.7"

[[deps.PolyesterWeave]]
deps = ["BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "Static", "ThreadingUtilities"]
git-tree-sha1 = "7e597df97e46ffb1c8adbaddfa56908a7a20194b"
uuid = "1d0040c9-8b98-4ee7-8388-3f51789ca0ad"
version = "0.1.5"

[[deps.PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

[[deps.PreallocationTools]]
deps = ["Adapt", "ArrayInterface", "ForwardDiff", "LabelledArrays"]
git-tree-sha1 = "6c138c8510111fa47b5d2ed8ada482d97e279bee"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.2.4"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "d3538e7f8a790dc8903519090857ef8e1283eecd"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.5"

[[deps.Primes]]
git-tree-sha1 = "984a3ee07d47d401e0b823b7d30546792439070a"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "0c03844e2231e12fda4d0086fd7cbe4098ee8dc5"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+2"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "78aadffb3efd2155af139781b8a8df1ef279ea39"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.4.2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Random123]]
deps = ["Random", "RandomNumbers"]
git-tree-sha1 = "afeacaecf4ed1649555a19cb2cad3c141bbc9474"
uuid = "74087812-796a-5b5d-8853-05524746bad3"
version = "1.5.0"

[[deps.RandomExtensions]]
deps = ["Random", "SparseArrays"]
git-tree-sha1 = "062986376ce6d394b23d5d90f01d81426113a3c9"
uuid = "fb686558-2515-59ef-acaa-46db3789a887"
version = "0.4.3"

[[deps.RandomNumbers]]
deps = ["Random", "Requires"]
git-tree-sha1 = "043da614cc7e95c703498a491e2c21f58a2b8111"
uuid = "e6cf234a-135c-5ec9-84dd-332b85af5143"
version = "1.5.3"

[[deps.RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "995a812c6f7edea7527bb570f0ac39d0fb15663c"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.5.1"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterface", "ChainRulesCore", "DocStringExtensions", "FillArrays", "LinearAlgebra", "RecipesBase", "Requires", "StaticArrays", "Statistics", "ZygoteRules"]
git-tree-sha1 = "f5dd036acee4462949cc10c55544cc2bee2545d6"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "2.25.1"

[[deps.RecursiveFactorization]]
deps = ["LinearAlgebra", "LoopVectorization", "Polyester", "StrideArraysCore", "TriangularSolve"]
git-tree-sha1 = "7ad4c2ef15b7aecd767b3921c0d255d39b3603ea"
uuid = "f2c3362d-daeb-58d1-803e-2bc74f2840b4"
version = "0.2.9"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Referenceables]]
deps = ["Adapt"]
git-tree-sha1 = "e681d3bfa49cd46c3c161505caddf20f0e62aaa9"
uuid = "42d2dcc6-99eb-4e98-b66c-637b7d73030e"
version = "0.1.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "cdbd3b1338c72ce29d9584fdbe9e9b70eeb5adca"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "0.1.3"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.ResettableStacks]]
deps = ["StaticArrays"]
git-tree-sha1 = "256eeeec186fa7f26f2801732774ccf277f05db9"
uuid = "ae5879a3-cd67-5da8-be7f-38c6eb64a37b"
version = "1.1.1"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "cdc1e4278e91a6ad530770ebb327f9ed83cf10c4"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.3"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMDDualNumbers]]
deps = ["ForwardDiff", "IfElse", "SLEEFPirates", "VectorizationBase"]
git-tree-sha1 = "62c2da6eb66de8bb88081d20528647140d4daa0e"
uuid = "3cdde19b-5bb0-4aaf-8931-af3e248e098b"
version = "0.1.0"

[[deps.SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[deps.SLEEFPirates]]
deps = ["IfElse", "Static", "VectorizationBase"]
git-tree-sha1 = "d4c366b135fc2e1af7a000473e08edc5afd94819"
uuid = "476501e8-09a2-5ece-8869-fb82de89a1fa"
version = "0.6.31"

[[deps.SafeTestsets]]
deps = ["Test"]
git-tree-sha1 = "36ebc5622c82eb9324005cc75e7e2cc51181d181"
uuid = "1bc83da4-3b8d-516f-aca4-4fe02f6d838f"
version = "0.0.1"

[[deps.SciMLBase]]
deps = ["ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "RecipesBase", "RecursiveArrayTools", "StaticArrays", "Statistics", "Tables", "TreeViews"]
git-tree-sha1 = "c086056df381502621dc6b5f1d1a0a1c2d0185e7"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "1.28.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "Requires"]
git-tree-sha1 = "38d88503f695eb0301479bc9b0d4320b378bafe5"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "0.8.2"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SparseDiffTools]]
deps = ["Adapt", "ArrayInterface", "Compat", "DataStructures", "FiniteDiff", "ForwardDiff", "Graphs", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays", "VertexSafeGraphs"]
git-tree-sha1 = "87efd1676d87706f4079e8e717a7a5f02b6ea1ad"
uuid = "47a9eef4-7e08-11e9-0b38-333d64bd3804"
version = "1.20.2"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "5ba658aeecaaf96923dce0da9e703bd1fe7666f9"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.4"

[[deps.SplittablesBase]]
deps = ["Setfield", "Test"]
git-tree-sha1 = "39c9f91521de844bad65049efd4f9223e7ed43f9"
uuid = "171d559e-b47b-412a-8079-5efa626c420e"
version = "0.1.14"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "7f5a513baec6f122401abfc8e9c074fdac54f6c1"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.4.1"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "6976fab022fea2ffea3d945159317556e5dad87c"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.4.2"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "c3d8ba7f3fa0625b062b82853a7d5229cb728b6b"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.2.1"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "8977b17906b0a1cc74ab2e3a05faa16cf08a8291"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.16"

[[deps.StatsFuns]]
deps = ["ChainRulesCore", "HypergeometricFunctions", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "25405d7016a47cf2bd6cd91e66f4de437fd54a07"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "0.9.16"

[[deps.SteadyStateDiffEq]]
deps = ["DiffEqBase", "DiffEqCallbacks", "LinearAlgebra", "NLsolve", "Reexport", "SciMLBase"]
git-tree-sha1 = "3e057e1f9f12d18cac32011aed9e61eef6c1c0ce"
uuid = "9672c7b4-1e72-59bd-8a11-6ac3964bc41f"
version = "1.6.6"

[[deps.StochasticDiffEq]]
deps = ["Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DiffEqJump", "DiffEqNoiseProcess", "DocStringExtensions", "FillArrays", "FiniteDiff", "ForwardDiff", "LevyArea", "LinearAlgebra", "Logging", "MuladdMacro", "NLsolve", "OrdinaryDiffEq", "Random", "RandomNumbers", "RecursiveArrayTools", "Reexport", "SparseArrays", "SparseDiffTools", "StaticArrays", "UnPack"]
git-tree-sha1 = "24d8b3ab7e91b351ccbed5e54499a1864a64a6c6"
uuid = "789caeaf-c7a9-5a7d-9973-96adeb23e2a0"
version = "6.45.0"

[[deps.StrideArraysCore]]
deps = ["ArrayInterface", "CloseOpenIntervals", "IfElse", "LayoutPointers", "ManualMemory", "Requires", "SIMDTypes", "Static", "ThreadingUtilities"]
git-tree-sha1 = "49d616ef230fec080d02ada0ca5639e652cca06b"
uuid = "7792a7ef-975c-4747-a70f-980b88e8d1da"
version = "0.2.13"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "57617b34fa34f91d536eb265df67c2d4519b8b98"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.5"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+0"

[[deps.Sundials]]
deps = ["CEnum", "DataStructures", "DiffEqBase", "Libdl", "LinearAlgebra", "Logging", "Reexport", "SparseArrays", "Sundials_jll"]
git-tree-sha1 = "e0805213754f0d871f9333eacd77862a44acb46d"
uuid = "c3572dad-4567-51f8-b174-8c6c989267f4"
version = "4.9.3"

[[deps.Sundials_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "OpenBLAS_jll", "Pkg", "SuiteSparse_jll"]
git-tree-sha1 = "04777432d74ec5bc91ca047c9e0e0fd7f81acdb6"
uuid = "fb77eaff-e24c-56d4-86b1-d163f2edb164"
version = "5.2.1+0"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "Bijections", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "IfElse", "LabelledArrays", "LinearAlgebra", "Metatheory", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "TermInterface", "TimerOutputs"]
git-tree-sha1 = "bfa211c9543f8c062143f2a48e5bcbb226fd790b"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "0.19.7"

[[deps.Symbolics]]
deps = ["ArrayInterface", "ConstructionBase", "DataStructures", "DiffRules", "Distributions", "DocStringExtensions", "DomainSets", "Groebner", "IfElse", "Latexify", "Libdl", "LinearAlgebra", "MacroTools", "Metatheory", "NaNMath", "RecipesBase", "Reexport", "Requires", "RuntimeGeneratedFunctions", "SciMLBase", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicUtils", "TermInterface", "TreeViews"]
git-tree-sha1 = "759d6102719068d95acae0b5480c157fa278ca82"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "4.3.1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "5ce79ce186cc678bbb5c5681ca3379d1ddae11a1"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.7.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.1"

[[deps.TermInterface]]
git-tree-sha1 = "7aa601f12708243987b88d1b453541a75e3d8c7a"
uuid = "8ea1fca8-c5ef-4a55-8b96-4e9afe9c9a3c"
version = "0.2.3"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.ThreadingUtilities]]
deps = ["ManualMemory"]
git-tree-sha1 = "f8629df51cab659d70d2e5618a430b4d3f37f2c3"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.5.0"

[[deps.ThreadsX]]
deps = ["ArgCheck", "BangBang", "ConstructionBase", "InitialValues", "MicroCollections", "Referenceables", "Setfield", "SplittablesBase", "Transducers"]
git-tree-sha1 = "d223de97c948636a4f34d1f84d92fd7602dc555b"
uuid = "ac1d9e8a-700a-412c-b207-f0111f4b6c0d"
version = "0.1.10"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "d60b0c96a16aaa42138d5d38ad386df672cb8bd8"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.16"

[[deps.Tokenize]]
git-tree-sha1 = "0952c9cee34988092d73a5708780b3917166a0dd"
uuid = "0796e94c-ce3b-5d07-9a54-7f471281c624"
version = "0.5.21"

[[deps.Transducers]]
deps = ["Adapt", "ArgCheck", "BangBang", "Baselet", "CompositionsBase", "DefineSingletons", "Distributed", "InitialValues", "Logging", "Markdown", "MicroCollections", "Requires", "Setfield", "SplittablesBase", "Tables"]
git-tree-sha1 = "c76399a3bbe6f5a88faa33c8f8a65aa631d95013"
uuid = "28d57a85-8fef-5791-bfe6-a80928e7c999"
version = "0.4.73"

[[deps.TreeViews]]
deps = ["Test"]
git-tree-sha1 = "8d0d7a3fe2f30d6a7f833a5f19f7c7a5b396eae6"
uuid = "a2a6695c-b41b-5b7d-aed9-dbfdeacea5d7"
version = "0.3.0"

[[deps.TriangularSolve]]
deps = ["CloseOpenIntervals", "IfElse", "LayoutPointers", "LinearAlgebra", "LoopVectorization", "Polyester", "Static", "VectorizationBase"]
git-tree-sha1 = "b8d08f55b02625770c09615d96927b3a8396925e"
uuid = "d5829a12-d9aa-46ab-831f-fb7c9ab06edf"
version = "0.1.11"

[[deps.URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["ConstructionBase", "Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "b649200e887a487468b71821e2644382699f1b0f"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.11.0"

[[deps.Unzip]]
git-tree-sha1 = "34db80951901073501137bdbc3d5a8e7bbd06670"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.1.2"

[[deps.VectorizationBase]]
deps = ["ArrayInterface", "CPUSummary", "HostCPUFeatures", "Hwloc", "IfElse", "LayoutPointers", "Libdl", "LinearAlgebra", "SIMDTypes", "Static"]
git-tree-sha1 = "1901efb08ce6c4526ddf7fdfa9181dc3593fe6a2"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.21.25"

[[deps.VertexSafeGraphs]]
deps = ["Graphs"]
git-tree-sha1 = "8351f8d73d7e880bfc042a8b6922684ebeafb35c"
uuid = "19fa3120-7c27-5ec5-8db8-b0b0aa330d6f"
version = "0.2.0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e45044cd873ded54b6a5bac0eb5c971392cf1927"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.2+0"

[[deps.ZygoteRules]]
deps = ["MacroTools"]
git-tree-sha1 = "8c1a8e4dfacb1fd631745552c8db35d0deb09ea0"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.2"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "ece2350174195bb31de1a63bea3a41ae1aa593b6"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "0.9.1+5"
"""

# ╔═╡ Cell order:
# ╟─507bfc5b-2f6b-43ea-be67-f224aa842afa
# ╟─8e41bbfa-7241-44c3-bc2d-1a018549693e
# ╟─9914fc0f-4fd3-416d-8c3f-9b1c26999b20
# ╟─057ba01c-7f78-456d-a92e-9c57bf183add
# ╟─5ed95079-b774-4cbe-8c00-1a3c1388539c
# ╟─e3979b13-75a1-4a6c-9e77-40bd01e89404
# ╟─fee3f3ef-a9a2-4bb7-965d-46c0ad7d289c
# ╟─96208e5e-2937-41b1-a270-f3f06e1e7ad9
# ╟─a40571e8-2112-420f-94a0-dc31706e3e08
# ╟─e02fd130-545e-4983-b12f-14eddb6f4bf3
# ╟─585ed61f-8c22-4deb-91b3-7dac9789fe25
# ╟─0984844f-868a-4c58-9650-cddb4bd75bba
# ╟─a28457ea-be99-4cb8-a87e-797f5fee7421
# ╟─cf486828-89ef-4698-b658-4205e847ede1
# ╟─2f66f76c-01fb-4176-93bc-3de8d5831e8f
# ╟─dac8f83a-bab8-4ead-9a64-066218e77d5b
# ╟─a1dab4c7-2629-4427-8bff-e380dfae5bc8
# ╟─3812d970-60fb-4fa7-810f-0f92dc4eb28b
# ╟─55434867-6774-4863-926c-990c9deac192
# ╟─20bbb617-1d2f-4210-8a3b-19a1631ba70d
# ╟─b87efce9-71a9-4407-9d77-bd029e6aa557
# ╟─ddcc3158-a0f7-4d56-80d2-1b6d53a68d96
# ╟─dc68a2ab-773b-40bd-a5d3-dfdd2c987489
# ╟─30afb560-a44a-11ec-1ebf-6fbfadad9353
# ╟─8fca7a96-b7d5-4bac-86f4-400f260db607
# ╟─3fe78776-a90a-457a-9176-7a495d0d3e30
# ╠═3f4053b5-f0ed-4b6b-8a13-eb7b039a8fad
# ╠═8801582a-dabd-4444-832c-8cc137225396
# ╟─435a8401-1c27-4a85-bb78-ae8a9c39eeb4
# ╟─a4034df8-7185-41dd-8c4e-a46e4ac6028f
# ╠═75703719-6fc1-487b-b0c9-14ba9e339239
# ╠═3021edfe-198c-49d0-8006-667b83c568c5
# ╟─a04787f9-8262-43dc-921f-be23863a019c
# ╟─7eb42169-410b-4f1a-bc4c-f9230cd3a904
# ╟─c4e8b093-aed0-4d49-8cbd-b57f950a29fa
# ╟─75e15198-e559-4201-a292-0d5e85898987
# ╟─feeb275b-f7f2-4d10-9449-48657dc3599a
# ╟─25aa3755-d698-4e01-beff-90567b54fa62
# ╟─3ae1813b-4741-44b9-9e96-0ddb8246d57e
# ╟─5847adf5-8444-4563-9562-81de6ee87d7b
# ╟─708b2261-19ca-4d90-839a-65eedf9458e3
# ╟─fac136db-3cfa-47d3-8329-b21c2cca74bb
# ╟─8f1a6377-1894-4c0b-8dfe-727c3d3c1c7e
# ╟─1215f027-d12a-4fc1-ab1c-82c27fd75d95
# ╟─78f83d23-863f-4113-a2ba-07282397e9f2
# ╟─f46dd83a-acb9-4e3f-88d9-5ab068478df7
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
