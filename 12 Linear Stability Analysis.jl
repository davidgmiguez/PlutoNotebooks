### A Pluto.jl notebook ###
# v0.19.24

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

# ╔═╡ 2e6401b9-9c89-4b4b-8492-d0a83003579b
using Plots, PlutoUI, ParameterizedFunctions, DifferentialEquations

# ╔═╡ 8ba695de-0650-407c-864a-c394c9c7d3ed
html"<button onclick='present()'>present</button>"

# ╔═╡ 0102b1ee-7f41-4400-aef3-c33e18b5b716
begin
	struct Foldable{C}
	    title::String
	    content::C
	end
	
	function Base.show(io, mime::MIME"text/html", fld::Foldable)
	    write(io,"<details><summary>$(fld.title)</summary><p>")
	    show(io, mime, fld.content)
	    write(io,"</p></details>")
	end
end

# ╔═╡ df36d23e-3760-48c1-8c0c-d52f15e7f520
md" ## Equilibrium And Stability 

Differential equations derived for these systems of interacting entities quickly become impossible to solve analitically. Therefore we have to use other mathematical tools to study them. The most common in the context of nonlineat dynamical systems is called 'Linear Stability Analysis'. 
"

# ╔═╡ 0b0440bf-4f75-47bd-bbd8-315066e79c75
Resource("https://www.researchgate.net/publication/346614443/figure/fig11/AS:964859927728140@1607051941592/Classification-of-fixed-points-for-two-dimensional-nonlinear-systems-As-a-function-of.png")

# ╔═╡ 279eaa4f-41ca-4168-b492-f36543bb8204
md"
##

Linear stability analysis is a method that allows us to study how a system behaves near an equilibrium point. It will help us to know if equilibrium is stable or unstable, and the bifurcations that occur in these equilibrium points, based only on a simplified linearized version of the system of study:

The model for Logistic Growth is simple enough that it can be solved analytically, but let's illustrate how linear stability analysis works this with the logistic model,.

```math
\begin{align*}
\frac{\mathrm{d} N}{\mathrm{d} t}=\mu N(1 -\frac{N}{K}) \tag{17}
\end{align*}
```
##

The main equation can be written in the following generic form:

```math
\begin{eqnarray}
\frac{\partial N}{\partial t} = f(N) 
\end{eqnarray}
```

where N is a variable and $f(N)$  is a functions which governs its temporal evolution. 
##
The first step is to calculate the fixed points of Eq. \ref{base1}. In continuous systems, the steady states occur when there is no change in the amount of our quantity $N$

```math
\begin{eqnarray}
\frac{\partial \overline{N}}{\partial t} = 0 = f(\overline{N})
\end{eqnarray}
```

Where we denote $\overline{N}$ as the value of our variable in steady state."

# ╔═╡ 00037d0a-8be0-4c75-94bc-23c2c639e6e8
md" 
##
Let's find the equilibrium points of our logistic model. 

```math
\begin{align*}
\frac{\mathrm{d} N}{\mathrm{d} t}= \mu N(1 -\frac{N}{K}) \tag{17}
\end{align*}
```

##
By taking the derivate to zero:

```math
\begin{align*}
\frac{\mathrm{d} \overline{N}}{\mathrm{d} t}=0= \mu \overline{N}(1 -\frac{\overline{N}}{K}) \tag{17}
\end{align*}
```
##
 we can easily see that we have two equilibrium states: 
```math
\begin{align*}
\overline{N}=0\\
1 -\frac{\overline{N}}{K} = 0 \Rightarrow \overline{N} = K
\end{align*}
```


"

# ╔═╡ 64dc539b-b268-43b6-a96e-2c4fa7b48474
md" The next step after finding the equilibrium points is to check if they are stable or unstable, or what type of equilibrium we have. We do that by introducing small perturbations.
##
Let’s take a fixed point $\overline{N}$ and perturb it by an infinitesimal amount n. We are interested in the dynamics of N = $\overline{N}$ + n and whether N will move away away or towards $\overline{N}$ as time progresses. 

```math
\begin{align*}
\frac{\mathrm{d} N}{\mathrm{d} t}= \frac{\mathrm{d} \overline{N} + n}{\mathrm{d} t}= \frac{\mathrm{d} \overline{N} }{\mathrm{d} t} + \frac{\mathrm{d} n}{\mathrm{d} t}
\end{align*}
```
##
and since 

```math
\begin{align*}
\frac{\mathrm{d} \overline{N} }{\mathrm{d} t} = 0 
\end{align*}
```

##
by definition, we have that the change in time of the population is equivalent to the change in the time of the small perturbation. 

```math
\begin{align*}
\frac{\mathrm{d} N}{\mathrm{d} t}= \frac{\mathrm{d} n}{\mathrm{d} t}
\end{align*}
```
##
Since $n$ is very small by definition, we can linearize the dynamics around the ficxed points $\overline{N}$, using a Taylor expansion
. 

```math
\begin{align*}
\frac{\mathrm{d} N}{\mathrm{d} t}= \frac{\mathrm{d} n}{\mathrm{d} t} = f(\overline{N}+n)= f(\overline{N})+ f'(\overline{N})·n +···
\end{align*}
```

being $f'(\overline{N})$ the value of derivative of function $f(N)$ in the point of equilibrium $\overline{N}$
##
Since $f(\overline{N})$ = 0 by definition, we get  

```math
\begin{align*}
 \frac{\mathrm{d} n}{\mathrm{d} t} = f'(\overline{N})·n +···
\end{align*}
```
"

# ╔═╡ f5faefd6-f496-40fb-8043-ef6b9d172ae0
md"

If the perturbation $n$ is very small, then linear and nonlinear evolution are in fact approximately the same. But as $n$ increases in size the nonlinear effects become increasingly more important and evolving y with the linearized dynamics or the full nonlinear dynamics is no longer equivalent. 

So, this means that if perturbations $n$ are small, we can forget about higher order terms and simply assume that this has a form of $\dot{n}= \lambda \cdot n$. 
##
So if we integrate we have a solution that is an exponential fucntion. 


```math
\begin{align*}
n(t)=n(0) e^{\lambda t}
\end{align*}
```
##
So, depending on the sign of $\lambda$, this perturbation $n$ may increase ($\lambda >0$) or decrease ($\lambda <0$). As an example, let's test this for the logistic system. 

```math
\begin{align*}
\frac{\mathrm{d} N}{\mathrm{d} t}=\mu N(1 -\frac{N}{K}) = \frac{\mathrm{d} (\mu \overline{N}(1 -\frac{\overline{N}}{K})}{\mathrm{d} t} n \tag{17}
\end{align*}
```
##
For simplicity, we will decompose the fucntion in to two terms:
```math
\begin{align*}
\mu N(1 -\frac{N}{K}) = g(N) \cdot h(N)
\end{align*}
```
##
being 


```math
\begin{align*}
g(N) &= \mu N \\
h(N) &= 1 -\frac{N}{K}
\end{align*}
```
##
now we need the detivative 
```math
\begin{align*}
\frac{\mathrm{d} [g(N) \cdot h(N)]}{\mathrm{d} t}= \frac{\mathrm{d} g(N) }{\mathrm{d} t}  h(N) + g(N) \frac{\mathrm{d} h(N) }{\mathrm{d} t} 
\end{align*}
```
##
since 
```math
\begin{align*}
\frac{\mathrm{d} g(N) }{\mathrm{d} t} &=\mu\\
\frac{\mathrm{d} h(N) }{\mathrm{d} t} &= -\frac{1}{K}
\end{align*}
```
##
so,

```math
\begin{align*}
\frac{\mathrm{d} (\mu N(1 -\frac{N}{K})}{\mathrm{d} t} =\mu (1 -\frac{N}{K}) - \frac{\mu \cdot N}{K}
\end{align*}
```

"

# ╔═╡ 32b3bc77-c4f4-40ab-a5a9-e035022189ac
md"
##
now we substitute our steady state values, for $\overline{N}=0$, we have 

```math
\begin{align*}
\frac{\mathrm{d} (\mu \overline{N}(1 -\frac{\overline{N}}{K})}{\mathrm{d} t}\Biggr\rvert_{\overline{N}=0} =\mu (1 -\frac{0}{K}) - \frac{\mu \cdot 0}{K} = \mu
\end{align*}
```
##
So, as long as $\mu > 0$ this fixed point is unestable. If we have $N=0$, the systems remains there. As soon as we perturb the number (and we can only perturb by slightly increasing, these perturbations grow exponentially. For the second steady state. 

```math
\begin{align*}
\frac{\mathrm{d} (\mu \overline{N}(1 -\frac{\overline{N}}{K})}{\mathrm{d} t}\Biggr\rvert_{\overline{N}=K} =\mu (1 -\frac{K}{K}) - \frac{\mu \cdot K}{K} = - \mu
\end{align*}
```
##
So, as long as $\mu > 0$ this fixed point is stable. Perturbations will allways decrease exponentially. 

To observe this grafically, lets plot the function $\mu N(1 -\frac{N}{K})$


"

# ╔═╡ c70c7a06-d1aa-40a9-902b-4bd4b23a6392
begin
	T_slide = @bind T html"<input type=range min=5 max=15 step=1>"
	K_slide = @bind K html"<input type=range min=100 max=500 step=10>"
	md"""
	**Set the Cell Cycle Length and the carrying capacity**
	
	value of T: $(T_slide) 
	value of K: $(K_slide)
	
	"""
end

# ╔═╡ 2cfd31c4-803a-4880-b2a4-5af6594c0975
begin
	N=collect(0:0.1:K);
	plot(N, log(2)/T .* N .* (1 .- (N ./ K)),label="T= $T, K= $K",seriestype=:line,xlabel=("N"),ylabel=("f(N)"),ylims = (0,20),xlims = (0,500))
end

# ╔═╡ 2170f883-1b01-4d16-907e-49c72118e224
md"
##
Equilibrium points are the values where the function is zero. We see the two of them, the unsable $\overline{N}=0$ and the stable $\overline{N}=K$. In this system, every perturbation moves away from $\overline{N}=0$ towards $\overline{N}=K$, which is the carrying capacity of the system.

The slope represents how fast the change occurs, so increasing $T$, means that we reach the carrying capacity of the system faster. 
##
So, this is a liner model, and we can solve it analitically, so performinng a perturbation analysis does not provide extra information (we now the full dynamics because we have an analitical solution). The advantage of this perturbation analysis (or linear stabilty analysis) is when we work with systems that cannot be solved anallyically, such as systems with multiple variables, and systems with nonlinearities. 

"

# ╔═╡ 7d89a7bd-baba-42c1-af4f-26f74138199a
md"## Stability analysis for the logistic growth with Allee effect

```math
\begin{align*}
\frac{\partial N}{\partial t}= \mu N \left(\frac{N}{A} -1\right) \left(1-\frac{N}{K}\right) 
\end{align*}
```



"

# ╔═╡ 9db950f3-6989-4225-845f-dae93d52a9b9
begin
	TT_slide = @bind TT html"<input type=range min=5 max=15 step=1>"
		KK_slide = @bind KK html"<input type=range min=100 max=500 step=10>"
		AA_slide = @bind AA html"<input type=range min=5 max=100 step=10>"
	
		md"""
		**Set the parameters**
		
		value of T: $(TT_slide) 
		value of K: $(KK_slide)
		value of A: $(AA_slide)
		
		"""
end

# ╔═╡ c620fb38-6e3d-4d54-b515-a9a5fc4360e4
begin
	#N=collect(0:0.1:K);
	plot(N, log(2)/TT  .* N .* (1 .- (N ./ KK)) .* ((N ./ AA) .-1 ),label="T= $TT, K= $KK, A = $AA",seriestype=:line,xlabel=("N"),ylabel=("f(N)"),ylims = (-5,10),xlims = (0,200))
end

# ╔═╡ 66ea1e28-b50a-43fc-acf5-9a146b51d505
md" 

```math
\begin{align*}
 \mu N \left(\frac{N}{A} -1 \right) \left(1-\frac{N}{K}\right) =0
\end{align*}
```

We have three fixed points, or points where the derivative is zero N=0, N=A and N=K. To calculate the stability we need the derivative  
```math
\begin{align*}
\frac{-A \cdot K + 2 \cdot A \cdot N + 2 \cdot K \cdot N - 3 \cdot N^2}{A \cdot K}
\end{align*}
```
"

# ╔═╡ b22fd4f1-ee14-4cb4-aff9-94153f470fe1
begin
		#N=collect(0:0.1:K);
		plot(N, 1 ./ 10 .* log(2)/TT  .* N .* (1 .- (N ./ KK)) .* ((N ./ AA) .-1 ),label="T= $T, K= $K",seriestype=:line,xlabel=("N"),ylabel=("f(N)"),ylims = (-5,10),xlims = (0,200))
		plot!(N, log(2)/TT  .* (-AA .* KK .+ 2 .* AA .* N .+ 2 .* KK .* N .- 3 .*N .^2) ./(AA .* KK),label="T= $TT, K= $KK",seriestype=:line,xlabel=("N"),ylabel=("f(N)"),ylims = (-1,1),xlims = (0,200))
end

# ╔═╡ 555aa5d3-9e29-459a-8798-51a576d252b4
md"At the firts fixed point $N=0$, wich means the extintion. so it is stable because the derivative is unstable.  At the secont fixed point, N=A is unstable, the derivative is positive at this point. The third one N=K is unstable again. "

# ╔═╡ 8c809d92-4edb-40ca-ac80-18b7e8c72fae
md"## 
> Exercise: Find the equilibrium and the stability of the logistic growth model with Allee effect and including a term that represents a constant loss of individuals due to harvesting. 
> 
>
>```math
>\begin{align*}
>\frac{\partial N}{\partial t}= \mu N \left(\frac{N}{A} -1\right) \left(1-\frac{N}{K}\right) - B
>\end{align*}
>```
> where B is the harvesting rate


"

# ╔═╡ 6c87979a-690a-400b-9ab8-7ef26b829195
md"## Stability analysis of nonlinear systems: Lotka-Volterra model

The Lotka-Volterra model, also known as Predator-Prey model describes the inteactions between two populations of species where one feeds into the other.  One can think of rabbits $x$ and foxes $y$, such that rabits multiply where there is no foxes (assuming an infinite amount of food for the rabbits). This is basically a first order production of rabbits, an autocalalitic system that results in exponetial growth of rabbits.
##

Next, foxes $y$ feed on the rabbits and multiply due to the good food conditions. Then also foxes die at a constant rate (rabbits also die, but the mdoel assumes that their rate of birth is much higher than the rate of death). The scheme of interactions for this very simple system is the following:

```math
\begin{align}
 x &\overset{k_1}{\longrightarrow} 2 x   \\
 x + y &\overset{k_2}{\longrightarrow} 2 y \\
 y &\overset{k_3}{\longrightarrow} 0 
 \end{align}
```
##
First, we find the differential equations that govern the dynamcis of the following system, and the steady state values. we start by writting the matrices of stoichiometric coefficients:

```math
A=\begin{bmatrix}
  1 & 0  \\  
  1 & 1   \\
  0 & 1       \end{bmatrix} ;
B=\begin{bmatrix}
  2 & 0   \\ 
  0 & 2   \\
  0 & 0       \end{bmatrix}; \tag{8}
```

"

# ╔═╡ 04c14475-24cd-49c3-8a7b-83a9425664be
A = [1 0;1 1;0 1];B = [2 0;0 2;0 0];(B-A)'

# ╔═╡ 919bf3e4-6aa0-440b-b76f-34a4812c6752
md"##
in this particular case

```math
K=\begin{pmatrix}
 k_1 & 0 & 0  \tag{9}\\ 
 0 &  k_2 & 0  \\ 
 0 &  0 & k_3 
\end{pmatrix}
```
##
and 


```math
X^A=\begin{pmatrix}
X_1^1\cdot X_2^0  \\
X_1^1\cdot X_2^1  \\
X_1^0\cdot X_2^1   
\end{pmatrix} = \begin{pmatrix}
 X_1 \\
 X_1 \cdot X_2\\
 X_2
\end{pmatrix} \tag{10}
```
##
so, the equations that define the system are

```math
\begin{align}
 \begin{bmatrix}
\frac{\mathrm{d} X_1}{\mathrm{d} t}\\ \frac{\mathrm{d} X_2}{\mathrm{d} t} \end{bmatrix}& 
=  \begin{bmatrix} 1 & -1 & 0  \\ 0 & 1 &  -1  \end{bmatrix}
\begin{pmatrix}
 k_1 & 0 & 0  \tag{9}\\ 
 0 &  k_2 & 0  \\ 
 0 &  0 & k_3 
\end{pmatrix}
 \begin{pmatrix}
 X_1 \\
 X_1 \cdot X_2\\
 X_2
\end{pmatrix} 
\end{align}
```
##
and after multiplying the matrices
```math
\begin{align}
 \begin{bmatrix}
\frac{\mathrm{d} X_1}{\mathrm{d} t}\\ \frac{\mathrm{d} X_2}{\mathrm{d} t} \end{bmatrix}& 
=  \begin{bmatrix} 1 & -1 & 0  \\ 0 & 1 &  -1  \end{bmatrix}
 \begin{pmatrix}
 k_1 \cdot X_1 \\
 k_2 \cdot X_1 \cdot X_2\\
 k_3 \cdot X_2
\end{pmatrix} 
\end{align}
```
##
so finally, 

```math
\begin{align}
 \begin{bmatrix}
\frac{\mathrm{d} X_1}{\mathrm{d} t}\\ \frac{\mathrm{d} X_2}{\mathrm{d} t} 
\end{bmatrix}&= 
 \begin{pmatrix}
   k_1 \cdot X_1 - k_2 \cdot X_1 \cdot X_2 \\
  k_2 \cdot X_1 \cdot X_2 - k_3 \cdot X_2
\end{pmatrix} \tag{11}
\end{align}
```
##
Therefore, the  equations for the evolution of `[x]` and `[y]` are as follows:

```math
\begin{align}       
            \frac{ dx }{dt} &=  k_1 \cdot x - k_2 \cdot x \cdot y  \tag{5}\\ 
            \frac{ dy }{dt} &= k_2 \cdot x \cdot y - k_3 \cdot y  \tag{6} 
            \end{align} 
```


"

# ╔═╡ 855b8530-861e-427e-b6cc-c1b0283568ca
md" 
##
which in general form, we can write as:

```math
\begin{eqnarray}
\frac{\partial x}{\partial t} = f(x, y) \\
\frac{\partial y}{\partial t} = g(x, y) 
\end{eqnarray}
```
##
where $f(x, y)$ and $g(x, y)$ are nonlinear equations that govern the temporal evolution and couple the behavior of the two variables $x$ and $y$:

Next, we need to calculate the fixed points:

```math
\begin{eqnarray}
\frac{\partial \overline{x}}{\partial t} =0= f(\overline{x},\overline{y})  \\
\frac{\partial \overline{y}}{\partial t}= 0 = g(\overline{x},\overline{y})
\end{eqnarray}
```
##
For the particular case of the Lotcka-Volterra, we just set eqs. 5 and 6 to zero



```math
\begin{align}       
            k_1 \cdot \overline{x} - k_2 \cdot \overline{x} \cdot \overline{y}  \tag{5} &= 0\\ 
            k_2 \cdot \overline{x} \cdot \overline{y} - k_3 \cdot \overline{y}  \tag{6} &= 0
            \end{align} 
```
##
and solve for `x` and `y`. We obtain two solutions,  $\overline{x}=\overline{y}=0$ and 


```math
\begin{align}       
            \overline{x} &= \frac{k_3}{k_2}\tag{5} \\ 
            \overline{y} &= \frac{k_1}{k_2} \tag{6}   
            \end{align}    
```


	"

# ╔═╡ c112e62c-9c38-48dc-8294-42911b02ccc5
md" ##
The next step is to find the characteristics of the steady states. We do that by following the same rationalle of the logistic model, i.e., to expand our equations as Taylor series around the steady state $(\overline{x},\overline{y})$., but now for multi-variable equations. 

```math
\begin{eqnarray}
\frac{\partial x}{\partial t} = M_{11} \cdot x + M_{12} \cdot y + ... \\
\frac{\partial y}{\partial t} = M_{21} \cdot x + M_{22} \cdot y + ... 
\end{eqnarray}
```
##            
Where $M_{ij}$ are the components of the Jacobian matrix, evaluated at the steady state $(\overline{x},\overline{y})$. 

```math
\begin{align}
 J=\begin{bmatrix} 
 M_{11} & M_{12} \\ 
 M_{21} & M_{22}
 \end{bmatrix}_{\overline{x},\overline{y}}= \begin{bmatrix} 
\frac{\partial  f(x,y)}{\partial x}\Biggr\rvert_{\overline{x},\overline{y}} & \frac{\partial  f(x,y)}{\partial y}\Biggr\rvert_{\overline{x},\overline{y}} \\ 
\frac{\partial  g(x,y)}{\partial x}\Biggr\rvert_{\overline{x},\overline{y}} & \frac{\partial  g(x,y)}{\partial y}\Biggr\rvert_{\overline{x},\overline{y}}
 \end{bmatrix} 
 \end{align} 
```
##
which for the particular case of the Lotka-Volterra is 


```math
\begin{align}
 J=\begin{bmatrix} 
 M_{11} & M_{12} \\ 
 M_{21} & M_{22}
 \end{bmatrix}_{\overline{x},\overline{y}} = \begin{bmatrix} 
 k_1 - k_2 \cdot \overline{y}  &  - k_2 \cdot \overline{x}   \\ 
k_2  \cdot \overline{y}  & k_2 \cdot \overline{x}  - k_3 
 \end{bmatrix}
 \end{align} 
```
##
Next, to investigate the stability, we check solutions in the form of small perturbations as follows:

```math
\begin{eqnarray}
(x,y) = (X_0,Y_0) e^{\lambda t}  
\end{eqnarray}
```
Here, $\lambda$ is the growth rate of the perturbations, also refered as eigenvalue. Each steady state will behave differently in terms of the dynamcis of the perturbations. Therefore, each steady state will have an associated eigenvalue. To find the eigen values, we solve the characteristic polynomial $det[J-\lambda I]=0$. 
##
```math
\begin{eqnarray}
Det \left(\begin{array}{cc}M_{11}-\lambda & M_{12} \\M_{21}& M_{22}-\lambda \end{array}\right) =0 
\end{eqnarray}
```
##
which gives us the corresponding equation:
```math
\begin{eqnarray}
\lambda^{2} -\lambda Tr(M) + Det(M)=0
\end{eqnarray}
```
##
where:
```math
\begin{eqnarray}
Tr(M)= M_{11}+M_{22} \\
Det(M)= M_{11}M_{22}-M_{12}M_{21} 
\end{eqnarray}
```

"

# ╔═╡ 6c96da07-b5ea-4096-a14c-5cd8f0f47c04
md"##
for the lotka-volterra case:

```math
\begin{eqnarray}
Tr(M)= k_{1} + k_2 (\overline{x}-\overline{y}) - k_3  \\
Det(M)= (k_1- k_2 \cdot \overline{y})(k_2 \cdot \overline{x} - k_3) - (k_2 \cdot\overline{y})(-k_2 \cdot \overline{x})  )
\end{eqnarray}
```
##
calculating 

```math
\begin{eqnarray}
Tr(M)= k_{1} + k_2 (\overline{x}-\overline{y}) - k_3  \\
Det(M)= - k_1 \cdot k_3 - k_2^2 \cdot \overline{y} \cdot \overline{x} + k_2  \cdot k_3 \cdot \overline{y} + k_1 \cdot k_2 \cdot \overline{x}+ k_2^2 \cdot\overline{y} \cdot \overline{x} 
\end{eqnarray}
```
##
and 

```math
\begin{eqnarray}
Tr(M)= k_{1} + k_2 (\overline{x}-\overline{y}) - k_3  \\
Det(M)= k_2  \cdot k_3 \cdot \overline{y} + k_1 \cdot k_2 \cdot \overline{x} - k_1 \cdot k_3 
\end{eqnarray}
```
##
so the chracteristic equation becomes:
```math
\begin{eqnarray}
\lambda^{2} -\lambda (k_{1} + k_2 (\overline{x}-\overline{y}) - k_3) + k_2  \cdot k_3 \cdot \overline{y} + k_1 \cdot k_2 \cdot \overline{x} - k_1 \cdot k_3=0
\end{eqnarray}
```
##
So, the polynomium for the fixed point $\overline{x},\overline{y}=[0,0]$ is 

```math
\begin{eqnarray}
\lambda^{2} +\lambda (k_3 - k_{1})  - k_1 \cdot k_3=0
\end{eqnarray}
```


"

# ╔═╡ 532acbb6-e504-4deb-88d5-e001a476af15
md"

```math
\begin{eqnarray}
\lambda &=& \frac{-(k_3 - k_{1})\pm \sqrt{(k_3 - k_{1})^2 + 4 \cdot k_1 \cdot k_3}}{2} \\
\lambda &=& \frac{-k_3 + k_{1}\pm \sqrt{k_3^2 + k_{1}^2 +2 \cdot k_1 \cdot k_3}}{2} \\ 
\lambda &=& \frac{-k_3 + k_{1}\pm \sqrt{(k_3 + k_{1})^2}}{2}  \\
\lambda &=& \frac{-k_3 + k_{1}\pm (k_3 + k_{1})}{2}  \\
\end{eqnarray}
```
so the two solutions are
```math
\begin{eqnarray}
\lambda_1 &=&  \frac{-k_3 + k_{1} + k_3 + k_{1}}{2} = k_{1} \\
\lambda_2 &=&  \frac{-k_3 + k_{1} - k_3 - k_{1}}{2} = -k_{3} \\
\end{eqnarray}
```


"

# ╔═╡ 906a0931-e86f-48f3-be11-eb692639e8cc
md" ##
Again, since $\lambda$ is the exponent that sets the dynamcis of the perturbations, depending on its value, the steady state is stable or unstable. 

since $k_1$ and $k_3$ are allways positive, we have one positive eigenvalue and one negative eigenvalue. 

For our Locka Volterra case, we can evaluate first the eigenvalues for the first steady state. It means that perturbations in one variable grow while perturbations in the other variable decay.  


"

# ╔═╡ c374e1bd-ec3d-4cf7-aff1-072d6ba60b27
function quadratic(a, b, c)
          discr = b^2 - 4*a*c
          discr >= 0 ?   ( (-b + sqrt(discr))/(2a), (-b - sqrt(discr))/(2a) ) : error("Only complex roots")
        end

# ╔═╡ 2b91a06c-3803-4877-ac9a-1669fbe58d30
begin
	k1_slide = @bind k1 html"<input type=range min=1 max=5 step=.1>"
	k2_slide = @bind k2 html"<input type=range min=1 max=5 step=.1>"
	k3_slide = @bind k3 html"<input type=range min=1 max=5 step=.1>"
	md"""
	**Set the values of the kinetic constants**
	
	value of k1: $(k1_slide)
	value of k2: $(k2_slide)
	value of k3: $(k3_slide)
	
	"""
end

# ╔═╡ e8ca5682-cce0-4b2c-ba18-47f2e38192c7
begin
	a= 1
	b= k3 - k1 
	c= - k1 * k3
	quadratic(a,b,c)
end

# ╔═╡ 17e34ce6-a1fa-4193-9201-3e61a188b48f
md"  
##
The stability of this fixed point [0,0] is of importance. If it both ewigenvalues are negative, the point will be stable, and non-zero populations might be attracted towards it, and as such the dynamics of the system might lead towards the extinction of both species for many cases of initial population levels.

##
However, as the steady state at the origin is unstable in one of the variables, we find that the extinction of both species is difficult in the model. In fact, in teh absence of foxes and rabits, a small increase in the foxes will lead to extintion (no food), while a small increase in the amount of rabbits will read to exponential increase (the unstable branch). 

These type of points are called a saddle node . 

##
For the other solution $[\overline{x},\overline{y}]=[\frac{k_3}{k_2},\frac{k_1}{k_2}]$

so the chracteristic equation becomes:
```math
\begin{eqnarray}
\lambda^{2} -\lambda (k_{1} + k_2 (\frac{k_3}{k_2}-\frac{k_1}{k_2}) - k_3) + k_2  \cdot k_3 \cdot \frac{k_1}{k_2} + k_1 \cdot k_2 \cdot \frac{k_3}{k_2} - k_1 \cdot k_3=0
\end{eqnarray}
```
##
```math
\begin{eqnarray}
\lambda^{2}  + k_3 \cdot k_1 =0
\end{eqnarray}
```

so, solving 

```math
\begin{eqnarray}
\lambda =  \pm \sqrt{- k_3 \cdot k_1} 
\end{eqnarray}
```

so, the two eigenvalues are:

```math
\begin{eqnarray}
\lambda_1 = + i \cdot k_3 \cdot k_1 \\
\lambda_2 = - i \cdot k_3 \cdot k_1
\end{eqnarray}
```

"

# ╔═╡ 87df7c88-7e16-4c50-bd80-df5bfea26ee8
begin
	
	aa=1
	bb= 0
	cc=- k1  + k2 *( k1 * k3)
	quadratic(aa,bb,cc)
end

# ╔═╡ 9433f4d1-3abf-4dea-8107-1070ee5feb2e
md"##
The two values are purely imaginary so we cannot say much about the stability. A small perturbation will not experience repulsion or atraction towards this steady state. There is no stable state (no atractor), and trajectories circulate about the fixed point in a stable orbit. This is called a _center_. 

The solutions travel periodically around the level sets in the counterclockwise direction
##
To test this, we solve numerically the system 

We assume as initial conditions:

```math
\begin{align}       
            x (0) &= 1 \tag{7} \\ 
            y (0) &= 1   \tag{8}   
\end{align}            
```
       
"

# ╔═╡ 2b9e43fd-d00d-40f1-ac98-07d48c53d861
md" ##
Similarly to what we did in the previous case, we would try to see the fixed points graphically. To do that in two dimensional systems, we find the functions where $\dot{x} = 0$ and $\dot{y} = 0$. These lines will represent the boundaries  between increase and decrease in $x$ and $y$. 
##
These curves are called the  nullclines. The method of nullclines is a technique for determining the global behavior of solutions of competing species models. This method provides an effective means of finding trapping regions for some differential equations. In a competition model, if a species population x is above a certain level, the fact of limited resources will cause x to decrease. 
##
Let's illustrate this again with the Lotka-Volterra. The functions that satisfy that the defivatives of $x$ and $y$ are zero are:

```math
\begin{align}       
            k_1 \cdot \overline{x} - k_2 \cdot \overline{x} \cdot \overline{y}  \tag{5} &= 0\\ 
            k_2 \cdot \overline{x} \cdot \overline{y} - k_3 \cdot \overline{y}  \tag{6} &= 0
            \end{align} 
```
##
In this particuular case, the lines are very simple, just constant values. 

```math
\begin{align}       
            \overline{x} &= \frac{k_3}{k_2}\tag{5} \\ 
            \overline{y} &= \frac{k_1}{k_2} \tag{6}   
            \end{align}    
```

"

# ╔═╡ 3349c942-15a8-4b1e-8b4a-7f5154f48b12
lv! = @ode_def LotkaVolterra begin
  dx = k1*x - k2*x*y
  dy =  k2*x*y - k3*y 
    end k1 k2 k3

# ╔═╡ c3828c90-0559-4686-9f56-41b5b6d48176
begin
	x₀=1
	y₀=1
	t₀=0.0
	final_time=10.0;
	prob = ODEProblem(lv!,[x₀,y₀],(t₀,final_time),(k1,k2,k3))
	sol = solve(prob)
	plot(sol,ylims = (0, 5))

	title!("Lotka-Volterra ")
	xlabel!("time [a.u.]")
    ylabel!("Amplitude [a.u.]")
	
end

# ╔═╡ ff0c4a06-de7b-4cf2-b920-cd84ac9927ea
[u[1] for (u,t) in tuples(sol)]


# ╔═╡ f8238a6e-1e5e-4629-ac20-2c9cb964c53e
tuples(sol)

# ╔═╡ df81eb5c-16ef-4bf6-9ad0-2cc899c44cb6
begin
	vline([k3/k2],ylims = (0, 5),xlims = (0, 10));
	hline!([k1/k2],ylims = (0, 5),xlims = (0, 10));
	title!("Null-Clines of the Lotka Volterra ")
	xlabel!("x [a.u.]")
    ylabel!("y [a.u.]")
	plot!([u[1] for (u,t) in tuples(sol)],[u[2] for (u,t) in tuples(sol)],ylims = (0, 15))

prob2 = ODEProblem(lv!,[x₀*2,y₀*2],(0.0,10.0),(k1,k2,k3))
	sol2 = solve(prob2)
	plot!([u[1] for (u,t) in tuples(sol2)],[u[2] for (u,t) in tuples(sol2)],ylims = (0, 15))
	
end

# ╔═╡ a5b38e90-cdbe-442f-a078-dd8598a0a4c2


md"
##
More concretely, if $Re(\lambda) < 0$, the perturbation decays in time and the steady state ($\overline{x},\overline{y}$) is stable. 

On the contrary, when $Re(\lambda) > 0$, the perturbation grows exponentially and the steady sate is unstable. 
##
More concretely, the steady state is stable if the following conditions are fulfilled:
```math
\begin{eqnarray}
Tr(M) < 0 \\
Det(M) >0
\end{eqnarray}
```
##
We can write the eigenvalue expression separating real and imaginary part:
```math
\begin{eqnarray}
\lambda=\mu \pm i \omega  
\end{eqnarray}
```
##
where

```math
\begin{eqnarray}
\mu=\frac{1}{2} Tr(M)\\
\omega=\sqrt{-\frac{1}{4} Tr(M)^2 + Det(M)}
\end{eqnarray}
```
##
The Hopf bifurcation takes place when $Det(M) > (1 / 4) Tr(M)^2 $ and $Tr(M) > 0$. In this case the eigenvalue has nonzero imaginary part and the solution of the system is oscillatory. In the Hopf threshold ($M_{11}=-M_{22}$) the complex part of the eigenvalue becomes:

```math
\begin{eqnarray}
\omega^2= \omega_{c}^{2}=-M_{11}^{2}-M_{12}M_{21}>0 \\
M_{12}M_{21}>M_{11}^{2}  (>0) 
\end{eqnarray}
```
##
One of the values must be positive, and the other negative. We choose $M_{11}>0$ $\longrightarrow$ $M_{22}<0$ and $M_{12}>0$ $\longrightarrow$ $M_{21}<0$. This way, we can write Eq. \ref{lineal1} and  Eq. \ref{lineal2} as follows:

```math
\begin{eqnarray}
\frac{\partial x}{\partial t} = M_{11} x - |M_{12}| y + ... \\
\frac{\partial y}{\partial t} = M_{21} x - |M_{22}| y + ... 
\end{eqnarray}
```
"

# ╔═╡ 3118e37d-f7bd-488e-ab3a-74e0194d3108
md"
## Stability and Eignevalues

Eigenvalues can be used to determine if a fixed point (also known as an equilibrium point) is stable or unstable. 

The eigenvalues of a system linearized around a fixed point can determine the stability behavior of a system around the fixed point. 

The particular stability behavior depends upon the existence of real and imaginary components of the eigenvalues, along with the signs of the real components and the distinctness of their values. 

We will examine each of the possible cases below.
##
### Complex Eigenvalues: 
-  If the real part is positive, the system is unstable and behaves as an unstable oscillator. This can be visualized as a vector tracing a spiral away from the fixed point.
- If the real part is zero, the system behaves as an undamped oscillator.
- If the real part is negative, then the system is stable and behaves as a damped oscillator. 
##
### Real Eigenvalues:
-  If equal to zero, the system will be unstable. This is just a trivial case of the complex eigenvalue that has a zero part.
- When all eigenvalues are real, positive, and distinct, the system is unstable. On a gradient field, a spot on the field with multiple vectors circularly surrounding and pointing out of the same spot (a node) signifies all positive eigenvalues. This is called a source node
- When all eigenvalues are real, negative, and distinct, the system is unstable. Graphically on a gradient field, there will be a node with vectors pointing toward the fixed point. This is called a sink node.
- If the set of eigenvalues for the system has both positive and negative eigenvalues, the fixed point is an unstable saddle point. A saddle point is a point where a series of minimum and maximum points converge at one area in a gradient field, without hitting the point. It is called a saddle point because in 3 dimensional surface plot the function looks like a saddle.

"

# ╔═╡ 48dbf14c-ab77-4e6a-9503-8dc53f5dfba8
Resource("https://eng.libretexts.org/@api/deki/files/18509/image-67.jpeg?revision=1")

# ╔═╡ cf72446a-a512-11ec-2b47-ef706c91c6a0
md" ## Cubic Autocatalator model

To study the behavior of nonlinear systems, a set of mathematical tools is  commonly used. Here, we will outline its main aspects from a simplified point of view, trying to introduce the reader to the mathematics inside the nonlinear pattern formation field. In addition we will try to illustrate the problem using a very simple autocatalitic model: the _Cubic Autocatalor Model_.


```math
\begin{align}
 2u + v &\overset{k_1=1}{\longrightarrow} 3u   \\
 u  &\overset{k_2=1}{\longrightarrow} 0 \\
 0 &\overset{k_3=\mu}{\longrightarrow} v
 \end{align}
```
  
##
The main equation which governs the aspects of pattern formation systems is the following nonlinear equations:

```math
\begin{eqnarray}
\frac{\partial u}{\partial t} = f(\mu, u, v) \\
\frac{\partial v}{\partial t} = g(\mu, u, v)
\end{eqnarray}
```
##
where $u$ and $v$ correspond to the concentration of activator and inhibitor.
Here, $f(\mu, u, v)$ and $g(\mu, u, v)$ are nonlinear functions which govern the temporal evolution of the variables. 
##
The first step is to calculate the fixed points of Eq. \ref{base1} and \ref{base2}, i.e., the values of the variables where the null-clines are in coincidence and equal to zero. This defines the steady state for the variables in a zero dimensional system.

```math
\begin{eqnarray}
f(\mu, u, v)=0  \\
g(\mu, u, v)=0
\end{eqnarray}
```
##
An example of the null-clines for the \textit{Cubic Autocatalor} Model can be seen in Fig. \ref{nullclines_cubic}. The equations for this specific model are:

```math
\begin{eqnarray}
\frac{\partial u}{\partial t} = u^2 v -u \\
\frac{\partial v}{\partial t} = \mu - u^2 v 
\end{eqnarray}
```

The steady state for this model is ($u_0,v_0$)= ($\mu, 1/\mu$)."

# ╔═╡ b67ed8c0-8ca7-4928-9ebc-98da513c432f
md"
##
The following step to study the evolution of the system is to linearize Eq. \ref{base1} and \ref{base2} around the steady state $(u_{0},v_{0})$. 

```math
\begin{eqnarray}
\frac{\partial u}{\partial t} = M_{11} u + M_{12} v + f_{2}(\mu, u, v) +   f_{3}(\mu, u, v) + ... \\
\frac{\partial v}{\partial t} = M_{21} u + M_{22} v + g_{2}(\mu, u, v) + g_{3}(\mu, u, v) + ... 
\end{eqnarray}
```
##
where $M_{ij}$ is calculated in the steady state ($u_{0},v_{0}$) as follows:
```math
\begin{eqnarray}
M_{11} = \frac{\partial  f(\mu, u, v)}{\partial u}\\
M_{12} = \frac{\partial  f(\mu, u, v)}{\partial v} \\
M_{21} = \frac{\partial  g(\mu, u, v)}{\partial u}\\
M_{22} = \frac{\partial  g(\mu, u, v)}{\partial v}
\end{eqnarray}
```

##
Where $M_{ij}$ are the components of the Jacobian matrix, evaluated at the steady state $(\overline{x},\overline{y})$. 

```math
\begin{align}
 J=\begin{bmatrix} 
 M_{11} & M_{12} \\ 
 M_{21} & M_{22}
 \end{bmatrix}_{\overline{x},\overline{y}}= \begin{bmatrix} 
\frac{\partial  f(x,y)}{\partial x}\Biggr\rvert_{\overline{x},\overline{y}} & \frac{\partial  f(x,y)}{\partial y}\Biggr\rvert_{\overline{x},\overline{y}} \\ 
\frac{\partial  g(x,y)}{\partial x}\Biggr\rvert_{\overline{x},\overline{y}} & \frac{\partial  g(x,y)}{\partial y}\Biggr\rvert_{\overline{x},\overline{y}}
 \end{bmatrix} 
 \end{align} 
```
##
which for the particular case of the Cubic autocatalator is 


```math
\begin{align}
 J=\begin{bmatrix} 
 M_{11} & M_{12} \\ 
 M_{21} & M_{22}
 \end{bmatrix}_{\overline{x},\overline{y}} = \begin{bmatrix} 
2 \overline{u}  \overline{v} - 1  &  \overline{u}^2   \\ 
-2  \overline{u}  \overline{v}  &  - \overline{u}^2
 \end{bmatrix}
 \end{align} 
```

##



To investigate the stability, we check solutions in the form of small perturbations as follows:
```math
\begin{eqnarray}
(u,v) = (U,V) e^{\lambda t}  
\end{eqnarray}
```
##
Here, $\lambda$ is the growth rate of the perturbations. The next step is to solve the eigenvalue problem, resulting of the introduction of Eq. \ref{solucion1} in the linearized system:
```math
\begin{eqnarray}
Det \left(\begin{array}{cc}M_{11}-\lambda & M_{12} \\M_{21}& M_{22}-\lambda \end{array}\right) =0 
\end{eqnarray}
```
##
and the corresponding equation:
```math
\begin{eqnarray}
\lambda^{2} -\lambda Tr(M) + Det(M)=0
\end{eqnarray}
```
##
where:
```math
\begin{eqnarray}
Tr(M)= M_{11}+M_{22} \\
Det(M)= M_{11}M_{22}-M_{12}M_{21} 
\end{eqnarray}
```
##
for this particular case

```math
\begin{eqnarray}
Tr(M)= 2 \overline{u}  \overline{v} - 1 - \overline{u}^2  \\
Det(M)= - \overline{u}^2 (2 \overline{u}  \overline{v} - 1) + 2   \overline{v} \overline{u}^3  = \overline{u}^2 
\end{eqnarray}
```
##
and the corresponding equation:
```math
\begin{eqnarray}
\lambda^{2} -\lambda (2 \overline{u}  \overline{v} - 1 - \overline{u}^2) +  \overline{u}^2=0
\end{eqnarray}


```
"




# ╔═╡ 05fe1502-dbac-4e9a-9944-7e5dcb090ae4
md"
##
We can sustitute now the value of the steady state

```math
\begin{eqnarray}
\lambda^{2} -\lambda (2 \frac{\mu}{\mu} - 1 - \mu^2) + \mu^2=0
\end{eqnarray}
```
##
so 


```math
\begin{eqnarray}
\lambda^{2} + \lambda (\mu^2 -1) + \mu^2=0
\end{eqnarray}
```

"

# ╔═╡ 8673be1a-122a-44e7-a9c7-45d12afb466b
md" ##
if $Re(\lambda) < 0$, the perturbation decays in time and the steady state ($\overline{x},\overline{y}$) is stable. 

On the contrary, when $Re(\lambda) > 0$, the perturbation grows exponentially and the steady sate is unstable. 

So, the steady state is stable if the following conditions are fulfilled:
```math
\begin{eqnarray}
Tr(M) < 0 \\
Det(M) >0
\end{eqnarray}
```
##


In this case, since $\mu$ is a reaction rate, it is allways positive, 

```math
\begin{eqnarray}
Det(M) = \mu^2 > 0
\end{eqnarray}
```

##
On the other hand, 

```math
\begin{eqnarray}
Tr(M) = 1 - \mu^2 < 0 \\ 
Tr(M) = \mu^2  > 1 \\ 
\end{eqnarray}
```

##
The Tr(M) is negative if $\mu > 1$, so the system is stable for values of $\mu>1$ and inestable for values of $\mu<1$. 

At this particular point $\mu<=1$, the polynomial becomes 

```math
\begin{eqnarray}
\lambda^{2} + 1 = 0\\
\lambda^{2} = -1\\
\lambda_{1,2} = \pm \sqrt{-1} = \pm i\\
\end{eqnarray}
```

The eigenvalue is purely imaginary, so perturbations do not decay or grow, but the system oscillates. Son at $\mu=1$ we have a stable center. It is a Hopf bifurcation. 


"

# ╔═╡ c36cfb05-9d05-42f0-87eb-a898b6c66317
md" ## 
if we solve the quadratic

```math
\begin{eqnarray}
\lambda = \frac{- (\mu^2 -1) + \sqrt{(\mu^2 -1)^2 - 4\mu^2 }}{2}\\
\lambda = \frac{- (\mu^2 -1) - \sqrt{(\mu^2 -1)^2 - 4\mu^2 }}{2}
\end{eqnarray}
```
we obtain 

```math
\begin{eqnarray}
\lambda = \frac{1 - \mu^2  + \sqrt{(\mu^2 -1)^2 - 4\mu^2 }}{2}\\
\lambda = \frac{1 - \mu^2 - \sqrt{(\mu^2 -1)^2 - 4\mu^2 }}{2}
\end{eqnarray}
```

we study when the number inside the square root changes it sign:
```math
\begin{eqnarray}
(\mu^2 -1)^2 - 4\mu^2 =0\\
(\mu^2 -1)^2 = 4\mu^2\\
\mu^2 -1 = 2\mu\\
\mu^2 +2 \mu -1 = 0
\end{eqnarray}
```

we solve 

```math
\begin{eqnarray}
\mu = \frac{-2 \pm \sqrt{2^2+4}}{2}\\
\mu = \frac{-2 \pm \sqrt{8}}{2}\\
\mu = \frac{-2 \pm 2\sqrt{2}}{2}\\
\mu = -1 \pm \sqrt{2}\\
\end{eqnarray}
```
so, between $ -1 - \sqrt{2} $ and $-1 + \sqrt{2}$ the value of $\mu < 0 $ so the square root is imaginary. This means we can get oscillations. Therefore, we have the following cases,

```math
\begin{eqnarray}
\mu  \in (0, -1 + \sqrt{2}) &\rightarrow& imaginary\\
\mu  \in (-1 + \sqrt{2},\infty) &\rightarrow& real\\
\end{eqnarray}
```

so, based on the value of the Trace, stable oscillations occur $\mu  \in (-1 + \sqrt{2},1)$, below this value $\mu  \in (0, -1 + \sqrt{2})$ we have a unstable node, if $\mu  \in (-1 + \sqrt{2},\infty)$ we have an stable spiral.

"

# ╔═╡ de338095-9ae6-4baa-bb2b-6a4b17f53e06
cubic! = @ode_def CubicAutocatalator begin
  dx = x^2*y - x
  dy =  µ - x^2*y 
end µ

# ╔═╡ 268c3afc-d66f-4419-8cc7-e0bf8ff000a0
begin
	µ_slide = @bind µ html"<input type=range min=0 max=2 step=.01>"
	md"""
	**Set the values of the kinetic constants**
	
	value of µ: $(µ_slide)
	
	"""
end

# ╔═╡ c8b683ca-9620-428d-98c4-2e5b18ec1720
begin
	aaa=1
	bbb= µ^2-1
	ccc=µ^2
	quadratic(aaa,bbb,ccc)
end

# ╔═╡ 2e76d7b1-b293-410c-9ab3-48bed3c2a388
prob3 = ODEProblem(cubic!,[x₀,y₀],(0.0,500.0),µ)

# ╔═╡ 7cb2a5e9-a1a3-438e-a44d-252b7009638d
sol3 = solve(prob3);

# ╔═╡ 1c798e63-9a9d-4e27-91a3-a2e0e1e89b75
begin
	plot(sol3,ylims = (0, 10))
	title!("Solution for for Cubic autocatalor µ = $µ ")
end

# ╔═╡ 88299113-dc0c-4956-9db8-783c763d957b
begin
	plot([u[1] for (u,t) in tuples(sol3)],[u[2] for (u,t) in tuples(sol3)],ylims = (0, 4),xlims = (0, 4))
	title!("Phase plane for Cubic autocatalor for µ = $µ ")
end

# ╔═╡ ec119596-0b9e-467f-b05b-abf77d4d0a62
begin
	vline([k3/k2],ylims = (0, 4),xlims = (0, 4));
	hline!([k1/k2],ylims = (0, 4),xlims = (0,4));
	title!("Null-Clines of the Cubic autocatalator model ")
	xlabel!("x [a.u.]")
    ylabel!("y [a.u.]")
	plot!([u[1] for (u,t) in tuples(sol)],[u[2] for (u,t) in tuples(sol)],)


	plot!([u[1] for (u,t) in tuples(sol2)],[u[2] for (u,t) in tuples(sol2)])
	
end

# ╔═╡ d533c313-2dcd-404d-8740-65294a096954
md"##

"

# ╔═╡ 67091261-3b6c-4a94-9a70-86519a1eed76
md"

### Stability of Spatial Systems

The next step is to to consider the spatial dimensions of the system in Eq. \ref{base1} and Eq. \ref{base2}.

```math
\begin{eqnarray}
\frac{\partial u}{\partial t} = f(\mu, u, v)  + D_u \frac{\partial^{2} u}{\partial \vec{r}^{2}}\\
\frac{\partial v}{\partial t} = g(\mu, u, v) + D_v \frac{\partial^{2} v}{\partial \vec{r}^{2}}
\end{eqnarray}
```
##
Here, $D_u$ and $D_v$ are the diffusion coefficients of activator and inhibitor and $ \vec{r}$ is the spatial coordinate. We will scale the diffusion coefficients in a way that we can reduce to a variable which only takes account of the ratio between the diffusion coefficients: $d=D_{v}/D_{u}$. 
##
Now we have to check solutions with the spatial part: 
```math
\begin{eqnarray}
(u,v) = (U,V) e^{\lambda t + i \vec{k}\vec{r}}  
\end{eqnarray}
```
##
The Jacobian matrix $M$ of the system is:
```math
\begin{eqnarray}
M=\left(\begin{array}{cc}M_{11} - k^{2} & M_{12} \\M_{21}& M_{22} - d k^{2} \end{array}\right) 
\end{eqnarray}
```
##
If we solve the eigenvalue problem ($Det (M-\lambda I)=0$), as in the previous case without spatial dimensions (unstable steady state) some other conditions are required to get positive eigenvalues. The equation is:
```math
\begin{eqnarray}
\lambda^2+\lambda(k^2(1+d)-Tr(M))+Det(M)=0  
\end{eqnarray}
```
##
The solution is in the form:
```math
\begin{eqnarray}
\lambda=\frac{1}{2}(-k^{2}(1+d)+Tr (M) \pm 
 \sqrt{(k^{2}(1+d)-Tr (M))^{2}-4 B} 
 \end{eqnarray}
```
##
 where

```math 
 \begin{eqnarray}
 d &=& \frac{D_{v}}{D_{u}} \\
 Tr (M) &=& M_{11} + M_{22} \\
 Det(M) &=& M_{11}  M_{22} - M_{21}  M_{12} \\
 B &=& d k^{4} - d k^{2} M_{11} - k^{2} M_{22} + Det (M) 
\end{eqnarray}
```
##
So, the system will be unstable if one of the following conditions is fulfilled: 
```math
\begin{eqnarray}
k^2(1+d)-Tr(M) < 0 \\
 Det(M) < 0  
\end{eqnarray}
```
##
In addition, if the eigenvalues are positive and real, which means that $k^2(1+d)-Tr(M))^{2}> 4 Det(M)$, the system will grow exponentially (Turing bifurcation) The system, now with spatial dimensions, develops steady periodic patterns. There is a window of unstable wavelengths which the system may exhibit ($k$ with Re[$\lambda_{1,2}] > 0$). 
##
But there is one with maximum growth rate, which can be easily calculated by solving Eq. \ref{eigen1}:
```math
\begin{eqnarray}
\frac{\partial \lambda}{\partial k}=0  
\end{eqnarray}
```
##

Fig. Re dispersion is a plot of the real part of one of the the eigenvalues $\lambda_{1}$ which has a region of positive growth for some wavenumbers in the  Lengyel-Epstein model (see Sec sec: LE model ). This means that a perturbation with a wavenumber with positive eigenvalue will grow exponentially in time. The other eigenvalue is negative, so it does not influence the behavior of the system. In addition Fig. Im_dispersion  shows the imaginary part of both eigenvalues. Positive imaginary values of the growth rate are outside of the regime of positive real values in Fig. Re_dispersion , so the periodic pattern (with wavenumber $k$) is steady in time."

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DifferentialEquations = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
ParameterizedFunctions = "65888b18-ceab-5e60-b2b9-181511a3b968"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
DifferentialEquations = "~7.7.0"
ParameterizedFunctions = "~5.15.0"
Plots = "~1.38.12"
PlutoUI = "~0.7.51"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.5"
manifest_format = "2.0"
project_hash = "11b8f1107256f4f23d70078cc8c5956c8274aa40"

[[deps.AbstractAlgebra]]
deps = ["GroupsCore", "InteractiveUtils", "LinearAlgebra", "MacroTools", "Random", "RandomExtensions", "SparseArrays", "Test"]
git-tree-sha1 = "3ee5c58774f4487a5bf2bb05e39d91ff5022b4cc"
uuid = "c3fe647b-3220-5bb0-a1ea-a7954cac585d"
version = "0.29.4"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.AbstractTrees]]
git-tree-sha1 = "faa260e4cb5aba097a73fab382dd4b5819d8ec8c"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.4"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "76289dc51920fdc6e0013c872ba9551d54961c24"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.6.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "62e51b39331de8911e4a7ff6f5aaf38a5f4cc0ae"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.2.0"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra", "Requires", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "c4d9efe93662757bca4cc24df50df5f75e659a2d"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.4.4"

[[deps.ArrayInterfaceCore]]
deps = ["LinearAlgebra", "SnoopPrecompile", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "e5f08b5689b1aad068e01751889f2f615c7db36d"
uuid = "30b0a656-2188-435a-8636-2ec0e6a096e2"
version = "0.1.29"

[[deps.ArrayLayouts]]
deps = ["FillArrays", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "4aff5fa660eb95c2e0deb6bcdabe4d9a96bc4667"
uuid = "4c555306-a7a7-4459-81d9-ec55ddd5c99a"
version = "0.8.18"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.BandedMatrices]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "SnoopPrecompile", "SparseArrays"]
git-tree-sha1 = "6ef8fc1d77b60f41041d59ce61ef9eb41ed97a83"
uuid = "aae01518-5342-5314-be14-df237901396f"
version = "0.17.18"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bijections]]
git-tree-sha1 = "fe4f8c5ee7f76f2198d5c2a06d3961c249cce7bd"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.1.4"

[[deps.BitFlags]]
git-tree-sha1 = "43b1a4a8f797c1cddadf60499a8a077d4af2cd2d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.7"

[[deps.BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "0c5f81f47bbbcf4aea7b2959135713459170798b"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.5"

[[deps.BoundaryValueDiffEq]]
deps = ["BandedMatrices", "DiffEqBase", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "NLsolve", "Reexport", "SciMLBase", "SparseArrays"]
git-tree-sha1 = "ed8e837bfb3d1e3157022c9636ec1c722b637318"
uuid = "764a87c0-6b3e-53db-9096-fe964310641d"
version = "2.11.0"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CEnum]]
git-tree-sha1 = "eb4cb44a499229b3b8426dcfb5dd85333951ff90"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.2"

[[deps.CPUSummary]]
deps = ["CpuId", "IfElse", "Static"]
git-tree-sha1 = "2c144ddb46b552f72d7eafe7cc2f50746e41ea21"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.2.2"

[[deps.CSTParser]]
deps = ["Tokenize"]
git-tree-sha1 = "3ddd48d200eb8ddf9cb3e0189fc059fd49b97c1f"
uuid = "00ebfdb7-1f24-5e51-bd34-a7502290713f"
version = "3.3.6"

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
git-tree-sha1 = "e30f2f4e20f7f186dc36529910beaedc60cfa644"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.16.0"

[[deps.ChangesOfVariables]]
deps = ["LinearAlgebra", "Test"]
git-tree-sha1 = "f84967c4497e0e1955f9a582c232b02847c5f589"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.7"

[[deps.CloseOpenIntervals]]
deps = ["Static", "StaticArrayInterface"]
git-tree-sha1 = "70232f82ffaab9dc52585e0dd043b5e0c6b714f1"
uuid = "fb6a15b2-703c-40df-9091-08a04967cfa9"
version = "0.1.12"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "9c209fb7536406834aa938fb149964b985de6c83"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.1"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "be6ab11021cd29f0344d5c4357b163af05a48cba"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.21.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "600cc5508d66b78aae350f7accdb58763ac18589"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.10"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonMark]]
deps = ["Crayons", "JSON", "PrecompileTools", "URIs"]
git-tree-sha1 = "532c4185d3c9037c0237546d817858b23cf9e071"
uuid = "a80b9123-70ca-4bc0-993e-6e3bcb318db6"
version = "0.8.12"

[[deps.CommonSolve]]
git-tree-sha1 = "9441451ee712d1aec22edad62db1a9af3dc8d852"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.3"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "7a60c856b9fa189eb34f5f8a6f6b5529b7942957"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.6.1"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.1+0"

[[deps.CompositeTypes]]
git-tree-sha1 = "02d2316b7ffceff992f3096ae48c7829a8aa0638"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.3"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "96d823b94ba8d187a6d8f0826e731195a74b90e9"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.2.0"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "738fec4d684a9a6ee9598a8bfee305b26831f28c"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.2"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.CpuId]]
deps = ["Markdown"]
git-tree-sha1 = "fcbb72b032692610bfbdb15018ac16a36cf2e406"
uuid = "adafc99b-e345-5852-983c-f28acb93d879"
version = "0.3.1"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "8da84edb865b0b5b0100c0666a9bc9a0b71c553c"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.15.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelayDiffEq]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "LinearAlgebra", "Logging", "OrdinaryDiffEq", "Printf", "RecursiveArrayTools", "Reexport", "SciMLBase", "SimpleNonlinearSolve", "SimpleUnPack"]
git-tree-sha1 = "89f3fbfe78f9d116d1ed0721d65b0b2cf9b36169"
uuid = "bcd4f6db-9728-5f36-b5f7-82caef46ccdb"
version = "5.42.0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[deps.DiffEqBase]]
deps = ["ArrayInterface", "ChainRulesCore", "DataStructures", "Distributions", "DocStringExtensions", "EnumX", "FastBroadcast", "ForwardDiff", "FunctionWrappers", "FunctionWrappersWrappers", "LinearAlgebra", "Logging", "Markdown", "MuladdMacro", "Parameters", "PreallocationTools", "Printf", "RecursiveArrayTools", "Reexport", "Requires", "SciMLBase", "Setfield", "SparseArrays", "Static", "StaticArraysCore", "Statistics", "Tricks", "TruncatedStacktraces", "ZygoteRules"]
git-tree-sha1 = "ed1108bd9a68977d5e0cbd8b2882293337c15f1c"
uuid = "2b5f629d-d688-5b77-993f-72d75c75574e"
version = "6.124.0"

[[deps.DiffEqCallbacks]]
deps = ["DataStructures", "DiffEqBase", "ForwardDiff", "LinearAlgebra", "Markdown", "NLsolve", "Parameters", "RecipesBase", "RecursiveArrayTools", "SciMLBase", "StaticArraysCore"]
git-tree-sha1 = "63b6be7b396ad395825f3cc48c56b53bfaf7e69d"
uuid = "459566f4-90b8-5000-8ac3-15dfb0a30def"
version = "2.26.1"

[[deps.DiffEqNoiseProcess]]
deps = ["DiffEqBase", "Distributions", "GPUArraysCore", "LinearAlgebra", "Markdown", "Optim", "PoissonRandom", "QuadGK", "Random", "Random123", "RandomNumbers", "RecipesBase", "RecursiveArrayTools", "Requires", "ResettableStacks", "SciMLBase", "StaticArrays", "Statistics"]
git-tree-sha1 = "2c4ed3eedb87579bfe9f20ecc2440de06b9f3b89"
uuid = "77a26b50-5914-5dd7-bc55-306e6241c503"
version = "5.16.0"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "a4ad7ef19d2cdc2eff57abbbe68032b1cd0bd8f8"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.13.0"

[[deps.DifferentialEquations]]
deps = ["BoundaryValueDiffEq", "DelayDiffEq", "DiffEqBase", "DiffEqCallbacks", "DiffEqNoiseProcess", "JumpProcesses", "LinearAlgebra", "LinearSolve", "NonlinearSolve", "OrdinaryDiffEq", "Random", "RecursiveArrayTools", "Reexport", "SciMLBase", "SteadyStateDiffEq", "StochasticDiffEq", "Sundials"]
git-tree-sha1 = "ac145e3d718157c679fc4febf2fcef73ec77b067"
uuid = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
version = "7.7.0"

[[deps.Distances]]
deps = ["LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "49eba9ad9f7ead780bfb7ee319f962c811c6d3b2"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.8"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "4f59fe4eb1308011bd33b390369cbad74e46eea4"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.92"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.DomainSets]]
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "Random", "StaticArrays", "Statistics"]
git-tree-sha1 = "698124109da77b6914f64edd696be8dccf90229e"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.6.6"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.DynamicPolynomials]]
deps = ["DataStructures", "Future", "LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Pkg", "Reexport", "Test"]
git-tree-sha1 = "8b84876e31fa39479050e2d3395c4b3b210db8b0"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.4.6"

[[deps.EnumX]]
git-tree-sha1 = "bdb1942cd4c45e3c678fd11569d5cccd80976237"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.4"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bad72f730e9e91c08d9427d5e8db95478a3c323d"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.8+0"

[[deps.ExponentialUtilities]]
deps = ["Adapt", "ArrayInterface", "GPUArraysCore", "GenericSchur", "LinearAlgebra", "Printf", "SnoopPrecompile", "SparseArrays", "libblastrampoline_jll"]
git-tree-sha1 = "fb7dbef7d2631e2d02c49e2750f7447648b0ec9b"
uuid = "d4d017d3-3776-5f7e-afef-a10c40355c18"
version = "1.24.0"

[[deps.ExprTools]]
git-tree-sha1 = "c1d06d129da9f55715c6c212866f5b1bddc5fa00"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.9"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "74faea50c1d007c85837327f6775bea60b5492dd"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+2"

[[deps.FastBroadcast]]
deps = ["ArrayInterface", "LinearAlgebra", "Polyester", "Static", "StaticArrayInterface", "StrideArraysCore"]
git-tree-sha1 = "d1248fceea0b26493fd33e8e9e8c553270da03bd"
uuid = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
version = "0.2.5"

[[deps.FastClosures]]
git-tree-sha1 = "acebe244d53ee1b461970f8910c235b259e772ef"
uuid = "9aa1b823-49e4-5ca5-8b0f-3971ec8bab6a"
version = "0.3.2"

[[deps.FastLapackInterface]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "c1293a93193f0ae94be7cf338d33e162c39d8788"
uuid = "29a986be-02c6-4525-aec4-84b980013641"
version = "1.2.9"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "7072f1e3e5a8be51d525d64f63d3ec1287ff2790"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.13.11"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "Setfield", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "6604e18a0220650dbbea7854938768f15955dd8e"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.20.0"

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
git-tree-sha1 = "00e252f4d706b3d55a8863432e742bf5717b498d"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.35"

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
git-tree-sha1 = "d62485945ce5ae9c0c48f124a84998d755bae00e"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.3"

[[deps.FunctionWrappersWrappers]]
deps = ["FunctionWrappers"]
git-tree-sha1 = "b104d487b34566608f8b4e1c39fb0b10aa279ff8"
uuid = "77dc65aa-8811-40c2-897b-53d922fa7daf"
version = "0.1.3"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "d972031d28c8c8d9d7b41a536ad7bb0c2579caca"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.8+0"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "1cd7f0af1aa58abc02ea1d872953a97359cb87fa"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.4"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Preferences", "Printf", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "UUIDs", "p7zip_jll"]
git-tree-sha1 = "efaac003187ccc71ace6c755b197284cd4811bfe"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.72.4"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4486ff47de4c18cb511a0da420efebb314556316"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.72.4+0"

[[deps.GenericSchur]]
deps = ["LinearAlgebra", "Printf"]
git-tree-sha1 = "fb69b2a645fa69ba5f474af09221b9308b160ce6"
uuid = "c145ed77-6b09-5dd9-b285-bf645a82121e"
version = "0.5.3"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "d3b3624125c1474292d0d8ed0f65554ac37ddb23"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.74.0+2"

[[deps.Glob]]
git-tree-sha1 = "97285bbd5230dd766e9ef6749b80fc617126d496"
uuid = "c27321d9-0574-5035-807b-f59d2c89b15c"
version = "1.3.1"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "1cf1d7dcb4bc32d7b4a5add4232db3750c27ecb4"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.8.0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.Groebner]]
deps = ["AbstractAlgebra", "Combinatorics", "Logging", "MultivariatePolynomials", "Primes", "Random", "SnoopPrecompile"]
git-tree-sha1 = "b6c3e9e1eb8dcc6fd9bc68fe08dcc7ab22710de6"
uuid = "0b43b601-686d-58a3-8a1c-6623616c7cd4"
version = "0.3.4"

[[deps.GroupsCore]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "9e1a5e9f3b81ad6a5c613d181664a0efc6fe6dd7"
uuid = "d5909c97-4eac-4ecc-a3dc-fdd0858a4120"
version = "0.4.0"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "41f7dfb2b20e7e8bf64f6b6fae98f4d2df027b06"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.9.4"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.HostCPUFeatures]]
deps = ["BitTwiddlingConvenienceFunctions", "IfElse", "Libdl", "Static"]
git-tree-sha1 = "734fd90dd2f920a2f1921d5388dcebe805b262dc"
uuid = "3e5b6fbb-0976-4d2c-9146-d79de83f2fb0"
version = "0.1.14"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "84204eae2dd237500835990bcade263e27674a93"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.16"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

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
git-tree-sha1 = "5cd07aab533df5170988219191dfad0519391428"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.3"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "f366daebdfb079fd1fe4e3d560f99a0c892e15bc"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IntervalSets]]
deps = ["Dates", "Random", "Statistics"]
git-tree-sha1 = "16c0cc91853084cb5f58a78bd209513900206ce6"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.4"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "6667aadd1cdee2c6cd068128b3d226ebc4fb0c67"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.9"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IterativeSolvers]]
deps = ["LinearAlgebra", "Printf", "Random", "RecipesBase", "SparseArrays"]
git-tree-sha1 = "1169632f425f79429f245113b775a0e3d121457c"
uuid = "42fd0dbc-a981-5370-80f2-aaf504508153"
version = "0.9.2"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "f377670cda23b6b7c1c0b3893e37451c5c1a2185"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.5"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6f2675ef130a300a112286de91973805fcc5ffbc"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.91+0"

[[deps.JuliaFormatter]]
deps = ["CSTParser", "CommonMark", "DataStructures", "Glob", "Pkg", "PrecompileTools", "Tokenize"]
git-tree-sha1 = "01ab91fcce19c965b3c68a82eb26260a3a7af271"
uuid = "98e50ef6-434e-11e9-1051-2b60c6c9e899"
version = "1.0.29"

[[deps.JumpProcesses]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "FunctionWrappers", "Graphs", "LinearAlgebra", "Markdown", "PoissonRandom", "Random", "RandomNumbers", "RecursiveArrayTools", "Reexport", "SciMLBase", "StaticArrays", "TreeViews", "UnPack"]
git-tree-sha1 = "50bd271af7f6cc23be7d24c8c4804809bb5d05ae"
uuid = "ccbc3e58-028d-4f4c-8cd5-9ae44345cda5"
version = "9.6.3"

[[deps.KLU]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse_jll"]
git-tree-sha1 = "764164ed65c30738750965d55652db9c94c59bfe"
uuid = "ef3ab10e-7fda-4108-b977-705223b18434"
version = "0.4.0"

[[deps.Krylov]]
deps = ["LinearAlgebra", "Printf", "SparseArrays"]
git-tree-sha1 = "0356a64062656b0cbb43c504ad5de338251f4bda"
uuid = "ba0b0d4f-ebba-5204-a429-3ac8c609bfb7"
version = "0.9.1"

[[deps.KrylovKit]]
deps = ["ChainRulesCore", "GPUArraysCore", "LinearAlgebra", "Printf"]
git-tree-sha1 = "1a5e1d9941c783b0119897d29f2eb665d876ecf3"
uuid = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"
version = "0.6.0"

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
deps = ["ArrayInterface", "ChainRulesCore", "ForwardDiff", "LinearAlgebra", "MacroTools", "PreallocationTools", "RecursiveArrayTools", "StaticArrays"]
git-tree-sha1 = "cd04158424635efd05ff38d5f55843397b7416a9"
uuid = "2ee39098-c373-598a-b85f-a56591580800"
version = "1.14.0"

[[deps.LambertW]]
git-tree-sha1 = "c5ffc834de5d61d00d2b0e18c96267cffc21f648"
uuid = "984bce1d-4616-540c-a9ee-88d1112d94c9"
version = "0.4.6"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "8c57307b5d9bb3be1ff2da469063628631d4d51e"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.21"

[[deps.LayoutPointers]]
deps = ["ArrayInterface", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "88b8f66b604da079a627b6fb2860d3704a6729a1"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.14"

[[deps.Lazy]]
deps = ["MacroTools"]
git-tree-sha1 = "1370f8202dac30758f3c345f9909b97f53d87d3f"
uuid = "50d2b5c4-7a5e-59d5-8109-a42b560f39c0"
version = "0.15.1"

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
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

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
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "7bbea35cec17305fc70a0e5b4641477dc0789d9d"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.2.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LinearSolve]]
deps = ["ArrayInterface", "DocStringExtensions", "EnumX", "FastLapackInterface", "GPUArraysCore", "IterativeSolvers", "KLU", "Krylov", "KrylovKit", "LinearAlgebra", "Preferences", "RecursiveFactorization", "Reexport", "SciMLBase", "SciMLOperators", "Setfield", "SnoopPrecompile", "SparseArrays", "Sparspak", "SuiteSparse", "UnPack"]
git-tree-sha1 = "4a4f8cc7a59fadbb02d1852d1e0cef5dca3a9460"
uuid = "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"
version = "1.42.0"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "0a1b7c2863e44523180fdb3146534e265a91870b"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.23"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "cedb76b37bc5a6c702ade66be44f831fa23c681e"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.0"

[[deps.LoopVectorization]]
deps = ["ArrayInterface", "ArrayInterfaceCore", "CPUSummary", "ChainRulesCore", "CloseOpenIntervals", "DocStringExtensions", "ForwardDiff", "HostCPUFeatures", "IfElse", "LayoutPointers", "LinearAlgebra", "OffsetArrays", "PolyesterWeave", "PrecompileTools", "SIMDTypes", "SLEEFPirates", "SpecialFunctions", "Static", "StaticArrayInterface", "ThreadingUtilities", "UnPack", "VectorizationBase"]
git-tree-sha1 = "3bb62b5003bc7d2d49f26663484267dc49fa1bf5"
uuid = "bdcacae8-1622-11e9-2a5c-532679323890"
version = "0.12.159"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.ManualMemory]]
git-tree-sha1 = "bcaef4fc7a0cfe2cba636d84cda54b5e4e4ca3cd"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.8"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "Random", "Sockets"]
git-tree-sha1 = "03a9b9718f5682ecb107ac9f7308991db4ce395b"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.7"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.ModelingToolkit]]
deps = ["AbstractTrees", "ArrayInterface", "Combinatorics", "Compat", "ConstructionBase", "DataStructures", "DiffEqBase", "DiffEqCallbacks", "DiffRules", "Distributed", "Distributions", "DocStringExtensions", "DomainSets", "ForwardDiff", "FunctionWrappersWrappers", "Graphs", "IfElse", "InteractiveUtils", "JuliaFormatter", "JumpProcesses", "LabelledArrays", "Latexify", "Libdl", "LinearAlgebra", "MacroTools", "NaNMath", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLBase", "Serialization", "Setfield", "SimpleNonlinearSolve", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicIndexingInterface", "SymbolicUtils", "Symbolics", "UnPack", "Unitful"]
git-tree-sha1 = "e52b0337d18fd6ca8b1e1e3d801636dba203b966"
uuid = "961ee093-0014-501f-94e3-6117800e7a78"
version = "8.55.1"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.MuladdMacro]]
git-tree-sha1 = "cac9cc5499c25554cba55cd3c30543cff5ca4fab"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.4"

[[deps.MultivariatePolynomials]]
deps = ["ChainRulesCore", "DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "eaa98afe2033ffc0629f9d0d83961d66a021dfcc"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.4.7"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "964cb1a7069723727025ae295408747a0b36a854"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.3.0"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "a0b464d183da839699f4c79e7606d9d186ec172c"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.3"

[[deps.NLsolve]]
deps = ["Distances", "LineSearches", "LinearAlgebra", "NLSolversBase", "Printf", "Reexport"]
git-tree-sha1 = "019f12e9a1a7880459d0173c182e6a99365d7ac1"
uuid = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
version = "4.5.1"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.NonlinearSolve]]
deps = ["ArrayInterface", "DiffEqBase", "EnumX", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "LinearSolve", "RecursiveArrayTools", "Reexport", "SciMLBase", "SimpleNonlinearSolve", "SnoopPrecompile", "SparseArrays", "SparseDiffTools", "StaticArraysCore", "UnPack"]
git-tree-sha1 = "a6000c813371cd3cd9cbbdf8a356fc3a97138d92"
uuid = "8913a72c-1f9b-4ce2-8d82-65094dcecaec"
version = "1.6.0"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "82d7c9e310fe55aa54996e6f7f94674e2a38fcb4"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.12.9"

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

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "51901a49222b09e3743c65b8847687ae5fc78eb2"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.1"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9ff31d101d987eb9d66bd8b176ac7c277beccd09"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.20+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Optim]]
deps = ["Compat", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "a89b11f0f354f06099e4001c151dffad7ebab015"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.7.5"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "d321bf2de576bf25ec4d3e4360faca399afca282"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.0"

[[deps.OrdinaryDiffEq]]
deps = ["Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "ExponentialUtilities", "FastBroadcast", "FastClosures", "FiniteDiff", "ForwardDiff", "FunctionWrappersWrappers", "IfElse", "LineSearches", "LinearAlgebra", "LinearSolve", "Logging", "LoopVectorization", "MacroTools", "MuladdMacro", "NLsolve", "NonlinearSolve", "Polyester", "PreallocationTools", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLNLSolve", "SimpleNonlinearSolve", "SimpleUnPack", "SparseArrays", "SparseDiffTools", "StaticArrayInterface", "StaticArrays", "TruncatedStacktraces"]
git-tree-sha1 = "47fc5cf4174a7d45fa541669abc5405d9ef6b8df"
uuid = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"
version = "6.51.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.40.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "67eae2738d63117a196f497d7db789821bce61d1"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.17"

[[deps.ParameterizedFunctions]]
deps = ["DataStructures", "DiffEqBase", "DocStringExtensions", "Latexify", "LinearAlgebra", "ModelingToolkit", "Reexport", "SciMLBase"]
git-tree-sha1 = "78ab7ecc18b307e00abba28bb29d7ed6bf11b9f7"
uuid = "65888b18-ceab-5e60-b2b9-181511a3b968"
version = "5.15.0"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7302075e5e06da7d000d9bfa055013e3e85578ca"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.9"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

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
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "1f03a2d339f42dca4a4da149c7e15e9b896ad899"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.1.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "f92e1315dadf8c46561fb9396e525f7200cdc227"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.5"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Preferences", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "d03ef538114b38f89d66776f2d8fdc0280f90621"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.38.12"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "b478a748be27bd2f2c73a7690da219d0844db305"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.51"

[[deps.PoissonRandom]]
deps = ["Random"]
git-tree-sha1 = "a0f1159c33f846aa77c3f30ebbc69795e5327152"
uuid = "e409e4f3-bfea-5376-8464-e040bb5c01ab"
version = "0.4.4"

[[deps.Polyester]]
deps = ["ArrayInterface", "BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "ManualMemory", "PolyesterWeave", "Requires", "Static", "StaticArrayInterface", "StrideArraysCore", "ThreadingUtilities"]
git-tree-sha1 = "0fe4e7c4d8ff4c70bfa507f0dd96fa161b115777"
uuid = "f517fe37-dbe3-4b94-8317-1923a5111588"
version = "0.7.3"

[[deps.PolyesterWeave]]
deps = ["BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "Static", "ThreadingUtilities"]
git-tree-sha1 = "240d7170f5ffdb285f9427b92333c3463bf65bf6"
uuid = "1d0040c9-8b98-4ee7-8388-3f51789ca0ad"
version = "0.2.1"

[[deps.PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

[[deps.PreallocationTools]]
deps = ["Adapt", "ArrayInterface", "ForwardDiff", "Requires"]
git-tree-sha1 = "f739b1b3cc7b9949af3b35089931f2b58c289163"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.4.12"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "259e206946c293698122f63e2b513a7c99a244e8"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.1.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "7eb1686b4f04b82f96ed7a4ea5890a4f0c7a09f1"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.0"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "311a2aa90a64076ea0fac2ad7492e914e6feeb81"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.3"

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
git-tree-sha1 = "6ec7ac8412e83d57e313393220879ede1740f9ee"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.8.2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Random123]]
deps = ["Random", "RandomNumbers"]
git-tree-sha1 = "552f30e847641591ba3f39fd1bed559b9deb0ef3"
uuid = "74087812-796a-5b5d-8853-05524746bad3"
version = "1.6.1"

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
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterface", "DocStringExtensions", "GPUArraysCore", "IteratorInterfaceExtensions", "LinearAlgebra", "RecipesBase", "Requires", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables"]
git-tree-sha1 = "02ef02926f30d53b94be443bfaea010c47f6b556"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "2.38.5"

[[deps.RecursiveFactorization]]
deps = ["LinearAlgebra", "LoopVectorization", "Polyester", "SnoopPrecompile", "StrideArraysCore", "TriangularSolve"]
git-tree-sha1 = "9088515ad915c99026beb5436d0a09cd8c18163e"
uuid = "f2c3362d-daeb-58d1-803e-2bc74f2840b4"
version = "0.2.18"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "90bc7a7c96410424509e4263e277e43250c05691"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.0"

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
git-tree-sha1 = "f65dcb5fa46aee0cf9ed6274ccbd597adc49aa7b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.1"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6ed52fdd3382cf21947b15e8870ac0ddbff736da"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.4.0+0"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "d7d9ebe28062161c1e314ed643097b0c6fe657d9"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.7"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[deps.SLEEFPirates]]
deps = ["IfElse", "Static", "VectorizationBase"]
git-tree-sha1 = "cda0aece8080e992f6370491b08ef3909d1c04e7"
uuid = "476501e8-09a2-5ece-8869-fb82de89a1fa"
version = "0.6.38"

[[deps.SciMLBase]]
deps = ["ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "PrecompileTools", "Preferences", "RecipesBase", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLOperators", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables", "TruncatedStacktraces"]
git-tree-sha1 = "e803672f8d58e9937f59923dd3b159c9b7e1838b"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "1.92.0"

[[deps.SciMLNLSolve]]
deps = ["DiffEqBase", "LineSearches", "NLsolve", "Reexport", "SciMLBase"]
git-tree-sha1 = "a8eb97c56cac50c21096582afb2a0110784dc36e"
uuid = "e9a6253c-8580-4d32-9898-8661bb511710"
version = "0.1.6"

[[deps.SciMLOperators]]
deps = ["ArrayInterface", "DocStringExtensions", "Lazy", "LinearAlgebra", "Setfield", "SparseArrays", "StaticArraysCore", "Tricks"]
git-tree-sha1 = "90163ebc767cba9f126ea00aeef1d75ed74fe7b0"
uuid = "c0aeaf25-5076-4817-a8d5-81caf7dfa961"
version = "0.2.8"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "30449ee12237627992a99d5e30ae63e4d78cd24a"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.SimpleNonlinearSolve]]
deps = ["ArrayInterface", "DiffEqBase", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "Reexport", "Requires", "SciMLBase", "SnoopPrecompile", "StaticArraysCore"]
git-tree-sha1 = "54c78ac3cc0343a16785adabe5bbf4063c737967"
uuid = "727e6d20-b764-4bd8-a329-72de5adea6c7"
version = "0.1.14"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.SimpleUnPack]]
git-tree-sha1 = "58e6353e72cde29b90a69527e56df1b5c3d8c437"
uuid = "ce78b400-467f-4804-87d8-8f486da07d0a"
version = "1.1.0"

[[deps.SnoopPrecompile]]
deps = ["Preferences"]
git-tree-sha1 = "e760a70afdcd461cf01a575947738d359234665c"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "a4ada03f999bd01b3a25dcaa30b2d929fe537e00"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.0"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SparseDiffTools]]
deps = ["Adapt", "ArrayInterface", "Compat", "DataStructures", "FiniteDiff", "ForwardDiff", "Graphs", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays", "VertexSafeGraphs"]
git-tree-sha1 = "e19ac47477c9a8fcca06dab5e5471417d5d9d723"
uuid = "47a9eef4-7e08-11e9-0b38-333d64bd3804"
version = "1.31.0"

[[deps.Sparspak]]
deps = ["Libdl", "LinearAlgebra", "Logging", "OffsetArrays", "Printf", "SparseArrays", "Test"]
git-tree-sha1 = "342cf4b449c299d8d1ceaf00b7a49f4fbc7940e7"
uuid = "e56a9233-b9d6-4f03-8d0f-1825330902ac"
version = "0.3.9"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "ef28127915f4229c971eb43f3fc075dd3fe91880"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.2.0"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "dbde6766fc677423598138a5951269432b0fcc90"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.8.7"

[[deps.StaticArrayInterface]]
deps = ["ArrayInterface", "Compat", "IfElse", "LinearAlgebra", "Requires", "SnoopPrecompile", "SparseArrays", "Static", "SuiteSparse"]
git-tree-sha1 = "33040351d2403b84afce74dae2e22d3f5b18edcb"
uuid = "0d7ed370-da01-4f52-bd93-41d350b8b718"
version = "1.4.0"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "8982b3607a212b070a5e46eea83eb62b4744ae12"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.25"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6b7ba252635a5eff6a0b0664a41ee140a1c9e72a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "45a7769a04a3cf80da1c1c7c60caf932e6f4c9f7"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.6.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.StatsFuns]]
deps = ["ChainRulesCore", "HypergeometricFunctions", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "f625d686d5a88bcd2b15cd81f18f98186fdc0c9a"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.0"

[[deps.SteadyStateDiffEq]]
deps = ["DiffEqBase", "DiffEqCallbacks", "LinearAlgebra", "NLsolve", "Reexport", "SciMLBase"]
git-tree-sha1 = "564451a262696334a3bab19108a99dd90d5a22c8"
uuid = "9672c7b4-1e72-59bd-8a11-6ac3964bc41f"
version = "1.15.0"

[[deps.StochasticDiffEq]]
deps = ["Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DiffEqNoiseProcess", "DocStringExtensions", "FillArrays", "FiniteDiff", "ForwardDiff", "JumpProcesses", "LevyArea", "LinearAlgebra", "Logging", "MuladdMacro", "NLsolve", "OrdinaryDiffEq", "Random", "RandomNumbers", "RecursiveArrayTools", "Reexport", "SciMLBase", "SparseArrays", "SparseDiffTools", "StaticArrays", "UnPack"]
git-tree-sha1 = "073da86200349ddf4ef8bc3e3f3acd62e1d554f7"
uuid = "789caeaf-c7a9-5a7d-9973-96adeb23e2a0"
version = "6.60.0"

[[deps.StrideArraysCore]]
deps = ["ArrayInterface", "CloseOpenIntervals", "IfElse", "LayoutPointers", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface", "ThreadingUtilities"]
git-tree-sha1 = "5ffcee1813efc849f188dce82ca1553bd5f3a476"
uuid = "7792a7ef-975c-4747-a70f-980b88e8d1da"
version = "0.4.14"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+0"

[[deps.Sundials]]
deps = ["CEnum", "DataStructures", "DiffEqBase", "Libdl", "LinearAlgebra", "Logging", "PrecompileTools", "Reexport", "SciMLBase", "SparseArrays", "Sundials_jll"]
git-tree-sha1 = "ace8080f882c5181d61c8dbb749ac9aa72a49bd0"
uuid = "c3572dad-4567-51f8-b174-8c6c989267f4"
version = "4.17.0"

[[deps.Sundials_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "OpenBLAS_jll", "Pkg", "SuiteSparse_jll"]
git-tree-sha1 = "04777432d74ec5bc91ca047c9e0e0fd7f81acdb6"
uuid = "fb77eaff-e24c-56d4-86b1-d163f2edb164"
version = "5.2.1+0"

[[deps.SymbolicIndexingInterface]]
deps = ["DocStringExtensions"]
git-tree-sha1 = "f8ab052bfcbdb9b48fad2c80c873aa0d0344dfe5"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.2.2"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "Bijections", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "IfElse", "LabelledArrays", "LinearAlgebra", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "TimerOutputs", "Unityper"]
git-tree-sha1 = "5cb1f963f82e7b81305102dd69472fcd3e0e1483"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "1.0.5"

[[deps.Symbolics]]
deps = ["ArrayInterface", "ConstructionBase", "DataStructures", "DiffRules", "Distributions", "DocStringExtensions", "DomainSets", "Groebner", "IfElse", "LaTeXStrings", "LambertW", "Latexify", "Libdl", "LinearAlgebra", "MacroTools", "Markdown", "NaNMath", "RecipesBase", "Reexport", "Requires", "RuntimeGeneratedFunctions", "SciMLBase", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicUtils", "TreeViews"]
git-tree-sha1 = "e23ec62c083ca8f15a4b7174331b3b8d1c511e47"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "5.3.1"

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
git-tree-sha1 = "1544b926975372da01227b382066ab70e574a3ec"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.1"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.ThreadingUtilities]]
deps = ["ManualMemory"]
git-tree-sha1 = "c97f60dd4f2331e1a495527f80d242501d2f9865"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.5.1"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "f548a9e9c490030e545f72074a41edfd0e5bcdd7"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.23"

[[deps.Tokenize]]
git-tree-sha1 = "90538bf898832b6ebd900fa40f223e695970e3a5"
uuid = "0796e94c-ce3b-5d07-9a54-7f471281c624"
version = "0.5.25"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "9a6ae7ed916312b41236fcef7e0af564ef934769"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.13"

[[deps.TreeViews]]
deps = ["Test"]
git-tree-sha1 = "8d0d7a3fe2f30d6a7f833a5f19f7c7a5b396eae6"
uuid = "a2a6695c-b41b-5b7d-aed9-dbfdeacea5d7"
version = "0.3.0"

[[deps.TriangularSolve]]
deps = ["CloseOpenIntervals", "IfElse", "LayoutPointers", "LinearAlgebra", "LoopVectorization", "Polyester", "Static", "VectorizationBase"]
git-tree-sha1 = "31eedbc0b6d07c08a700e26d31298ac27ef330eb"
uuid = "d5829a12-d9aa-46ab-831f-fb7c9ab06edf"
version = "0.1.19"

[[deps.Tricks]]
git-tree-sha1 = "aadb748be58b492045b4f56166b5188aa63ce549"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.7"

[[deps.TruncatedStacktraces]]
deps = ["InteractiveUtils", "MacroTools", "Preferences"]
git-tree-sha1 = "7bc1632a4eafbe9bd94cf1a784a9a4eb5e040a91"
uuid = "781d530d-4396-4725-bb49-402e4bee1e77"
version = "1.3.0"

[[deps.URIs]]
git-tree-sha1 = "074f993b0ca030848b897beff716d93aca60f06a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.2"

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
git-tree-sha1 = "ba4aa36b2d5c98d6ed1f149da916b3ba46527b2b"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.14.0"

[[deps.Unityper]]
deps = ["ConstructionBase"]
git-tree-sha1 = "d5f4ec8c22db63bd3ccb239f640e895cfde145aa"
uuid = "a7c27f48-0311-42f6-a7f8-2c11e75eb415"
version = "0.1.2"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.VectorizationBase]]
deps = ["ArrayInterface", "CPUSummary", "HostCPUFeatures", "IfElse", "LayoutPointers", "Libdl", "LinearAlgebra", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "b182207d4af54ac64cbc71797765068fdeff475d"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.21.64"

[[deps.VertexSafeGraphs]]
deps = ["Graphs"]
git-tree-sha1 = "8351f8d73d7e880bfc042a8b6922684ebeafb35c"
uuid = "19fa3120-7c27-5ec5-8db8-b0b0aa330d6f"
version = "0.2.0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "ed8d92d9774b077c53e1da50fd81a36af3744c1c"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+0"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "93c41695bc1c08c46c5899f4fe06d6ead504bb73"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.10.3+0"

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
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "49ce682769cd5de6c72dcf1b94ed7790cd08974c"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.5+0"

[[deps.ZygoteRules]]
deps = ["ChainRulesCore", "MacroTools"]
git-tree-sha1 = "977aed5d006b840e2e40c0b48984f7463109046d"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.3"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "868e669ccb12ba16eaf50cb2957ee2ff61261c56"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.29.0+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

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
git-tree-sha1 = "9ebfc140cc56e8c2156a15ceac2f0302e327ac0a"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+0"
"""

# ╔═╡ Cell order:
# ╟─2e6401b9-9c89-4b4b-8492-d0a83003579b
# ╟─8ba695de-0650-407c-864a-c394c9c7d3ed
# ╟─0102b1ee-7f41-4400-aef3-c33e18b5b716
# ╟─df36d23e-3760-48c1-8c0c-d52f15e7f520
# ╟─0b0440bf-4f75-47bd-bbd8-315066e79c75
# ╟─279eaa4f-41ca-4168-b492-f36543bb8204
# ╟─00037d0a-8be0-4c75-94bc-23c2c639e6e8
# ╟─64dc539b-b268-43b6-a96e-2c4fa7b48474
# ╟─f5faefd6-f496-40fb-8043-ef6b9d172ae0
# ╟─32b3bc77-c4f4-40ab-a5a9-e035022189ac
# ╟─c70c7a06-d1aa-40a9-902b-4bd4b23a6392
# ╟─2cfd31c4-803a-4880-b2a4-5af6594c0975
# ╟─2170f883-1b01-4d16-907e-49c72118e224
# ╟─7d89a7bd-baba-42c1-af4f-26f74138199a
# ╟─9db950f3-6989-4225-845f-dae93d52a9b9
# ╟─c620fb38-6e3d-4d54-b515-a9a5fc4360e4
# ╟─66ea1e28-b50a-43fc-acf5-9a146b51d505
# ╟─b22fd4f1-ee14-4cb4-aff9-94153f470fe1
# ╟─555aa5d3-9e29-459a-8798-51a576d252b4
# ╟─8c809d92-4edb-40ca-ac80-18b7e8c72fae
# ╟─6c87979a-690a-400b-9ab8-7ef26b829195
# ╟─04c14475-24cd-49c3-8a7b-83a9425664be
# ╟─919bf3e4-6aa0-440b-b76f-34a4812c6752
# ╟─855b8530-861e-427e-b6cc-c1b0283568ca
# ╟─c112e62c-9c38-48dc-8294-42911b02ccc5
# ╟─6c96da07-b5ea-4096-a14c-5cd8f0f47c04
# ╟─532acbb6-e504-4deb-88d5-e001a476af15
# ╟─906a0931-e86f-48f3-be11-eb692639e8cc
# ╟─c374e1bd-ec3d-4cf7-aff1-072d6ba60b27
# ╟─2b91a06c-3803-4877-ac9a-1669fbe58d30
# ╟─e8ca5682-cce0-4b2c-ba18-47f2e38192c7
# ╟─17e34ce6-a1fa-4193-9201-3e61a188b48f
# ╠═87df7c88-7e16-4c50-bd80-df5bfea26ee8
# ╟─9433f4d1-3abf-4dea-8107-1070ee5feb2e
# ╟─c3828c90-0559-4686-9f56-41b5b6d48176
# ╟─2b9e43fd-d00d-40f1-ac98-07d48c53d861
# ╠═ff0c4a06-de7b-4cf2-b920-cd84ac9927ea
# ╠═f8238a6e-1e5e-4629-ac20-2c9cb964c53e
# ╠═3349c942-15a8-4b1e-8b4a-7f5154f48b12
# ╟─df81eb5c-16ef-4bf6-9ad0-2cc899c44cb6
# ╟─a5b38e90-cdbe-442f-a078-dd8598a0a4c2
# ╟─3118e37d-f7bd-488e-ab3a-74e0194d3108
# ╟─48dbf14c-ab77-4e6a-9503-8dc53f5dfba8
# ╟─cf72446a-a512-11ec-2b47-ef706c91c6a0
# ╟─b67ed8c0-8ca7-4928-9ebc-98da513c432f
# ╟─05fe1502-dbac-4e9a-9944-7e5dcb090ae4
# ╟─8673be1a-122a-44e7-a9c7-45d12afb466b
# ╟─c36cfb05-9d05-42f0-87eb-a898b6c66317
# ╠═de338095-9ae6-4baa-bb2b-6a4b17f53e06
# ╟─c8b683ca-9620-428d-98c4-2e5b18ec1720
# ╠═2e76d7b1-b293-410c-9ab3-48bed3c2a388
# ╠═7cb2a5e9-a1a3-438e-a44d-252b7009638d
# ╟─268c3afc-d66f-4419-8cc7-e0bf8ff000a0
# ╠═1c798e63-9a9d-4e27-91a3-a2e0e1e89b75
# ╟─88299113-dc0c-4956-9db8-783c763d957b
# ╠═ec119596-0b9e-467f-b05b-abf77d4d0a62
# ╠═d533c313-2dcd-404d-8740-65294a096954
# ╟─67091261-3b6c-4a94-9a70-86519a1eed76
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
