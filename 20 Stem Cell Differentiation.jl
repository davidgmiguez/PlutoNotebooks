### A Pluto.jl notebook ###
# v0.19.20

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
using Plots , DifferentialEquations, HypertextLiteral,PlutoUI

# ╔═╡ 7f30ad4d-9964-4140-8003-051653d5f1e4
md" # Differentiation dynamics

what happens when we analyze the dynamcis of a population of cells that not only proliferates, but also differentiates into another type of cell? For instance, a population of stem cells, in a developing organ. We assume a developing organ as a population of cycling progenitors 'P' that cycle with an average cell cycle $T$. Some of these cells terminally differentiate, exit the cell cycle and acquire a given specialized phenotype 'D'. We start from a initial population of progenitors $P_0$ and differentiated $D_0$ cells. A common approac is top characterize the dynamics of the population focusing on the outcome of the cell division of the progenitors. In principle, each division of a 'P' cell can give two progenitors ('pp' division), two differentiated cells ('dd' division) and also an assymetric mode of divisin where a progenitor and a differentiated cell is generated ('pd' division). If we calculate the average amount of 'P' and 'D' generated after a single cell cycle (n=1) we can write the number of progenitor and differentiated cells as:

```math
\begin{eqnarray}
P_1&=&P_{0}(2pp+pd) \tag{35}\\
D_1&=&D_{0}+P_{0}(2dd+pd)\tag{36}
\end{eqnarray}
```

where using the condition $pp+pd+dd=1$, 

```math
\begin{eqnarray}
P_1&=&P_{0} (1+pp-dd)\tag{37}\\
D_1&=&D_{0}+P_{0}(1+dd-pp)\tag{38}
\end{eqnarray}
```

Therfore, for n=2,

```math
\begin{eqnarray}
P_2&=&P_{1} (1+pp-dd)\tag{39}\\
D_2&=&D_{1}+P_{1}(1+dd-pp)\tag{40}
\end{eqnarray}
```

applying eqs. 37 and 38, we obtain

```math
\begin{eqnarray}
P_2&=&P_{0}(1+pp-dd)(1+pp-dd)=P_{0} (1+pp-dd)^2\tag{41}\\
D_2&=&D_{0}+P_{0}(1+dd-pp)+P_{0}(1+pp-dd)(1+dd-pp)\tag{42}
\end{eqnarray}
```

and rearranging terms in eq. 42:


```math
\begin{eqnarray}
D_2&=&D_{0}+P_{0}(1+dd-pp)(1+(1+pp-dd))\tag{43}\\
\end{eqnarray}
```

Subsequently, for n=3

```math
\begin{eqnarray}
P_3&=&P_{2} (1+pp-dd)\tag{44}\\
D_3&=&D_{2}+P_{2}(1+dd-pp) \tag{45}
\end{eqnarray}
```


applying eqs. 41 and 42, we obtain

```math
\begin{eqnarray}
P_3&=&P_{0} (1+pp-dd)^2 (1+pp-dd)= P_{0} (1+pp-dd)^3 \tag{46}\\
D_3&=&D_{0}+P_{0}(1+dd-pp)(1+(1+pp-dd) + (1+pp-dd)^2) \tag{47}
\end{eqnarray}
```

therefore, for $n$ steps, we obtain, 

```math
\begin{eqnarray}
P_n&=&P_{0} (1+pp-dd)^n \tag{48}\\
D_n&=&D_{0}+P_{0}(1+dd-pp)(1+(1+pp-dd)+(1+pp-dd)^2+...+(1+pp-dd)^{n-1})\tag{49}
\end{eqnarray}
```

where the second term in eq. 49 can be written as

```math
\begin{eqnarray}
1+(1+pp-dd)+(1+pp-dd)^2+...+(pp-dd)^{n-1}=\displaystyle\sum_{i=0}^{n-1} (1+pp-dd)^i\tag{50}
\end{eqnarray}
```

which renaming $r=1+pp-dd$ is equivalent to 

```math
\begin{eqnarray}
\displaystyle\sum_{i=0}^{n-1} r^i=\frac{1-r^n}{1-r}\tag{51}
\end{eqnarray}
```

therefore, eq. 50 can be written as

```math
\begin{eqnarray}
D_n&=&D_{0}+P_{0}(1+dd-pp)\frac{1-(1+pp-dd))^n}{1-(1+pp-dd)} \tag{51}
\end{eqnarray}
```

which, after simplifying terms, can be rewritten as

```math
\begin{eqnarray}
D_n&=&D_{0}+P_{0}(1-(1+pp-dd)^n)\frac{1+dd-pp}{dd-pp}\tag{52}
\end{eqnarray}
```

and taking into account eq. 48, we obtain the final equation for the number of progenitors 'P' and differentiated 'D' cells in a stem cell population that is growing and differentiating:

```math
\begin{eqnarray}
D_n&=&D_{0}+(P_{0}-P_{n})\frac{1+dd-pp}{dd-pp}=D_{0}+(P_{n}-P_{0})\frac{1+dd-pp}{pp-dd} \tag{53}
\end{eqnarray}
```

Interestingly, both equations depend on the iteration step $n$ only via the number of progenitors at a given time in the system $P_n$. 

"


# ╔═╡ ede201b5-8f82-4cb0-be2a-be03b9140c50
md"For simplicity, the system of equations has been derived for a situation of discrete $n=\Delta t/T$, (n=1,2,3...), i.e., with the time step equal to the average cell cycle $\Delta t=T$. If we instead consider the time step as half of the cell cycle ($\Delta t=T/2$), then $n=2 \Delta t/T$, (n=1,2,3...), and eq. 37 is now:

```math
\begin{eqnarray}
P_1&=&P_{0} (1+pp-dd)^{\frac{1}{2}}\tag{54}
\end{eqnarray}
```

and following identical iteration steps we arrive at

```math
\begin{eqnarray}
P_n&=&P_{0} (1+pp-dd)^{\frac{n}{2}}\tag{55}
\end{eqnarray}
```

while eqs. for 'D' cells remain the same, since it does not depend explicitly on the iteration step. This way, for $\delta t=1$, and following the same process we obtain

```math
\begin{eqnarray}
P_n&=&P_{0} (1+pp-dd)^{\frac{n}{T}}\tag{56}
\end{eqnarray}
```

which can be generalize for the continuum limit $n=t$,

```math
\begin{eqnarray}
P_{t}&=&P_0 (1+pp-dd)^{\frac{t}{T}}\tag{57}\\
D_{t}&=&D_{0}+(P_{t}-P_{0})\frac{1+dd-pp}{pp-dd}\tag{58}
\end{eqnarray}
```

Finally, f we rewrite  $\Delta$$P=P_{t}-P_{0}$, we obtain the expression. 

```math
\begin{eqnarray}
P_{t}&=&P_0 (1+pp-dd)^{\frac{t}{T}}\tag{59}\\
D_{t}&=&D_{0}+\Delta P\frac{1+dd-pp}{pp-dd}\tag{60}
\end{eqnarray}
```
"

# ╔═╡ fffe908a-a8ff-404b-9138-bb23210a5fae
begin
	pp_slide = @bind pp html"<input type=range min=0.0 max=1.0 step=0.1>"
	dd_slide = @bind dd html"<input type=range min=0.0 max=1.0 step=0.1>"
	
	md"""
	**Set the proliferation and differentiation dynamics**
	
	value of pp: $(pp_slide)
	
	value of dd: $(dd_slide)
	
	"""
end

# ╔═╡ 470fe58c-4020-4333-861f-a6afe54a9e53
begin
		P₀=100
		D₀=50
		#T=24
		t=collect(0:10)
		plot(t,t->P₀*(1+pp-dd)^t,label="P",seriestype=:line)
		plot!(t,t->D₀+P₀*(((1+pp-dd)^t)-1)*((1-pp+dd)/(pp-dd)),label="D",seriestype=:line,ylims = (0,400))
		plot!(t,t->P₀*(1+pp-dd)^t+D₀+P₀*(((1+pp-dd)^t)-1)*((1-pp+dd)/(pp-dd)),label="T",seriestype=:line,ylims = (0,400))
end

# ╔═╡ 79536c27-2f74-40a2-ae9f-b44ed021f208
begin
	T_slide = @bind T html"<input type=range min=1.0 max=5.0 step=0.1>"
	
	md"""
	**Set the value of the cell cycle**
	
	value of T: $(T_slide)

	
	"""
end

# ╔═╡ 20980ae1-0efa-4836-96e8-8ce71f513fae
begin
		plot(t,t->P₀*(1+pp-dd)^(t./T),label="P",seriestype=:line)
		plot!(t,t->D₀+P₀*(((1+pp-dd)^(t./T))-1)*((1-pp+dd)/(pp-dd)),label="D",seriestype=:line,ylims = (0,400))
		plot!(t,t->P₀*(1+pp-dd)^(t./T)+D₀+P₀*(((1+pp-dd)^(t./T))-1)*((1-pp+dd)/(pp-dd)),label="T",seriestype=:line,ylims = (0,400))
end

# ╔═╡ aecc1472-9167-4b09-a08c-9e701def7d54
md"The true power of these equations is that they are analytical, in the sense that now we turn the equations down to see if we can predict the correct value of pp-dd and T"

# ╔═╡ 497ebc18-b94a-420d-a21a-708d699dec5c
begin
	
	P=P₀.*(1+pp-dd).^(t./T)
	D=D₀.+P₀.*(((1 .+pp-dd).^(t./T)).-1).*((1-pp+dd)/(pp-dd))
	#P4=plot(t,P,label="P",seriestype=:line,ylims = (0,4300))
    #plot!(t,D,label="D",seriestype=:line,ylims = (0,400))
	P_=P[2:end]
	D_=D[2:end]
	t_=t[2:end]
	P__=P[1:end-1]
	D__=D[1:end-1];
	t__=t[1:end-1];

	gamma=1
	pp_dd=(P_ .-P__) ./(P_ .-P__ .+D_ .-D__);
	T_=(t_ .-t__) .* log.(1 .+(gamma.*(abs.(pp_dd)))) ./ abs.(log.(P_ ./P__));

	P2=scatter(t_,pp_dd, ylims = (-1,1),label=" Estimated pp-dd")
	hline!([pp-dd],label=" true pp-dd")
	
	P3=scatter(t_,T_,ylims = (0,2),label=" Estimated T")
	hline!([T], label=" True T")

	plot(P2,P3,layout=(1,2),legend=true,size = (800, 500))
	
end

# ╔═╡ ae184b0d-f61d-4aea-853b-17cb597d1087
md"
to do:
    
    it will be cool to include a term of saturations as the logistic equation does. 
    
"

# ╔═╡ ca6a3577-ac75-4993-b2e5-4975ef4aaf9f
md"The equations for the cell cycle start to fail whene we go to values lower than pp-dd=0"

# ╔═╡ 3e8066bd-d8dc-47d9-98e5-0b8b0759e179
md" ### Including Apoptosis

what if we include apoptosis"

# ╔═╡ 2e9621c6-50c4-4d30-84f5-f26ea707a808
pp=0.6
dd=0.2
ø=0.1
P₀=100
D₀=50
T=24
t=collect(0:0.1:100)
P1=plot(t,t->P₀*(1+pp-dd-ø)^(t/T),label="P",seriestype=:line,ylims = (0,4300))
plot!(t,t->D₀+P₀*(((1+pp-dd-ø)^(t/T))-1)*((1-pp+dd-ø)/(pp-dd-ø)),label="D",seriestype=:line,ylims = (0,400))
plot!(t,t->P₀*(((1+pp-dd-ø)^(t/T))-1)*(ø/(pp-dd-ø)),label="ø",seriestype=:line,ylims = (0,400))



pp=0.0001
dd=0.00001
ø=0.1
P2=plot(t,t->P₀*(1+pp-dd-ø)^(t/T),label="P",seriestype=:line,ylims = (0,4300))
plot!(t,t->D₀+P₀*(((1+pp-dd-ø)^(t/T))-1)*((1-pp+dd-ø)/(pp-dd-ø)),label="D",seriestype=:line,ylims = (0,400))
plot!(t,t->P₀*(((1+pp-dd-ø)^(t/T))-1)*(ø/(pp-dd-ø)),label="ø",seriestype=:line,ylims = (0,400))


pp=0.0001
dd=0.3
ø=0.1
P3=plot(t,t->P₀*(1+pp-dd-ø)^(t/T),label="P",seriestype=:line,ylims = (0,4300))
plot!(t,t->D₀+P₀*(((1+pp-dd-ø)^(t/T))-1)*((1-pp+dd-ø)/(pp-dd-ø)),label="D",seriestype=:line,ylims = (0,400))
plot!(t,t->P₀*(((1+pp-dd-ø)^(t/T))-1)*(ø/(pp-dd-ø)),label="ø",seriestype=:line,ylims = (0,400))

plot(P1,P2,P3,layout=(1,3),legend=true,size = (800, 500))

# ╔═╡ f0413f68-1eb4-4a3d-bdcf-61e8aa96c0e7
md"### testing robustness

Developing systems need to be robust, we can check how sensitive are this type of differentiation systems to perturbations


"

# ╔═╡ 54dd5337-c76f-4a23-8bec-a0cd57b433dd
md"compared with exponential growth"

# ╔═╡ 276f3941-8637-4274-b98f-6e86a4d1329f
begin
	plot(x -> log2(1+x), collect(-1:0.1:1),seriestype=:line,xlabel=("amount of perturbation"),ylabel=("change in terms of cell cycle units"))
	plot!(x -> log(1+x)/log(1+pp-dd), collect(-1:0.1:1),seriestype=:line,xlabel=("amount of perturbation"),ylabel=("change in terms of cell cycle units"),xlims = (-0.25,0.25),ylims = (-1,1))
end

# ╔═╡ 98fd21ec-3654-44a8-b1fa-dbf74c9bab6e
md"let's find out what is teh cell cycle to reach in 10 hours from 10 to 100 cells" 

# ╔═╡ 662d5c9d-bfb7-4ad9-a497-e3d6344b02af
1*log(1+pp-dd)/log(20/10)

# ╔═╡ Cell order:
# ╠═d096a6be-65a1-428d-9bfb-da7fe89f4c19
# ╟─7f30ad4d-9964-4140-8003-051653d5f1e4
# ╠═470fe58c-4020-4333-861f-a6afe54a9e53
# ╟─ede201b5-8f82-4cb0-be2a-be03b9140c50
# ╟─fffe908a-a8ff-404b-9138-bb23210a5fae
# ╟─79536c27-2f74-40a2-ae9f-b44ed021f208
# ╟─20980ae1-0efa-4836-96e8-8ce71f513fae
# ╟─aecc1472-9167-4b09-a08c-9e701def7d54
# ╠═497ebc18-b94a-420d-a21a-708d699dec5c
# ╠═ae184b0d-f61d-4aea-853b-17cb597d1087
# ╟─ca6a3577-ac75-4993-b2e5-4975ef4aaf9f
# ╟─3e8066bd-d8dc-47d9-98e5-0b8b0759e179
# ╠═2e9621c6-50c4-4d30-84f5-f26ea707a808
# ╠═f0413f68-1eb4-4a3d-bdcf-61e8aa96c0e7
# ╠═54dd5337-c76f-4a23-8bec-a0cd57b433dd
# ╠═276f3941-8637-4274-b98f-6e86a4d1329f
# ╠═98fd21ec-3654-44a8-b1fa-dbf74c9bab6e
# ╠═662d5c9d-bfb7-4ad9-a497-e3d6344b02af
