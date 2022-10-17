### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 0c768bc9-4148-48e5-b6bf-637d25c57fd3
begin
	using Images, Plots, DifferentialEquations
	img_water=load("/home/user/img/water.jpg");
	img_solar = load("https://media.wired.com/photos/5934776af061de0423ccdf98/191:100/pass/solar_system.jpg");
	img_pathway=load("https://www.spandidos-publications.com/article_images/or/29/1/OR-29-01-0003-g01.jpg");
	img_brain=load("/home/user/img/brain.jpg");
	img_clock=load("/home/user/img/Clockwork.jpg");
	img_dices=load("/home/user/img/dices.jpg");
	img_plant=load("/home/user/img/plant.jpg");
	img_complexity=load("/home/user/img/complexity.jpg");
	
end

# ╔═╡ 37367e24-963e-11ec-0d47-ff893b72b17b
md" # 1. Complex Systems

## 1.1 What is a Complex System?

Before defining what do we understand by __complex system__, let’s define first what do we understand as a __system__.

 $img_solar

*System*: A group of interacting, interrelated, or interdependent elements forming a complex whole.

The parts that form a system do not need to be equivalent. i. e. , it can be diverse.

$img_water

a set of water molecules is not diverse

$img_pathway
 a set of of inteeracting proteins in a signaling pathway is diverse

Now that we know what is a system, we are ready to define what is a __Complex System__

__Complex System__: are a specific type of systems that behave in a way that cannot be inferred from the properties of their individual pants. 

$img_brain The brain is a perfect example of a complex system, For instance, conciusness cannot be inferred from the behaviors of isolated neurons

"

# ╔═╡ 1c107775-c064-4be0-a683-8410cd290086
md" In this sence, it is quite dificuly to define what makes a system __complex__. 
Lets start by clarifying the difference between Complicated and complex: a machine can be complicated, but its function and properties can be fully predicted by the properties of their parts:

- It has been rationally designed and built to do a particular task. 
- It is not robust, if one part fails the system fails. 
- It does not adapt to changes. it is not __fluid__.

$img_clock *A clock is not a complex system, it is just a complicated system*
"

# ╔═╡ 8b5e6428-736a-4c59-8df7-538b8c4e6ec2
md" ## 1.2 Characteristics of complex systems

### i). A complex system is _not normal_: 

In the sense that they do not fit normal gaussian distributions. 

Lets illustrate this concept of normal, versus not-normal system. Let's for instance, play a game where we throw two dices and compute the sum. Repeat the experiments an number of times and plot the results in a histogram, 
$img_dices "

# ╔═╡ 881f031e-7fc8-47e5-a357-dca3db43f538
begin
	🎲_slide = @bind 🎲 html"<input type=range min=100 max=10000 step=100>"

	
	md"""
	Move the slider to increase the number of times you throw the 🎲
	
	Number of attempts: $(🎲_slide)

	"""
end

# ╔═╡ 499cdd4a-7e24-4120-be8c-a40956b7de4d
begin
	total_normal = Array{Float64, 2}(undef, 1, 🎲);
	for i=1:🎲
	    value1=rand((1,2,3,4,5,6))
	    value2=rand((1,2,3,4,5,6))
	    total_normal[i]=value1+value2
	end
end

# ╔═╡ e0e2382f-fbfc-4aad-87cc-89617531381e
begin
	p1=histogram(total_normal[:],bins = range(2,12, step = 1));
	title!("Normal")
	xlabel!("value")
	ylabel!("counts");
end

# ╔═╡ 41c92a9a-0cd3-45d0-8c0e-4b64578f621a
md"Now, lets repeat the same experiment, but set the rule that if two dices are similar, do it again. Thsi a sort of simple feedback loop, where the output affect the input. "

# ╔═╡ a2057b67-3e43-4ef8-81e2-fe89f2b7c57c
begin
	total_not_normal = Array{Float64, 2}(undef, 1, 🎲);
for i=1:🎲
    value1=rand((1,2,3,4,5,6))
    value2=rand((1,2,3,4,5,6))
    if  value1==value2
        value1=rand((1,2,3,4,5,6))
        value2=rand((1,2,3,4,5,6))
    end  
    total_not_normal[i]=value1+value2
end
end

# ╔═╡ 8fab198c-496b-41a9-95fc-d2a886679f54
begin
	p2=histogram(total_not_normal[:],bins = range(2,12, step = 1));
	title!("not Normal")
	xlabel!("value")
	ylabel!("counts");
	plot(p1,p2,layout=(1,2),legend=false)
end

# ╔═╡ f9cba1c0-ced9-406f-a9d2-1e43097e58c0
md" we can see that the two distributions are clearly different, the left one is more similar to a typical gaussian distribution, while the right one, despite having the same average value, has different  tails.

### ii). Their response  can be difficult to predict 

For instance, one of the most ssimplest quantilites of a complex system is teh dependence on inital conditions

"

# ╔═╡ a68b0146-6bcf-4b1b-9804-2c5dbefdcf0a
NoFeedback! = @ode_def ab2 begin
   dM = -γ_M*M+α_M*T^n/(K^n +T^n)
   dP =   α_P * M - γ_P * P
    end α_M γ_M T n α_P γ_P K

# ╔═╡ 010f87b7-a625-466c-b11f-bd11e4454f6a
PositiveFeedback! = @ode_def ab begin
   dM = -γ_M*M+α_M*P^n/(K^n +P^n)
   dP =   α_P * M - γ_P * P
    end α_M γ_M T n α_P γ_P K

# ╔═╡ 06f53c91-8e22-4f73-a6ae-1d44de7c9308
begin
	dog_slide = @bind 🐶 html"<input type=range min=0.01 max=0.5 step=0.1>"
	cat_slide = @bind 🐱 html"<input type=range min=0.01 max=0.5 step=0.1>"
	
	md"""
	**How many molecules do you have?**
	
	Initial concentration of a: $(dog_slide)
	
	Initial concentration of b: $(cat_slide)
	"""
end

# ╔═╡ aa37be5e-507f-4269-8235-accd44e55093
begin
	u₀ = [🐶,🐱]
	tspan = (0.0,50.0)
	n=3
	K=1
	k_m=1
	D=1
	T=3
	α_M=k_m*D
	γ_M=0.1
	α_P=0.6
	γ_P=0.5
	p=[α_M,γ_M,T,n,α_P,γ_P,K];
	
prob1 = ODEProblem(PositiveFeedback!,u₀,tspan,p)
prob2 = ODEProblem(NoFeedback!,u₀,tspan,p)

sol1 = solve(prob1)
P1=plot(sol1,label=["mRNA" "Protein"],ylims = (0,12))
title!("Positive feedback")
xlabel!("Time [s]")
ylabel!("Concentration [M]")

sol2 = solve(prob2)
P2=plot(sol2,label=["mRNA" "Protein"],ylims = (0,12))
title!("Linear")
xlabel!("Time [s]")
ylabel!("Concentration [M]")

plot(P1,P2,layout=(1,2),legend=true)
end

# ╔═╡ 6a3a5f82-df35-4d53-b0bd-1d6c45798ef4
md" ### iii) Complex systems can be __robust__
Complex systems can have the property of beiong robust, in the sense that they are less affected by changes in the environment, compared to simple systems. A perfect example of thsi robustness is the property called __adaptation__"

# ╔═╡ 19aa372b-52ef-430e-9ee5-f03e0ce63cf0
FeedForward2! = @ode_def ab3 begin
   dM1 = -γ_M*M1+α_M*T^n/(K^n +T^n)
   dP1 =   α_P * M1 - γ_P * P1 
   dM2 = -γ_M*M2+α_M*T^n/(K^n +T^n)*K2^n/(K2^n +P1^n)
   dP2 =   α_P * M2 - γ_P * P2 
    end α_M γ_M T n α_P γ_P K K2

# ╔═╡ 04ad7a1e-e546-4c81-b65a-09ddaa85aee5
begin
	TT = @bind 😺 html"<input type=range min=1.0 max=5.0 step=1.0>"
	
	md"""
	**slide to change amount of transcription factor T**
	
	Initial concentration of transcription factor : $(TT)
	
	
	"""
end

# ╔═╡ beb934fd-8daf-4b10-8712-a5f9e7fe8ba2
prob5 = ODEProblem(FeedForward2!,[0.01,0.00001,0.01,0.00001],tspan,[α_M,γ_M,😺,n,α_P,γ_P,K,0.1]);

# ╔═╡ afcbb8d7-30b6-4229-8c58-b5a321e818db
begin
	sol5 = solve(prob5)
	P5=plot(sol5,label=["mRNA1" "Protein1" "mRNA2" "Protein2"],ylims = (0,12))
	title!("Linear")
	xlabel!("Time [s]")
	ylabel!("Concentration [M]")
end

# ╔═╡ 76340e6c-a438-4c2e-b536-0e9f3da22193
md"We can see that higher values of the amount of trasncription factor in thsi network do not affect the dynamics. The system is almost insensitive and very robust aginast changes in this parameter"

# ╔═╡ f2503077-3550-472d-a58f-84d2ee4a5f09
md" ### iv) A complex system is characterized by __global__ scale properties 

A complex system, despite being formed by many interacting parts,  seems to behave as one single entity. The typical example is the flock of birds. The 2021 nobel price in Physics was awarded to Giorgio Parisi, for his contribution the study of thesy type of coordinated dynamics. _They may seem very far from spin glasses, but there is something in common, What they share, and what is very interesting, is how complex behaviors arise. This is a theme recurrent in physics and biology, and most of the research that I have done is to get at this thing: how complex collective behavior may arise from elements that each have a simple behavior_. "

# ╔═╡ 2b5d86fa-2c9a-4ef3-a3d6-24f2fb59c213
html"""
<iframe width="700" height="400" src="https://www.youtube.com/embed/V4f_1_r80RY" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
"""

# ╔═╡ f111e935-b711-439c-bec4-5e53d224b571
md" ### v) A complex system is characterized by __emergence__ of properties 

The basic idea of emergence is that there are properties at the upper hierarchical levels of nature that are not derivable from or reducible to the properties and laws of the lower levels. Emergence is the opposte concept or approach to Reductionism, that by contrast, argues that everything can be explained by (reduced to) the basic laws of physics. 

For instance, information is neither matter nor energy, although it needs matter to be embodied and energy to be communicated. How information travels inside a system is key to understantd the system. But you cannot study how information propagartes by looking at a part of a system in isolation. So, you need to study the interactions. And from the interactions emerge new properties. One of teh manin ontributors to thsi idea of emergence is Illa Prigogine (recipient of the Nobel Price in 1977)

```math
\begin{align*}
System > \sum_{i} part_i \tag{16} \\
\end{align*}
```
in fact 

```math
\begin{align*}
System = \sum_{i} part_i +  interactions \tag{17} \\
\end{align*}
```

A System is more than the sum of the parts, a system is parts + interactions 
$img_plant Life is a good example of an emerging property"

# ╔═╡ Cell order:
# ╟─0c768bc9-4148-48e5-b6bf-637d25c57fd3
# ╟─37367e24-963e-11ec-0d47-ff893b72b17b
# ╟─1c107775-c064-4be0-a683-8410cd290086
# ╟─8b5e6428-736a-4c59-8df7-538b8c4e6ec2
# ╟─881f031e-7fc8-47e5-a357-dca3db43f538
# ╠═499cdd4a-7e24-4120-be8c-a40956b7de4d
# ╠═e0e2382f-fbfc-4aad-87cc-89617531381e
# ╟─41c92a9a-0cd3-45d0-8c0e-4b64578f621a
# ╠═a2057b67-3e43-4ef8-81e2-fe89f2b7c57c
# ╠═8fab198c-496b-41a9-95fc-d2a886679f54
# ╟─f9cba1c0-ced9-406f-a9d2-1e43097e58c0
# ╠═a68b0146-6bcf-4b1b-9804-2c5dbefdcf0a
# ╠═010f87b7-a625-466c-b11f-bd11e4454f6a
# ╟─06f53c91-8e22-4f73-a6ae-1d44de7c9308
# ╠═aa37be5e-507f-4269-8235-accd44e55093
# ╟─6a3a5f82-df35-4d53-b0bd-1d6c45798ef4
# ╠═19aa372b-52ef-430e-9ee5-f03e0ce63cf0
# ╟─04ad7a1e-e546-4c81-b65a-09ddaa85aee5
# ╠═beb934fd-8daf-4b10-8712-a5f9e7fe8ba2
# ╠═afcbb8d7-30b6-4229-8c58-b5a321e818db
# ╟─76340e6c-a438-4c2e-b536-0e9f3da22193
# ╟─f2503077-3550-472d-a58f-84d2ee4a5f09
# ╟─2b5d86fa-2c9a-4ef3-a3d6-24f2fb59c213
# ╟─f111e935-b711-439c-bec4-5e53d224b571
