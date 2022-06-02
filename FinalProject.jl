### A Pluto.jl notebook ###
# v0.19.5

using Markdown
using InteractiveUtils

# ╔═╡ 1371b6b0-d960-11ec-229a-8113497269c9
begin
	using CoolProp
	using Unitful
end

# ╔═╡ 6532de65-0460-480f-92c8-9c34de2d49c1
# Assumptions
begin
	
	# Reactions is complete: CH4(g) --> 2H2(g) + C(s)
	# Isenthalpic throttling valve
	
	# Given Constants
	P_in = 1u"atm"
	T_in = 298u"K"
	P_react = 10u"kPa"
	T_tank = 300u"K"
	P_tank_min = 200u"bar"
	P_tank_max = 220u"bar"
	T_pyro = 1600u"K"

	Q_solar = 15u"kW"
	H_rxn = 74.6u"kJ/mol"
	Hf_meth = -74850u"kJ/kmol"
	M_h2 = 2*1.00794u"g/mol"
	M_meth = 16.04u"g/mol"
	
	## Initial Calcuations
	n_dot_meth = uconvert(u"kmol/s", 0.1 * Q_solar / H_rxn)
	n_dot_h2 = 2 * n_dot_meth
	n_dot_c = n_dot_meth
	
end


# ╔═╡ 7a998866-7b3d-4379-b24f-84eef4ca2c61
## Analysis of Valve Exergy
let
	P1 = P_in
	P2 = P_react
	T1 = T_in
	T2 = 297.56999u"K"

	S1 = PropsSI("Smolar", "T", T1, "P", P1, "CH4")
	S2 = PropsSI("Smolar", "T", T2, "P", P2, "CH4")

	σ = uconvert(u"kW/K", n_dot_meth*(S2 - S1))
end

# ╔═╡ 163b31bf-8a16-4f55-aa13-2b3c18b23df6
function molar_enthalpy(T, P, name::String, Hf=0u"kJ/kmol")
	Href = PropsSI("Hmolar", "T", 298.15u"K", "P", 1u"atm", name)
	H =  PropsSI("Hmolar", "T", T, "P", P, name)
	return uconvert(u"kJ/kmol", Hf + H - Href)
end

# ╔═╡ 07058eee-9cd4-47f9-a216-a26f2b3c1dc7
## Analysis of valve
# Need to 
let
	P1 = P_in
	P2 = P_react
	T1 = T_in
	H1 = molar_enthalpy(T1, P1, "CH4")
	println(H1)
	println(P1)
	println(P2)
	δ = 0.001u"K"
	T2 = 297u"K"
	while true
		H2 = molar_enthalpy(T2, P2, "CH4")
		if abs(H1 - H2) ≤ 0.1u"kJ/kmol"
			println(H2)
			break
		end
		T2 = T2 + δ
	end
	print(T2)
end

# ╔═╡ ec499bec-41df-4565-ae28-79d90e1bcd78
## Analysis of reactor CV

let
	T1 = T_in
	P1 = P_in
	T3 = 1473u"K" ##assume pyrolysis temp
	P3 = P_react
	Q_in = Q_solar
	H2 = molar_enthalpy(T1, P1, "CH4", Hf_meth)

	H3 = molar_enthalpy(T3, P3, "H2")
	Q_net = H3*n_dot_h2 - H2*n_dot_meth

	Q_loss = uconvert(u"kW", Q_in - Q_net)
end

# ╔═╡ 48911248-52e5-4f69-84d0-9041f9cf0921
## Analysis of reactor CV (not legit)
let
	T1 = T_in
	P1 = P_in
	P3 = P_react
	#using energy balance of reactor CV, no work done in reactor
	H2 = molar_enthalpy(T1, P1, "CH4", Hf_meth)
	LHS = uconvert(u"kW", Q_solar + H2*n_dot_meth)
	println(LHS)
	T3 = 9000u"K"
	δ = 1u"K"
	while true
		RHS = molar_enthalpy(T3, P3, "H2") * n_dot_h2
		if abs(LHS - RHS) ≤ 0.1u"kW"
			println(RHS)
			break
		end
		T3 = T3 + δ
	end
	T3
end

# ╔═╡ dfc30cf6-097f-44ae-865f-e8f5e8199e02
function molar_entropy(T, P, name::String, Sa=0u"kJ/kmol/K")
	S = PropsSI("Smolar", "T", T, "P", P, name)
	
end

# ╔═╡ a307c64b-0cfd-4b6e-87d9-9b06ec1f2dd7
## Analysis of Reactor Exergy

let
	## Thermomechanical Exergy
	Q_net = 2.929811315030972u"kW"
	T2 = T_in
	P2 = P_react
	T3 = 1600u"K" ##assume pyrolysis temp
	P3 = P_react
	Tₒ= 298u"K"
	T_film = (T2 + T3) / 2
	T_film = Tₒ

	H2 = molar_enthalpy(T2, P2, "CH4", Hf_meth)
	H3 = molar_enthalpy(T3, P3, "H2")
	S2 = molar_entropy(T2, P2, "CH4", Hf_meth)
	S3 = molar_entropy(T3, P3, "H2")
	σ = uconvert(u"kW/K", S3*n_dot_h2 - S2*n_dot_meth - Q_net/T_film)
	σ |> println
	ΔEx = uconvert(u"kW", n_dot_meth*(H3-2*H2 - Tₒ*(S3- 2*S3)))
	ΔEx |> println
	
	## Chemical Exergy
	# Gibbs function of formation
	Gf_h2 = 0u"kJ/kmol"
	Gf_meth = -50.7u"kJ/mol"
	Gf_c = 0u"kJ/kmol"

	# Using model 1
	e_ch_c = 404_580u"kJ/kmol"
	e_ch_h2 = 235_250u"kJ/kmol"

	chemical_ex_1 = Gf_meth + 2*Gf_h2 + Gf_c + e_ch_c + e_ch_h2

	# Using model 2 
	e_ch_c = 410_260u"kJ/kmol"
	e_ch_h2 = 236_100u"kJ/kmol"

	chemical_ex_2 = Gf_meth + 2*Gf_h2 + Gf_c + e_ch_c + e_ch_h2
	

end

# ╔═╡ c1c55222-8d78-4428-98b3-c09da7696405
function chemical_exergy(T, P, names::Vector{String}, coeff::Vector{Float64})
	
end

# ╔═╡ f7767aa1-70b9-4ef7-b068-ebf52857e980
function print_values(items...)
	"Printing: " |> println
	for item in items
		item |> println
	end
end

# ╔═╡ b10ca75e-2184-4c21-bce3-8b6359c80d9c
## Analysis of Compressor Work
let
	T3 = T_pyro
	P3 = P_react
	T4 = 2288.778u"K"
	P4 = P_tank_min
	
	H3 = PropsSI("Hmolar", "T", T3, "P", P3, "H2")
	H4 = PropsSI("Hmolar", "T", T4, "P", P4, "H2")
	print_values(H3, H4)

	
	W = uconvert(u"kW", n_dot_h2*(H4 - H3))
end

# ╔═╡ 8e9a5c43-efda-4d72-8b14-dca45800f579
function state_after_pump(η, T1, P1, P2)
	H1 = PropsSI("H", "T", T1, "P", P1, "H2")
	S1 = PropsSI("S", "T", T1, "P", P1, "H2")
	H2s = PropsSI("H", "S", S1, "P", P2, "H2")

	H2r = H1 + (H2s-H1)/η
	return uconvert(u"K", PropsSI("T", "H", H2r, "P", P2, "H2"))
end

# ╔═╡ 60659234-4ec4-4753-8591-e3e1abd01c7f
## Analysis of Compressor Isentropic Efficieny
let
	η = 0.6 #Assume isentropic efficiency is about 60-80%
	T3 = T_pyro
	P3 = P_react
	P4 = 1u"bar" # should be a little bigger than P_tank
	state_after_pump(0.8, 300u"K", 1u"atm", 300u"atm")
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CoolProp = "e084ae63-2819-5025-826e-f8e611a84251"
Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[compat]
CoolProp = "~0.1.0"
Unitful = "~1.11.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f74e9d5388b8620b4cee35d4c5a618dd4dc547f4"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.3.0"

[[CoolProp]]
deps = ["CoolProp_jll", "Markdown", "Unitful"]
git-tree-sha1 = "94062163b5656b1351f7f7a784341b8fe13c1ca1"
uuid = "e084ae63-2819-5025-826e-f8e611a84251"
version = "0.1.0"

[[CoolProp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "1166ddeb518f330ea1708922fd8c5d1a04ec7a73"
uuid = "3351c21f-4feb-5f29-afb9-f4fcb0e27549"
version = "6.4.1+1"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[Unitful]]
deps = ["ConstructionBase", "Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "b649200e887a487468b71821e2644382699f1b0f"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.11.0"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╠═1371b6b0-d960-11ec-229a-8113497269c9
# ╠═6532de65-0460-480f-92c8-9c34de2d49c1
# ╠═07058eee-9cd4-47f9-a216-a26f2b3c1dc7
# ╠═7a998866-7b3d-4379-b24f-84eef4ca2c61
# ╠═ec499bec-41df-4565-ae28-79d90e1bcd78
# ╠═a307c64b-0cfd-4b6e-87d9-9b06ec1f2dd7
# ╠═48911248-52e5-4f69-84d0-9041f9cf0921
# ╠═b10ca75e-2184-4c21-bce3-8b6359c80d9c
# ╠═60659234-4ec4-4753-8591-e3e1abd01c7f
# ╠═163b31bf-8a16-4f55-aa13-2b3c18b23df6
# ╠═dfc30cf6-097f-44ae-865f-e8f5e8199e02
# ╠═c1c55222-8d78-4428-98b3-c09da7696405
# ╠═f7767aa1-70b9-4ef7-b068-ebf52857e980
# ╠═8e9a5c43-efda-4d72-8b14-dca45800f579
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
