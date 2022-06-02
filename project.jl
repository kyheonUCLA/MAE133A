### A Pluto.jl notebook ###
# v0.19.5

using Markdown
using InteractiveUtils

# ╔═╡ 3b5aaae0-e20e-11ec-0ea1-f39b52d3602e
begin
	using Unitful
	using CoolProp
end

# ╔═╡ 94fe6dee-89b4-4e6a-8af7-f35b049066d4
# Adiabatic H2 Compressor
begin
	mutable struct Pump 
		W::Any
		ratio::AbstractFloat
		Pin::Any
		Tin::Any
	end
end

# ╔═╡ 3a036fda-e16f-4cbf-9bc5-372f575af9a9
# Isenthalpic CH4 Valve 
struct Valve
	Tin::Any
	Pin::Any
	Pex::Any
end

# ╔═╡ 719eef9e-ac7a-4b2d-bb05-19740d43c9ba
# Heat Exhanger does no work
struct Cooler
	Tin::Any
	Tex::Any
	Pin::Any
	Pex::Any
end

# ╔═╡ a7f4fe29-6485-438b-94d4-835ff2e97750
# CH4 decomposition reactor
# Reactor does no work
struct Reactor
	Tin::Any
	Tex::Any
	P::Any
end

# ╔═╡ e3a33cab-f8f3-4434-b087-6321698e975b
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
	V_tank = 1u"m^3"
	T_pyro = 1600u"K"

	Q_solar = 15u"kW"
	H_rxn = 74.6u"kJ/mol"
	Hf_meth = -74850u"kJ/kmol"
	M_h2 = 2*1.00794u"g/mol"
	M_meth = 16.04u"g/mol"
	
	## Initial Calcuations
	n_dot_meth = uconvert(u"kmol/s", 0.1 * Q_solar / H_rxn)
	n_dot_meth |> println
	n_dot_h2 = 2 * n_dot_meth
	n_dot_c = n_dot_meth
	
end


# ╔═╡ 5dd18145-5491-4137-9194-d4fe10b76816
ans3 = let
	x=1u"K"
end

# ╔═╡ 7841815b-c3a5-4c5c-8970-cf18cfba676d
md"""
3) For your compressor configuration in part A, estimate the overall isentropic efficiency of the compression process.

Answer: $(round(typeof(ans3), ans3,digits=2))
"""

# ╔═╡ 172aa52c-1a62-435a-876f-91a6bf543cb5
ans4 = let
	R = 8.3145u"J/mol/K"
	n_H2 = 1u"kg"/M_h2
	P_H2 = n_H2 * R * T_tank / V_tank
	P_final = uconvert(u"bar", P_tank_min + P_H2)
end

# ╔═╡ c5e700ce-7b4b-4b9a-987c-f1b98616efde
md"""
4) Calculate the final pressure of the tube cylinder at the end of one day.

Answer: $(round(typeof(ans4), ans4,digits=2))
"""

# ╔═╡ bb15cd28-6812-4402-bffe-a0c849528643
ans10 = let
	x = 5u"K"
end

# ╔═╡ 0965a736-3af7-478d-b07c-e42d9b939473
md"""
10) Calculate the total exergy destroyed as the pressurized hydrogen loses heat from the tank to the surroundings

Answer: $(round(typeof(ans10),ans10,digits=2))
"""

# ╔═╡ 1d65754b-57ec-420d-a5be-adea52476969
function pumpTemperature(pump::Pump)
	Hin = uconvert(u"kJ/kg", PropsSI("H", "T", pump.Tin, "P", pump.Pin, "H2"))
	Hex = uconvert(u"kJ/kg", Hin - pump.W/(n_dot_h2*M_h2))	
	Pex = pump.ratio*pump.Pin
	return uconvert(u"K", PropsSI("T", "H", Hex, "P", Pex, "H2"))
end

# ╔═╡ 5dd20bf6-6416-4214-bbae-a3ef287a47b5
function pumpEfficiency(pump::Pump)
	Hin = PropsSI("H", "T", pump.Tin, "P", pump.Pin, "H2")
	Hex = Hin - pump.W/(n_dot_h2*M_h2)
	S = PropsSI("S", "T", pump.Tin, "P", pump.Pin, "H2")
	Hs = PropsSI("H", "S", S, "P", pump.Pex, "H2")
	η = uconvert(u"K/K", (Hs - Hin)/(Hex - Hin))
	return η
end

# ╔═╡ 2ff23a31-b092-4841-bf5c-d92076c61517
function pumpExergyDestroyed(pump::Pump, Tₒ)
	Tex = pumpTemperature(pump)
	Pex = pump.Pin*pump.ratio 
	Sin = PropsSI("Smolar", "T", pump.Tin, "P", pump.Pin, "H2")
	Sex = PropsSI("Smolar", "T", Tex, "P", Pex, "H2")
	σ = n_dot_h2*(Sex - Sin)
	return uconvert(u"kW", Tₒ*σ)
end

# ╔═╡ d79ef2bf-9cba-4f52-9a68-15ac05c68d3d
function pumpExergy(pump::Pump, Tₒ)
	return uconvert(u"kW", -pump.W - pumpExergyDestroyed(pump, Tₒ))
end

# ╔═╡ a5099c86-0ad1-4b6d-9574-8a4b2aec492d
ans9 = let
	pump = Pump(-0.5u"kW", 10, 0.1u"bar", 300u"K")
	pumpTemperature(pump)
	# pumpEfficiency(pump)
	pumpExergyDestroyed(pump, 298u"K")
	pumpExergy(pump, 298u"K")
end

# ╔═╡ 62634867-2e70-43bb-a0cd-b9629723f244
md"""
9) Calculate the rate of exergy destruction in the compression process.

Answer: $(round(typeof(ans9),ans9,digits=2))
"""

# ╔═╡ 7dc1dcf2-9cf2-4c34-b10d-0900937e4d2e
function valveTemperature(valve::Valve)
	Hin = PropsSI("H", "T", valve.Tin, "P", valve.Pin, "CH4")
	return uconvert(u"K", PropsSI("T", "H", Hin, "P", valve.Pex, "CH4"))
end

# ╔═╡ b91c94fb-8d7e-41ab-b642-0bbeca3d1ce4
function valveExergyDestroyed(valve::Valve, Tₒ)
	Tex = valveTemperature(valve)
	Sin = PropsSI("Smolar", "T", valve.Tin, "P", valve.Pin, "CH4")
	Sex = PropsSI("Smolar", "T", Tex, "P", valve.Pex, "CH4")
	σ = n_dot_meth*(Sex-Sin)
	destroyed = Tₒ*σ
	return uconvert(u"kW", destroyed)
end

# ╔═╡ 2041825c-155c-4ad5-9360-706ecf96ea46
ans5 = let
	Tₒ = 298u"K"
	valve = Valve(T_in, P_in, P_react)
	valveExergyDestroyed(valve, Tₒ)
end

# ╔═╡ c8ac1365-8599-4ab6-b4be-8626fa63b225
md"""
5) Calculate the rate of exergy destruction at the mass flow controller at nominal efficiency

Answer: $(round(typeof(ans5), ans5,digits=2))

"""

# ╔═╡ 9d5f5420-b801-4e14-b061-02a234598653
function valveExergy(valve::Valve, Tₒ)
	Tex = valveTemperature(valve)
	Hin = PropsSI("Hmolar", "T", valve.Tin, "P", valve.Pin, "CH4")
	Hex = PropsSI("Hmolar", "T", Tex, "P", valve.Pex, "CH4")
	return uconvert(u"kW", n_dot_meth*(Hex-Hin)-valveExergyDestroyed(valve, Tₒ))
end

# ╔═╡ 96d3e8a4-9a50-4e17-bddd-4ba7a4714d1c
function reactorChemicalExergy(react::Reactor, Tₒ)
	G_meth = PropsSI("Gmolar", "T", react.Tin, "P", react.P, "CH4")
	G_h2 = PropsSI("Gmolar", "T", react.Tex, "P", react.P, "H2")
	R = 8.3145u"J/mol/K"
	value = -R*Tₒ*log((0.666)^2)
	Ex_chem = n_dot_meth*G_meth - n_dot_h2*G_h2 + n_dot_meth*value
	uconvert(u"kW", Ex_chem)
end

# ╔═╡ 4c0a7f58-a659-4cda-9037-335bbc56d953
function coolerHeat(hx::Cooler)
	Hin = PropsSI("Hmolar", "T", hx.Tin, "P", hx.Pin, "H2")
	Hex = PropsSI("Hmolar", "T", hx.Tex, "P", hx.Pex, "H2")
	return uconvert(u"kW", n_dot_h2*(Hex-Hin))
end

# ╔═╡ 554e0545-78ae-4fb3-9bbc-d00edcc91045
function molar_enthalpy(T, P, name::String, Hf=0u"kJ/kmol")
	Href = PropsSI("Hmolar", "T", 298.15u"K", "P", 1u"atm", name)
	H =  PropsSI("Hmolar", "T", T, "P", P, name)
	return uconvert(u"kJ/kmol", Hf + H - Href)
end

# ╔═╡ 46ba6b8a-c7e6-4609-a301-a05eaae599f0
function reactorHeat(react::Reactor) 
	Hin = molar_enthalpy(react.Tin, react.P, "CH4", Hf_meth)
	Hex = molar_enthalpy(react.Tex, react.P, "H2")
	Q =  n_dot_h2*Hex-n_dot_meth*Hin
	return uconvert(u"kW", Q)
end

# ╔═╡ 1ab27cb2-34db-4501-89d0-1256d6060705
function reactorExergyDestroyed(react::Reactor, Tₒ)
	Sin = PropsSI("Smolar", "T", react.Tin, "P", react.P, "CH4")
	Sex = PropsSI("Smolar", "T", react.Tex, "P", react.P, "H2")
	Tb = 298u"K"
	Q = reactorHeat(react)
	σ = n_dot_h2*Sex-n_dot_meth*Sin-Q/Tb
	return uconvert(u"kW", Tₒ*σ)
end

# ╔═╡ 8a6a970b-491c-41b3-ac76-8c82fa87077e
ans6 = let
	# Obviously wrong
	valve = Valve(T_in, P_in, P_react)
	T_in_to_react = valveTemperature(valve)
	react = Reactor(T_in_to_react, T_pyro, P_react)
	reactorExergyDestroyed(react, 298u"K")
	
end

# ╔═╡ 3a582d65-a2f9-4ccf-aab0-e6d38acbbd50
md"""
6) Calculate the rate of exergy destruction in the reactor at nominal efficiency

Answer: $(round(typeof(ans6), ans6,digits=2))
"""

# ╔═╡ 31e31fa0-bd10-4ee2-9222-be7283884dc7
ans7 = let
	Ed = ans6 #getting reactor exergy
	Q_c = 1u"kg"
	# need to use mc delta T
end

# ╔═╡ 4a1b2c62-0048-4f58-a597-c0e728c9bdd0
md"""
7) What fraction of this reactor exergy destruction is consumed by cooling the graphite from 1600 K to room temperature?

Answer: $(round(typeof(ans7), ans7,digits=2))
"""

# ╔═╡ 45e0cbfe-06a3-4504-8ebe-79b7915e7e22
ans8 = let
	Ed = ans6
	hx = Cooler(T_pyro, 298u"K", P_react, P_react) #isobaric cooling
	Q_h2 = coolerHeat(hx)
	fraction = Q_h2/Ed
end

# ╔═╡ 61691472-8282-4588-8ea1-daac04a472d2
md"""
8) What fraction of this reactor exergy destruction is consumed by cooling the product hydrogen from 1600 K to room temperature?

Answer: $(round(ans8,digits=2))
"""

# ╔═╡ f4a60594-ebec-4c1f-8c60-66a1215c4784
function reactorExergy(react::Reactor, Tₒ)
	Hin = molar_enthalpy(react.Tin, react.P, "CH4", Hf_meth)
	Hex = molar_enthalpy(react.Tex, react.P, "H2")
	Ex_chem = reactorChemicalExergy(valve, Tₒ)
	Ex_d = reactorExergyDestroyed(valve, Tₒ) 
	return uconvert(u"kW", n_dot_h2*Hex-n_dot_meth*Hin + Ex_chem - Ex_d)
end

# ╔═╡ a5f8c399-a773-4868-b91b-3c8c7ac35103
function get_flow_rate(m_dot, T=273.15u"K", P=1u"atm")
	ρ = PropsSI("D", "T", T, "P", P, "CH4")
	return uconvert(u"ft^3/minute", m_dot/ρ)
end

# ╔═╡ ee854bbd-10bf-465c-8346-f6725aebdc5b
ans2 = let
	scfm = get_flow_rate(n_dot_meth*M_meth)
	uconvert(u"cm^3/minute", scfm)
end

# ╔═╡ 7a13c45f-bc2e-452a-9c5c-2f455aaf5597
md"""
2) What size of methane mass flow controller (in sccm) is required at nominal efficiency? Find a suitable product. 

Answer: $(round(typeof(ans2), ans2,digits=2))
"""

# ╔═╡ 47f72e30-757e-4806-b189-2d79fdaa96b2
md"""

Resources: https://www.hydrogen.energy.gov/pdfs/9013_energy_requirements_for_hydrogen_gas_compression.pdf

"""

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
# ╠═3b5aaae0-e20e-11ec-0ea1-f39b52d3602e
# ╠═94fe6dee-89b4-4e6a-8af7-f35b049066d4
# ╠═3a036fda-e16f-4cbf-9bc5-372f575af9a9
# ╠═719eef9e-ac7a-4b2d-bb05-19740d43c9ba
# ╠═a7f4fe29-6485-438b-94d4-835ff2e97750
# ╠═e3a33cab-f8f3-4434-b087-6321698e975b
# ╠═7a13c45f-bc2e-452a-9c5c-2f455aaf5597
# ╠═ee854bbd-10bf-465c-8346-f6725aebdc5b
# ╟─7841815b-c3a5-4c5c-8970-cf18cfba676d
# ╠═5dd18145-5491-4137-9194-d4fe10b76816
# ╟─c5e700ce-7b4b-4b9a-987c-f1b98616efde
# ╠═172aa52c-1a62-435a-876f-91a6bf543cb5
# ╟─c8ac1365-8599-4ab6-b4be-8626fa63b225
# ╠═2041825c-155c-4ad5-9360-706ecf96ea46
# ╟─3a582d65-a2f9-4ccf-aab0-e6d38acbbd50
# ╠═8a6a970b-491c-41b3-ac76-8c82fa87077e
# ╟─4a1b2c62-0048-4f58-a597-c0e728c9bdd0
# ╟─31e31fa0-bd10-4ee2-9222-be7283884dc7
# ╟─61691472-8282-4588-8ea1-daac04a472d2
# ╠═45e0cbfe-06a3-4504-8ebe-79b7915e7e22
# ╟─62634867-2e70-43bb-a0cd-b9629723f244
# ╠═a5099c86-0ad1-4b6d-9574-8a4b2aec492d
# ╟─0965a736-3af7-478d-b07c-e42d9b939473
# ╠═bb15cd28-6812-4402-bffe-a0c849528643
# ╠═1d65754b-57ec-420d-a5be-adea52476969
# ╠═5dd20bf6-6416-4214-bbae-a3ef287a47b5
# ╠═d79ef2bf-9cba-4f52-9a68-15ac05c68d3d
# ╠═2ff23a31-b092-4841-bf5c-d92076c61517
# ╠═7dc1dcf2-9cf2-4c34-b10d-0900937e4d2e
# ╠═b91c94fb-8d7e-41ab-b642-0bbeca3d1ce4
# ╠═9d5f5420-b801-4e14-b061-02a234598653
# ╠═46ba6b8a-c7e6-4609-a301-a05eaae599f0
# ╠═96d3e8a4-9a50-4e17-bddd-4ba7a4714d1c
# ╠═f4a60594-ebec-4c1f-8c60-66a1215c4784
# ╠═1ab27cb2-34db-4501-89d0-1256d6060705
# ╠═4c0a7f58-a659-4cda-9037-335bbc56d953
# ╠═554e0545-78ae-4fb3-9bbc-d00edcc91045
# ╠═a5f8c399-a773-4868-b91b-3c8c7ac35103
# ╠═47f72e30-757e-4806-b189-2d79fdaa96b2
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
