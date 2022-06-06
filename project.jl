### A Pluto.jl notebook ###
# v0.19.8

using Markdown
using InteractiveUtils

# ╔═╡ 3b5aaae0-e20e-11ec-0ea1-f39b52d3602e
begin
	using Unitful
	using CoolProp
end

# ╔═╡ 715ed4ad-d7c7-4698-9df1-61b92d8e39f2
md"""
## Group Members

Gianna Braga, Brian Burrous, Ky Heon
"""

# ╔═╡ 22a88a37-1952-4435-a5d6-ec3d6d1f5a13
md"""
## Project setup
We're using a series of compressors with heat exchangers to cool the gas back to room temperature after each compression. The heat exchanger process is assumed to be isobaric. 

Diagram 1
![Hello](https://github.com/kyheonUCLA/MAE133A/blob/main/diagram1.jpeg?raw=true)
Diagram 2

![](https://github.com/kyheonUCLA/MAE133A/blob/main/diagram2.jpeg?raw=true)
"""

# ╔═╡ 94fe6dee-89b4-4e6a-8af7-f35b049066d4
# Adiabatic H2 Compressor
struct Pump 
	W::Any
	ratio::AbstractFloat
	P::Any
	T::Any
	name::AbstractString
	η::Number

end

# ╔═╡ 3a036fda-e16f-4cbf-9bc5-372f575af9a9
# Isenthalpic CH4 Valve 
struct Valve
	Tin::Any
	Pin::Any
	Pex::Any
	name::AbstractString
end

# ╔═╡ 719eef9e-ac7a-4b2d-bb05-19740d43c9ba
# ideal Heat Exhanger, approximatly isobaric cooling
struct Cooler
	Tin::Any
	Tex::Any
	P::Any
	Tb::Any
	name::AbstractString
end

# ╔═╡ a7f4fe29-6485-438b-94d4-835ff2e97750
# CH4 decomposition reactor
# Reactor does no work
struct Reactor
	Tin::Any
	Tex::Any
	P::Any
	Tb::Any
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
	M_carbon = 12u"g/mol"
	C_carbon = 0.71u"J/g/K"
	R = 8.3145u"J/mol/K"
	
	MolarMass = Dict("C"=>12.01u"g/mol", "H2"=>2.02u"g/mol", "CH4"=>16.04u"g/mol")
	MODEL2 = Dict("C"=>410_260u"kJ/kmol", "H2"=>236_100u"kJ/kmol", "CH4"=>831_650u"kJ/kmol")

	## Initial Calcuations
	n_dot_meth = uconvert(u"kmol/s", 0.1 * Q_solar / H_rxn)
	n_dot_meth |> println
	n_dot_h2 = 2 * n_dot_meth
	n_dot_c = n_dot_meth
	n_dot_meth = 0.00688973u"mol/s"
	n_dot_h2 = 2 * n_dot_meth
	
end


# ╔═╡ addb44fe-d353-42f4-94ba-0a8a26bacdc6
begin
	massH2 = 1u"kg"
	massC = 3u"kg"
	
	day = 10u"hr"
	molH2 = massH2/MolarMass["H2"]
	molC = massC/MolarMass["C"]
	molCH4 = molC

	MolarFlow = Dict("C"=>molC/day, "H2"=>molH2/day, "CH4"=>molC/day)
	MolarFlow["C"] = uconvert(u"mol/s", MolarFlow["C"])
	MolarFlow["H2"] = uconvert(u"mol/s", MolarFlow["H2"])
	MolarFlow["CH4"] = uconvert(u"mol/s", MolarFlow["CH4"])
	
end

# ╔═╡ e5903da9-977b-4a30-a51a-2f84b9bb7a6d
md"""
1) Conduct a search of available hydrogen compressors. Does a single compressor exist that serves the required function (gas flow, pressure ratio)? If not, what multiple compressor configuration would satisfy the specifications? Explain and include links where appropriate.
"""

# ╔═╡ 62f35fcc-0333-40ec-952d-f199031db50c
md"""
We will use a series of 4 compressors, assuming compressor ratios of 10, 10, 10, and 2, to get to the desired tank pressure of 200 bar. The first three stanges can be accomplished using the following compressor

https://www.directindustry.com/prod/fornovo-gas-spa/product-223666-2423388.html

The final compressor stage exceeds the maximum pressure of this compressor, so we will also be using this compressor for the final compression stage. 

https://www.directindustry.com/prod/jp-sauer-sohn-maschinenbau-gmbh/product-5707-1982415.html 
"""

# ╔═╡ 66ffbd89-03a3-4931-ab2b-05f828a2973c


# ╔═╡ 7a13c45f-bc2e-452a-9c5c-2f455aaf5597
md"""
2) What size of methane mass flow controller (in sccm) is required at nominal efficiency? Find a suitable product. 

"""

# ╔═╡ ee854bbd-10bf-465c-8346-f6725aebdc5b
ans2 = let
	ρ = PropsSI("D", "T", 273u"K", "P", 1u"atm", "CH4")
	scfm = MolarFlow["CH4"]*MolarMass["CH4"] / ρ
	uconvert(u"cm^3/minute", scfm)
end

# ╔═╡ 7841815b-c3a5-4c5c-8970-cf18cfba676d
md"""
3) For your compressor configuration in part A, estimate the overall isentropic efficiency of the compression process.

"""

# ╔═╡ 54fb5386-ead8-4be0-b539-58c9c959f894
md"""
We couldn't find any good data for our specific hydrogen pump, but various papers suggest an isentropic effeciency of 65-85%. 
"""

# ╔═╡ 5dd18145-5491-4137-9194-d4fe10b76816
# ans3 = let
# 	pump_size = -1u"kW"
# 	pump = Pump(pump_size/MolarFlow["H2"], 10, P_react, 35u"°C", "H2", 0.8)
# 	@show temperature(pump)
# end

# ╔═╡ c5e700ce-7b4b-4b9a-987c-f1b98616efde
md"""
4) Calculate the final pressure of the tube cylinder at the end of one day.

"""

# ╔═╡ 172aa52c-1a62-435a-876f-91a6bf543cb5
ans4 = let
	molH2 = 1u"kg"/M_h2
	P_H2 = molH2 * R * T_tank / V_tank
	P_final = uconvert(u"bar", P_tank_min + P_H2)
end

# ╔═╡ c8ac1365-8599-4ab6-b4be-8626fa63b225
md"""
5) Calculate the rate of exergy destruction at the mass flow controller at nominal efficiency


"""

# ╔═╡ 3a582d65-a2f9-4ccf-aab0-e6d38acbbd50
md"""
6) Calculate the rate of exergy destruction in the reactor at nominal efficiency

"""

# ╔═╡ 4a1b2c62-0048-4f58-a597-c0e728c9bdd0
md"""
7) What fraction of this reactor exergy destruction is consumed by cooling the graphite from 1600 K to room temperature?

"""

# ╔═╡ 31e31fa0-bd10-4ee2-9222-be7283884dc7
ED_C = let
	Tₒ = 298u"K"
	Pₒ = 1u"atm"
	Tb = T_pyro
	massCarbon = 3u"kg"
	ΔT = (T_pyro - Tₒ)
	Q = massCarbon*C_carbon*ΔT/10u"hr"
	uconvert(u"kW", (1-Tₒ/Tb)*Q)
	
	
end

# ╔═╡ 3f75aca3-7c5c-48e6-800d-2b6411171d49
md"""
### Answer: 
3% of exergy destruction consumed by cooling graphite from 1600K to room temperature. 
"""

# ╔═╡ 61691472-8282-4588-8ea1-daac04a472d2
md"""
8) What fraction of this reactor exergy destruction is consumed by cooling the product hydrogen from 1600 K to room temperature?
"""

# ╔═╡ 1b6134f6-1643-4f90-831d-3a25c5c76253
md"""
### Answer
34% of the exergy destruction of the reactor is destroyed cooling H2. 
"""

# ╔═╡ 62634867-2e70-43bb-a0cd-b9629723f244
md"""
9) Calculate the rate of exergy destruction in the compression process.

"""

# ╔═╡ 0965a736-3af7-478d-b07c-e42d9b939473
md"""
10) Calculate the total exergy destroyed as the pressurized hydrogen loses heat from the tank to the surroundings

"""

# ╔═╡ bb15cd28-6812-4402-bffe-a0c849528643
md"""
Our design uses heat exchangers before and after every compressor stage. This means that the gas enters the tank at ambient temperature and tank pressure. Therefore there is not exergy change/no exergy destruction at this point. 

Answer: 0
"""

# ╔═╡ 0fd0f1e2-f6d4-4ed0-b88e-006742e4e40d
md"""
### 11:
If the heat lost by the hydrogen stream between the reaction zone and the compressor inlet could be sent completely back into the reaction zone, then by what percentage could the methane mass flow rate increase? What would be the exergy destruction? Would it increase or decrease? By how much?
"""

# ╔═╡ 09f63e0d-ed35-4a76-9cc0-f27b7a07f532
let
	h1 = PropsSI("Hmolar", "T", T_pyro, "P", P_react, "H2")
	h2 = PropsSI("Hmolar", "T", 298u"K", "P", 200u"bar", "H2") 

	Q_cooling = (h1-h2)*n_dot_h2

	Q_inital = uconvert(u"kW", MolarFlow["CH4"]*H_rxn/0.1)
	Q_final = uconvert(u"kW", Q_inital + Q_cooling)

	n_final_new = uconvert(u"mol/s", Q_final*0.1/H_rxn)
	uconvert(u"K/K", (n_final_new-MolarFlow["CH4"])/MolarFlow["CH4"])*100
end

# ╔═╡ 2374b686-40a3-498f-b044-72248177966e
md"""
Answer: The molar (and mass) flow rate of methane would increase by 10.5%. 
"""

# ╔═╡ be0ad564-d61b-41de-b607-00ac5731a7f1
md"""
### 12.
What are the technical and/or societal benefits and drawbacks of using solar heating in this design, as compared to burning methane or hydrogen for pyrolysis?
"""

# ╔═╡ 7720cdad-3c6b-4f4b-a4fe-691f6546b944
md"""
It is assumed that a solar concentrator is being used for this system. A solar concentrator works by having a curved surface or a series of flat surfaces all angled toward one focus, so that the sunlight bounces of the surfaces and reflects all to the same point. At that point a receiver is placed to collect the energy from the sun. This set up is fairly cheap to manufacture and set up. However, there is no way to store the energy so the solar concentrator can only be used when the sun is out and shining. It would be ideal to position this somewhere that the efficiency of this process is not decreased significantly.

If solar panels were being used, a battery can be attached to store the energy that isn't being used during the day, allowing the process to continue, even when the sun it not present. While having a battery solves the issue of the sun's limited availablility as an energy source, a new issue is created, cost. The website for the Office of Energy Efficiency and Renewable Energy states that solar batteries can range from \$12,000 to $22,000. In general, there is a high upfront cost with solar panels due to the price of purchasing the panels and their installation if that is required. However, this cost could be made up through the money saved when not paying for other energy sources. The US Government is also starting to offer incentives for using solar panels so those may aid in the initial cost of the solar energy set up.

If the hydrogen and graphite being produced for the assumed 10 hours a day is suffcient, and an irregular production rate due to the sun exposure fluctuation is no issue, than a solar concentrator is a great option. The most significant benefit of solar power over other energy sources, such as burning methane, is that the sun is a renewable energy source. It will never run out of heat to provide, at least not for roughly 5 billion more years. It is also a green energy source. The use of solar energy decreases greenhouse gas production which is vital to slowing the effects of climate change.   

Sources:

https://lightningsolar.com.au/comparative-guide-advantages-disadvantages-of-solar-panels/

https://powersolarphoenix.com/commercial-solar-panels-cost/

https://www.energy.gov/eere/solar/articles/should-i-get-battery-storage-my-solar-energy-system

"""

# ╔═╡ 029cc470-a3d6-439f-addf-79ed481a774f
md"""
### 13. 
More broadly, what are the technical and societal benefits and drawbacks of using methane (from natural gas) as a source for graphite and hydrogen?

"""

# ╔═╡ a5ef9be5-3c63-444c-9233-653d06fd0a68
md"""
There is no shortage of methane production in our society. The farming and agriculture industries are huge contributors to the methane content in the atmosphere. When manure is produced by livestock, such as pigs or cattle, methane is released. When the manure is spread onto crops for fertilization, methane is released. A fight has begun to reduce meat and dairy consumption in the US has begun. Alternative milks and vegetarian/vegan options are becoming more widespread than ever before. It is apparent however that much more time and effort is needed to change the general population’s eating habits so methane production continues.

Another significant source of methane in our society is landfills. As organic materials decompose under piles of waste, there are not significant sources of oxygen to facilitate aerobic decomposition. Bacteria that perform anaerobic decomposition come into play. These bacteria are huge sources of methane through their decomposition processes. Composting can help decrease these methane emissions. By collecting organic materials and allowing them to decompose in environments that provide aeration, organisms that do not produce methane can process the materials. Composting is not a widespread habit. Many communities do not have access to composting programs and especially programs that accept food waste, so much of the organic waste still ends up in landfills, producing methane.

Methane is a significant cause for the trapping of heat within Earth’s atmosphere that is leading to climate change. According to UNECE, methane has a warming potential roughly 28-34% higher than that of CO2. Methane capture is a solution being presented to the high emissions of methane from the above sources. In that set up, perforated tubes are arranged through the methane producing areas to collect the methane and transport it to facilities that will burn it. This will prevent much of the methane being produced in industrial settings from reaching the atmosphere. 
 
A downside of burning methane is that CO2 is produced by the process. CO2 in the atmosphere is another cause to the greenhouse effect. However, if the methane being burned is harvested from the atmosphere then the benefits of destroying methane outweighs the drawbacks of CO2 production. The burning of this collected methane is also beneficial because it eliminates reliance on other energy sources that also produce greenhouse gasses. As outlined above, due to the processes used for food production and waste management that this society relies upon, there is an abundance of methane to be utilized for burning and heating. If it is not utilized, it is only making the climate change issue worse.

Sources: 

https://www.epa.gov/lmop/basic-information-about-landfill-gas

https://unece.org/challenge
"""

# ╔═╡ 1d65754b-57ec-420d-a5be-adea52476969
begin
	function temperature(pump::Pump)
		k = 1.4
		polytropicExp = (k-1)/k
		Ts = uconvert(u"K", pump.T)*(pump.ratio)^polytropicExp
		Hs = PropsSI("Hmolar", "T", Ts, "P", pump.P*pump.ratio, pump.name)
		H = PropsSI("Hmolar", "T", pump.T, "P", pump.P, pump.name)
		Hr = H + (Hs-H)/pump.η
		P2 = pump.P*pump.ratio
		T = PropsSI("T", "P", P2, "Hmolar", Hr, pump.name)
		return uconvert(u"K", T) 
	end
	
	function temperature(valve::Valve)
		Hin = PropsSI("H", "T", valve.Tin, "P", valve.Pin, valve.name)
		return uconvert(u"K", PropsSI("T", "H", Hin, "P", valve.Pex, valve.name))
	end
end

# ╔═╡ 43a8e290-19c9-4150-9c56-a1dad4f1338c
begin
	function pressure(hx::Cooler)
		k = 1.4
		polytropicExp = k/(k-1)
		P = hx.P*(hx.Tex/hx.Tin)^polytropicExp
		return uconvert(u"bar", P)
	end
end

# ╔═╡ 451a4e43-8a9e-46a3-8919-750648c8289e
let
	#Tₒ = 298u"K"
	#hx = Cooler(T_pyro, Tₒ, P_react, T_pyro, "H2")
	#uconvert(u"bar/bar", P_react/pressure(hx))
end

# ╔═╡ 554e0545-78ae-4fb3-9bbc-d00edcc91045
function molar_enthalpy(T, P, name::String, Hf=0u"kJ/kmol")
 	Href = PropsSI("Hmolar", "T", 298.15u"K", "P", 1u"atm", name)
 	H =  PropsSI("Hmolar", "T", T, "P", P, name)
	return uconvert(u"kJ/kmol", Hf + H - Href)
end

# ╔═╡ 46ba6b8a-c7e6-4609-a301-a05eaae599f0
begin
	function heat(react::Reactor) 
		Hin = molar_enthalpy(react.Tin, react.P, "CH4", Hf_meth)
		Hex = molar_enthalpy(react.Tex, react.P, "H2")
		Q =  n_dot_h2*Hex-MolarFlow["CH4"]*Hin
		return uconvert(u"kW", Q)
	end

	function heat(hx::Cooler)
		Hin = PropsSI("Hmolar", "T", hx.Tin, "P", hx.P, hx.name)
		Hex = PropsSI("Hmolar", "T", hx.Tex, "P", hx.P, hx.name)
		return uconvert(u"kJ/mol", (Hex-Hin))
	end
end

# ╔═╡ a5f8c399-a773-4868-b91b-3c8c7ac35103
function get_flow_rate(m_dot, T=273.15u"K", P=1u"atm")
	ρ = PropsSI("D", "T", T, "P", P, "CH4")
	return uconvert(u"ft^3/minute", m_dot/ρ)
end

# ╔═╡ fc17e89c-8495-42c7-902a-b5b005a4ed7e
function get_flow_exergy(T, P, Tₒ, Pₒ, name::String)
	Hₒ = PropsSI("Hmolar", "T", Tₒ, "P", Pₒ, name)
	Sₒ = PropsSI("Smolar", "T", Tₒ, "P", Pₒ, name)
	H = PropsSI("Hmolar", "T", T, "P", P, name)
	S = PropsSI("Smolar", "T", T, "P", P, name)
	flow = H - Hₒ - Tₒ*(S - Sₒ)
	return uconvert(u"kJ/kmol", +flow)
end

# ╔═╡ 2ff23a31-b092-4841-bf5c-d92076c61517
begin
	function exergyDestroyed(pump::Pump, Tₒ, Pₒ)
		T = temperature(pump)
		P = pump.ratio*pump.P
		e_in = get_flow_exergy(pump.T, pump.P, Tₒ, Pₒ, pump.name)
		e_ex = get_flow_exergy(T, P, Tₒ, Pₒ, pump.name)

		Hex = PropsSI("Hmolar", "T", T, "P", P, pump.name)
		Hin = PropsSI("Hmolar", "T", pump.T, "P", pump.P, pump.name) 
		return uconvert(u"kJ/kmol", Hex-Hin + (e_in - e_ex))
	end

	function exergyDestroyed(valve::Valve, Tₒ, Pₒ)
		Tex = temperature(valve)
		e_in = get_flow_exergy(valve.Tin, valve.Pin, Tₒ, Pₒ, valve.name)
		e_ex = get_flow_exergy(Tex, valve.Pex, Tₒ, Pₒ, valve.name)
		return uconvert(u"kJ/mol", (e_in - e_ex))
	end

	function exergyDestroyed(react::Reactor, Tₒ, Pₒ)
		Tb = (react.Tex + react.Tin)/2
		e_in = get_flow_exergy(react.Tin,react.P,Tₒ,Pₒ,"CH4")
		e_ex = get_flow_exergy(react.Tex,react.P,Tₒ,Pₒ,"H2")
		chem_Ex = MODEL2["CH4"] - 2*MODEL2["H2"] - MODEL2["C"]
		Q = -heat(react)
		Ed = (1-Tₒ/Tb)*Q + MolarFlow["CH4"]*e_in-MolarFlow["H2"]*e_ex + MolarFlow["CH4"]*chem_Ex
		return uconvert(u"kW", Ed)
	end

	function exergyDestroyed(hx::Cooler, Tₒ, Pₒ)
		e_in = get_flow_exergy(hx.Tin, hx.P, Tₒ, Pₒ, hx.name)
		e_ex = get_flow_exergy(hx.Tex, hx.P, Tₒ, Pₒ, hx.name)
		Tb = (hx.Tin + hx.Tex)/2
		Q = -heat(hx)
		Ex = (1-Tₒ/Tb)*Q + e_in - e_ex
		return uconvert(u"kJ/mol", Ex)
	end
end

# ╔═╡ 2041825c-155c-4ad5-9360-706ecf96ea46
ans5 = let
	Tₒ = 298u"K"
	Pₒ = 1u"atm"
	valve = Valve(T_in, P_in, P_react, "CH4")
	ED = exergyDestroyed(valve, Tₒ, Pₒ)*MolarFlow["CH4"]
	uconvert(u"kW", ED)
end

# ╔═╡ 378cd3e8-c855-4e7e-a996-31785cebf62d
ED_reactor = let
	Tₒ = 298u"K"
	Pₒ = 1u"atm"
	T_react = 298u"K"
	Tpump = 300u"K"
	Tb = 300u"K"
	#Tb = T_pyro
	rxn = Reactor(T_react, T_pyro, P_react, Tb)
	Ed_rxn = abs(exergyDestroyed(rxn, Tₒ, Pₒ))

	#Hydrogen Cooling
	hx = Cooler(T_pyro, Tpump, P_react, Tb, "H2")
	Ed_H2 = exergyDestroyed(hx, Tₒ, Pₒ)*MolarFlow["H2"]

	#Carbon Cooling
	Q = massC*C_carbon*(T_pyro - Tₒ)/day
	Ed_C = Q*(1-Tₒ/Tb)
	
	Ed_reactor = Ed_rxn + Ed_H2 + Ed_C
	uconvert(u"kW", Ed_reactor)
end

# ╔═╡ aaafbbca-052b-4993-931e-ac5c9a827809
# Final ratio (percentage)
(ED_C/ED_reactor)*100

# ╔═╡ cc1d1854-0bde-4102-8ecb-262e4b52ac35
ED_H2 = let
	Tₒ = 298u"K"
	Pₒ = 1u"atm"
	Tpump = 300u"K"
	Tb = 300u"K"
	#Tb = T_pyro
	#Hydrogen Cooling
	hx = Cooler(T_pyro, Tpump, P_react, Tb, "H2")
	Ed_H2 = exergyDestroyed(hx, Tₒ, Pₒ)*MolarFlow["H2"]

	uconvert(u"kW", Ed_H2)
end

# ╔═╡ 45e0cbfe-06a3-4504-8ebe-79b7915e7e22
ans8 = let
	uconvert(u"K/K", ED_H2/ED_reactor)*100
end

# ╔═╡ 0e0f87f3-c47e-493e-ba21-21f00f3df827
ans9 = let
	Tₒ = 298u"K"
	Pₒ = 1u"atm"
	Tb = 300u"K"
	Tp = 300u"K"
	P = P_react
	T = Tp		
	power = -5u"kW"	
	pumps = []
	hxs = []
	P = P_react
	T = Tp
	for i = 1:3
		pump = Pump(power/MolarFlow["H2"], 10, P, Tp, "H2", 0.65)
		push!(pumps, pump)
		P = pump.P * pump.ratio
		T = temperature(pump)
		push!(hxs, Cooler(T, Tp, P, Tb, "H2"))
	end
	finalPump = Pump(power/MolarFlow["H2"], 2, P, Tp, "H2", 0.65)
	push!(pumps, finalPump)
	T = temperature(finalPump)
	push!(hxs, Cooler(T, Tp, finalPump.ratio*P, Tb, "H2"))

	Ed_pumps = sum([exergyDestroyed(pump, Tₒ, Pₒ) for pump in pumps])
	Ed_hxs = sum([exergyDestroyed(hx, Tₒ, Pₒ) for hx in hxs])
	Ed_compression = MolarFlow["H2"] * (Ed_pumps + Ed_hxs)
	uconvert(u"kW", Ed_compression)
end

# ╔═╡ 21029ebb-c76a-46ec-a29a-a16b3037622c
function info(η) 
		Tₒ = 298u"K"
		Pₒ = 1u"atm"
		
		power = -5u"kW"	
		pumps = []
		P = P_react
		for i = 1:3
			pump = Pump(power/n_dot_h2, 10, P, 35u"°C", "H2", η)
			push!(pumps, pump)
			P = pump.P * pump.ratio
		end
		push!(pumps, Pump(power/n_dot_h2, 2, P, 35u"°C", "H2", η))
		
		for pump in pumps
			pump |> println
			@show temperature(pump)
			@show exergyDestroyed(pump, Tₒ, Pₒ)
			"\n" |> print
		end
end

# ╔═╡ 47f72e30-757e-4806-b189-2d79fdaa96b2
md"""

Resources: [here] https://www.hydrogen.energy.gov/pdfs/9013_energy_requirements_for_hydrogen_gas_compression.pdf

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

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

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
deps = ["Libdl", "libblastrampoline_jll"]
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

[[OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

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
deps = ["SHA", "Serialization"]
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

[[libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╠═3b5aaae0-e20e-11ec-0ea1-f39b52d3602e
# ╟─715ed4ad-d7c7-4698-9df1-61b92d8e39f2
# ╟─22a88a37-1952-4435-a5d6-ec3d6d1f5a13
# ╠═94fe6dee-89b4-4e6a-8af7-f35b049066d4
# ╠═3a036fda-e16f-4cbf-9bc5-372f575af9a9
# ╠═719eef9e-ac7a-4b2d-bb05-19740d43c9ba
# ╠═a7f4fe29-6485-438b-94d4-835ff2e97750
# ╠═e3a33cab-f8f3-4434-b087-6321698e975b
# ╠═addb44fe-d353-42f4-94ba-0a8a26bacdc6
# ╠═e5903da9-977b-4a30-a51a-2f84b9bb7a6d
# ╟─62f35fcc-0333-40ec-952d-f199031db50c
# ╠═66ffbd89-03a3-4931-ab2b-05f828a2973c
# ╟─7a13c45f-bc2e-452a-9c5c-2f455aaf5597
# ╠═ee854bbd-10bf-465c-8346-f6725aebdc5b
# ╟─7841815b-c3a5-4c5c-8970-cf18cfba676d
# ╟─54fb5386-ead8-4be0-b539-58c9c959f894
# ╠═5dd18145-5491-4137-9194-d4fe10b76816
# ╟─c5e700ce-7b4b-4b9a-987c-f1b98616efde
# ╠═172aa52c-1a62-435a-876f-91a6bf543cb5
# ╟─c8ac1365-8599-4ab6-b4be-8626fa63b225
# ╠═2041825c-155c-4ad5-9360-706ecf96ea46
# ╟─3a582d65-a2f9-4ccf-aab0-e6d38acbbd50
# ╠═378cd3e8-c855-4e7e-a996-31785cebf62d
# ╟─4a1b2c62-0048-4f58-a597-c0e728c9bdd0
# ╠═31e31fa0-bd10-4ee2-9222-be7283884dc7
# ╠═aaafbbca-052b-4993-931e-ac5c9a827809
# ╟─3f75aca3-7c5c-48e6-800d-2b6411171d49
# ╟─61691472-8282-4588-8ea1-daac04a472d2
# ╠═cc1d1854-0bde-4102-8ecb-262e4b52ac35
# ╠═45e0cbfe-06a3-4504-8ebe-79b7915e7e22
# ╟─1b6134f6-1643-4f90-831d-3a25c5c76253
# ╟─62634867-2e70-43bb-a0cd-b9629723f244
# ╠═0e0f87f3-c47e-493e-ba21-21f00f3df827
# ╟─0965a736-3af7-478d-b07c-e42d9b939473
# ╟─bb15cd28-6812-4402-bffe-a0c849528643
# ╟─0fd0f1e2-f6d4-4ed0-b88e-006742e4e40d
# ╠═09f63e0d-ed35-4a76-9cc0-f27b7a07f532
# ╟─2374b686-40a3-498f-b044-72248177966e
# ╠═be0ad564-d61b-41de-b607-00ac5731a7f1
# ╟─7720cdad-3c6b-4f4b-a4fe-691f6546b944
# ╟─029cc470-a3d6-439f-addf-79ed481a774f
# ╟─a5ef9be5-3c63-444c-9233-653d06fd0a68
# ╠═1d65754b-57ec-420d-a5be-adea52476969
# ╠═43a8e290-19c9-4150-9c56-a1dad4f1338c
# ╠═2ff23a31-b092-4841-bf5c-d92076c61517
# ╠═46ba6b8a-c7e6-4609-a301-a05eaae599f0
# ╠═451a4e43-8a9e-46a3-8919-750648c8289e
# ╠═554e0545-78ae-4fb3-9bbc-d00edcc91045
# ╠═a5f8c399-a773-4868-b91b-3c8c7ac35103
# ╠═fc17e89c-8495-42c7-902a-b5b005a4ed7e
# ╠═21029ebb-c76a-46ec-a29a-a16b3037622c
# ╠═47f72e30-757e-4806-b189-2d79fdaa96b2
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
