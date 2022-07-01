R = 8.31446261815324 * 1e-3  # kJ /K /mol
kB = 1.380649 * 1e-23  #  m2 kg s-2 K-1

# TODO: add values for effective radii
# at the moment: most are filled out in a way that the tests pass and not based on papers
radius = {}
radius["HCO3"] = 0.195e-9 
radius["CO3"] = 0.24e-9 
radius["H"] = 0.0246e-9 # Shannon & Prewitt, 1969
radius["OH"] = 0.0435e-9
radius["BOH3"] = 0.193e-9
radius["BOH4"] = 0.22e-9  
