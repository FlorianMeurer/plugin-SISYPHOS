import os

elements = ["Au"]

dat_names = ["mu", 
             "wavelength", 
             "F000", 
             "tot_reflIns", 
             "goof", 
             "R_all", 
             "R1", 
             "wR2", 
             "last Shift"]

corr_filts = ["exptl_absorpt_coefficient_mu", 
              "diffrn_radiation_wavelength", 
              "exptl_crystal_F_000",
              "diffrn_reflns_number",
              "refine_ls_goodness_of_fit_ref",
              "refine_ls_R_factor_all",
              "refine_ls_R_factor_gt",
              "refine_ls_wR_factor_ref",
              "REM Shift_max"]

def parse_cif(loc):
    out = {}
    disp_dict = {}

    with open(loc, "r") as incif:
        for line in incif:
            for i,filter in enumerate(corr_filts):
                if filter in line:
                    out[f"{dat_names[i]}"] = float(line.split()[-1])

    with open(loc, "r") as incif:
        switch2 = False
        for line in incif:
              for elem in elements:
                switch = True                
                if line.startswith(f" {elem}  ") & switch:
                    print(line)
                    if " . . " in line:
                      continue

    with open(loc, "r") as incif:
        switch3 = False
        for line in incif:
              if switch3 and line in ["\n", "\r\n"]:
                  switch3 = False
              if switch3:
                  adr = line.split()
                  print(adr)
                  disp_dict[adr[0]] = (adr[1],adr[2])
              if line.startswith("  _atom_site_dispersion_imag"):
                  switch3 = True
    
    return out, disp_dict

data = ["MoC6O6_20001.cif", "MoC6O6_20100.cif"]

for dat in data:
    df1,df2 = parse_cif(dat)
    print(df1,df2)            
