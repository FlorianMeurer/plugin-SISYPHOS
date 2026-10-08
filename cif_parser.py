import os

elements = ["Mo"]
disp = True

def parse_cif( loc: str) -> dict:
      """Parses the cif given by loc and returns a dictionary of parsed information

      Args:
          loc (str): Path to the .cif file to be analyzed

      Returns:
          dict: Result dictionary from cif
      """
      print("loc", loc)
      dat_names = ["mu", 
        "wavelength", 
        "F000", 
        "tot_reflIns", 
        "goof", 
        "R_all", 
        "R1", 
        "wR2", 
        "last Shift"]

      corr_filts = ["_diffrn_radiation_wavelength",  
                    "_exptl_crystal_F_000",
                    "_diffrn_reflns_number",
                    "_refine_ls_goodness_of_fit_ref",
                    "_refine_ls_R_factor_all",
                    "_refine_ls_R_factor_gt",
                    "_refine_ls_wR_factor_ref",
                    "REM Shift_max"]
      out_dict = {}
      disp_dict = {}
      try:
        with open(loc, "r") as incif:
            lines = incif.readlines()
            for idx, line in enumerate(lines):
                for i, filter_str in enumerate(corr_filts):
                    if filter_str in line:
                        parts = line.split()
                        if len(parts) > 1:
                            try:
                                value = float(parts[-1])
                                out_dict[f"{dat_names[i]}"] = value
                            except (ValueError, IndexError) as e:
                                print(f"Failed to parse value from {filter_str}: line='{line.strip()}', error={e}")
                        else:
                            if idx + 1 < len(lines):
                                try:
                                    value = float(lines[idx + 1].strip())
                                    out_dict[f"{dat_names[i]}"] = value
                                except (ValueError, IndexError) as e:
                                    print(f"Failed to parse value from next line for {filter_str}: line='{lines[idx + 1].strip()}', error={e}")
      except Exception as e:
        print("fail 1", e)
      try:
        with open(loc, "r") as incif:
          print("outdict",out_dict)
          switch2 = False
          for line in incif:
            if line.startswith("  _atom_site_refinement_flags_occupancy"):
              switch2 = True
              continue
            if switch2:
              if line.startswith("\n"):
                switch2 = False
              else:
                lin = line.split(" ")
                atom = lin[1]
                ueq = lin[6].split("(")[0]
                ueq_delta = lin[6].split("(")[1][:-1]
                out_dict[f"{atom}_ueq"] = (float(ueq), int(ueq_delta))
        with open(loc, "r") as incif:
          if disp:
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
      except Exception as e:
        print("fail 2", e)
      return out_dict, disp_dict
  
  
  
data = ["MoC6O6_20001.cif", "MoC6O6_20100.cif"]

for dat in data:
    df1,df2 = parse_cif(dat)
    print(df1,df2)            
