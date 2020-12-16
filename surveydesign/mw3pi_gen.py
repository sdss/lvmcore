""" Generate MW3Pi targets """

import sys
import yaml
import numpy as np
import copy
# Template
# GUM_ext4:
#     coords: [288.19, -6.84]
#     region_type: rectangle
#     frame: galactic
#     region_params:
#        width: 160.0
#        height: 0.10
#        pa: 90
#     priority: 9
#     observatory: BOTH
#     telescope: LVM-160
#     max_airmass: 1.75
#     exptime: 900
#     n_exposures: 1
#     min_exposures: 1
#     min_moon_dist: 60
#     max_lunation: 1.0
#     overhead: 1.1
#     tiling_strategy: center_first
#     group: ["GUM"]

def help():
    print("Usage: python3 mw3pi_gen.py")
    print("arguments are key:value pairs joined by ':'")
    print("Supported arguments:")
    print("Nx: int number of columns")
    print("Ny: int number of rows")
    print("coords: str 'gal' or 'radec': coordinate frame")
    print("xmin/max: float min/max of x")
    print("ymin/max: float min/max of y")

class read_grid_parameters(object):
    """
    docstring
    """
    pass
    def __init__(self):
        self.params = {
        "Nx": 100,
        "Ny": 0,
        "frame":'galactic',
        "xmin":-90,
        "xmax":90,
        "ymin":-30,
        "ymax":30}
    
        self.parse_args()

    def parse_args(self):
        def is_float(s):
            try:
                float(s)
                return True
            except ValueError:
                return False

        def is_int(s):
            try:
                float(s)
                return True
            except ValueError:
                return False

        if len(sys.argv) > 1:
            for arg_val in sys.argv[1:]:
                arg, val = arg.val.split(":")

                if arg in ["Nx","Ny"]:
                    try:
                        val = int(val)
                    except ValueError:
                        sys.exit("N[x|y] must be an integer")

                elif arg in ["xmin","xmax","ymin","ymax"]:
                    try:
                        val = float(val)
                    except ValueError:
                        sys.exit("[x|y][min|max] must be a float")

                self.params[arg] = val

def make_grid():
    grid = read_grid_parameters()
    targets = {}

    template = {
        "coords": [0, 0],
        "region_type": "rectangle",
        "frame": "galactic",
        "region_params":{
            "width": 0.0,
            "height": 0.10,
            "pa": 90},
        "priority": 1,
        "observatory": "BOTH",
        "telescope": "LVM-160",
        "max_airmass": 1.75,
        "min_shadowheight": 1000.0,
        "exptime": 900,
        "n_exposures": 1,
        "min_exposures": 1,
        "min_moon_dist": 60,
        "max_lunation": 1.0,
        "overhead": 1.1,
        "tiling_strategy": "center_first",
        "group": ["MW3PI"]
    }

    # Target number for both row and column
    target_i = 0

    if grid.params["Nx"] > 0:
        dx = (grid.params["xmax"] - grid.params["xmin"])/(grid.params["Nx"]+1)
        # calculation of region size is done outside of loop because there is no dec dependence.
        l = abs(grid.params["ymax"] - grid.params["ymin"])
        y = (grid.params["ymax"] - grid.params["ymin"])/2.0

        for n in range(grid.params["Nx"]):
            name = "MW3PI_%i"%target_i
            targets[name] = copy.deepcopy(template)
            targets[name]["frame"] = grid.params["frame"]

            # Add +1 when calculating dx to prevent overlap when going from 0-360
            x = grid.params["xmin"] + dx * n
            targets[name]["coords"] = [x,y]

            targets[name]["region_params"] = {"width":l, "height":0.10, "pa":90}
        
            target_i += 1

    if grid.params["Ny"] > 0:
        dy = (grid.params["ymax"] - grid.params["ymin"])/(grid.params["Ny"]+1)
        x = (grid.params["xmax"] - grid.params["xmin"])/2.0

        for n in range(grid.params["Nx"]):
            name = "MW3PI_%i"%target_i
            targets[name] = copy.deepcopy(template)
            targets[name]["frame"] = grid.params["frame"]

            # Add +1 when calculating dx to prevent overlap when going from 0-360
            y = grid.params["ymin"] + dy * n
            targets[name]["coords"] = [x,y]

            # calculation of region size is done inside of loop because there is dec dependence.
            l = abs(grid.params["xmax"] - grid.params["xmin"])*np.cos(y * np.pi/180) * 0.9
            targets[name]["region_params"] = {"width":l, "height":0.10, "pa":0.0}

            target_i += 1

    with open('MW3Pi_lvm_target.yaml', 'w') as outfile:
        yaml.dump(targets, outfile,default_flow_style=None)

if __name__ == "__main__":
    make_grid()


