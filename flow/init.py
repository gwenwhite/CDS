#!/usr/bin/env python
"""Initialize the project's data space.

Iterates over all defined state points and initializes
the associated job workspace directories.
The result of running this file is the creation of a signac workspace:
    - signac.rc file containing the project name
    - signac_statepoints.json summary for the entire workspace
    - workspace/ directory that contains a sub-directory of every individual statepoint
    - signac_statepoints.json within each individual statepoint sub-directory.

"""

import signac
import flow
import logging
from collections import OrderedDict
from itertools import product


def get_parameters():
    ''''''
    parameters = OrderedDict()
    # System and model params:
    # MN is the combo of number of chains (M) and chain lengths (N)
    parameters["x"] = [15]
    parameters["y"] = [15]
    parameters["z"] = [15]
    parameters["php"] = [1/3]
    parameters["kT"] = [1.0]
    parameters["n_steps"] = [5e6]
    parameters["n_pillars"] = [2]
    parameters["gsd_period"] = [1e4]
    parameters["r_cut"] = [9.037500000000001]
    parameters["dt"] = [0.0005]
    parameters["tau"] = [100]
    parameters["sim_seed"] = [42]
    return list(parameters.keys()), list(product(*parameters.values()))


def main():
    project = signac.init_project() # Set the signac project name
    param_names, param_combinations = get_parameters()
    # Create the generate jobs
    for params in param_combinations:
        statepoint = dict(zip(param_names, params))
        job = project.open_job(statepoint)
        job.init()
        job.doc.setdefault("equilibrated", False)
        job.doc.setdefault("runs", 0)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
