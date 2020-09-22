from fenics import *

def build_xdmf_file(fname):
    xdmf = XDMFFile(fname)
    xdmf.parameters["flush_output"] = True
    xdmf.parameters["functions_share_mesh"] = True
    xdmf.parameters["rewrite_function_mesh"] = False

    return xdmf

def write_soln_to_file(xdmf, solns, times, save_freq=1):
    for j in range(0, len(times), save_freq):
        xdmf.write(solns[j], times[j])

    return xdmf
