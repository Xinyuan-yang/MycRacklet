import numpy as _aka_np
import py11_cracklet

private_keys = set(dir(py11_cracklet)) - set(dir())

for k in private_keys:
    globals()[k] = getattr(py11_cracklet,k)

def initialize(*args, **kwargs):
    raise RuntimeError("No need to call initialize,"
                       " use parseInput to read an input file")


def finalize(*args, **kwargs):
    raise RuntimeError("No need to call finalize")
