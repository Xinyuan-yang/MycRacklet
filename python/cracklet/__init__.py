import py11_cracklet

private_keys = set(dir(py11_cracklet)) - set(dir())

for k in private_keys:
    globals()[k] = getattr(py11_cracklet,k)

def initialize(*args, **kwargs):
    raise RuntimeError("No need to call initialize")


def finalize(*args, **kwargs):
    raise RuntimeError("No need to call finalize")
