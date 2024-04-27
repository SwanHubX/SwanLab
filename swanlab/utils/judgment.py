def in_jupyter():
    try:
        from IPython import get_ipython

        if "IPKernelApp" not in get_ipython().config:  # Kernel does not appear to be running
            return False
    except:
        return False
    return True
