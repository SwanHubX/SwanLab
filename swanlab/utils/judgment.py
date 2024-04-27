def in_jupyter():
    try:
        _ = __IPYTHON__
        return True
    except NameError:
        return False
