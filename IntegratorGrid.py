def int1(sysODE, sysInit, untilTime, ODEpars):

    import heyoka as hy
    from numpy import linspace
    ta = hy.taylor_adaptive(sysODE,sysInit, pars = ODEpars)

    grid = linspace(0, untilTime, 1000)
    #intOut = ta.propagate_until(t = untilTime)
    output = ta.propagate_grid(grid)

    return(output[4])