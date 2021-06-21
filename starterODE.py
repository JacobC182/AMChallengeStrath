def ODE():

    import heyoka as hk
    import numpy as np


    GMe = 3.986004407799724e+5 # [km^3/sec^2]
    GMo = 1.32712440018e+11 #[km^3/sec^2]
    GMm = 4.9028e+3 #[km^3/sec^2]
    Re = 6378.1363 #[km]
    C20 = -4.84165371736e-4
    C22 = 2.43914352398e-6
    S22 = -1.40016683654e-6
    theta_g = (np.pi/180)*280.4606 #[rad]
    nu_e = (np.pi/180)*(4.178074622024230e-3) #[rad/sec]
    nu_o = (np.pi/180)*(1.1407410259335311e-5) #[rad/sec]
    nu_ma = (np.pi/180)*(1.512151961904581e-4) #[rad/sec]
    nu_mp = (np.pi/180)*(1.2893925235125941e-6) #[rad/sec]
    nu_ms = (np.pi/180)*(6.128913003523574e-7) #[rad/sec]
    alpha_o = 1.49619e+8 #[km]
    epsilon = (np.pi/180)*23.4392911 #[rad]
    phi_o = (np.pi/180)*357.5256 #[rad]
    Omega_plus_w = (np.pi/180)*282.94 #[rad]
    PSRP = 4.56e-3 #[kg/(km*sec^2)]

    X,Y,Z = hk.make_vars("X","Y","Z")
    VX,VY,VZ = hk.make_vars("VX","VY","VZ")

    #Sun's position
    lo = phi_o + nu_o*hk.time
    lambda_o = Omega_plus_w + lo + (np.pi/180)*( (6892/3600)*hk.sin(lo) + (72/3600)*hk.sin(2*lo) )
    ro = (149.619 - 2.499*hk.cos(lo) - 0.021*hk.cos(2*lo))*(10**6)

    Xo = ro*hk.cos(lambda_o)
    Yo = ro*hk.sin(lambda_o)*np.cos(epsilon)
    Zo = ro*hk.sin(lambda_o)*np.sin(epsilon)

    #Moon's position
    phi_m = nu_o*hk.time
    phi_ma = nu_ma*hk.time
    phi_mp = nu_mp*hk.time
    phi_ms = nu_ms*hk.time
    L0 = phi_mp + phi_ma + (np.pi/180)*218.31617
    lm = phi_ma + (np.pi/180)*134.96292
    llm = phi_m + (np.pi/180)*357.5256
    Fm = phi_mp + phi_ma + phi_ms + (np.pi/180)*93.27283
    Dm = phi_mp + phi_ma - phi_m  + (np.pi/180)*297.85027

    rm = 385000 - 20905*hk.cos(lm) - 3699*hk.cos(2*Dm - lm) - 2956*hk.cos(2*Dm) - \
        570*hk.cos(2*lm) + 246*hk.cos(2*lm - 2*Dm) - 205*hk.cos(llm - 2*Dm) - \
        171*hk.cos(lm + 2*Dm) - 152*hk.cos(lm + llm - 2*Dm)
        
    lambda_m = L0 + (np.pi/180)*( (22640/3600)*hk.sin(lm) + (769/3600)*hk.sin(2*lm) - (4856/3600)*hk.sin(lm - 2*Dm) + \
        (2370/3600)*hk.sin(2*Dm) - (668/3600)*hk.sin(llm) - (412/3600)*hk.sin(2*Fm) - \
        (212/3600)*hk.sin(2*lm - 2*Dm) - (206/3600)*hk.sin(lm + llm - 2*Dm) + \
        (192/3600)*hk.sin(lm + 2*Dm) - (165/3600)*hk.sin(llm - 2*Dm) + \
        (148/3600)*hk.sin(lm - llm) - (125/3600)*hk.sin(Dm) - (110/3600)*hk.sin(lm + llm) - \
        (55/3600)*hk.sin(2*Fm - 2*Dm) )
        
    βm = (np.pi/180)*( (18520/3600)*hk.sin(Fm + lambda_m - L0 + (np.pi/180)*((412/3600)*hk.sin(2*Fm) + (541/3600)*hk.sin(llm)) ) - \
        (526/3600)*hk.sin(Fm - 2*Dm) + (44/3600)*hk.sin(lm + Fm - 2*Dm) - (31/3600)*hk.sin(-lm + Fm -2*Dm) - \
        (25/3600)*hk.sin(-2*lm + Fm) - (23/3600)*hk.sin(llm + Fm - 2*Dm) + (21/3600)*hk.sin(-lm + Fm) + \
        (11/3600)*hk.sin(-llm + Fm - 2*Dm) )
        
    Xm =  hk.cos(βm)*hk.cos(lambda_m)*rm
    Ym = -np.sin(epsilon)*hk.sin(βm)*rm + np.cos(epsilon)*hk.cos(βm)*hk.sin(lambda_m)*rm
    Zm =  np.cos(epsilon)*hk.sin(βm)*rm + hk.cos(βm)*np.sin(epsilon)*hk.sin(lambda_m)*rm

    #Earth's Keplerian terms
    magR2 = X**2 + Y**2 + Z**2
    fKepX = -GMe*X/(magR2**(3./2))
    fKepY = -GMe*Y/(magR2**(3./2))
    fKepZ = -GMe*Z/(magR2**(3./2))

    #Earth's J2 terms
    J2term1 = GMe*(Re**2)*np.sqrt(5)*C20/(2*magR2**(1./2))
    J2term2 = 3/(magR2**2)
    J2term3 = 15*(Z**2)/(magR2**3)
    fJ2X = J2term1*X*(J2term2 - J2term3)
    fJ2Y = J2term1*Y*(J2term2 - J2term3)
    fJ2Z = J2term1*Z*(3*J2term2 - J2term3)

    #Earth's C22 and S22 terms
    x =  X*hk.cos(theta_g + nu_e*hk.time) + Y*hk.sin(theta_g + nu_e*hk.time)
    y = -X*hk.sin(theta_g + nu_e*hk.time) + Y*hk.cos(theta_g + nu_e*hk.time)
    z = Z
    magr2 = x**2 + y**2 + z**2

    C22term1 = 5*GMe*(Re**2)*np.sqrt(15)*C22/(2*magr2**(7./2))
    C22term2 = GMe*(Re**2)*np.sqrt(15)*C22/(magr2**(5./2))
    fC22x = C22term1*x*(y**2 - x**2) + C22term2*x
    fC22y = C22term1*y*(y**2 - x**2) - C22term2*y
    fC22z = C22term1*z*(y**2 - x**2)

    S22term1 = 5*GMe*(Re**2)*np.sqrt(15)*S22/(magr2**(7./2))
    S22term2 = GMe*(Re**2)*np.sqrt(15)*S22/(magr2**(5./2))
    fS22x = -S22term1*(x**2)*y + S22term2*y
    fS22y = -S22term1*x*(y**2) + S22term2*x
    fS22z = -S22term1*x*y*z

    fC22X = fC22x*hk.cos(theta_g + nu_e*hk.time) - fC22y*hk.sin(theta_g + nu_e*hk.time)
    fC22Y = fC22x*hk.sin(theta_g + nu_e*hk.time) + fC22y*hk.cos(theta_g + nu_e*hk.time)
    fC22Z = fC22z

    fS22X = fS22x*hk.cos(theta_g + nu_e*hk.time) - fS22y*hk.sin(theta_g + nu_e*hk.time)
    fS22Y = fS22x*hk.sin(theta_g + nu_e*hk.time) + fS22y*hk.cos(theta_g + nu_e*hk.time)
    fS22Z = fS22z

    #Sun's gravity
    magRo2 = Xo**2 + Yo**2 + Zo**2
    magRRo2 = (X - Xo)**2 + (Y - Yo)**2 + (Z - Zo)**2
    fSunX = -GMo*( (X - Xo)/(magRRo2**(3./2)) + Xo/(magRo2**(3./2)) )
    fSunY = -GMo*( (Y - Yo)/(magRRo2**(3./2)) + Yo/(magRo2**(3./2)) )
    fSunZ = -GMo*( (Z - Zo)/(magRRo2**(3./2)) + Zo/(magRo2**(3./2)) )

    #Moon's gravity 
    magRm2 = Xm**2 + Ym**2 + Zm**2
    magRRm2 = (X - Xm)**2 + (Y - Ym)**2 + (Z - Zm)**2
    fMoonX = -GMm*( (X - Xm)/(magRRm2**(3./2)) + Xm/(magRm2**(3./2)) )
    fMoonY = -GMm*( (Y - Ym)/(magRRm2**(3./2)) + Ym/(magRm2**(3./2)) )
    fMoonZ = -GMm*( (Z - Zm)/(magRRm2**(3./2)) + Zm/(magRm2**(3./2)) )

    #Sun's radiation pressure (AOM is a heyoka parameter hy.par[0]. We
    #will be able to set it later without recompiling the integartor)
    SRPterm = hk.par[0]*PSRP*(alpha_o**2)/(magRRo2**(3./2))
    fSRPX = SRPterm*(X - Xo)
    fSRPY = SRPterm*(Y - Yo)
    fSRPZ = SRPterm*(Z - Zo)

    dXdt = VX
    dYdt = VY
    dZdt = VZ
    dVXdt = fKepX + fJ2X + fC22X + fS22X + fSunX + fMoonX + fSRPX
    dVYdt = fKepY + fJ2Y + fC22Y + fS22Y + fSunY + fMoonY + fSRPY
    dVZdt = fKepZ + fJ2Z + fC22Z + fS22Z + fSunZ + fMoonZ + fSRPZ

    systemODE = [(X,dXdt),(Y,dYdt),(Z,dZdt),(VX,dVXdt),(VY,dVYdt),(VZ,dVZdt)]

    return (systemODE)