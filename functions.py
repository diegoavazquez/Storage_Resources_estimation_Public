def Temp(z,thermal=0.03, Ts=5):
    '''Calculates Temperature
    Inputs: z[m], Thermal gradient in C/m (default 0.03), Ts Temp @ surface in C (default 5)'''
    T = z * thermal + Ts
    return T
    
    
def Press_ini(z, depl, hydro= 0.112):
    '''Calculates Initial Pressure Pini. 
    Inputs: z[m], depl[bar], hydrostatic gradient in bar/m (default 1.12/10)'''
    Pini = z * hydro - depl
    return Pini

  
def compress_r(por, alpha=2.7122106e-05, beta=-0.415):
    '''Calculates Rock compressibility in bar-1 - ref. Hall (1953)
    Inputs: por, alpha (default 2.7122057006e-05), beta (default -0.415)'''
    comp_r = alpha*(por**beta)
    return comp_r

  
def compress_w(z, press, thermal=0.03, Ts=5):
    '''Calculates Water Compressibility (ref. Brill and Beggs, 1978), which depends on Pini. 
    Inputs: z[m], Thermal gradient in C/m (default 0.03), Ts Temp @ surface in C (default 5),
    hydrostatic gradient in bar/m (default 1.12/10)'''
    #por = np.array(por)
    T = z * thermal + Ts
    c1 = 3.8546 - 0.000134 * press *14.5038
    c2 = -0.01052 + 0.000000477 * press * 14.5038
    c3 = 0.000039267 - 0.00000000088 * press * 14.5038
    comp_w = (c1 + c2 * (T * 9/5+ 32) + c3* (T* 9/5+ 32)**2) *0.000001 * 14.5038
    return comp_w

  
def dyn_Pmax(z, crest, SF=0.15, hydro= 0.112, Shmin= None):
    '''Depth dependant relationship to estimate Maximum Allowable Pressure in the "tank" 
    Lower Bound Shmin - SF % safety factor at crest + hydrostatic head
    Inputs: z [m], crest[m], SF (default 0.15), hydrostatic gradient in bar/m (default 1.12/10), Shmin in bar/m. 
    If Shmin not provided, it is calculated according to SNS depth relationship'''
    if Shmin is None:
        Pmax = (1-SF)* (0.76*0.00000000263*crest**3 + 0.15*0.00000802 * crest**2 + 1.15*0.12*crest) + (z-crest)*hydro
    else:
        Pmax = (1-SF)* Shmin * z + (z-crest)*hydro
    return Pmax
  

def den_CO2(P,T):
    '''Calculation of CO2 density in Kg/m3 = gr/l depending on Pressure in bar and Temperature in C.
    CO2 behaves as a supercritical fluid above its critical temperature (304.13 K, 31.0 °C, 87.8 °F)
    and critical pressure (7.3773 MPa, 72.8 atm, 1,070 psi, 73.8 bar), expanding to fill its container 
    like a gas but with a density like that of a liquid. 
    This fucntion calculates the RHOCO2 dependant on P and T based on Ouyang (2011)
    https://benthamopen.com/contents/pdf/TOPEJ/TOPEJ-4-13.pdf
    Inputs: Pressure in psia, and the correlation coefficients A0, A1 – A4 are solely associated 
    with Temperature in Celsius.'''
    P = P * 14.5038 # transform P from bar to psi
    if P < 3000:
        b00 = -2.148322085348E+05
        b01 = 1.168116599408E+04
        b02 = -2.302236659392E+02
        b03 = 1.967428940167E+00
        b04 = -6.184842764145E-03 
        b10 = 4.757146002428E+02
        b11 = -2.619250287624E+01
        b12 = 5.215134206837E-01
        b13 = -4.494511089838E-03
        b14 = 1.423058795982E-05 
        b20 = -3.713900186613E-01
        b21 = 2.072488876536E-02
        b22 = -4.169082831078E-04
        b23 = 3.622975674137E-06
        b24 = -1.155050860329E-08
        b30 = 1.228907393482E-04
        b31 = -6.930063746226E-06
        b32 = 1.406317206628E-07
        b33 = -1.230995287169E-09
        b34 = 3.948417428040E-12 
        b40 = -1.466408011784E-08
        b41 = 8.338008651366E-10
        b42 = -1.704242447194E-11
        b43 = 1.500878861807E-13
        b44 = -4.838826574173E-16 
    elif P >= 3000:
        b00 = 6.897382693936E+02
        b01 = 2.730479206931E+00
        b02 = -2.254102364542E-02
        b03 =  -4.651196146917E-03
        b04 = 3.439702234956E-05 
        b10 = 2.213692462613E-01
        b11 = -6.547268255814E-03
        b12 = 5.982258882656E-05
        b13 = 2.274997412526E-06
        b14 = -1.888361337660E-08 
        b20 = -5.118724890479E-05
        b21 = 2.019697017603E-06
        b22 = -2.311332097185E-08
        b23 = -4.079557404679E-10
        b24 = 3.893599641874E-12
        b30 = 5.517971126745E-09
        b31 = -2.415814703211E-10
        b32 = 3.121603486524E-12
        b33 = 3.171271084870E-14
        b34 = -3.560785550401E-16       
        b40 = -2.184152941323E-13
        b41 = 1.010703706059E-14
        b42 = -1.406620681883E-16
        b43 = -8.957731136447E-19
        b44 = 1.215810469539E-20 
    
    A0 = b00 + b01*T + b02*T**2  + b03*T**3 + b04*T**4 
    A1 = b10 + b11*T + b12*T**2  + b13*T**3 + b14*T**4 
    A2 = b20 + b21*T + b22*T**2  + b23*T**3 + b24*T**4 
    A3 = b30 + b31*T + b32*T**2  + b33*T**3 + b34*T**4 
    A4 = b40 + b41*T + b42*T**2  + b43*T**3 + b44*T**4 
        
    rho = A0 + A1*P + A2*P**2 + A3*P**3 + A4*P**4
    
    if P < 1070:
        rho = P * 44 / (0.0821 * (T + 273.15))  
        '''Calculation of CO2 density in Kg/m3 = gr/l depending on Pressure in bar and Temperature in C
    The molar mass of CO2 = atomic mass of C + 2×atomic mass of O =12+2×16 = 44g/mol
    The gas constant is R = 0.0821 Latmmol−1K−1 '''
    return rho
