import numpy as np
from scipy.optimize import root
import math

print("&     ELEMENTARY VALUES    \n")
R_earth = 6371  # km
G = 6.67430e-11
M_earth = 5.9722e24
mu = 398600  # Standard Gravitational Parameter (G * M)

h = float(input("Enter Altitude (km): "))
#For LOE = 400 TO 2000 KM since that is for sample of satellites of observation and iss

print(f"Radius of Earth: {R_earth} km")
print(f"Altitude: {h} km")
print(f"Mass of Earth: {M_earth} kg")
print(f"Standard Gravitational Parameter (μ): {mu} km³/s²\n")

print("&     SHAPE AND SIZE PARAMETERS    \n")
rp = float(input("Perigee(km)"))
#For perigee = 6784 km of iss
ra = float(input("Apogee(km)"))
#For apogee = 6793 km of iss
e = (ra - rp)/(ra + rp)
a = (rp + ra)/2  # Semi-major axis
print(f"Semi-major Axis (a): {a} km")
print(f"Eccentricity (e): {e}\n")

print("&     ORIENTATION PARAMETERS    \n")
i = float(input("Inclination (deg): "))
i_r= i*0.0174533
omega = float(input("Argument of Periapsis (ω)(deg): "))
#For orbital inclination = 51.64° of iss
RAAN = float(input("RAAN (Ω) (deg): "))
RAAN_D = RAAN * 0.0174533
print(f"Inclination: {i}°")
print(f"RAAN: {RAAN}°")

print("&     POSITION AND TIMING PARAMETERS    \n")
#for loop
tp = float ( input("Time at which satellite passes periapses(in seconds)"))
t0 = float (input ( "Epoch / reference time(in seconds)"))
# Calculations
# Period T = 2 * pi * sqrt(a^3 / mu)
T = 2 * np.pi * np.sqrt(a**3 / mu)
# Mean Motion n = 2 * pi / T
n = (2 * np.pi) / T
Time_Interval = [0, 60 , 120, 180 , 240, 300, 360, 420 , 480 , 540 , 600 ]
print("                 ")
for t in Time_Interval:
    print(f"-----------TIME INTERVAL = {t} s-----------\n")
    # Mean Anomaly M ( with reference time epoch) = M0 + n (t - t0)
    M0 = n * (t - tp)
    MA = M0 + n * (t - t0)

    print(f"Orbital Period (T): {round(T, 2)} seconds (approx {round(T/3600, 2)} hours)\n")
    print(f"Mean Anomaly (M): {round(MA, 4)} radians\n")

    #Kepler's equation
    #MA = E - e sinE
    #We solve for E using MA and e
    #MA - E + esinE = 0
    def kepler_eq(E, e, M):
        return E - e*np.sin(E) - M

    E0 = MA 
    E_solution = root(kepler_eq, E0, args=(e, MA))
    E = E_solution.x[0]
    print(f"Eccentric Anomaly (E): {E:.4f} radians\n")


    print(f"Root Solver Result: {E_solution.success}") 
    #Function to caluculatw the orbital True Anamoly:
    def True_Anamoly(upsilon , E , e):
        return np.tan(upsilon)-(np.sqrt((1+e)/(1-e))*np.tan(E/2))

    u0=0
    TA_solution = root(True_Anamoly, u0 , args=(E,e))
    
    upsilon = TA_solution.x[0]
    print(f"True Anamoly (TA) : {upsilon:.4f} radians")
    
    TA_deg = np.degrees(upsilon) % 360
    print(f"True Anomaly: {TA_deg:.2f} degrees\n")
    
    print(f"Root solver Result: {TA_solution.success}")
    
    #Function to calculate Orbital radius:
    def Orbital_radius(R_orbit ,a, e , E):
        return R_orbit - a*(1-e*np.cos(E))
    
    r0 = 0
    R_O_solution =root(Orbital_radius , r0 , args = (a,e,E))
    
    R_orbit = R_O_solution.x[0]
    print(f"Orbital Radius: {R_orbit:.4f} m")
    
    print(f"Root solver Result: {R_O_solution.success}\n")

    #POSITION IN ORBITAL PLANE
    def Orbital_plane_x(xp , R_orbit , upsilon):
        return xp - R_orbit*np.cos(upsilon)
    
    px0 = 0
    Orbital_p_solution = root(Orbital_plane_x, px0 , args = (R_orbit , upsilon))
    
    xp = Orbital_p_solution.x[0]
    print(f"Xp coordinate : {xp:.1f} km ")
    
    print(f"Root solver Result: {Orbital_p_solution.success}\n")#?
    
    def Orbital_plane_x(yp , R_orbit , upsilon):
        return yp - R_orbit*np.sin(upsilon)
    
    px0 = 0
    Orbital_p_solution = root(Orbital_plane_x, px0 , args = (R_orbit , upsilon))

    yp = Orbital_p_solution.x[0]
    print(f"Yp coordinate : {yp:.1f} km ")
    
    print(f"Root solver Result: {Orbital_p_solution.success}\n")
    
    Zp = 0
    
    print(f"Zp coordinate: {Zp:.1f} km \n")
    
    #STEP 7
    cos_RAAN = math.cos(RAAN_D)
    sin_RAAN = math.sin(RAAN_D)
    cos_i=math.cos(i_r)
    sin_i=math.sin(i_r)
    cos_o = math.cos(omega)
    sin_o = math.sin(omega)
    sin_j = math.sin(upsilon)
    cos_j= math.cos(upsilon)
    R3_Raan = np.array([[cos_RAAN , -sin_RAAN , 0] ,
                        [sin_RAAN , cos_RAAN, 0],
                        [0, 0 ,1]])
    R3_omega = np.array([[cos_o , -sin_o , 0],
                        [sin_o , cos_o, 0],
                        [0, 0 ,1]])
    R1_inclination = np.array([[1,0,0] , 
                               [ 0 , cos_i , -sin_i] , 
                               [0 ,sin_i , cos_i]])
    rj = ([[R_orbit*cos_j] ,
           [R_orbit*sin_j], 
           [0]])
    R_ECI = R3_Raan@R1_inclination@R3_omega@rj
    
    #Expanded
    cos_ou = math.cos(upsilon)
    sin_ou = math.sin(upsilon)
    x = R_orbit*((cos_RAAN * cos_ou) - (sin_RAAN*sin_ou*cos_i))
    y = R_orbit*((sin_RAAN * cos_ou) - (cos_RAAN*sin_ou*cos_i))
    z = R_orbit*(sin_ou*sin_i)
    print(f"x = {x:.2f} km")
    print(f"y = {y:.2f} km")
    print(f"z = {z:.2f} km\n")
