import matplotlib.pyplot as plt
import scipy as sp
from scipy.integrate import quad
import numpy as np
import os
from scipy.interpolate import interp1d
import math

file_path= os.path.join(os.getcwd(), 'raw_testG9.txt')
map_path= os.path.join(os.getcwd(), 'mapping.txt') #r'C:\Users\frede\Downloads\raw_testG9.txt'
path_upper = os.path.join(os.getcwd(), 'Airfoil Upper.txt') #r'C:\Users\sande\OneDrive\Bureaublad\Airfoil Upper.txt'
path_lower = os.path.join(os.getcwd(), 'Airfoil Lower.txt') #r'C:\Users\sande\OneDrive\Bureaublad\Airfoil Lower.txt'

points_upper = np.flip(np.genfromtxt(path_upper, delimiter = ''),0)
points_lower = np.genfromtxt(path_lower, delimiter = '')
slopes_upper = np.empty(len(points_upper)-1)
slopes_lower = np.empty(len(points_lower)-1)

for i in range(0, len(slopes_upper)):
    dy = points_upper[i+1][1]-points_upper[i][1]
    dx = points_upper[i+1][0]-points_upper[i][0]
    slope = dy/dx
    slopes_upper[i] = slope
    
for i in range(0, len(slopes_lower)):
    dy = points_lower[i+1][1]-points_lower[i][1]
    dx = points_lower[i+1][0]-points_lower[i][0]
    slope = dy/dx
    slopes_lower[i] = slope
    
def axialComponent(pressure, pressX, isUpper):
    slopes = 0
    slopesX = 0
    liftPress = 0
    dragPress = 0
    slope = None
    if isUpper:
        slopes = slopes_upper
        slopesX = points_upper[:,0]
    else:
        slopes = slopes_lower
        slopesX = points_lower[:,0]
    
    
    i = 0
    while slope == None: # If not already solved by edge cases
        if i < len(slopes)-1:
            if pressX < slopesX[i+1]:
                
                slope = slopes[i]
        else:
            slope = slopes[-1]
        i +=1
    Axial_pressure = -slope*pressure
    
    return(Axial_pressure)
List=np.genfromtxt(file_path,delimiter='')
v=22
chord = 0.16
with open(map_path, 'r') as file:
    locs = [float(line.strip()) for line in file]
airfoilUpperPerc = np.array(locs[:25])/100*chord
airfoilLowerPerc = np.array(locs[25:49])/100*chord

totalWakemm = np.array(locs[49:96])/1000
staticWakemm = np.array(locs[96:109])/1000
    
def ambient(Line,List):
    #rho=List[Line, 7]
    pt = List[Line, 104]
    p_bar = List[Line, 4]
    T = List[Line, 5] +273.15
    pb= List[Line, 3]

    #rho = 28.96 * p_bar /(T* 8.314)
    rho = p_bar*100/(287*T)

    S = 110.4
    T0 = 273.15
    viscosity = 1.716* 10**(-5) * (T/T0)**(3/2) * ((T0+S)/(T+S))

    q_inf = 0.211804 + 1.928442 * pb + 1.879374*10**(-4) * pb**2
    ps_inf = p_bar*100 - q_inf#(pt-q_inf)/2p_bar*100 - q_inf

    U_inf = np.sqrt(2*q_inf/rho)
    Rey = rho*U_inf*chord /viscosity

    return rho,viscosity,U_inf,Rey,ps_inf,p_bar



def pressure(Line,p_inf):
    pos = []
    p_top = []
    p_bottom = []
    p_backtotal = []
    p_backstatic = []
    pos2 = []
    x = 33
    b = 0
    for i in range(8,x):
        p_top.append(List[Line,i]+p_inf)
        pos.append(b)
        b+=chord/(25-1)#0.667
    b=0
    for i in range(x,57):
        pos2.append(b)
        p_bottom.append(List[Line,i]+p_inf)
        b+=chord/(24-1)
    
    backtotal = np.linspace(0,1,len(range(57,104)))
    backstatic = np.linspace(0,1,len(range(105,117)))
    for i in range(57,104):
        p_backtotal.append(List[Line,i]+p_inf)
    for i in range(105, 117):
        p_backstatic.append(List[Line,i]+p_inf)

    return p_top,pos,p_bottom,pos2,p_backtotal,backtotal,p_backstatic,backstatic



def cn(upper, lower, upperPos, lowerPos):
    lower_interp = interp1d(lowerPos, lower, bounds_error=False, fill_value="extrapolate")
    upper_interp = interp1d(upperPos, upper, bounds_error=False, fill_value="extrapolate")
    cn, _ = quad(lambda x: (lower_interp(x)-upper_interp(x))/chord, 0, chord)
    return cn

def cm(upper, lower, upperPos, lowerPos):
    lower_interp = interp1d(lowerPos, lower, bounds_error=False, fill_value="extrapolate")
    upper_interp = interp1d(upperPos, upper, bounds_error=False, fill_value="extrapolate")
    cm, _ = quad(lambda x: (upper_interp(x)-lower_interp(x))*x/chord**2, 0, chord)
    return cm

def dragRake(Uinf, Pinf,rho, Uwake, Pwake, Upos, pPos):
    Uw_interpd = interp1d(Upos, Uwake, bounds_error=False, fill_value="extrapolate")
    Pw_interpd = interp1d(pPos, Pwake, bounds_error=False, fill_value="extrapolate")
    firstTerm, _ = quad(lambda y: (Uinf-Uw_interpd(y))*Uw_interpd(y),Upos[0], Upos[-1])
    thirdTerm, _ = quad(lambda y: Pinf-Pw_interpd(y),Upos[0], Upos[-1])
    cd = (rho*firstTerm+thirdTerm)/(0.5*rho*U_inf**2*chord)
    return cd

if __name__ == '__main__':
    Line= 21
    single = False
    if single == False:
        AOA = []
        CL_Wake_List = []
        CL_List = []
        CD_Wake_List = []
        CD_List = []
        CM_List = []
        for j in range(3,22): #22 for no hysteresis, 35 for all data

            rho, viscous, U_inf, Rey, p_inf, pATM= ambient(j,List)
            p_top, pos, p_bottom, pos2, p_backtotal, backtotal, p_backstatic, backstatic = pressure(j,p_inf)


            p_axial_top = np.empty_like(p_top)
            p_axial_bottom = np.empty_like(p_bottom)

            for i in range(len(pos)):
                p_top[i]= (p_top[i]-p_inf)/(0.5 * rho *U_inf**2)
                p_axial_top[i] = axialComponent(p_top[i], airfoilUpperPerc[i]/.16, True)
            for i in range(len(pos2)):
                p_bottom[i]= (p_bottom[i]-p_inf)/(0.5 * rho *U_inf**2)
                p_axial_bottom[i] = axialComponent(p_bottom[i], airfoilLowerPerc[i]/.16, False)
            for i in range(len(backtotal)):
                cpt = (p_backtotal[i]-p_inf)/(0.5*rho*U_inf**2)
                if cpt > 1: cpt = 1
                p_backtotal[i]=cpt*0.5*rho*U_inf**2+p_inf
                
            

            backstatic_interpd = interp1d(staticWakemm, p_backstatic, bounds_error=False, fill_value="extrapolate")
            velocity_wake = []
            velocity_pos = []
            for i,x in enumerate(totalWakemm):
                if x < staticWakemm[0] or x > staticWakemm[-1]: continue
                vw = (2*(p_backtotal[i]-backstatic_interpd(x))/rho)**0.5
                velocity_wake.append(vw)
                velocity_pos.append(x)

            CN = cn(p_top, p_bottom, airfoilUpperPerc, airfoilLowerPerc)
            CA = cn(p_axial_top, p_axial_bottom, airfoilUpperPerc, airfoilLowerPerc)
            CM = cm(p_top, p_bottom, airfoilUpperPerc, airfoilLowerPerc)
            CM25 = CM + 0.25*CN
            xCP = -CM/CN
            cdRake = dragRake(U_inf, p_inf, rho, velocity_wake, p_backstatic, velocity_pos, staticWakemm)
            alpha = List[j, 2]
            ar = np.deg2rad(alpha)
            CL_Wake = CN*(np.cos(ar)+(np.sin(ar)**2)/np.cos(ar))-cdRake*np.tan(ar)
            
            Ax = math.cos(ar)*CA
            Ay = math.sin(ar)*CA
            
            Nx = math.sin(ar)*CN
            Ny = math.cos(ar)*CN
            
            CL = -Ay + Ny
            CD = Ax + Nx
            
            AOA.append(alpha)
            CL_Wake_List.append(CL_Wake)
            CD_Wake_List.append(cdRake)
            
            CL_List.append(CL)
            CD_List.append(CD)
            CM_List.append(CM25)
        plt.rcParams["figure.figsize"] = (8,4.5)
        plt.subplot(2, 2, 1)
        plt.plot(AOA,CL_Wake_List,'o--', label='Wake Rake Pressure Data', markersize=4, linewidth=1, color = 'orange', markerfacecolor='orange', markeredgecolor='black')
        plt.plot(AOA,CL_List,'^-', label='Surface Pressure Data', markersize=5, linewidth=1, color = 'cornflowerblue', markerfacecolor='cornflowerblue', markeredgecolor='black')
        plt.legend(loc = 'lower right')
        plt.grid()
        plt.xlabel(r'$\mathregular{\alpha}$ $[^\circ]$')
        plt.ylabel("$\mathregular{C_L}$ $[-]$")
        plt.subplot(2,2,2)
        plt.plot(AOA,CM_List, '^-', label='Surface Pressure Data', markersize=6, linewidth=1, color = 'cornflowerblue', markerfacecolor='cornflowerblue', markeredgecolor='black')
        plt.legend(loc = 'lower left')
        plt.grid()
        plt.xlabel(r'$\mathregular{\alpha}$ $[^\circ]$')
        plt.ylabel("$\mathregular{C_M}$ $[-]$")
        plt.subplot(2,2,3)
        plt.plot(CD_Wake_List,CL_Wake_List, 'o--', label='Wake Rake Pressure Data', markersize=5, linewidth=1, color = 'orange', markerfacecolor='orange', markeredgecolor='black')
        plt.plot(CD_List, CL_List,'^-', label='Surface Pressure Data', markersize=6, linewidth=1, color = 'cornflowerblue', markerfacecolor='cornflowerblue', markeredgecolor='black')
        plt.legend(loc = 'lower right')
        plt.xlabel("$\mathregular{C_D}$ $[-]$")
        plt.ylabel("$\mathregular{C_L}$ $[-]$")
        plt.grid()
        plt.subplot(2,2,4)
        plt.plot(AOA,CD_Wake_List, 'o--', label='Wake Rake Pressure Data', markersize=5, linewidth=1, color = 'orange', markerfacecolor='orange', markeredgecolor='black')
        plt.plot(AOA, CD_List, '^-', label='Surface Pressure Data', markersize=6, linewidth=1, color = 'cornflowerblue', markerfacecolor='cornflowerblue', markeredgecolor='black')
        plt.legend(loc = 'upper left')
        plt.xlabel(r'$\mathregular{\alpha}$ $[^\circ]$')
        plt.ylabel("$\mathregular{C_D}$ $[-]$")
        plt.grid()
        plt.show()

    if single == True:

        rho, viscous, U_inf, Rey, p_inf, pATM= ambient(Line,List)
        p_top, pos, p_bottom, pos2, p_backtotal, backtotal, p_backstatic, backstatic = pressure(Line,p_inf)

        p_axial_top = np.empty_like(p_top)
        p_axial_bottom = np.empty_like(p_bottom)

        for i in range(len(pos)):
            p_top[i]= (p_top[i]-p_inf)/(0.5 * rho *U_inf**2)
            p_axial_top[i] = axialComponent(p_top[i], airfoilUpperPerc[i]/.16, True)
        for i in range(len(pos2)):
            p_bottom[i]= (p_bottom[i]-p_inf)/(0.5 * rho *U_inf**2)
            p_axial_bottom[i] = axialComponent(p_bottom[i], airfoilLowerPerc[i]/.16, False)
        for i in range(len(backtotal)):
            cpt = (p_backtotal[i]-p_inf)/(0.5*rho*U_inf**2)
            if cpt > 1: cpt = 1
            p_backtotal[i]=cpt*0.5*rho*U_inf**2+p_inf
        backstatic_interpd = interp1d(staticWakemm, p_backstatic, bounds_error=False, fill_value="extrapolate")
        velocity_wake = []
        velocity_pos = []
        for i,x in enumerate(totalWakemm):
            if x < staticWakemm[0] or x > staticWakemm[-1]: continue
            vw = (2*(p_backtotal[i]-backstatic_interpd(x))/rho)**0.5
            velocity_wake.append(vw)
            velocity_pos.append(x)

        CN = cn(p_top, p_bottom, airfoilUpperPerc, airfoilLowerPerc)
        CA = cn(p_axial_top, p_axial_bottom, airfoilUpperPerc, airfoilLowerPerc)
        CM = cm(p_top, p_bottom, airfoilUpperPerc, airfoilLowerPerc)
        CM25 = CM + 0.25*CN
        xCP = -CM/CN
        cdRake = dragRake(U_inf, p_inf, rho, velocity_wake, p_backstatic, velocity_pos, staticWakemm)
        alpha = List[Line, 2]
        ar = np.deg2rad(alpha)
        CLWake = CN*(np.cos(ar)+(np.sin(ar)**2)/np.cos(ar))-cdRake*np.tan(ar)
        
        Ax = math.cos(ar)*CA
        Ay = math.sin(ar)*CA
        
        Nx = math.sin(ar)*CN
        Ny = math.cos(ar)*CN
        
        CL = -Ay + Ny
        CD = Ax + Nx
            
        print(f"AOA = {alpha}")
        print(f"normal coefficient = {CN}")
        print(f"axial coefficient = {CA}")
        print(f"lift coefficient = {CL}")
        print(f"drag coefficient = {CD}")
        print(f"moment coefficient = {CM25}")
        print(f"wake lift coefficient = {CLWake}")
        print(f"center of pressure = {xCP}")
        print(f"wake drag coefficient = {cdRake}")
        print(f"Reynolds numver = {Rey:.1e}")

        velocity = True
        if velocity:
            #plt.scatter(velocity_pos, velocity_wake, color='red', s=10, zorder=5)
            plt.rcParams.update({
                'font.size': 14,        # Set global font size
                'font.weight': 'medium'   # Set global font weight to bold
            })
            plt.plot(np.array(velocity_pos)*1000,velocity_wake, 'o-', label='Experimental data', markersize=4, linewidth=1, color = 'black', markerfacecolor='orange', markeredgecolor='black')
            plt.xlabel('Total pressure probe position [mm]', fontweight='medium')
            plt.ylabel('Wake velocity [m/s]', fontweight='medium')
            plt.ylim(16, 23)
            plt.grid()
            plt.legend(loc = 'lower right')
        else:
            plt.plot(pos,p_top)
            plt.plot(pos2,p_bottom)
            plt.gca().invert_yaxis()
            plt.grid()
        plt.show()
