import numpy as np
from scipy.integrate import simps


class ABCDModel():
    
    def __init__(self):
        self.elements = []
        
        
    def add_qubit(self, Lq,Cq, name = 'name'):
        self.elements.append({'type': 'qubit', name: name, 'Lq': Lq, 'Cq':Cq})
        
    def add_capacitor(self, Cg, name = 'name'):
        self.elements.append({'type': 'capacitor', name: name, 'Cg':Cg})
        
    def add_TL(self, l, Z0, name = 'name'):
        self.elements.append({'type': 'TL', name: name, 'l': l, 'Z0': Z0})
    

    def qubit(self, omega, element):
        Lq =   element['Lq']/np.abs(np.cos(self.flux))
        return np.matrix([[1, 0], [(1j*omega*Lq+1/(1j*omega*element['Cq']))/(1j*omega*Lq*1/(1j*omega*element['Cq'])), 1]])
    
    def capacitor(self, omega, element):
        return np.matrix([[1, 1/(1j*omega*element['Cg'])], [0, 1]])
        
    def TL(self,omega, element):
        L = 405e-9
        C = 171e-12
        #Z0 = np.sqrt(L/C)
        Y0 = 1 / element['Z0']
        beta = omega*np.sqrt(L*C)
        return np.matrix([[np.cos(beta*element['l']), 1j*element['Z0']*np.sin(beta*element['l'])], [1j*Y0*np.sin(beta*element['l']), np.cos(beta*element['l'])]])
    
    
    def compute_ABCD(self, omegas, fluxes):
        
        self.omegas = omegas   
        self.fluxes = fluxes
        self.ABCD = np.empty((len(self.fluxes), len(self.omegas), 2, 2), dtype = np.complex64) 
        
        for j in range(0,len(fluxes)):
            self.flux = fluxes[j]

            for i in range (0,len(omegas)):
                product = np.ones((2,2))
                for element in self.elements:
                    elem = getattr(self, element['type'])
                    m = elem(omegas[i], element)
                    product = product@m
                    
                self.ABCD[j][i] = product
        
        
            
    def calculate_S21(self):
        
        self.S21 = np.empty((len(self.fluxes), len(self.omegas)), dtype = np.float32)
        for i in range(0,len(self.ABCD)):
            for j in range(0, len(self.omegas)):
                M = self.ABCD[i][j]
                R = 12
                A = M[0,0]
                B = M[0,1]
                C = M[1,0]
                D = M[1,1]

                self.S21[i][j] = (np.abs(2/(A+B/R+R*C+D)))
        
       
    
    
    def calculate_photonic_power(self, Ts, Td):
        k =  1.38065*1e-23
        e=1.60e-19
        h=6.63e-34
     
        P = []
        for i in range(0,len(self.S21)):
            P.append(simps(h*self.omegas*(((1/(np.exp((h*self.omegas)/(k*Ts))-1))) - ((1/(np.exp((h*self.omegas)/(k*Td))-1))))*(np.power(self.S21[i],2)),self.omegas))
        return P

    
    
    def heat_flux(self,frequency,fluxes):
        
        self.compute_ABCD(frequency, fluxes)
        self.calculate_S21()
        P = self.calculate_photonic_power(0.35,0.12)
        
        return P

 