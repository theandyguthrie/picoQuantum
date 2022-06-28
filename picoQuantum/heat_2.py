
import numpy as np
from scipy.integrate import simps
import multiprocessing
import time

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
        Lq = (element['Lq']/np.abs(np.cos(self.delta))*(np.sqrt(1+(0.1**2)*(np.tan(self.delta))**2)))
        return np.matrix([[1, 0], [(1j*omega*Lq+1/(1j*omega*element['Cq']))/(1j*omega*Lq*1/(1j*omega*element['Cq'])), 1]])
    
    def capacitor(self, omega, element):
        return np.matrix([[1, 1/(1j*omega*element['Cg'])], [0, 1]])
        
    def TL(self,omega, element):
        L = 405e-9
        C = 171e-12
        Y0 = 1 / element['Z0']
        beta = omega*np.sqrt(L*C)
        return np.matrix([[np.cos(beta*element['l']), 1j*element['Z0']*np.sin(beta*element['l'])], [1j*Y0*np.sin(beta*element['l']), np.cos(beta*element['l'])]])
    
    
    
    def do_matrix_multiplication(self, j):
        

        ABCD = np.empty((len(self.omegas), 2, 2), dtype = np.complex64) 
        self.delta = self.deltas[j]
        
        for i in range(0,len(self.omegas)):
            product = [[1, 0], [0, 1]]
            for element in self.elements:
                elem = getattr(self, element['type'])
                m = elem(self.omegas[i], element)
                product = product@m
            
            ABCD[i] = product
            
            
        return ABCD
            
        
    
    def compute_ABCD(self, omegas, deltas):
        
        self.omegas = omegas   
        self.deltas = deltas
        self.ABCD = np.empty((len(self.deltas), len(self.omegas), 2, 2), dtype = 'c') 
        
        a_pool = multiprocessing.Pool(processes=3)
        
        self.ABCD = a_pool.map(self.do_matrix_multiplication, range(0,len(deltas)))
        

        return self.ABCD
    
    
    
    
    def calculate_S21(self, omegas, deltas):
        t0 = time.time()
        self.compute_ABCD(omegas,deltas)
        print(time.time() - t0)
        self.S21 = np.empty((len(self.deltas), len(self.omegas)), dtype = np.float32)
        for i in range(0,len(self.ABCD)):
            for j in range(0, len(self.omegas)):
                M = self.ABCD[i][j]
                R = 12
                A = M[0][0]
                B = M[0][1]
                C = M[1][0]
                D = M[1][1]
                self.S21[i][j] = np.abs(2/(A+B/R+R*C+D))
        return self.S21
        
       

    
    
    def calculate_photonic_power(self, Ts, Td):
        k = 1.38065*1e-23
        e = 1.60e-19
        h = 6.63e-34
        self.frequency = self.omegas / (2*np.pi)
        P = []
        for i in range(0,len(self.S21)):
            P.append(np.trapz(h*self.frequency*(((1/(np.exp((h*self.frequency)/(k*Ts))-1))) - ((1/(np.exp((h*self.frequency)/(k*Td))-1))))*(np.power(self.S21[i],2)),self.frequency))
        return P

    
    
    def heat_flux(self,omegas,deltas):
        
        self.calculate_S21(omegas, deltas)
        P = self.calculate_photonic_power(0.35,0.12)
        
        return P

 