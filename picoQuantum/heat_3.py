import numpy as np
from scipy.integrate import simps
import time
from progress.bar import Bar


class ABCDModel():
    
    def __init__(self, port_resistance = 12):
        self.elements = []
        self.port_resistance = port_resistance
        
        
        
    def add_qubit(self, Lq,Cq, name = 'name'):
        self.elements.append({'type': 'qubit', 'name': name, 'Lq': Lq, 'Cq':Cq})
        
    def add_capacitor(self, Cg, name = 'name'):
        self.elements.append({'type': 'capacitor', 'name': name, 'Cg':Cg})
        
    def add_TL(self, l, Z0, name = 'name'):
        self.elements.append({'type': 'TL', 'name': name, 'l': l, 'Z0': Z0})
    

    def qubit(self, omega, element):
        return np.matrix([[1, 0], [(1j*omega*element['Lq']+1/(1j*omega*element['Cq']))/(1j*omega*element['Lq']*1/(1j*omega*element['Cq'])), 1]])
    
    def capacitor(self, omega, element):
        return np.matrix([[1, 1/(1j*omega*element['Cg'])], [0, 1]])
        
    def TL(self,omega, element):
        L = 405e-9
        C = 171e-12
        Y0 = 1 / element['Z0']
        beta = omega*np.sqrt(L*C)
        return np.matrix([[np.cos(beta*element['l']), 1j*element['Z0']*np.sin(beta*element['l'])], [1j*Y0*np.sin(beta*element['l']), np.cos(beta*element['l'])]])
    
    
    def compute_ABCD(self, omegas):
        
        self.omegas = omegas   
        self.ABCD = np.empty((len(self.omegas), 2, 2), dtype = np.complex64) 
        

        for i in range (0,len(omegas)):
            product = [[1, 0], [0, 1]]
            for element in self.elements:
                elem = getattr(self, element['type'])
                m = elem(omegas[i], element)
                product = product@m
                
            self.ABCD[i] = product
    
        
        return self.ABCD
    
    
    
    def calculate_S21(self, omegas):

        self.compute_ABCD(omegas)
        self.S21 = np.empty((len(self.omegas)), dtype = np.float64)

        for j in range(0, len(self.omegas)):
            M = self.ABCD[j]
            R = self.port_resistance
            A = M[0,0]
            B = M[0,1]
            C = M[1,0]
            D = M[1,1]

            self.S21[j] = np.abs(2/(A+B/R+R*C+D))
        return self.S21
        
       

    
    
    def calculate_photonic_power(self, Ts, Td):
        k = 1.38065*1e-23
        e = 1.60e-19
        h = 6.63e-34
        self.frequency = self.omegas / (2*np.pi)
        
        P = np.trapz(h*self.frequency*(((1/(np.exp((h*self.frequency)/(k*Ts))-1))) - ((1/(np.exp((h*self.frequency)/(k*Td))-1))))*(np.power(self.S21,2)),self.frequency)
        return P

    
    
    def calculate_heat_flux(self):
    
        self.calculate_S21(self.omegas)
        P = self.calculate_photonic_power(0.35,0.12)
        
        return P
    
    def add_sweep_parameter(self, component, parameter, sweep_values):
        
        self.sweep_component_name = component
        self.sweep_parameter = parameter
        self.sweep_values = sweep_values
        
        
        
        
    def set_parameter(self, component, parameter, value):
        for index in range(len(self.elements)):
                if  self.elements[index]['name'] == component:
                    self.elements[index][parameter] =  value
        
        print('Setting ' + component + ' ' + parameter + ' to ' + str(value))
        
        
    def do_sweep(self, omegas):
        self.omegas = omegas
        self.heat_flux = np.empty((len(self.sweep_values)), dtype = np.float64)
        
        with Bar('Sweeping ' + self.sweep_component_name + ' '+  self.sweep_parameter,max = len(self.sweep_values)) as bar:
            for i in range(0,len(self.sweep_values)):
                for index in range(len(self.elements)):
                    if  self.elements[index]['name'] == self.sweep_component_name:
                        self.elements[index][self.sweep_parameter] =  self.sweep_values[i]
                    
                self.heat_flux[i] = self.calculate_heat_flux()        
                bar.next()
        return self.heat_flux
            
            
        
        
        

 