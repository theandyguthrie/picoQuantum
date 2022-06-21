from picoQuantum.heat import ABCDModel

import sys


import numpy as np
from matplotlib import pyplot as plt

model = ABCDModel()
model.add_TL(5.2e-3, 50, name = 'TL1')
model.add_capacitor(10.25e-15, name = 'cap1')
model.add_qubit(4.54e-9, 110e-15, name = 'qubit1')
model.add_capacitor(10.25e-15, name = 'cap2')
model.add_TL(5.2e-3, 50, name = 'TL2')

omegas = 2*np.pi*np.linspace(1e9,10e9,300)
deltas = np.linspace(0, 2*np.pi, 100)

#model.compute_ABCD(omegas,[0])
#plt.plot(1e-9*omegas/(2*np.pi), model.calculate_S21()[0])

plt.plot(np.divide(deltas,np.pi),np.multiply(1e15,model.heat_flux(omegas,deltas)))
plt.show()