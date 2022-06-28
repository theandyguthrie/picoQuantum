#from picoQuantum.heat import ABCDModel
if __name__ ==  '__main__':
    from heat_3 import ABCDModel
    import sys
    import numpy as np
    from matplotlib import pyplot as plt

    model = ABCDModel(port_resistance = 12)
    model.add_TL(5.2e-3, 50, name = 'TL1')
    model.add_capacitor(10.25e-15, name = 'cap1')
    model.add_qubit(4.54e-9, 50e-15, name = 'qubit1')
    model.add_capacitor(10e-15, name = 'cap3')
    model.add_qubit(4.54e-9, 50e-15, name = 'qubit2')
    model.add_capacitor(10.25e-15, name = 'cap2')
    model.add_TL(5.2e-3, 50, name = 'TL2')

    omegas = 2*np.pi*np.linspace(1e9,10e9,600)
    deltas = np.linspace(-0.5*np.pi, 0.5*np.pi, 100)


    Lq = 4.54e-9/np.abs(np.cos(deltas))

    model.add_sweep_parameter('qubit1', 'Lq', Lq)
    
    model.set_parameter('qubit2', 'Lq', 4.54e-9/np.abs(np.cos(0)))
    heatflux1 = model.do_sweep(omegas)
    
    model.set_parameter('qubit2', 'Lq', 4.54e-9/np.abs(np.cos(0.4)))
    heatflux2 = model.do_sweep(omegas)
    
    model.set_parameter('qubit2', 'Lq', 4.54e-9/np.abs(np.cos(0.5)))
    heatflux3 = model.do_sweep(omegas)

    plt.plot(np.divide(deltas,np.pi),np.multiply(1e15,heatflux1))
    plt.plot(np.divide(deltas,np.pi),np.multiply(1e15,heatflux2))
    plt.plot(np.divide(deltas,np.pi),np.multiply(1e15,heatflux3))
    plt.show()