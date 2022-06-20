#input: Ts(T source), Td(T drain)
#output: Power at different delta(phase)
def calculate_photonic_power(data, Ts, Td):
    P = []
    for sParams in data:
         P.append(np.trapz(1e9*h*sParams['Frequency [GHz]'][1:]*(((1/(np.exp((1e9*h*sParams['Frequency [GHz]'][1:])/(k*Ts))-1))) - ((1/(np.exp((1e9*h*sParams['Frequency [GHz]'][1:])/(k*Td))-1))))*((sParams['S21'][1:])**2),1e9*sParams['Frequency [GHz]'][1:]))
    return P

 
def readCSV(path, N, ind):
   data_raw = np.array(pd.read_csv(path,skiprows = 3, delimiter=',', usecols=[0,ind], names=['Frequency', 'S21']))

   #Find number of frequency points
   try:
      for ind in range(4,3000): np.double(data_raw[ind])
   except ValueError:  numberofrequencypoints = ind 
   print('Number of frequency points: ',  numberofrequencypoints)
   
  


   #Loop over values of delta to turn 1D data into 2D array
   output = []
   for ind in range(0,N):
      skip = (ind)*3 + (ind*numberofrequencypoints)
      data = data_raw[skip: skip + numberofrequencypoints]
      

      output.append({'Frequency [GHz]' : np.double(data.T[0]), 'S21': np.double(data.T[1])})
      
   return output


def calculate_modulation_ratio(power):
   return (np.max(power) - np.min(power))/ np.max(power)

def pi_phase_shift(power):
   pi_ind = int(len(power)/4)
   return np.concatenate((power[pi_ind:],power[:pi_ind]))