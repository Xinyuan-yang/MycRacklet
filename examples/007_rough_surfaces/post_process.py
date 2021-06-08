import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys
import os
plt.style.use('fabstyle')

def extractData(output_dir,filename,nb_elem):
    time = np.loadtxt(output_dir+'Timer_'+filename)
    nb_t_step = time.shape[0]
    bin_data = np.memmap(output_dir+filename, dtype='float', mode='r')
    new_t_step = int(np.floor(bin_data.shape[0]/(nb_elem*nb_elem)))
    return np.reshape(bin_data[0:nb_elem*nb_elem*new_t_step],(new_t_step,nb_elem,nb_elem))

def loadSimuParameters(output_dir):

    param_file = open(output_dir+'Parameters.cra')
    dico = {}
    
    line = param_file.readline()
    values = line.split()	
    while (len(values) > 0):
        dico[values[0]] = float(values[1])	
        line = param_file.readline()
        values = line.split()

    return dico

assert(len(sys.argv)>1), 'python post_process.py <output_dir> '

output_dir = sys.argv[1]

try:
    os.mkdir(output_dir+'fig')
except:
    pass

parameters = loadSimuParameters(output_dir)
print(parameters.keys())
dom_size = float(parameters['dom_size_x'])
nu = float(parameters['nu_top'])
cs = float(parameters['cs_top'])
nb_elem = int(parameters['nb_ele_x'])
E = float(parameters['E_top'])

shear_velo = extractData(output_dir,'ST_Diagram_shear_velocity_jump.cra',nb_elem)
shear_strength = extractData(output_dir,'ST_Diagram_shear_strength.cra',nb_elem)

nb_steps = shear_velo.shape[0]

cmax = np.nanmax(np.log10(shear_velo[::,:,:]))
a = np.log10(shear_velo[::,:,:])
cmin = np.nanmin(a[a != -np.inf])

print(cmax,cmin)

figr, axr = plt.subplots(1, 2, figsize=(12, 4))
for tt in range(0,nb_steps):
    
    axr[0].imshow(np.log10(shear_velo[tt,:,:]), cmap='magma', origin='lower', extent=(0,1,0,1),clim=(cmin,cmax))
    axr[1].imshow(shear_velo[tt,:,:],origin='lower', extent=(0,1,0,1))
    axr[0].set_title('Slip velocity (log-scale)')
    axr[1].set_title('Slip velocity')
    figr.savefig(output_dir+"fig/rough"+str(10*nb_steps+tt)+".png",dpi=200)
    plt.close()
    
