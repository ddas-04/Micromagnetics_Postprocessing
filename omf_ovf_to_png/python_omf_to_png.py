import numpy as np
import os
import fnmatch
import re
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size':32})
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams.update({'font.weight':'bold'})
plt.rcParams["font.family"] = "Times New Roman"
rc = {"font.family" : "serif", 
      "mathtext.fontset" : "stix"}
plt.rcParams.update(rc)
plt.rcParams["font.serif"] = ["Times New Roman"] + plt.rcParams["font.serif"]
############### Simulation parameter values used ############################
Ms=1
############ Convert list to string #########################################
def listToString(s):     
    # initialize an empty string
    str1 = " "     
    # return string  
    return (str1.join(s))
############ Get numerical data for the variable ##############################
def findnumericaldata(list_var,lines):
	matching = [s for s in lines if list_var in s]
	match_string=listToString(matching)
	s = [float(s) for s in re.findall(r'-?\d+\.?\d*', match_string)]
	s=np.asarray(s)
	if len(s)==2:
		s_num=s[0]*10**(s[1])
		s=[]
		s=s_num
	return s 
############# OMF file data extraction #####################################
def omfdataextraction(filename):
	data=open(filename)
	lines=data.readlines()
	data.close()
	var=['xmin', 'xmax','xnodes', 'ymin','ymax','ynodes', 'zmin','zmax', 'znodes', 'Total simulation time']
	numerical_data=np.zeros(len(var))
	for i in range(len(var)):
		var_i=var[i]
		#print(findnumericaldata(var_i))
		numerical_data[i]=findnumericaldata(var_i,lines)
	# create x, y, z array
	x=np.linspace(numerical_data[0],numerical_data[1],int(numerical_data[2]))
	y=np.linspace(numerical_data[3],numerical_data[4],int(numerical_data[5]))
	z=np.linspace(numerical_data[6],numerical_data[7],int(numerical_data[8]))
	# create meshgrid for 2d data
	[X,Y,Z]=np.meshgrid(x,y,z)
	# Obtain Ms values begining and end index
	matching = [s for s in lines if "Data Text" in s]
	Ms_begin_index=lines.index(matching[0])+1
	Ms_end_index=lines.index(matching[1])
	Ms_val=np.loadtxt(lines[Ms_begin_index:Ms_end_index])
	time_val=numerical_data[9]
	return x, X, y, Y, z, Z, Ms_val,time_val
############ Get the number and list of omf files in the directory #################
N=0
omf_files_list=[]
for file in sorted(os.listdir('.')):
	if fnmatch.fnmatch(file, '*.omf'):
		omf_files_list.append(file)
		N=N+1

Sim_time=np.zeros(N)
filename=omf_files_list[-1]
[x,X,y,Y,z,Z,S,t_i]=omfdataextraction(filename)
Total_time=t_i*1e9
Pulse_signal=[]
time_i=[]
tp=1
tw=0.5
for i in range(N):
	filename=omf_files_list[i]
	print(filename)
	[x,X,y,Y,z,Z,S,t_i]=omfdataextraction(filename)
	save_filename=filename.replace('.omf', '')
	
	Sim_time[i]=t_i*1e9
	time_i=np.append(time_i,Sim_time[i])
	n=int(Sim_time[i]/tp)
	
	time=t_i*1e9
	max_x = max(x)*1e9
	max_y = max(y) * 1e9
	
	
	X_mat=X[:,:,0]*1e9 # converted into nm
	Y_mat=Y[:,:,0]*1e9 # converted into nm
	#print(X_mat.shape)
	S=S/Ms    # normalize the magnetization vector
	mz_1d=S[:,2]
	mz_reshaped=np.reshape(mz_1d,(len(z),len(y),len(x)))
	
	l=len(z)-1
	mz_2d_top=mz_reshaped[l,:,:]
	
	
	plt.figure(figsize=(12,3))
	
	plt.pcolor(X_mat, Y_mat, mz_2d_top,cmap='bwr') #,shading='auto')
	
	plt.title('Time = %.2f ns' %(time))
	plt.xticks([0, max_x])
	plt.yticks([0,max_y])
	plt.tick_params(which='both', direction='in', length=6, width=2, colors='k')
	plt.savefig(save_filename+'.png', bbox_inches='tight', pad_inches=0.0)
	plt.close('all')

