import numpy as np
import os
import fnmatch
import re
import matplotlib.pyplot as plt

import matplotlib as mpl
mpl.rcParams['text.usetex'] = False

plt.rcParams.update({'font.size':24})
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams.update({'font.weight':'bold'})
plt.rcParams["font.family"] = "Times New Roman"
rc = {"font.family" : "serif", 
      "mathtext.fontset" : "stix"}
plt.rcParams.update(rc)
plt.rcParams["font.serif"] = ["Times New Roman"] + plt.rcParams["font.serif"]
############### Simulation parameter values used ############################
Ms=580e3
p=0.4
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
	# list of variable values that needs to be extracted
	var=['xmin', 'xmax','xnodes', 'ymin','ymax','ynodes', 'zmin','zmax', 'znodes', 'Desc:  Total simulation time:']
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
################ Normalize the weights from G ##############################	
def diff_norm(array):
	max_a=max(array)
	min_a=min(array)
	diff=array-min_a
	diff_norm=diff/(max_a-min_a)
	return diff_norm
################ Pulse number and Weights at pulse extractor ###################
def pulse_extractor(Sim_time, weights):	
	nt=0
	pulse_number=[]
	weight_pulse=[]
	for ti in range(len(Sim_time)):
		if Sim_time[ti]>=nt:
			pulse_number=np.append(pulse_number,nt+1)
			weight_pulse=np.append(weight_pulse,weights[ti])
			nt=nt+1
			
	return pulse_number, weight_pulse
################# function for sign change location for DW #####################
def sign_change_position(a):
	asign = np.sign(a)
	signchange = ((np.roll(asign, 1) - asign) != 0).astype(int)
	pos=np.asarray(np.where(signchange==1.0)).flatten()
	return pos
############ Get the number and list of omf files in the directory #################
N=0
omf_files_list=[]
for file in sorted(os.listdir('.')):
	if fnmatch.fnmatch(file, '*.omf'):
		omf_files_list.append(file)
		N=N+1

#G=np.zeros(N)
Sim_time=np.zeros(N)
Pulse_signal=[]
time_i=[]
tp=1.0
tw=0.5

filename=omf_files_list[-1]
print(filename)
[x,X,y,Y,z,Z,S,t_i]=omfdataextraction(filename)
Total_time=t_i*1e9
print(Total_time)

#exit()
DW_left_pos_array=[]
DW_right_pos_array=[]
DW_mid_pos_array=[]


for i in range(N):
	print('-------------------------------------------------------------------')
	filename=omf_files_list[i]
	print(filename)
	[x,X,y,Y,z,Z,S,t_i]=omfdataextraction(filename)
	save_filename=filename.replace('.omf', '')
	Sim_time[i]=t_i*1e9
	
	
	# for top layer 
	l=len(z)-1
	X_mat=X[:,:,l]*1e9 # converted into nm
	Y_mat=Y[:,:,l]*1e9 # converted into nm
	
	mz_1d=S[:,2]
	mz_reshaped=np.reshape(mz_1d,(len(z),len(y),len(x)))
	mz_2d_TL=mz_reshaped[l,:,:]
	[ny,nx]=mz_2d_TL.shape

	mid=int(ny*0.5)
	mz_along_x=mz_2d_TL[mid,:]
	pos_index=sign_change_position(mz_along_x)
	left_DW_pos_index=pos_index[0]
	right_DW_pos_index=pos_index[1]
	x=x*1e9	
	left_DW_pos=x[left_DW_pos_index]
	right_DW_pos=x[right_DW_pos_index]
	DW_mid_pos=0.5*(left_DW_pos+right_DW_pos)
	
	DW_left_pos_array=np.append(DW_left_pos_array,left_DW_pos)
	DW_right_pos_array=np.append(DW_right_pos_array,right_DW_pos)
	DW_mid_pos_array=np.append(DW_mid_pos_array,DW_mid_pos)
	
np.save('Sim_time.npy',Sim_time)
np.save('DW_left_pos_array.npy',DW_left_pos_array)
np.save('DW_right_pos_array.npy',DW_right_pos_array)
np.save('DW_mid_pos_array.npy',DW_mid_pos_array)

plt.figure(figsize=(10,7))
plt.plot(Sim_time,DW_left_pos_array, linewidth=3.0, label='left side')
plt.plot(Sim_time,DW_right_pos_array, linewidth=3.0, label='right side')
plt.plot(Sim_time,DW_mid_pos_array, linewidth=3.0, label='mid pos')
plt.ylim([-10, 520])
plt.legend()
#plt.yticks([0,256,512])
plt.grid()
plt.savefig('Domain_Step_2_LIF.png', bbox_inches='tight', pad_inches=0.1)
plt.show()

	
'''	
	
	
	
	
	
	#fig, ax = plt.subplots(2, 1,figsize=(21.33,11.25), gridspec_kw={'height_ratios': [1.5, 1]})
	fig, ax = plt.subplots(3, 1,figsize=(15,8), gridspec_kw={'height_ratios': [2.5, 1, 1]})
	
	ax[0].plot(time_i,Pulse_signal, linewidth=2.5)
	ax[0].set_xlabel('Time(ns)') #,fontsize=18)
	ax[0].set_ylabel(r"$J\mathrm{(GA/m^2)}$")
	ax[0].set_xlim([0, Total_time])
	ax[0].set_ylim([-0.4, 4.4])
	ax[0].set_yticks([0, 4])
	ax[0].set_title('Time = %.4f ns' %Sim_time[i])

	v1_x=370*np.array([1,1])
	v1_y=np.array([0, 32])
	v2_x=390*np.array([1,1])
	v2_y=np.array([0, 32])
	v3_x=430*np.array([1,1])
	v3_y=np.array([0, 32])
	v4_x=450*np.array([1,1])
	v4_y=np.array([0, 32])
	ax[1].pcolor(X_mat, Y_mat, mz_2d_TL,cmap='bwr') #,shading='auto')

	ax[1].plot(v1_x,v1_y,'k', linewidth=3.0)
	ax[1].plot(v2_x,v2_y,'k--', linewidth=3.0)
	ax[1].plot(v3_x,v3_y,'k--', linewidth=3.0)
	ax[1].plot(v4_x,v4_y,'k', linewidth=3.0)
	#ax[0].plot(len(x)*np.ones(len(y)),y*1e9, 'k', linewidth=4.0)
	ax[1].set_xticks([])#0, 256, 512,768,1024])
	ax[1].set_yticks([0, 32])

	ax[2].pcolor(X_mat, Y_mat, mz_2d_BL,cmap='bwr') #,shading='auto')
	ax[2].plot(v1_x,v1_y,'k', linewidth=3.0)
	ax[2].plot(v2_x,v2_y,'k--', linewidth=3.0)
	ax[2].plot(v3_x,v3_y,'k--', linewidth=3.0)
	ax[2].plot(v4_x,v4_y,'k', linewidth=3.0)
	ax[2].set_xticks([0, 256, 512])
	ax[2].set_yticks([0, 32])

	plt.savefig(save_filename+'.png', bbox_inches='tight', pad_inches=0.1)
	plt.close('all')
	
'''
