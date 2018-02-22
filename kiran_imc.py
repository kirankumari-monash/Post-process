#!/usr/bin/puthon 
import os
import numpy as np
from scipy.io import netcdf
#from numpy.linalg import inv
iteration = 0
ite_done = 0
iteration_read = 1  # 1 for reading the number of iteration from the file
#iteration_read = 0
phi_new = 1 # 1 for reading the cp ref and calculating the new phi
cp_ens =0
b_obs_ens =0
cp_ens_1d = 0
cp_error =0
r_cum = 0
r_error = 0
rg_cum = 0
rg_error = 0
b_obs_ens_1d = 0
b_obs_error = 0
covar_1_ens = 0
covar_1_error = 0
cpr =1.82
al = 1.5306333
be = 1.213115524
epsilonr =0.1
k_bT = 1 
slicer =5 #Initial samples to be discarded
trajt = 50
for traj in range(trajt):
  f = netcdf.NetCDFFile('net_%s.nc' %(traj+101),'r',mmap=False)
#########################################################
#####   READING FROM NetCDF
#########################################################

  Nsample = f.dimensions['No_of_samples']
  print Nsample
  Ndim = f.dimensions['Ndim']
  NBeads = f.dimensions['NBeads']
  timenet = f.variables['Time']
  time = timenet[:]*1
  confinet = f.variables['configuration']
  confi_net = confinet[:]*1
  confi = confi_net[:,:,slicer:Nsample]
  Nsample = Nsample - slicer
#  gradnet = f.variables['Gradient']
  f.close() 

##########################################################
#####    PROPERTY CALCULATION
##########################################################
# COnfiguration is shifted according to the centre of mass 
# All the calculations will be done based of the shifted confi
  cmass = np.zeros((Ndim,Nsample))
  cmass = np.sum(confi[:],axis=0)
  cmass = cmass/(NBeads)
  confi = confi -cmass # shifted confi
##############  Re

  dist= np.zeros((Ndim, Nsample)) 
  r_sqs= np.zeros((Nsample)) 
  dist = confi[NBeads-1,:,:]-confi[0,:,:]
  dist = np.square(dist)
  r_sqs = np.sum(dist[:],axis=0)

###############   Rg #####################################
  dist_rg = np.zeros((NBeads,Ndim, Nsample))
  rg_sqs = np. zeros((Nsample))
  cmass = np.zeros((Ndim,Nsample))
  rg_sqs = np.zeros((NBeads,Nsample))
  rg = np.zeros((Nsample))
  dist_rg = np.square(confi)
  rg_sqs = np.sum(dist_rg[:],axis=1)
  rg = np.sum(rg_sqs[:],axis=0)
  rg = rg/(NBeads)
#  print 'rg', rg



###########  CONTACT PROBABILITY  ########################
  deltaR =np.zeros((NBeads,NBeads,Nsample))
  cp = np.zeros((NBeads,NBeads,Nsample))
  b2b = np.zeros((NBeads,NBeads,Ndim,Nsample))
  for nu in range (1,NBeads):
    for mu in range (nu):
      b2b[nu,mu,:,:] = confi[nu,:,:] - confi[mu,:,:]
  b2b = np.square(b2b[:]) 
  deltaR = np.sum(b2b[:],axis=2)
  deltaR = np.sqrt(deltaR[:])
  for Ns in range (Nsample):
    for nu in range (1,NBeads):
      for mu in range (nu):
        if deltaR[nu,mu,Ns] < cpr:
           cp[nu,mu,Ns] = cp[nu,mu,Ns]+1
           cp[mu,nu,Ns] = cp[nu,mu,Ns]
###############   1D 
  arg = ((NBeads-1)*NBeads)/2
  argr=0
  cp_1d = np.zeros((arg,Nsample))
  for nu in range (NBeads):
    for mu in range (nu):
      cp_1d[argr,:] = cp[nu,mu,:]
      argr = argr +1
  
#  print 'cp',cp_1d 
#########    Observable B_i     ############################
  b_obs = np.zeros((NBeads,NBeads,Nsample))
  for Ns in range (Nsample):
    for nu in range (NBeads):
      for mu in range (nu):
        if deltaR[nu,mu,Ns] <= 1.122462048:
           b_obs[nu,mu,Ns] = 1
        if deltaR[nu,mu,Ns] > 1.122462048 and deltaR[nu,mu,Ns] < 1.5:
           a = deltaR[nu,mu,Ns]
           a = (al*a)+be           
           b_obs[nu,mu,Ns] = -0.5*(np.cos(a)+1.0)
           b_obs[mu,nu,Ns] = b_obs[nu,mu,Ns]
#    print 'deltar',deltaR
########### 1D
  argr=0
  b_obs_1d = np.zeros((arg,Nsample))
  for nu in range (NBeads):
    for mu in range (nu):
      b_obs_1d[argr,:] = b_obs[nu,mu,:]
      argr = argr +1
#  print 'b_obs1d', b_obs_1d
#  print b_obs_1d.shape
#  print cp_1d.shape
##########   COvariance matrix   for <BA>  ##########################
  covar_1time = np.zeros((arg,arg))
  for Ns in range (Nsample):
    covar_1 = np.zeros((arg,arg))
    B = np.zeros((arg,1))
    A = np.zeros((arg,1))
    B[:,0] = b_obs_1d[:,Ns]
    BT = np.transpose(B)
    A[:,0] = cp_1d[:,Ns]
    AT = np.transpose(A)
    covar_1[:,:] = np.matmul(B,AT)
    print 'cov', covar_1
    covar_1time = covar_1time + covar_1
  print 'cp', cp_1d
  print 'b_obs', b_obs_1d
  print 'covar_1', covar_1time
  
#    print 'sample completed', Ns    
#  print 'Completed traj', traj
################# global averages
  r_ens = np.sum(r_sqs[:],axis=0)
  r_ens = r_ens/(Nsample)
  r_cum = r_cum + r_ens
  r_ens_sq = np.square(r_ens)
  r_error = r_error + r_ens_sq

  rg_ens = np.sum(rg[:],axis=0)
  rg_ens = rg_ens/(Nsample)
  rg_cum = rg_cum + rg_ens
  rg_ens_sq = np.square(rg_ens)
  rg_error = rg_error + rg_ens_sq
 
  cp_ens = np.sum(cp_1d[:],axis=1)
  cp_ens = cp_ens/(Nsample)
#  print 'cp_ens', cp_ens
  cp_ens_1d = cp_ens_1d + cp_ens
  cp_ens_sq = np.square(cp_ens) 
#  print 'cp_ens_sq', cp_ens_sq
  cp_error = cp_error + cp_ens_sq
  
  b_obs_ens = np.sum(b_obs_1d[:],axis=1)
  b_obs_ens = b_obs_ens/(Nsample)
  b_obs_ens_1d = b_obs_ens_1d + b_obs_ens
  b_obs_sq = np.square(b_obs_ens)
  b_obs_error = b_obs_error + b_obs_sq

  covar_1time = covar_1time/(Nsample)
  print 'covar_1_time', covar_1time
  covar_1sq = np.square(covar_1time)
  covar_1_ens = covar_1_ens + covar_1time
  print 'covar_1_ens', covar_1_ens 
#  covar_1_ens_sq = np.square(covar_1_ens)
  print 'sq', covar_1_ens
  covar_1_error = covar_1_error + covar_1sq
  print 'Completed traj', traj, covar_1_error
  print '-----------------------------------'
print '##################################################################'
#  print cp
##########################################################
######### Ensemble averaging 
##########################################################
# Re and Re error
r_cum = r_cum/(trajt)
print 're2', r_cum
r_cum_sq = np.square(r_cum)
r_error = r_error/(trajt)
r_error = np.subtract(r_error,r_cum_sq) 
r_error = r_error/(trajt-1)
r_error = np.sqrt(r_error)
print 're2_error',r_error

rg_cum = rg_cum/(trajt)
print 'rg2', rg_cum
rg_cum_sq = np.square(rg_cum)
rg_error = rg_error/(trajt)
rg_error = np.subtract(rg_error,rg_cum_sq)
rg_error = rg_error/(trajt-1)
rg_error = np.sqrt(rg_error)
print 'rg2_error',rg_error

# cp and CP error
cp_ens_1d = cp_ens_1d/(trajt)
cp_ens_1d_sq = np.square(cp_ens_1d)
cp_error = cp_error/(trajt)
cp_error = np.subtract(cp_error,cp_ens_1d_sq)
cp_error = cp_error/(trajt-1)
cp_error = np.sqrt(cp_error)
print 'cp', cp_ens_1d
print 'cp_error', cp_error
b_obs_ens_1d = b_obs_ens_1d/(trajt)
b_obs_ens_1d_sq = np.square(b_obs_ens_1d)
b_obs_error = b_obs_error/(trajt)
b_obs_error = np.subtract(b_obs_error,b_obs_ens_1d_sq)
b_obs_error = b_obs_error/(trajt-1)
b_obs_error = np.sqrt(b_obs_error)
print 'b', b_obs_ens_1d
print 'b_error', b_obs_error
covar_1_ens = covar_1_ens/(trajt)
print 'covar_1ens', covar_1_ens
covar_1_ens_sq = np.square(covar_1_ens)
print 'covar_ens_sq', covar_1_ens_sq
covar_1_error = covar_1_error/(trajt)
print 'sq', covar_1_error, trajt
covar_1_error = np.subtract(covar_1_error,covar_1_ens_sq)
print 'pata nai', covar_1_error
covar_1_error = covar_1_error/(trajt-1)
print 'co error', covar_1_error
covar_1_error = np.sqrt(covar_1_error)
print 'covar1_ens', covar_1_ens
print 'covar1_error', covar_1_error
#########################################################
################## cp to 2d
argr=0
cp_check = np.zeros((NBeads,NBeads))
cp_check_error = np.zeros((NBeads,NBeads))
for nu in range (NBeads):
  for mu in range (nu):
    cp_check[nu,mu] = cp_ens_1d[argr]
    cp_check[mu,nu] = cp_check[nu,mu]
    cp_check_error[nu,mu] = cp_error[argr]
    cp_check_error[mu,nu] = cp_check_error[nu,mu]
    argr = argr +1

#####################################################
##### covariave 2nd part
covar_2 = np.zeros((arg,arg))
covar_2_error = np.zeros((arg,arg))
covar_2_mag = np.zeros((arg,arg))
for nu in range (arg):
  for mu in range (arg):
    covar_2[nu,mu] = (b_obs_ens_1d[nu]) *( cp_ens_1d[mu])
    print 'covar2', covar_2[nu,mu], cp_ens_1d[mu], b_obs_ens_1d[nu]
    covar_2_error1 = b_obs_error[nu]/b_obs_ens_1d[nu]
    covar_2_error1 = np.square(covar_2_error1)
    covar_2_error2 = cp_error[mu]/cp_ens_1d[mu]
    covar_2_error2 = np.square(covar_2_error2)
    covar_2_error[nu,mu] = covar_2_error1 + covar_2_error2
    covar_2_error[nu,mu] = np.sqrt(covar_2_error[nu,mu])
covar_2_mag = np.square(covar_2)
covar_2_mag = np.sqrt(covar_2_mag)
covar_2_error = covar_2_mag * covar_2_error     
print 'co2', covar_2
print 'covar_2 error', covar_2_error
covar_final = k_bT*np.subtract(covar_1_ens[:],covar_2[:])
covar_1_error = np.square(covar_1_error)
covar_2_error = np.square(covar_2_error)
covar_error = k_bT*(covar_1_error + covar_2_error)
covar_error = np.sqrt(covar_error)
print 'covar_final', covar_final
print 'covar_error', covar_error
#print 'final', covar_final

#### Inverse of covariance matrix
#covar_inv = np.linalg.inv(covar_final)
#print 'covar_inv', covar_inv
#print 'cond_num', np.linalg.cond(covar_final)
#covar_check_1 = np.matmul(covar_inv,covar_final)
#covar_check_2 = np.dot(covar_inv,covar_final)
#print 'check_1', covar_check_1
#print 'check_2', covar_check_2

f_1d = open("cp_1d.txt","w")
#f_1d.write("Pairs\tCP\terror\n")
sno =1
for argr in range(arg):
  f_1d.write("%i\t%f\t%f\n"%(argr+1,cp_ens_1d[argr],cp_error[argr]))
f_1d.close()
if iteration_read == 1:
  f_ite = open("iteration_number.txt","r")
  line = f_ite.readline()
  ite_done =int(line)
  f_ite.close()
  iteration = ite_done +1
f_ite = open("iteration_number.txt","w")
f_ite.write("%i"%iteration)
f_ite.close()
fcpi = open("contprob_itera.txt","a+")
fcpi.write("Iteration number %i\n"%iteration)
fcpi.write("Contact Probability")
fcpi.write("\n")
for nu in range (NBeads):
  for mu in range (NBeads):
    fcpi.write("%f\t" %(cp_check[nu,mu]))
  fcpi.write("\n")
fcpi.write("\n")
fcpi.write("Error in Contact Probability")
fcpi.write("\n")
for nu in range (NBeads):
  for mu in range (NBeads):
    fcpi.write("%f\t" %(cp_check_error[nu,mu]))
  fcpi.write("\n")
fcpi.write("\n")
fcpi.write("\n")
fcpi.close()

fcpi = open("contprob_itera_1d.txt","a+")
fcpi.write("Iteration number %i\n"%iteration)
for argr in range (arg):
  fcpi.write("%i\t%f\t%f\n"%(argr+1,cp_ens_1d[argr],cp_error[argr]))
fcpi.close()


fcp = open("contprob_new.txt","w")
fcp.write("Contact Probability")
fcp.write("\n")
for nu in range (NBeads):
  for mu in range (NBeads):
    fcp.write("%f\t" %(cp_check[nu,mu]))
  fcp.write("\n")
fcp.write("\n")
fcp.write("Error in Contact Probability")
fcp.write("\n")
for nu in range (NBeads):
  for mu in range (NBeads):
    fcp.write("%f\t" %(cp_check_error[nu,mu]))
  fcp.write("\n")
fcp.write("\n")
fcp.write("\n")
fcp.close()

#############################################################
#######         Refining \phi values              ###########
#############################################################
########        Reference Contact probability    ############
if phi_new == 1:
  cp_ref = np.zeros((NBeads,NBeads))
  cp_ref_error  = np.zeros((NBeads,NBeads))
  f = open("cp_ref.txt","r")
  lines =  f.readlines()
  lineno = 1
  cp_lineno = 0
  cp_err_lineno = 0
  for line in lines:
    line = line.strip()
    data= line.split()
    if lineno >1 and lineno < NBeads +2:
      for mu in range (NBeads):
        cp_ref[cp_lineno,mu] = float(data[mu])
      cp_lineno+=1
    if lineno > NBeads +3 and lineno < 2*NBeads + 4:
      for mu in range (NBeads):
        cp_ref_error[cp_err_lineno,mu] = float(data[mu])
      cp_err_lineno +=1 
    lineno+=1  
  print 'cpref', cp_ref
  print 'cp_ref_err', cp_ref_error
  argr = 0
  cp_ref_1d = np.zeros((arg))
  cp_ref_error_1d = np.zeros((arg))
  for nu in range (NBeads):
    for mu in range (nu):
      cp_ref_1d[argr] = cp_ref[nu,mu]
      cp_ref_error_1d[argr] = cp_ref_error[nu,mu]
      argr = argr +1
#  print 'old', cp_ref_1d
################################################
############# RMSD
  rmsd_sq = np.subtract(cp_ref_1d[:],cp_ens_1d[:])
  rmsd_denom = np.square(cp_ref_1d[:])
  rmsd_de = np.sum(rmsd_denom[:],axis=0)
#  print 'rmsd_denom',rmsd_de
  rmsd_sq = np.square(rmsd_sq)
#  print 'rmsd', rmsd_sq
  rmsd = np.sum(rmsd_sq[:],axis=0)
#  print 'rmsd_up',rmsd
#  rmsd = rmsd/rmsd_de
#  print 'rmsd',rmsd
#  rmsd =2* rmsd/(NBeads*(NBeads-1))
  rmsd = rmsd/rmsd_de
  rmsd = np.sqrt(rmsd)
  rmsd = rmsd * 100
  print 'rmsd', rmsd
  f_rmsd = open("error_imc.txt","a+")
  f_rmsd.write("Iteration number %i\n"%iteration)
  f_rmsd.write("Percentage Error = %f\n"%rmsd )
  f_rmsd.write("\n")
  f_rmsd.close()


########## 1d

  phi_old = np.zeros((NBeads,NBeads))
  f = open("phi.txt","r")
  lines =  f.readlines()
  nu=0
  for line in lines:
    line = line.strip()
    data= line.split()
#    ptin1 = float(data[0])
#    ptin2 = float(data[1])
#    ptin3 = float(data[2])
    for mu in range (NBeads):
      phi_old[nu,mu] = float(data[mu])
    nu=nu+1
#  print 'phi_old', phi_old


#  f_ite = open("phi_ite.txt","a+")
#  f_ite.write("Iteration number = %i\n"%iteration)
#  for nu in range (NBeads):
#    for mu in range (NBeads):
#      f_ite.write("%f\t" %(phi_old[nu,mu]))
#    f_ite.write("\n")
#  f_ite.write("\n")
#  f_ite.close()

#########  converting phi to 1d
  argr = 0
  phi_old_1d = np.zeros((arg))
  for nu in range (NBeads):
    for mu in range (nu):
      phi_old_1d[argr] = phi_old[nu,mu]
      argr = argr +1
#  print 'phi_old', phi_old_1d
#  print 'cp', cp_time_1d 
#########
  delta_cp =np.subtract(cp_ref_1d[:],cp_ens_1d[:])
  print 'cp_new', cp_ens_1d
  delta_cp_error1 = np.square(cp_error)
  delta_cp_error2 = np.square(cp_ref_error_1d)
  delta_cp_error = delta_cp_error1 + delta_cp_error2
  delta_cp_error = np.sqrt(delta_cp_error)
  delta_cp_error = epsilonr* delta_cp_error
#  print 'delta cp', delta_cp
  delta_cp = epsilonr*delta_cp
#  print 'cp_ref-cp', delta_cp
#  phi_dummy = np.matmul(covar_inv,delta_cp)
#  print 'phi_dummy',phi_dummy
#  print 'covar',covar_final
  phi_new_1d = np.linalg.solve(covar_final,delta_cp)
  print 'phi_new',phi_new_1d
  phi_error1 =  np.matmul(covar_error,phi_new_1d)
  print 'phi_err1', phi_error1
#  print 'phi_without', phi_new_1d
  phi_error1 = np.subtract(delta_cp_error,phi_error1)
  phi_error_1d = np.linalg.solve(covar_final,phi_error1)
  print 'phi_error', phi_error_1d
  phi_new_1d = phi_new_1d + phi_old_1d
  print 'phi', phi_new_1d
#  phi_new_1d = np.dot(covar_inv,delta_cp)
#  phi_new_1d = epsilonr*phi_new_1d
#  phi_new_1d = phi_old_1d + phi_new_1d
#  print 'phi_new',phi_new_1d
########
#############################################################
#converting phi back to matrix
  argr=0
  phi_new = np.zeros((NBeads,NBeads))
  phi_error = np.zeros((NBeads,NBeads))
  for nu in range (NBeads):
    for mu in range (nu):
      phi_new[nu,mu] = phi_new_1d[argr]
      phi_new[mu,nu] = phi_new[nu,mu]
      phi_error[nu,mu] = phi_error_1d[argr]
      phi_error[mu,nu] = phi_error[nu,mu]
      argr = argr +1
  f_ite = open("phi_ite.txt","a+")
  f_ite.write("Iteration number = %i\n"%iteration)
  for nu in range (NBeads):
    for mu in range (NBeads):
      f_ite.write("%f\t" %(phi_new[nu,mu]))
    f_ite.write("\n")
  f_ite.write("\n")
  f_ite.write("Error in phi\n")
  for nu in range (NBeads):
    for mu in range (NBeads):
      f_ite.write("%f\t" %(phi_error[nu,mu]))
    f_ite.write("\n")
  f_ite.write("\n")
  f_ite.close()

  f_1d = open("phi_ite_1d.txt","a+")
  for argr in range(arg):
    f_1d.write("Iteration number = %i\n"%iteration)
    f_1d.write("%i\t%f\t%f\n"%(argr+1,phi_new_1d[argr],phi_error_1d[argr]))
  f_1d.write("\n")
  f_1d.close()


  print 'phi', phi_new
  file = open("phi.txt","w")
  for nu in range (NBeads):
    for mu in range (NBeads):
#      if phi_new[nu,mu] <0:
#        phi_new[nu,mu] = 0
#      if phi_new[nu,mu] >1:
#        phi_new[nu,mu] = 1
      file.write("%f\t" %(phi_new[nu,mu]))
    file.write("\n")
  file.close()
#  iteration = iteration +1 
#  f_ite = open("phi_ite.txt","a+")
#  f_ite.write("Iteration number = %i\n"%iteration)
#  for nu in range (NBeads):
#    for mu in range (NBeads):
#      if phi_new[nu,mu] <0:
#        phi_new[nu,mu] = 0
#      if phi_new[nu,mu] >1:
#        phi_new[nu,mu] = 1
#      f_ite.write("%f\t" %(phi_new[nu,mu]))
#    f_ite.write("\n")
#  f_ite.write("\n")
#  f_ite.close()

##############################################################
############           Trial for writting   
######
#   x = ptInfo[0]
#   print 'x',x
#print file_phi.read()
#for line in range (NBeads+1):
#  test = f.readlines()
#  print 'line', test
#cp_ref = file_phi.read()
#print 'ref',cp_ref[1,1]
#input = np.loadtxt("phi.txt", dtype='f', delimiter='\t')
#print(input)
#l = []
#with open('phi.txt', 'r') as f:
#  for line in f:
#    line = line.strip()
#    if len(line) > 0:
#      l.append(map(float, line.split('\t')))
#print l

