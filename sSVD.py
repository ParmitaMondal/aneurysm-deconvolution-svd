import os
import numpy as np
import pandas as pd
import pydicom
import matplotlib.pyplot as plt
from scipy.signal import medfilt
from scipy.signal import correlate
from PIL import Image
from scipy import interpolate
from scipy.linalg import toeplitz, svd, pinv
from openpyxl import load_workbook
from scipy.stats import gamma
from scipy import signal
from scipy.optimize import curve_fit 
import math

def API_parameters(x, y):
    j = 0
    BAT_index = 0
    while j < len(x):
        if (y[j] > np.max(y)*0.1):
            BAT_index = j-1
            break
        j += 1
    PH=np.max(y)  
    j = 0
    PH_index = 0
    while j < len(x):
        if (y[j] == PH):
            PH_index = j
            break
        j += 1
    TTP=x[PH_index]-x[BAT_index]
    AUC=np.trapz(y, x)
    MTT=np.trapz(y*x, x)/AUC
    derivative=np.gradient(y[BAT_index:PH_index])
    max_Df=derivative.max()
    return x[BAT_index], MTT, TTP, PH, AUC, max_Df

def API_SVD(x_inl, y_inl, x_ANY, y_ANY): #p
      
    j = 0
    BAT_index = 0
    while j < len(x_inl):
        if (y_inl[j] > np.max(y_inl)*0.1):
            BAT_index = j-1
            break
        j += 1
   
    P_SVD= 0.04 #change 1
    
    inlet_zeros = np.zeros([len(y_inl)])
    inlet_toeplitz = toeplitz(y_inl, inlet_zeros)
    V, S, U = svd(np.matrix.transpose(inlet_toeplitz))
    S = np.diag(S)
    truncate=P_SVD*np.max(S)
    lam=0.1
    S[S<truncate]=0
    S[S>truncate]= S[S>truncate]/(S[S>truncate]**2 + lam**2)
     
     # Compute impulse response function
    IRF = abs(np.matmul(np.matmul(V, S), np.matmul(U, y_ANY)))
    IRF1 = IRF/np.max(IRF)  
    RBF = np.max(IRF)
    RBV = np.trapz(IRF, x_inl)
    MTT_SVD = RBV / RBF 
    
    return RBF, RBV, MTT_SVD, IRF

mask_aneurysm_path = "path..." 
mask_inlet_path="path..." 
excel_file = 'path...'
sheet = 'Sheet1'
save_file = 'path...'
df=load_workbook(save_file)
ws=df.worksheets[0]
caseInfo = pd.read_excel(excel_file)
cases_mask = caseInfo.iloc[:,0]
cases_dicom=caseInfo.iloc[:,1]

for x in range(0,1):  

    import re
    import glob
    proj = np.empty((511,205,1578), dtype=np.float32) 
    numbers = re.compile(r'(\d+)')
    def numericalSort(value):
      parts = numbers.split(value)
      parts[1::2] = map(int, parts[1::2])
      return parts

    filelist = sorted(glob.glob('path...'+"*.raw"), key=numericalSort) 
    for file, idx in enumerate(filelist):
       proj[:,:,file] = np.fromfile(open(filelist[file], 'rb'), dtype=np.float32).reshape((511,205)) 
    plt.imshow(proj[:,:,900])
    ANY_case=cases_dicom[x]
    ANY_case_aneurysm_mask=cases_mask[x]
    print(ANY_case)
    dicom_temp=proj
    dicom_image=dicom_temp
    dicom_image[dicom_image==-np.inf]=0
    dicom_mean=np.mean(dicom_image[:,:,400:np.shape(dicom_image)[2]-100], axis=2)
    dimension=dicom_temp.shape[0]
    TDC_average=np.zeros(shape=np.shape(dicom_image)[2], dtype=np.float32)
    TDC_inlet_average=np.zeros(shape=np.shape(dicom_image)[2], dtype=np.float32)
    time_vector=np.arange(0,1.578,0.001) 
    mask_aneurysm=np.zeros(shape=(dimension,dimension), dtype=np.float32)
    mask_inlet=np.zeros(shape=(dimension,dimension), dtype=np.float32)
    mask_aneurysm = Image.open(os.path.join(mask_aneurysm_path,  ANY_case_aneurysm_mask +'_0.tif')) 
    mask_inlet = Image.open(os.path.join(mask_inlet_path,  ANY_case_aneurysm_mask +'_0_inl.tif'))  
    im=np.asarray(mask_aneurysm, dtype=np.float32)
    im2=np.asarray(mask_inlet, dtype=np.float32)
    mask_display=np.asarray(mask_aneurysm)+np.asarray(mask_inlet)
    im=abs((im-255)/255)
    im2=abs((im2-255)/255)
    ind = np.transpose(np.nonzero(im))
    ind_inlet = np.transpose(np.nonzero(im2)) 
    print(im.sum(), im2.sum())
    for mm in range(0, len(ind)):
        TDC_average = TDC_average+dicom_image[ind[mm][0],ind[mm][1] ,: ]/im.sum()
       
    for nn in range(0, len(ind_inlet)):
        TDC_inlet_average = TDC_inlet_average+dicom_image[ind_inlet[nn][0], ind_inlet[nn][1], :]/im2.sum()
 
    TDC_average=medfilt(TDC_average, kernel_size=1)
    TDC_average[TDC_average<0]=0 
    
   
    ANY_results=np.array(API_parameters(time_vector, TDC_average))
    
    TDC_inlet_average=medfilt(TDC_inlet_average, kernel_size=1)
    TDC_inlet_average[TDC_inlet_average<0]=0 
    
    I_inlet_results=np.array(API_parameters(time_vector, TDC_inlet_average)) #without gamma-fitting
   
    I_inlet_results[0]=1
    cor_pre=correlate(TDC_inlet_average-np.mean(TDC_inlet_average), TDC_average-np.mean(TDC_average))/(len(TDC_average) * np.std(TDC_inlet_average) * np.std(TDC_average))
    ANY_results[0]=cor_pre.max()
    
    plt.plot(time_vector, TDC_average, time_vector, TDC_inlet_average)
    plt.xlabel('Time')
    plt.ylabel('Contrast density')
    plt.gca().legend(('TDC of the aneurysm dome ','TDC of the inlet '))
    plt.show()
  
    residual_func_=API_SVD(time_vector, TDC_inlet_average, time_vector, TDC_average) 
    residual_function_= residual_func_[3]
    residual_function1_ =residual_function_/np.max(residual_function_) #impulse response function
    k=signal.fftconvolve(TDC_inlet_average, residual_function_) #convolve the inlet always with original residue function to get the aneurysm TDC
    k1=k[0:len(time_vector),]  
    plt.plot(time_vector, k1)
    plt.show()
    res=residual_func_[0:3]
    
    plt.plot(time_vector, k1, time_vector, TDC_average) 
    plt.xlabel(('Time'),  fontsize="15")
    plt.ylabel(('Contrast density'), fontsize="15")
    plt.gca().legend(('$Q_{new}$', 'Q'), fontsize="12", loc="upper left")
    plt.show()
    plt.plot(time_vector, residual_function1_, color='red', label='Plot 1') #plotting impulse response function 
    plt.xlabel(('Time'),  fontsize="15")
    plt.ylabel(('Impulse Response Function'), fontsize="15")
    plt.show()
    
    
   
