"""
Usage:
python KG_processing.py <Gaussian Width in eV> <Temperature in eV> <frequency grid spacing in eV>
 <Gaussian Width in eV> should be greater or equal to that specified in SHRED input file
 <Temperature in eV> should be the same as in SHRED input file
 <frequency grid spacing in eV> sets cutoff of time signal or zero-padding to artificially add points

script reads Kubo_Greenwood_Cmn_t.data output from SHRED with run_Kubo_Greenwood=.true.
script outputs stoc_KG.txt with AC electrical conductivity, thermal conductivity and thermopower
script outputs KG_Lmn.txt with AC Onsager Coefficients
"""

import os
import sys
import scipy.interpolate
import numpy as np

def do_fft(filename2, filename3, width, temp, dw):
        data=np.loadtxt("Kubo_Greenwood_Cmn_t.data", skiprows=1)
        time=data[:,0]
        dt=time[1]-time[0]
        nt_of_dw=int(1.0/dw/dt)
        y=np.fft.fftfreq(nt_of_dw,d=dt)
        w=(y[:int(len(y)/2-1)]+1E-10)*(2.0*np.pi)
        print(dt)
        print(y)
        print(w*27.211)
        max_t=np.argmax(time)
        print("length of the signal: "+ str(max_t) + " max: " + str(time[max_t]), "dw: "+ str(dw*27.211*2*np.pi)   )
        print("max w "+ str(int(len(y)/2-1)*dw*27.211*2*np.pi))
        output=open(filename2,'w')
        output2=open(filename3,'w')
        sig_fft=np.zeros((9,int(len(y)/2-1)),dtype=complex)
        for i in range(1,5):
                signal=0.0+1j*data[:,i]
                smooth_sig=signal[:]*np.exp(-0.25*width*width*time[:]*time[:])*dt
                sig_fft[i,:]=np.fft.fft(smooth_sig[:],n=nt_of_dw)[:int(len(y)/2-1)]
        for i in range(1,5):
                signal=0.0+1j*data[:,i+4]
                smooth_sig=signal[:]*dt
                sig_fft[i+4,:]=np.fft.fft(smooth_sig[:],n=nt_of_dw)[:int(len(y)/2-1)]
        L11=sig_fft[1,:]/w[:] + 1E-10
        L12=sig_fft[2,:]/w[:]
        L21=sig_fft[3,:]/w[:]
        L22=sig_fft[4,:]/w[:]
        L11_=sig_fft[5,:]
        L12_=sig_fft[6,:]
        L21_=sig_fft[7,:]
        L22_=sig_fft[8,:]
        output.write("1:Freq.\t2:E.cond.\t3:T.cond\t4:T.powr")
        output.write("1:Freq.\t2:L11_R\t3:L12_R\t4:L21_R\t5:L22_R\t6:L11_I\t7:L12_I\t8:L21_I\t9:L22_I")
        output.write("\n")
        for j in range(1,int(len(y)/2-1)):
                output.write("{}\t".format(w[j]))
                output.write("{}\t{}\t{}\t".format(np.real(L11[j]),(np.real(L22[j])-np.real(L21[j])*np.real(L12[j])/np.real(L11[j]))/temp,np.real(L12[j])/np.real(L11[j])/temp))
                output.write("\n")
                output2.write("{}\t".format(w[j]))
                output2.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t".format(np.real(L11[j]),np.real(L12[j]),np.real(L21[j]),np.real(L22[j]),np.imag(L11[j]),np.imag(L12[j]),np.imag(L21[j]),np.imag(L22[j])))
                output2.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t".format(np.real(L11_[j]),np.real(L12_[j]),np.real(L21_[j]),np.real(L22_[j]),np.imag(L11_[j]),np.imag(L12_[j]),np.imag(L21_[j]),np.imag(L22_[j])))
                output2.write("\n")
        output.close
        output2.close
#Begin Program
if __name__ == "__main__":
        width=float(sys.argv[1])
        width=width/27.211
        temp=float(sys.argv[2])
        temp=temp/27.211
        dw=float(sys.argv[3])
        dw=dw/27.211/2.0/np.pi
        filename2 = "stoc_KG.txt"
        filename3 = "KG_Lmn.txt"
        do_fft(filename2, filename3, width, temp, dw)
