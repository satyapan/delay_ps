# Delay power spectra plotting task.
# Author: S. Munshi

from nenucal import msutils
import numpy as np
import matplotlib.pyplot as plt
from ps_eor import psutil
from scipy.signal import blackmanharris
from casacore.tables import *
from matplotlib.colors import LogNorm
import matplotlib as mpl
import imageio
import os

mpl.rcParams['image.cmap'] = 'Spectral_r'
mpl.rcParams['image.interpolation'] = 'nearest'
mpl.rcParams['axes.grid'] = True
mpl.rcParams.update({'font.size': 20})

class delay_ps:
    def __init__(self, delay_ps_arr, ms_bu_sort, w_terms_sort, delay_list, timeavg):
        self.ms_bu_sort = ms_bu_sort
        self.horizon_list = self.ms_bu_sort/3.0e8
        self.delay_ps_arr = delay_ps_arr
        self.w_terms_sort = w_terms_sort
        self.delay_list = delay_list
        self.flipped = False
        self.timeavg = timeavg
    def flip_w(self):
        delay_ps_flipped = np.copy(self.delay_ps_arr)
        w_terms_sort_flipped = np.copy(self.w_terms_sort)
        if self.timeavg:
            for i in range(np.shape(self.delay_ps_arr)[1]):
                if self.w_terms_sort[i] < 0:
                    delay_ps_flipped[:,i] = np.flip(delay_ps_flipped[:,i])
        else:
            for i in range(np.shape(self.delay_ps_arr)[1]):
                if self.w_terms_sort[i] < 0:
                    delay_ps_flipped[:,i,:] = np.flip(delay_ps_flipped[:,i,:], axis=0)
        self.delay_ps_arr = delay_ps_flipped
        self.flipped = not self.flipped
        
    def plot(self, ax=None, plot_w_line=True, flip_w=False, fig_name='delay_ps', **kargs):
        if flip_w != self.flipped:
            self.flip_w()          
        if self.timeavg:
            if ax == None:
                fig,ax = plt.subplots(figsize=(12,8))
            X,Y = np.meshgrid(self.ms_bu_sort,self.delay_list)
            im = ax.pcolormesh(X,Y,self.delay_ps_arr, norm=LogNorm(), **kargs)
            ax.plot(self.ms_bu_sort,self.horizon_list,color='black', linewidth=2, linestyle='dotted')
            ax.plot(self.ms_bu_sort,-self.horizon_list,color='black', linewidth=2, linestyle='dotted')
            if plot_w_line:
                if flip_w:
                    ax.plot(self.ms_bu_sort,(self.ms_bu_sort-abs(self.w_terms_sort))/3.0e8,color='black', linewidth=2, alpha=0.6)
                    ax.plot(self.ms_bu_sort,-(self.ms_bu_sort+abs(self.w_terms_sort))/3.0e8,color='black', linewidth=2, alpha=0.6)
                else:
                    ax.plot(self.ms_bu_sort,(self.ms_bu_sort-self.w_terms_sort)/3.0e8,color='black', linewidth=2, alpha=0.6)
                    ax.plot(self.ms_bu_sort,-(self.ms_bu_sort+self.w_terms_sort)/3.0e8,color='black', linewidth=2, alpha=0.6)
            fig = plt.gcf()
            cb = fig.colorbar(im, aspect=50)
            cb.ax.set_ylabel(r'Power ($\mathrm{Jy}^2$)')
            ax.set_xlabel('Baseline (m)')
            ax.set_ylabel('Delay (s)')
            fig.tight_layout()
            if ax == None:
                fig.tight_layout()
                fig.savefig(fig_name+'.png', dpi=100)
        else:
            X,Y = np.meshgrid(self.ms_bu_sort,self.delay_list)
            images = []
            for j in range(np.shape(self.delay_ps_arr)[2]):
                fig, ax = plt.subplots(figsize=(12,8))
                im = ax.pcolormesh(X,Y,self.delay_ps_arr[:,:,j],norm=LogNorm(), **kargs)
                ax.plot(self.ms_bu_sort,self.horizon_list,color='black', linewidth=2, linestyle='dotted')
                ax.plot(self.ms_bu_sort,-self.horizon_list,color='black', linewidth=2, linestyle='dotted')
                if plot_w_line:
                    if flip_w:
                        ax.plot(self.ms_bu_sort,(self.ms_bu_sort-abs(self.w_terms_sort))/3.0e8,color='black', linewidth=2, alpha=0.6)
                        ax.plot(self.ms_bu_sort,-(self.ms_bu_sort+abs(self.w_terms_sort))/3.0e8,color='black', linewidth=2, alpha=0.6)
                    else:
                        ax.plot(self.ms_bu_sort,(self.ms_bu_sort-self.w_terms_sort)/3.0e8,color='black', linewidth=2, alpha=0.6)
                        ax.plot(self.ms_bu_sort,-(self.ms_bu_sort+self.w_terms_sort)/3.0e8,color='black', linewidth=2, alpha=0.6)
                fig.colorbar(im, aspect=50)
                ax.set_xlabel('Baseline (m)')
                ax.set_ylabel('Delay (s)')
                ax.set_title('Segment %s out of %s'%(j+1,np.shape(self.delay_ps_arr)[2]))
                fig.tight_layout()
                fig.savefig(fig_name+str(j)+'.png')
                images.append(imageio.imread(fig_name+str(j)+'.png'))
                os.remove(fig_name+str(j)+'.png')
            imageio.mimsave(fig_name+'.gif', images, fps=1)
            
class delay_ps_gen:
    def __init__(self, ms_file, data_col='SIM_DATA', timeavg=True, n_timeavg=10, stokes='I'):
        self.ms_file = ms_file
        self.data_col = data_col
        self.timeavg = timeavg
        self.n_timeavg = n_timeavg
        self.stokes = stokes
        ms = msutils.MsDataCube.load(self.ms_file, 0, 2000, data_col=self.data_col, n_time_avg=self.n_timeavg)
        self.data = ms.data
        self.ant1 = ms.ant1
        self.ant2 = ms.ant2
        self.shape = ms.data.shape
        self.freq_list = ms.freq
        self.df = self.freq_list[1]-self.freq_list[0]
        self.delay_list = np.fft.fftshift(np.fft.fftfreq(self.shape[0], d=self.df))
        self.ms_bu_ind = np.argsort(ms.bu)
        self.ms_bu_sort = ms.bu[self.ms_bu_ind]
        self.w_terms_sort = self.get_w_terms()
    def get_w_terms(self):
        t = table(self.ms_file, readonly=True)
        ant1 = t.getcol('ANTENNA1')
        ant2 = t.getcol('ANTENNA2')
        uvw = t.getcol('UVW')
        t.close()
        ms_ant1_sort = self.ant1[self.ms_bu_ind]
        ms_ant2_sort = self.ant2[self.ms_bu_ind]
        w_terms_sort = []
        for i in range(self.shape[2]):
            ant1_val = ms_ant1_sort[i]
            ant2_val = ms_ant2_sort[i]
            for j in range(uvw.shape[0]):
                if ant1[j] == ant1_val and ant2[j] == ant2_val:
                    w_terms_sort.append(uvw[j,2])
                    break
        return(np.array(w_terms_sort))
    def get_delay_ps(self):
        if self.stokes == 'I':
            data_stokes = 0.5*(self.data[:,:,:,0]+self.data[:,:,:,3])
        elif self.stokes == 'Q':
            data_stokes = 0.5*(self.data[:,:,:,0]-self.data[:,:,:,3])
        elif self.stokes == 'U':
            data_stokes = 0.5*(self.data[:,:,:,1]+self.data[:,:,:,2])
        elif self.stokes == 'V':
            data_stokes = 0.5*(-1j)*(self.data[:,:,:,1]-self.data[:,:,:,2])
        bh_taper = blackmanharris(self.shape[0])
        if self.timeavg:
            data = np.nanmean(data_stokes, axis=1)
            ps_timeavg = np.empty((self.shape[0],self.shape[2]))
            for i in range(self.shape[2]):
                ps_timeavg[:,i] = abs(psutil.nudft(self.freq_list,data[:,i],w=bh_taper)[1])**2
            delay_ps_arr = ps_timeavg[:,self.ms_bu_ind]
            return delay_ps(delay_ps_arr,self.ms_bu_sort,self.w_terms_sort,self.delay_list,self.timeavg)
        else:
            data = data_stokes
            ps = np.empty((self.shape[0],self.shape[2],self.shape[1]))
            for i in range(self.shape[2]):
                for j in range(self.shape[1]):
                    ps[:,i,j] = abs(psutil.nudft(self.freq_list,data[:,j,i],w=bh_taper)[1])**2
            delay_ps_arr = ps[:,self.ms_bu_ind,:]
            return delay_ps(delay_ps_arr,self.ms_bu_sort,self.w_terms_sort,self.delay_list,self.timeavg)