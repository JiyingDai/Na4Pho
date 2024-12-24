# -*- coding: utf-8 -*-
"""
Created on Tue May  9 16:22:00 2023

@author: jdai
"""
%matplotlib qt5
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

#%% Signal modeling

# Define the signal
n = 20
a = np.zeros(n)

a[0:2] = 1
a[2:16] = 0.02
a[16:18] = 0
a[18:20] = 0.02

spatial_window = np.zeros(n)
spatial_window_size = 1
spatial_window[0:np.int(np.rint(n*(1-spatial_window_size)/2))] = 0
spatial_window[np.int(np.rint(n*(1-spatial_window_size)/2)):np.int(np.rint(n*(1+spatial_window_size)/2))] = 1
spatial_window[np.int(np.rint(n*(1+spatial_window_size)/2)):n] = 0

a_fft = np.fft.fft(a * spatial_window)
freq = np.fft.fftfreq(a.shape[-1])

# Define the readout bandwidth (truncation)
window_size = np.int(n*0.9)
trunc_window = np.zeros(n)

trunc_window[0:np.int(window_size/2)] = 1
trunc_window[np.int(n-window_size/2):n] = 1

a_fft_trunc = a_fft * trunc_window

# Reconstruction
a_rec = np.fft.ifft(a_fft)
a_rec_trunc = np.fft.ifft(a_fft_trunc)


#%% Coil combination modeling

# Coil sensitivity modeling
mu = 0
variance = 1
sigma = np.sqrt(variance)
x_g = np.linspace(mu - 3*sigma, mu + 3*sigma, 2*n+1)
G1 = stats.norm.pdf(x_g, mu, sigma)

B1 = G1[n+1:]/max(G1)
B2 = np.flip(B1)

# Received signal modeling (aliasing not taken into account)
S1 = a * B1
S2 = a * B2

S1_fft = np.fft.fft(S1)
S2_fft = np.fft.fft(S2)

S1_fft_trunc = S1_fft * trunc_window
S2_fft_trunc = S2_fft * trunc_window


S1_rec_proc = np.fft.ifft(S1_fft_trunc)
S2_rec_proc = np.fft.ifft(S2_fft_trunc)

fig1, axs = plt.subplots(1,4,figsize=(7.5, 3.5))
x = np.arange(n)
axs[0].plot(x,B1,linewidth = 3),axs[0].set_title('Coil file 1', fontsize=10), axs[0].set_ylim((-0.1,1.2))
axs[1].plot(x,B2,linewidth = 3),axs[1].set_title('Coil file 2', fontsize=10), axs[1].set_ylim((-0.1,1.2)), axs[1].set_yticks([])
axs[2].plot(S1_rec_proc,linewidth = 3),axs[2].set_title('Reconstructed S1(x)', fontsize=10), axs[2].set_ylim((-0.1,1.2)), axs[2].set_yticks([])
axs[3].plot(S2_rec_proc,linewidth = 3),axs[3].set_title('Reconstructed S2(x)', fontsize=10), axs[3].set_ylim((-0.1,1.2)), axs[3].set_yticks([])

#%%
# Coil combination
S_comb = (S1_rec_proc*B1 + S2_rec_proc*B2)/np.sqrt(np.square(B1)+np.square(B2))
S_comb_self = np.sqrt(np.square(S1_rec_proc)+np.square(S2_rec_proc))
S_comb_other = (S1_rec_proc*B1 + S2_rec_proc*B2)/(np.square(B1)+np.square(B2))
S_comb_sw = (S1_rec_proc*S1_rec_proc + S2_rec_proc*S2_rec_proc)/np.sqrt(np.square(S1_rec_proc)+np.square(S2_rec_proc))

fig2, axs = plt.subplots(1,4,figsize=(8.5, 3.5))
axs[0].plot(a_rec,linewidth = 3), axs[0].set_ylim((-0.1,1.2)),axs[0].set_title('True signal', fontsize=10)
axs[1].plot(S_comb_self,linewidth = 3),axs[1].set_title('Simple square root', fontsize=10), axs[1].set_ylim((-0.1,1.2)), axs[1].set_yticks([])
axs[2].plot(S_comb,linewidth = 3),axs[2].set_title('Roemer equal noise', fontsize=10), axs[2].set_ylim((-0.1,1.2)), axs[2].set_yticks([])
axs[3].plot(S_comb_other,linewidth = 3),axs[3].set_title('Reomer uniform sensitivity', fontsize=10), axs[3].set_ylim((-0.1,1.2)), axs[3].set_yticks([])



