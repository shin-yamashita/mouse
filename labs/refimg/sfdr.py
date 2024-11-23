#!/usr/bin/env python3

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as sig
import math
import tkinter as tk

def dB(x):
  return 10*np.log10(np.abs(x))

def INL(s, K2, K3, ls):
   return -s**3 * 10 ** -K3 - s**2 * 10**-K2 + s
   
def INLm(s, K2, K3, ls, slice, ms):
   m = ((s + 1.0) % 0.25) * ls / 512
   if slice == int((s + 1.0) / 0.25):
      m += ms / 512 *0.25
   return -s**3 * 10 ** -K3 - s**2 * 10**-K2 + s + m
   
def spectra(fc, fs, N, Nb, K2, K3, ls, amp, slice, msl):
  y = []
  ms = []
  lmsb = 4
  for i in range(N):
    t = i / fs
    
    s = math.sin(math.pi * 2 * fc * t) * 10 ** (amp/20) # -1.0~1.0 FS 
    msb = int((s + 1.0) * 4)  
    ulsb = int(((s + 1.0) * 4 - msb) * 4)
    lsb = s - ((msb - 4) / 4.0 + ulsb / 16.0)
    s = ((msb - 4) / 4.0 + ulsb / 16.0) + lsb
    dmsb = msb - lmsb
    lmsb = msb
    #s = dmsb
    s = s - 0.01 * dmsb / 4
    s = INLm(s, K2, K3, ls, slice, msl)
    for b in range(Nb):
        y.append(s)
        ms.append(msb)

  fr, Gxx = sig.welch(y, fs=fs*Nb, window='flattop', nperseg=N, scaling='spectrum')
  inl = []
  for i in range(4096):
    s = i / 2048 - 1.0
    sd = INLm(s, K2, K3, ls, slice, msl)
    inl.append(sd-s)
  inl = np.array(inl) - ((inl[-1]-inl[0]) * (np.arange(4096)/4096) + inl[0])
  return fr, Gxx, inl, y, ms

class Application(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title('SFDR')
        self.create_widgets()
        #self.ax = self.fig.add_subplots(2,1)
        self.draw_plot()

    def create_widgets(self):
        '''ウィジェットの配置'''
        self.canvas_frame = tk.Frame(self)
        self.canvas_frame.pack(side=tk.LEFT, expand=True, fill=tk.BOTH)
        self.control_frame = tk.Frame(self)
        self.control_frame.pack(side=tk.RIGHT)

        # figureの配置
        self.fig, self.ax = plt.subplots(3,1)
        #self.fig = plt.figure() #figsize=(5, 5))
        self.canvas = FigureCanvasTkAgg(self.fig, self.canvas_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=tk.TOP, expand=True, fill=tk.BOTH)
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.canvas_frame)
        self.toolbar.pack()
        # スライダーの配置
        self.fc = self.slider("fc", 0.9, 1.0, 16.0, 0.05)
        self.k3 = self.slider("k3", 3.5, 1.0, 8.0, 0.5)
        self.k2 = self.slider("k2", 3.5, 1.0, 8.0, 0.5)
        self.ls = self.slider("ls", 0, -5.0, 5.0, 0.5)
        self.amp = self.slider("amp", 0, 0.0, -30.0, 0.5)
        self.slice = self.slider("slice", 0, 0, 7, 1)
        self.ms = self.slider("ms", 0, -5, 5, 0.5)
        self.scan = self.button("scan", self.fc_scan)

    def button(self, name='k', command=None):
        frm = tk.Frame(self)
        label = tk.Label(frm, text=name)
        label.pack(side=tk.LEFT)
        ks= tk.Button(frm, command=command)
        ks.pack(side=tk.RIGHT)
        frm.pack()

    def slider(self, name='k', init=2, from_=1.0, to=8.0, resolution=0.5):
        frm = tk.Frame(self)
        label = tk.Label(frm, text=name)
        label.pack(side=tk.LEFT)
        k = tk.DoubleVar()
        ks= tk.Scale(frm, variable=k, from_=from_, to=to, resolution=resolution, orient=tk.HORIZONTAL, command=self.draw_plot)
        k.set(init)
        ks.pack(side=tk.RIGHT)
        frm.pack()
        return k

    def calc_wfm(self, fc):
        Nb = 8
        fs = 32
        N = 8192
        fr, gxx, inl, y, msb = spectra(fc, fs, N, Nb, 
                    self.k2.get(), self.k3.get(), self.ls.get(), self.amp.get(), self.slice.get(), self.ms.get())
        fr = fr[:512]
        gxx = gxx[:512]
        #ifc = fc/32*8192/8
        #print(fc, ifc)
        peaks, properties = sig.find_peaks(dB(gxx), prominence=10, height=-80)
        pklst = dB(gxx[peaks])
        #print(peaks, pklst)
        pk2 = sorted(pklst, reverse=True)[:3]
        sfdr = pk2[0] - pk2[1]
        ifs = peaks[list(pklst).index(pk2[1])]
        print(f" SFDR : {sfdr:.2f} dBc @ {fr[ifs]} {ifs}  pk2[0]{pk2[0]}")
        f3 = fc * 3
        f3 = f3 if f3 < 16.0 else abs(32.0 - f3)
        f2 = fc * 2
        f2 = f2 if f2 < 16.0 else abs(32.0 - f2)
        if1 = round(fc * 1024 / 32) % 512
        if3 = round(f3 * 1024 / 32) % 512
        if2 = round(f2 * 1024 / 32) % 512
        print(f" HD3 : {fc * 3:.2f} {f3:.2f} {if3} {fr[if3]} : {dB(gxx[if1]) - dB(gxx[if3]):.2f} (dBc)")
        print(f" HD2 : {fc * 2:.2f} {f2:.2f} {if2} {fr[if2]} : {dB(gxx[if1]) - dB(gxx[if2]):.2f} (dBc)")
        msbsl = [[]] * 7
        corr = np.zeros(len(msb))
        for b in range(7):
          msbsl[b] = (np.array(msb) > b)
          corr += ((msbsl[b]) & ~np.roll(msbsl[b], Nb) & ~np.roll(msbsl[b], 2*Nb))
          corr += -1 * (~(msbsl[b]) & np.roll(msbsl[b], Nb) & np.roll(msbsl[b], 2*Nb))
        _, Gxx = sig.welch(corr/8, fs=fs*Nb, window='flattop', nperseg=N, scaling='spectrum')

        return fr, gxx,inl, y, corr, Gxx, peaks, if1, if2, if3, ifs, msbsl

    def fc_scan(self):
      x = []
      hd1 = []
      hd2 = []
      hd3 = []
      sfdr = []
      c3 = []
      for ifc in range(9, 126):
        fc = ifc / 10
        fr, gxx,inl, y, corr, Gxx, peaks, if1, if2, if3, ifs, msbsl  = self.calc_wfm(fc)
        x.append(fr[if1])
        hd1.append(dB(gxx[if1]))
        hd2.append(dB(gxx[if2]))
        hd3.append(dB(gxx[if3]))
        sfdr.append(dB(gxx[ifs]))
        c3.append(dB(Gxx[if3]))
      self.ax[0].clear()
      self.ax[0].plot(x, hd1, 'o-', label='hd1')
      self.ax[0].plot(x, hd2, 'o-', label='hd2', color='red')
      self.ax[0].plot(x, hd3, 'o-', label='hd3', color='orange') 
      self.ax[0].plot(x, sfdr, 'o-', label='sfdr', color='green')
      self.ax[0].plot(x, c3, 'o-', label='c3')
      self.ax[0].grid()
      self.ax[0].set_ylim(-90,0)
      self.ax[0].legend()
      self.fig.tight_layout()
      self.canvas.draw()

    def draw_plot(self, *args):
        '''グラフ更新関数'''
        fc = self.fc.get()

        fr, gxx,inl, y, corr, Gxx, peaks, if1, if2, if3, ifs, msbsl  = self.calc_wfm(fc)

        self.ax[0].clear()
        self.ax[0].plot(fr, dB(gxx))
        self.ax[0].plot(fr[peaks], dB(gxx[peaks]), "x", color="green")
        self.ax[0].plot(fr[ifs], dB(gxx[ifs]), "o", color="yellow")
        self.ax[0].plot(fr[if3], dB(gxx[if3]), "v", color="orange")
        self.ax[0].plot(fr[if2], dB(gxx[if2]), "*", color="red")
        self.ax[0].set_ylim(-100,0)
        self.ax[0].set_xlim(0,16)
        self.ax[0].set_ylabel("spectrum(dB)")
        self.ax[0].set_xlabel("freq(GHz)")
        self.ax[0].grid()
        self.ax[1].clear()
        self.ax[1].plot(np.arange(4096)-2048, np.array(inl)*2048)
        self.ax[1].set_ylim(-10,10)
        self.ax[1].set_xlim(-2048,2048)
        self.ax[1].set_ylabel("INL(LSB)")
        self.ax[1].set_xlabel("code")
        self.ax[1].grid()
        self.ax[2].clear()
        self.ax[2].plot(np.arange(1024), np.array(y[:1024])*2048)

        for b in range(7):
          self.ax[2].plot(np.arange(1024), msbsl[b][:1024]*100 + b*200)

        self.ax[0].plot(fr[:512], dB(Gxx[:512]), "--")
        self.ax[0].plot(fr[if3], dB(Gxx[if3]), "v", color="orange")
        self.ax[0].plot(fr[if2], dB(Gxx[if2]), "*", color="red")        

        self.ax[2].plot(np.arange(1024), corr[:1024]*100 -1000)
        self.ax[2].set_xlim(0,400)
        self.ax[2].grid()
        self.fig.tight_layout()
        self.canvas.draw()




app = Application()
app.mainloop()
