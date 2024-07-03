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
   
def spectra(fc, fs, N, Nb, K2, K3, ls, amp, slice, ms):
  y = []
  for i in range(N):
    t = i / fs
    s = math.sin(math.pi * 2 * fc * t) * 10 ** (amp/20)
    s = INLm(s, K2, K3, ls, slice, ms)
    for b in range(Nb):
        y.append(s)

  fr, Gxx = sig.welch(y, fs=fs*Nb, window='flattop', nperseg=N, scaling='spectrum')
  inl = []
  for i in range(4096):
    s = i / 2048 - 1.0
    sd = INLm(s, K2, K3, ls, slice, ms)
    inl.append(sd-s)
  inl = np.array(inl) - ((inl[-1]-inl[0]) * (np.arange(4096)/4096) + inl[0])
  return fr, Gxx, inl

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
        self.fig, self.ax = plt.subplots(2,1)
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

    def draw_plot(self, *args):
        '''グラフ更新関数'''
        fc = self.fc.get()
        fr, gxx, inl = spectra(fc, 32, 8192, 8, 
                    self.k2.get(), self.k3.get(), self.ls.get(), self.amp.get(), self.slice.get(), self.ms.get())
        fr = fr[:512]
        gxx = gxx[:512]
        #ifc = fc/32*8192/8
        #print(fc, ifc)
        peaks, properties = sig.find_peaks(dB(gxx), prominence=10, height=-90)
        print(peaks, dB(gxx[peaks]))

        self.ax[0].clear()
        self.ax[0].plot(fr, dB(gxx))
        self.ax[0].plot(fr[peaks], dB(gxx[peaks]), "x")
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
        self.fig.tight_layout()
        self.canvas.draw()



app = Application()
app.mainloop()
