#!/usr/bin/env python3

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as sig
import math
import tkinter as tk

def dB(x):
  return 10*np.log10(np.abs(x))

def spectra(fc, fs, N, Nb, K2, K3):
  y = []
  for i in range(N):
    t = i / fs
    s = math.sin(math.pi * 2 * fc * t)
    s = s**3 * 10 ** -K3 + s**2 * 10**-K2 + s
    for b in range(Nb):
        y.append(s)

  fr, Gxx = sig.welch(y, fs=fs*Nb, window='blackman', nperseg=N)
  return fr, Gxx

class Application(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title('SFDR')
        self.create_widgets()
        self.ax = self.fig.add_subplot(111)
        self.draw_plot()

    def create_widgets(self):
        '''ウィジェットの配置'''
        self.canvas_frame = tk.Frame(self)
        self.canvas_frame.pack(side=tk.LEFT)
        self.control_frame = tk.Frame(self)
        self.control_frame.pack(side=tk.RIGHT)

        # figureの配置
        self.fig = plt.figure() #figsize=(5, 5))
        self.canvas = FigureCanvasTkAgg(self.fig, self.canvas_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        # スライダーの配置
        self.control_label = tk.Label(self.control_frame, text='fc:')
        self.control_label.pack(anchor=tk.NW)
        self.fc = tk.DoubleVar()
        self.x_scale = tk.Scale(self.control_frame,
            variable=self.fc, from_=1.0, to=16.0, resolution=0.1, orient=tk.HORIZONTAL)
        self.fc.set(8)
        self.x_scale.pack(anchor=tk.NW)
        label = tk.Label(self.control_frame, text='k3:')
        label.pack(anchor=tk.NW)
        self.k3 = tk.DoubleVar()
        k3= tk.Scale(self.control_frame,
            variable=self.k3, from_=1.0, to=8.0, resolution=0.5, orient=tk.HORIZONTAL)
        self.k3.set(2)
        k3.pack(anchor=tk.NW)
        label = tk.Label(self.control_frame, text='k2:')
        label.pack(anchor=tk.NW)
        self.k2 = tk.DoubleVar()
        k2= tk.Scale(self.control_frame,
            variable=self.k2, from_=1.0, to=8.0, resolution=0.5, orient=tk.HORIZONTAL)
        self.k2.set(2)
        k2.pack(anchor=tk.NW)

    def update_anim(self, dt):
        '''グラフ更新関数'''
        fr, gxx = spectra(self.fc.get(), 32, 8192, 8, self.k2.get(), self.k3.get())
        self.ax.clear()
        self.ax.plot(fr, dB(gxx))
        self.ax.set_ylim(-100,20)
        self.ax.set_xlim(0,16)
        self.ax.grid()

    def draw_plot(self, event=None):

        self.ani = FuncAnimation(
                  self.fig,  # Figureオブジェクト
                  self.update_anim,  # グラフ更新関数
                  init_func=None,  # 初期化関数
                  interval = 100,  # 更新間隔(ms)
                  blit = False,
                  )
        self.canvas.draw()


app = Application()
app.mainloop()
