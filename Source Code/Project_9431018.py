#By: Negin Firouzian
#9431018


import numpy as np
import scipy
from sympy import DiracDelta
import scipy.integrate
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
from mpl_toolkits.axes_grid1 import Size
import matplotlib.pyplot as plt
from PIL import ImageTk
from PIL import ImageTk
from PIL import Image
import tkinter
import _tkinter


def resource_path(relative_path):
    if hasattr(sys, '_MEIPASS'):
        return os.path.join(sys._MEIPASS, relative_path)
    return os.path.join(os.path.abspath("."), relative_path)


global f1
global f2
global s
f1 = lambda t: np.maximum(0, 1-abs(t))
f2 = lambda t: 1*(abs(t-0)<1).astype(float)




fig = plt.figure(figsize=(14, 6))
fig.subplots_adjust(left=0.25, bottom=0.25)


def showConvolution(f1, f2, t0):
    Fs = 50
    T = 5
    t = np.arange(-T, T, 1 / Fs)

    convolution = np.zeros(len(t))
    for n, t_ in enumerate(t):
        prod = lambda tau: f1(tau) * f2(t_ - tau)
        convolution[n] = scipy.integrate.simps(prod(t), t)

    f_shift = lambda t: f2(t0 - t)

    Heavy = lambda t: (t < t0) * 1
    convolution=convolution*Heavy(t)
    prod1 = lambda tau: f1(t) * f_shift(t)


    plt.subplot(311)
    plt.plot(t, f1(t), label=r'$f_1(\tau)$')
    plt.plot(t, f_shift(t), label=r'$f_2(t_0-\tau)$')


    plt.subplot(312)
    plt.plot(t, prod1(t), 'r-', label=r'$f_1(\tau)f_2(t_0-\tau)$')

    plt.subplot(313)
    plt.plot(t, convolution, 'g-', label='$(f_1*f_2)(t)$')


axis_color = 'lightgoldenrodyellow'
T_0 = 0

axT = plt.axes([0.35, 0.1, 0.45, 0.03], facecolor=axis_color)
sT = Slider(axT, 'T', -10, 10.0, valinit=T_0)

showConvolution(f1, f2, 4)

def sliders_on_changed(val):
    plt.subplot(312).cla()
    plt.subplot(311).cla()
    plt.subplot(313).cla()
    showConvolution(f1, f2,sT.val)
    fig.canvas.draw_idle()

sT.on_changed(sliders_on_changed)


#--------------------PULSE---------
PULSE = plt.imread("Rec.gif")
Pulse_button_ax = fig.add_axes([0.01, 0.25, 0.1, 0.1])
Pulse_button = Button(Pulse_button_ax, '', color=axis_color, hovercolor='0.975',image=PULSE)
def Pulse_button_on_clicked(mouse_event):
    global f1
    f1 = lambda t: 1 * (abs(t - 0) < 3).astype(float)
    plt.subplot(312).cla()
    plt.subplot(311).cla()
    plt.subplot(313).cla()
    showConvolution(f1, f2,sT.val)
    fig.canvas.draw_idle()

Pulse_button.on_clicked(Pulse_button_on_clicked)

PULSE2 = plt.imread("Rec2.gif")
Pulse2_button_ax = fig.add_axes([0.12, 0.25, 0.1, 0.1])
Pulse2_button = Button(Pulse2_button_ax, '', color=axis_color, hovercolor='0.975',image=PULSE2)
def Pulse2_button_on_clicked(mouse_event):
    global f2
    f2 = lambda t: 1 * (abs(t - 0) < 1).astype(float)
    plt.subplot(312).cla()
    plt.subplot(311).cla()
    plt.subplot(313).cla()
    showConvolution(f1, f2,sT.val)
    fig.canvas.draw_idle()

Pulse2_button.on_clicked(Pulse2_button_on_clicked)

#--------------------Sin---------
Sin = plt.imread("Sin.gif")
Sin_button_ax = fig.add_axes([0.01, 0.37, 0.1, 0.1])
Sin_button = Button(Sin_button_ax, '', color=axis_color, hovercolor='0.975',image=Sin)
def Sin_button_on_clicked(mouse_event):
    global f1
    f1 = lambda t: (t<4) * (t>-4) * np.sin(-2*t)
    plt.subplot(312).cla()
    plt.subplot(311).cla()
    plt.subplot(313).cla()
    showConvolution(f1, f2,sT.val)
    fig.canvas.draw_idle()

Sin_button.on_clicked(Sin_button_on_clicked)

Sin2 = plt.imread("Sin2.gif")
Sin2_button_ax = fig.add_axes([0.12, 0.37, 0.1, 0.1])
Sin2_button = Button(Sin2_button_ax, '', color=axis_color, hovercolor='0.975',image=Sin2)
def Sin2_button_on_clicked(mouse_event):
    global f2
    f2 = lambda t: (t<4) * (t>-4) * np.sin(-2*t)
    plt.subplot(312).cla()
    plt.subplot(311).cla()
    plt.subplot(313).cla()
    showConvolution(f1, f2,sT.val)
    fig.canvas.draw_idle()

Sin2_button.on_clicked(Sin2_button_on_clicked)

#--------------------Sinc---------
Sinc = plt.imread("Sinc.gif")
Sinc_button_ax = fig.add_axes([0.01, 0.49, 0.1, 0.1])
Sinc_button = Button(Sinc_button_ax, '', color=axis_color, hovercolor='0.975',image=Sinc)
def Sinc_button_on_clicked(mouse_event):
    global f1
    f1 = lambda t: (t<4) * (t>-4) * np.sinc(t)
    plt.subplot(312).cla()
    plt.subplot(311).cla()
    plt.subplot(313).cla()
    showConvolution(f1, f2,sT.val)
    fig.canvas.draw_idle()

Sinc_button.on_clicked(Sinc_button_on_clicked)

Sinc2 = plt.imread("Sinc2.gif")
Sinc2_button_ax = fig.add_axes([0.12, 0.49, 0.1, 0.1])
Sinc2_button = Button(Sinc2_button_ax, '', color=axis_color, hovercolor='0.975',image=Sinc2)
def Sinc2_button_on_clicked(mouse_event):
    global f2
    f2 = lambda t: (t<4) * (t>-4) * np.sinc(t)
    plt.subplot(312).cla()
    plt.subplot(311).cla()
    plt.subplot(313).cla()
    showConvolution(f1, f2,sT.val)
    fig.canvas.draw_idle()

Sinc2_button.on_clicked(Sinc2_button_on_clicked)

#--------------------Two---------
Two= plt.imread("Two.gif")
Two_button_ax = fig.add_axes([0.01, 0.61, 0.1, 0.1])
Two_button = Button(Two_button_ax, '', color=axis_color, hovercolor='0.975',image=Two)
def Two_button_on_clicked(mouse_event):
    global f1
    f1 = lambda t: (t<3) * (t>0) * 1 + (t>-3) * (t<0) * -1
    plt.subplot(312).cla()
    plt.subplot(311).cla()
    plt.subplot(313).cla()
    showConvolution(f1, f2,sT.val)
    fig.canvas.draw_idle()

Two_button.on_clicked(Two_button_on_clicked)

Two2= plt.imread("Two2.gif")
Two2_button_ax = fig.add_axes([0.12, 0.61, 0.1, 0.1])
Two2_button = Button(Two2_button_ax, '', color=axis_color, hovercolor='0.975',image=Two2)
def Two2_button_on_clicked(mouse_event):
    global f2
    f2 = lambda t: (t<3) * (t>0) * 1 + (t>-3) * (t<0) * -1
    plt.subplot(312).cla()
    plt.subplot(311).cla()
    plt.subplot(313).cla()
    showConvolution(f1, f2,sT.val)
    fig.canvas.draw_idle()

Two2_button.on_clicked(Two2_button_on_clicked)

#--------------------Exp---------
Exp= plt.imread("Exp.gif")
Exp_button_ax = fig.add_axes([0.01, 0.73, 0.1, 0.1])
Exp_button = Button(Exp_button_ax, '', color=axis_color, hovercolor='0.975',image=Exp)
def Exp_button_on_clicked(mouse_event):
    global f1
    f1 = lambda t: (t>0) * np.exp(-t)
    plt.subplot(312).cla()
    plt.subplot(311).cla()
    plt.subplot(313).cla()
    showConvolution(f1, f2,sT.val)
    fig.canvas.draw_idle()

Exp_button.on_clicked(Exp_button_on_clicked)

Exp2= plt.imread("Exp2.gif")
Exp2_button_ax = fig.add_axes([0.12, 0.73, 0.1, 0.1])
Exp2_button = Button(Exp2_button_ax, '', color=axis_color, hovercolor='0.975',image=Exp2)
def Exp2_button_on_clicked(mouse_event):
    global f2
    f2 = lambda t: (t>0) * np.exp(-t)
    plt.subplot(312).cla()
    plt.subplot(311).cla()
    plt.subplot(313).cla()
    showConvolution(f1, f2,sT.val)
    fig.canvas.draw_idle()

Exp2_button.on_clicked(Exp2_button_on_clicked)

#--------------------Max---------
Max= plt.imread("Max.gif")
Max_button_ax = fig.add_axes([0.01, 0.85, 0.1, 0.1])
Max_button = Button(Max_button_ax, '', color=axis_color, hovercolor='0.975',image=Max)
def Max_button_on_clicked(mouse_event):
    global f1
    f1 = lambda t: np.maximum(0, 1-abs(t))
    plt.subplot(312).cla()
    plt.subplot(311).cla()
    plt.subplot(313).cla()
    showConvolution(f1, f2,sT.val)
    fig.canvas.draw_idle()

Max_button.on_clicked(Max_button_on_clicked)

Max2= plt.imread("Max2.gif")
Max2_button_ax = fig.add_axes([0.12, 0.85, 0.1, 0.1])
Max2_button = Button(Max2_button_ax, '', color=axis_color, hovercolor='0.975',image=Max2)
def Max2_button_on_clicked(mouse_event):
    global f2
    f2 = lambda t: np.maximum(0, 1-abs(t))
    plt.subplot(312).cla()
    plt.subplot(311).cla()
    plt.subplot(313).cla()
    showConvolution(f1, f2,sT.val)
    fig.canvas.draw_idle()

Max2_button.on_clicked(Max2_button_on_clicked)

plt.show()