import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal
import cmath

def shift_fourier_spectrum(f_arr, y_arr, freq_shift):
  'función que dado un espectro en frecuencia dado, lo shiftea freq_shift hz'
  new_y = np.zeros(y_arr.size, dtype=np.complex_)
  for index, element in enumerate(y_arr):
    if index+freq_shift>=0 and index+freq_shift< new_y.size:
      new_y[index+freq_shift] = element

  return f_arr, new_y

def create_ideal_fourier(f_arr, y_arr, fs):

  n_max = math.floor((f_arr.size-1)/(2*fs))
  new_y = y_arr
  
  for shift_freq in range(1, n_max+2):

    new_y = new_y + shift_fourier_spectrum(f_arr, y_arr, shift_freq*fs)[1]
    new_y = new_y + shift_fourier_spectrum(f_arr, y_arr, -shift_freq*fs)[1]
  
  return f_arr, new_y

def create_instantaneous_fourier(f_arr, y_arr, tau, fs):
  '''recibe el espectro de una señal y lo transforma
  de modo que quede el espectro de sampleo instantáneo'''

  x_fourier = f_arr
  y_arr_ideal = create_ideal_fourier(f_arr, y_arr, fs)[1]
  vec_exp = np.vectorize(cmath.exp)
  y_insta = tau*fs*np.sinc(tau*f_arr)*y_arr_ideal*vec_exp(-1j*2*np.pi*(1/(2*fs))*(f_arr))
  return x_fourier, y_insta

def create_natural_fourier(f_arr, y_arr, tau, fs):
  '''recibe el espectro de una señal y lo transforma
  de modo que quede el espectro de sampleo instantáneo'''

  x_fourier = f_arr

  vec_exp = np.vectorize(cmath.exp)


  n_max = math.floor((f_arr.size-1)/(2*fs))


  new_y = tau*fs*np.sinc(0*tau*np.pi*fs)*shift_fourier_spectrum(f_arr, y_arr, 0*fs)[1]
  
  for n in range(1, n_max+2):

    new_y = new_y + tau*fs*np.sinc(n*tau*np.pi*fs)*shift_fourier_spectrum(f_arr, y_arr, n*fs)[1]
    new_y = new_y + tau*fs*np.sinc(n*tau*np.pi*fs)*shift_fourier_spectrum(f_arr, y_arr, -n*fs)[1]

  new_y = new_y*vec_exp(-1j*2*np.pi*(1/(2*fs))*(f_arr))

  return x_fourier, new_y

def create_syh_fourier(f_arr, y_arr, fs):
  x, y = create_instantaneous_fourier(f_arr, y_arr, 1/fs, fs)
  return x, y
  
def signal_from_fourier(spectrum, period, period_quantity):
    ''' Antitransformamos un espectro consistende de deltas, period y 
    period_quantity se utilizan para la construccion en el tiempo (cual
    es el periodo esperado y cuantos periodos queremos representar?)
    '''
    cant_muestras = 250*period_quantity
    t = np.linspace(0, period_quantity*period, cant_muestras, endpoint=False)
    y = np.zeros(cant_muestras)

    vec_exp = np.vectorize(cmath.exp)

    span = (spectrum.size-1)/2
    for idx, x in enumerate(spectrum):
      f = idx-span
      if abs(x)>0.0001:   
        y = y + x* vec_exp(2j*np.pi*f*t)
    return t, y

# para crear funciones en el tiempo (más que nada para testear las series de 
# fourier que usamos en las funciones de arriba)

def create_sawtooth(f):
  t = np.linspace(0, 2/f, num = 100)
  y = np.zeros(100)
  for k in range(1, 201):
    y = (-2/np.pi)*((-1)**k)/k*np.sin(k*2*np.pi*f*t) + y
  
  return t, y

def create_sine(f):
    t = np.linspace(0, 2/f, num = 100)
    y = np.sin(2*np.pi*f*t)

    return t,y

def create_cosine(f):
    t = np.linspace(0, 2/f, num = 100)
    y = np.cos(2*np.pi*f*t)

    return t,y

def create_square(f):
  t = np.linspace(0, 2/f, num = 100)
  y = np.zeros(100)
  for i, el in enumerate(y):
    if int(i/25)%2!=0:
      y[i]=1
  return t, y

def graph_instant(x,y,f,percent):
  Ts=1/(100*f)

  new_x = np.arange(x[0], x[-1], Ts)
  new_y = np.interp(new_x, x, y)
  syh_y = np.zeros(new_y.size)

 
  for index, element in enumerate(new_y):
    if index%100==0:
      curr_ind = index
      for n in range(percent):
        i = int(index+50-np.floor(percent/2)+n)
        if i<new_y.size:
            syh_y[i]=new_y[curr_ind]

  return new_x, syh_y

def graph_natural(x,y,f,percent):
  Ts=1/(100*f)

  new_x = np.arange(x[0], x[-1], Ts)
  new_y = np.interp(new_x, x, y)
  syh_y = np.zeros(new_y.size)

 
  for index, element in enumerate(new_y):
    if index%100==0:
      curr_ind = index
      for n in range(percent):
        i = int(index+50-np.floor(percent/2)+n)
        if i<new_y.size:
            syh_y[i]=new_y[i]

  return new_x, syh_y

# en estas funciones "span" se refiere a la frecuencia máxima representada en  
# el espectro de fourier generado.
# ej. span=100 significa que 100hz es la máxima frecuencia representada y que
# los arreglos en frecuencia tendrán 201 elementos (frec ϵ Z c [-100, 100]).

import math

def create_sine_fourier(f, span):
  

  if f!=0:

    f_delta_derecha = (span+f)
    f_delta_izquierda = (span-f)

    if f_delta_derecha<2*span+1 and f_delta_derecha>=0:
      delta_1 = signal.unit_impulse(2*span+1, f_delta_derecha)
    else:
      delta_1 = np.zeros(2*span+1, dtype=np.complex_)
    if f_delta_izquierda<2*span+1 and f_delta_izquierda>=0:
      delta_2 = signal.unit_impulse(2*span+1, f_delta_izquierda)
    else:
      delta_2 = np.zeros(2*span+1, dtype=np.complex_)
  
    delta_vector =  (delta_1 - delta_2)/(2j)

  else:
    delta_vector = signal.unit_impulse(2*span+1, span)

  x_fourier = np.array(range(-span, span+1))

  return x_fourier, delta_vector

def create_cosine_fourier(f, span):
  f_delta_derecha = (span+f)
  f_delta_izquierda = (span-f)

  if f_delta_derecha<2*span+1 and f_delta_derecha>=0:
    delta_1 = signal.unit_impulse(2*span+1, f_delta_derecha)
  else:
    delta_1 = np.zeros(2*span+1, dtype=np.complex_)
  if f_delta_izquierda<2*span+1 and f_delta_izquierda>=0:
    delta_2 = signal.unit_impulse(2*span+1, f_delta_izquierda)
  else:
    delta_2 = np.zeros(2*span+1, dtype=np.complex_)
  
  delta_vector =  (delta_1 + delta_2)/(2)

  x_fourier = np.array(range(-span, span+1))

  return x_fourier, delta_vector

def create_sawtooth_fourier(f, span):

  y_fourier = np.zeros(2*span+1)
  
  m=1
  while m*f<span:
    m+=1

  for k in range(1, m):
    y_fourier = y_fourier + (-2/np.pi)*((-1)**k)/k*create_sine_fourier(k*f, span)[1]
  
  x_fourier = np.array(range(-span, span+1))

  return x_fourier, y_fourier

def create_square_fourier(f, span):

  y_fourier = np.zeros(2*span+1)
  
  m=1
  while m*f<span:
    m+=1

  for k in range(1, m):
    if k%2!=0:
      y_fourier = y_fourier + 4/(np.pi*k) * create_sine_fourier(k*f, span)[1]
  
  x_fourier = np.array(range(-span, span+1))

  return x_fourier, y_fourier

def create_arb_square_fourier(f, duty_cycle, span):

  y_fourier = np.zeros(2*span+1)
  
  m=1
  while m*f<span:
    m+=1

  for k in range(1, m):
      y_fourier = y_fourier + np.sin(2*np.pi*k*duty_cycle)/(k*np.pi) * create_cosine_fourier(k*f, span)[1]
      y_fourier = y_fourier + 2*(np.sin(np.pi*k*duty_cycle)**2)/(k*np.pi) * create_sine_fourier(k*f, span)[1]

  y_fourier = y_fourier + duty_cycle*create_sine_fourier(0, span)[1]
  
  x_fourier = np.array(range(-span, span+1))

  return x_fourier, y_fourier

      # Filtro pasa bajos

def lp_filter_butter(farr_original, h_original, fpass, fstop, dbpass, dbstop):
  ord, wn = signal.buttord(2*np.pi*fpass, 
                           2*np.pi*fstop, 
                           dbpass, 
                           dbstop, 
                           analog=True)
  b, a  = signal.butter(ord, wn, 'low', analog=True)
  warr, h = signal.freqs(b, a)
  farr = warr/(2*np.pi)
  farr = np.concatenate((-np.flip(farr),np.array([-1, 0, 1]),farr))
  h = np.concatenate((np.conjugate(np.flip(h)),np.array([1,1,1]), h))
  filt_interp = np.interp(farr_original, farr, h)

  return farr_original, filt_interp * h_original

def lp_filter_cheby1(farr_original, h_original, fpass, fstop, dbpass, dbstop):
  ord, wn = signal.cheb1ord(2*np.pi*fpass, 
                            2*np.pi*fstop, 
                            dbpass, 
                            dbstop, 
                            analog=True)
  b, a  = signal.cheby1(ord, dbpass, wn, 'low', analog=True)
  warr, h = signal.freqs(b, a)
  farr = warr/(2*np.pi)
  farr = np.concatenate((-np.flip(farr),np.array([-1, 0, 1]),farr))
  h = np.concatenate((np.conjugate(np.flip(h)),np.array([1,1,1]), h))
  filt_interp = np.interp(farr_original, farr, h)

  return farr_original, filt_interp *  h_original
  

def lp_filter_cheby2(farr_original, h_original, fpass, fstop, dbpass, dbstop):
  ord, wn = signal.cheb2ord(2*np.pi*fpass, 
                            2*np.pi*fstop, 
                            dbpass, 
                            dbstop, 
                            analog=True)
  b, a  = signal.cheby2(ord, dbstop, wn, 'low', analog=True)
  warr, h = signal.freqs(b, a)
  farr = warr/(2*np.pi)
  farr = np.concatenate((-np.flip(farr),np.array([-1, 0, 1]),farr))
  h = np.concatenate((np.conjugate(np.flip(h)),np.array([1,1,1]), h))
  filt_interp = np.interp(farr_original, farr, h)

  print('filter order: ', ord)
  return farr_original, filt_interp * h_original