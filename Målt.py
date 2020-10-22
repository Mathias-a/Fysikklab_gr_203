
# TFY4104/4107/4115 Fysikk høsten 2020.
#
# Programmet bestemmer høyden til de 8 festepunktene ved å trekke
# tilfeldige heltall mellom 50 og 300 (mm). Når festepunktene er
# akseptable, beregnes baneformen med mm som enhet. 
# Etter dette regnes høydeverdiene om til meter som enhet.
# Hele banens form y(x) beregnes
# ved hjelp av 7 ulike tredjegradspolynomer, på en slik måte
# at både banen y, dens stigningstall dy/dx og dens andrederiverte
# d2y/dx2 er kontinuerlige i alle 6 indre festepunkter.
# I tillegg velges null krumning (andrederivert) 
# i banens to ytterste festepunkter (med bc_type='natural' nedenfor).
# Dette gir i alt 28 ligninger som fastlegger de 28 koeffisientene
# i de i alt 7 tredjegradspolynomene.

# Programmet aksepterer de 8 festepunktene først når
# følgende betingelser er oppfylt:
#
# 1. Starthøyden er stor nok og en del høyere enn i de øvrige 7 festepunktene.
# 2. Helningsvinkelen i startposisjonen er ikke for liten.
# 3. Banens maksimale helningsvinkel er ikke for stor.
#
# Med disse betingelsene oppfylt vil 
# (1) objektet (kula/skiva/ringen) fullføre hele banen selv om det taper noe 
#     mekanisk energi underveis;
# (2) objektet få en fin start, uten å bruke for lang tid i nærheten av
#     startposisjonen; 
# (3) objektet forhåpentlig rulle rent, uten å gli/slure.

# Vi importerer nødvendige biblioteker:
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import CubicSpline
import csv

#FRAM TIL LINJE 108 REGNES TEORETISKE VERDIER UT

# Horisontal avstand mellom festepunktene er 200 mm
h = 200
xfast=np.asarray([0,h,2*h,3*h,4*h,5*h,6*h,7*h])
# Vi begrenser starthøyden (og samtidig den maksimale høyden) til
# å ligge mellom 250 og 300 mm
ymax = 300
# yfast: tabell med 8 heltall mellom 50 og 300 (mm); representerer
# høyden i de 8 festepunktene

#yfast=np.asarray(np.random.randint(50, ymax, size=8))
yfast = np.asarray([252,173,143,75,27,58,65,24])
 
# Omregning fra mm til m:
xfast = xfast/1000
yfast = yfast/1000

#Programmet beregner deretter de 7 tredjegradspolynomene, et
#for hvert intervall mellom to nabofestepunkter.
#Med scipy.interpolate-funksjonen CubicSpline:
cs = CubicSpline(xfast, yfast, bc_type='natural')
xmin = 0.000
xmax = 1.401
dx = 0.001
x = np.arange(xmin, xmax, dx)
Nx = len(x)
y = cs(x)
dy = cs(x,1)
d2y = cs(x,2)
B = np.arctan(dy)/np.pi*180
B_r = np.arctan(dy)
K = d2y/((1+dy**2)**(3/2))
g = 9.81
v_x = np.sqrt(2*g*(yfast[0]-y)/(1+0.4))
a_s = v_x**2*K
N = 0.1*(g*np.cos(B_r)+a_s)
f = 0.1*0.4*g*np.sin(B_r)/(1.4)
v_x_n = v_x*np.cos(B_r)
t_x = [0]
u= 0

for i in range(0,len(x)-1):
    dt = dx/((v_x_n[i]+v_x_n[i+1])/2)
    u += dt
    t_x.append(u)

data = []
filename = "Video_02.csv"
with open(filename) as csvfile:
    csvreader = csv.reader(csvfile)

    for datapoint in csvreader:
        values = [float(value) for value in datapoint]
        data.append(values)


målt_t = ([p[0] for p in data])
#målt_x = ([p[1] for p in data])
#målt_y = ([p[2] for p in data])
målt_v = ([p[1] for p in data])

målt_v_x = (målt_x[len(målt_x)-1]-målt_x[len(målt_x)-2])/(målt_t[len(målt_t)-1]-målt_t[len(målt_t)-2])
målt_v_y = (målt_y[len(målt_y)-1]-målt_y[len(målt_y)-2])/(målt_t[len(målt_t)-1]-målt_t[len(målt_t)-2])
målt_v = np.sqrt(målt_v_x**2+målt_v_y**2)

print('Fart = ', målt_v)





#Plotting

#Bane

baneform = plt.figure('y(x)',figsize=(12,3))
plt.plot(x,y,xfast,yfast,'*')
plt.title('Banens form')
plt.xlabel('$x$ (m)',fontsize=20)
plt.ylabel('$y(x)$ (m)',fontsize=20)
plt.ylim(0,0.350)
plt.grid()
plt.show()
#baneform.savefig("baneform.pdf", bbox_inches='tight')
#baneform.savefig("baneform.png", bbox_inches='tight')


#Vinkel

vinkel = plt.figure('y(x)',figsize=(12,3))
plt.plot(x,B)
plt.title('Vinkel')
plt.xlabel('$x$ (m)',fontsize=20)
plt.ylabel('$B$ (grader)',fontsize=20)
#plt.ylim(0,0.350)
plt.grid()
plt.show()

#krumning

krumning = plt.figure('y(x)',figsize=(12,3))
plt.plot(x,K)
plt.title('Krumning')
plt.xlabel('$x$ (m)',fontsize=20)
plt.ylabel('$K$ (1/m)',fontsize=20)
#plt.ylim(0,0.350)
plt.grid()
plt.show()

#Fart

fart = plt.figure('y(x)',figsize=(12,3))
plt.plot(x,v_x)
plt.title('Fart')
plt.xlabel('$x$ (m)',fontsize=20)
plt.ylabel('$v$ (m/s)',fontsize=20)
#plt.ylim(0,0.350)
plt.grid()
plt.show()

#Normalkraft

Normalkraft = plt.figure('y(x)',figsize=(12,3))
plt.plot(x,N)
plt.title('Normalkraft')
plt.xlabel('$x$ (m)',fontsize=20)
plt.ylabel('N (N)',fontsize=20)
#plt.ylim(0,0.350)
plt.grid()
plt.show()

#friksjon

frikson  = plt.figure('y(x)',figsize=(12,3))
plt.plot(x,np.abs(f/N))
plt.title('Friksjon')
plt.xlabel('$x$ (m)',fontsize=20)
plt.ylabel('|f/N|',fontsize=20)
#plt.ylim(0,0.350)
plt.grid()
plt.show()


#x(t) 

Tid  = plt.figure('y(x)',figsize=(12,3))
plt.plot(t_x,x, label = "Teoretisk")
plt.plot(målt_t,målt_x, label = 'Målt')
plt.legend()
plt.title('pos av tid')
plt.xlabel('tid (s)',fontsize=20)
plt.ylabel('x (m)',fontsize=20)
#plt.ylim(0,0.350)
plt.grid()
plt.show()


#v(t)

Avstand  = plt.figure('y(x)',figsize=(12,3))
plt.plot(t_x,v_x,målt_t,målt_v)
plt.title('Fart av tid')
plt.xlabel('tid (s)',fontsize=20)
plt.ylabel('v (s)',fontsize=20)
#plt.ylim(0,0.350)
plt.grid()
plt.show()

