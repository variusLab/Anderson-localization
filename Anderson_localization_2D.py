# Anderson localization of cold atoms in a two-dimensionnal laser speckle disordered potential
# in the presence of Rashba-type spin-orbit interaction

import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft2, ifft2, fftshift
from scipy.integrate import simps
import matplotlib.animation as animation # requiert ffmpeg

# Definition des fonctions

def TF2_transform(f, a):
    return (a**2)*fftshift(fft2(fftshift(f)))

def TF2_inverse(f, a):
    return fftshift(ifft2(fftshift(f)))/(a**2)

def norm2(f):
    return np.real(f*np.conjugate(f))

def Heavside(X, Y, a):
    res0 = X**2 + Y**2
    n = np.size(X[0]);
    res = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            if res0[i,j] <= a**2:
                res[i,j] = 1.
    return res

def Gauss2D(x1, x2, a):
    return np.exp(-(x1*x1 + x2*x2)/(2*a*a))/(2*np.pi*a*a)

# Mise en place d'une grille spatio-temporelle
# et d'une grille dans l'espace des {k}

M = 2**8       # nombre de points dans la grille
dx = 0.2*np.pi # pas
dy = dx
L = M*dx       # longueur de l'intervalle selon x (= celui selon y)
n = np.arange(-M/2, M/2, 1)
x = dx*n
y = x
X,Y = np.meshgrid(x, y, indexing = 'ij')
xmax = x.max()
xmin = x.min()

dk = 2*np.pi/L # pas dans l'espace de Fourier
kx = dk*n      # grille dans l'espace de Fourier
ky = kx
Kx, Ky = np.meshgrid(kx, ky, indexing = 'ij')    
K2 = Kx*Kx + Ky*Ky
K4 = K2*K2

dt = 0.05      # pas temporel (dt<=0.1 : ok)  
tmax = 7.
t = np.arange(0., tmax, dt)
Q = np.size(t)

# Parametres
lambda_R = 3.    # constante de couplage spin-orbite
V0 = 0.2         # intensite du desordre
sigma_Gauss = 5. # largeur du paquet d'onde initial

# Generation d'un potentiel desordonne (champ electrique A(x,y) ss dim)
A = np.random.normal(0, 1., (M,M)) + 1j*np.random.normal(0, 1., (M,M))
k0 = 1.
sigma = 1/k0     # longueur de correlation
tfA = TF2_transform(A, dx) # transformee de Fourier de A

# redefinition de A avec correlation induite
A = TF2_inverse(tfA*Heavside(Kx, Ky, k0), dx)
V = norm2(A)
meanV = np.mean(V)
V = (V/meanV - 1)*V0

# Angle theta et phi
theta = np.zeros((M,M))
for j in range(M):
    for i in range(int(M/2)):
        theta[i,j] = np.arctan(ky[j]/kx[i])
    theta[int(M/2), j] = np.pi*np.sign(ky[j])/2.
    for i in np.arange(int(M/2) + 1, M, 1):
        theta[i,j] = np.arctan(ky[j]/kx[i])
        
phi = theta - np.pi/2
# Deux vecteurs propres de T (en tout point (kx, ky))
lambda1 = K2 + lambda_R*np.sqrt(K2)
lambda2 = K2 - lambda_R*np.sqrt(K2)

# Definition de l'etat initial
psi0 = Gauss2D(X, Y, sigma_Gauss)
Iy0 = simps(psi0**2, x=x, axis = 0) # integration selon l'axe des x
I0 = simps(Iy0, x = y)              # integration selon l'axe des y
psi0 = psi0/np.sqrt(2*I0)           # renormalisation
PsiUp_0 = psi0
PsiDown_0 = psi0

# Calcul de Psi (par la methode split-step FFT symetrique)
# procedure "V/2 - T - V/2", version non optimisee
PsiUp = np.zeros((Q, M, M), 'complex') # suite de Q matrices MxM
PsiDown = np.zeros((Q, M, M), 'complex')

# Initialisation
PsiUp[0, :, :] = PsiUp_0
PsiDown[0, :, :] = PsiDown_0

for q in range(Q-1):
    PsiUp_a = np.exp(-1j*V*dt/2.)*PsiUp[q,:,:]
    PsiDown_a = np.exp(-1j*V*dt/2.)*PsiDown[q,:,:]
    PsiUp_b_0 = TF2_transform(PsiUp_a, dx)
    PsiDown_b_0 = TF2_transform(PsiDown_a, dx)
    alpha_q = (PsiUp_b_0*np.exp(1j*phi) + PsiDown_b_0)/np.sqrt(2.)
    beta_q = (-PsiUp_b_0*np.exp(1j*phi) + PsiDown_b_0)/np.sqrt(2.)
    PsiUp_b = (alpha_q*np.exp(-1j*lambda1*dt)*np.exp(-1j*phi)-
               beta_q*np.exp(-1j*lambda2*dt)*np.exp(-1j*phi))/np.sqrt(2.)
    PsiDown_b = (alpha_q*np.exp(-1j*lambda1*dt)+
               beta_q*np.exp(-1j*lambda2*dt))/np.sqrt(2.)
    PsiUp_c = TF2_inverse(PsiUp_b, dx)
    PsiDown_c = TF2_inverse(PsiDown_b, dx)
    PsiUp[q+1,:,:] = np.exp(-1j*V*dt/2.)*PsiUp_c
    PsiDown[q+1,:,:] = np.exp(-1j*V*dt/2.)*PsiDown_c
    
Psi2 = norm2(PsiUp) + norm2(PsiDown)  # suite de Q matrices MxM : |psi|^2(x,y, t_q)
NormePsi2x = simps(Psi2, x=y, axis=2) # integration pour tout t selon y
NormePsi2 = simps(NormePsi2x, x=x, axis=1) # on obtient <Psi|Psi>(t)
Cmx = simps((X*X + Y*Y)*Psi2, x=y, axis=2)
Cm = simps(Cmx, x=x, axis=1) # correspond a l^2(t):=<psi(r,t)|r^2|psi(r,t)>

# Animation 2D
fig = plt.figure(1)
ims = []
for q in np.linspace(0, Q-1, 400):
    im = plt.imshow(Psi2[int(q),:,:], cmap=plt.get_cmap('jet'), animated=True,
                    aspect='equal', origin='lower', extent=(xmin,xmax,xmin,xmax))
    ims.append([im])
plt.colorbar()
plt.xlabel(r'$x/\sigma$', fontsize=18)
plt.ylabel(r'$y/\sigma$', fontsize=18)
plt.title(r'Space-time evolution of $|\psi|^2$', fontsize=18, pad = 15)
ani = animation.ArtistAnimation(fig, ims, interval=50, blit=True, repeat_delay=1000)
Writer = animation.writers['ffmpeg']
writer = Writer(fps=15, bitrate=1800)
filename = r'E://Python projetcs/Anderson localization/M %s-dt %s-tmax_%s-sigma %s-V0 %s-lambdaR_%s.mp4' %(M, dt, tmax, sigma_Gauss, V0, lambda_R) # e.g. gif, mp4
ani.save(filename, writer=writer)
plt.show()


plt.figure(2)
plt.imshow(Psi2[int(Q/2),:,:], cmap=plt.get_cmap('jet'), animated=True, aspect='equal',
           origin='lower', extent=(xmin,xmax,xmin,xmax))
plt.colorbar()
plt.xlabel(r'$x/\sigma$', fontsize=18)
plt.ylabel(r'$y/\sigma$', fontsize=18)
plt.title(r'Capture of $|\psi|^2$ profile at $t=%sdt$' %(int(Q/2)), fontsize=18, pad = 15)
plt.text(0.05*L/2, -0.97*L/2, r'$\sigma_G=%s \quad V_0=%s \quad \lambda_R=%s$' %(sigma_Gauss, V0, lambda_R), fontsize='small', color='white')
plt.text(-0.95*L/2, 0.91*L/2, r'$M=%s \quad t_{max}=%s \quad dt=%s$' %(M, tmax, dt), fontsize='small', color='white')


plt.figure(3) # plot de l(t) : taille moyenne du paquet d'onde l(t) = sqrt( <psi(r,t)|r^2|psi(r,t)> )
plt.plot(t, np.sqrt(Cm), 'b')
plt.xlabel(r'$t E_{\sigma}/\hbar$', fontsize=18)
plt.ylabel(r'$\ell(t)/\sigma$', fontsize=18)
plt.title(f'Time evolution of the root-mean-square\ndisplacement of $\psi$', fontsize=14, pad = 8)
plt.xlim(0, tmax)

plt.figure(4) # verification de la conservation de la norme <psi|psi>(t) pour tout t
plt.plot(t, NormePsi2)
plt.xlabel(r'$t E_{\sigma}/\hbar$', fontsize=18)
plt.ylabel(r'$<\psi|\psi>$', fontsize=18)
plt.title(r'Time evolution of $|\psi|^2$', fontsize=18, pad = 15)
plt.xlim(0, tmax)
plt.ylim(0, 1.5)
plt.show()


