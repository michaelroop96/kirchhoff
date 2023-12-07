#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 11:24:25 2023

@author: michaelroop
"""

import math
import numpy as np
import matplotlib.pyplot as plt

def initial_values_so(NMAT):
    W0 = np.random.rand(NMAT,NMAT)
    W0 = (W0-W0.T)/2
    return W0
def generic(NMAT):
    a = np.random.rand(NMAT,1)
    b = np.random.rand(NMAT,NMAT)
    c = np.random.rand(NMAT,NMAT)
    return a,b,c
def kir(N):
    a = np.random.rand(N)
    a[1] = a[0]
    c11 = np.random.rand()
    c22 = c11
    b11 = np.random.rand()
    b22 = b11
    print(c11)
    b = np.diag([b11,b22,np.random.rand()],k=0)
    c = np.diag([c11,c22,np.random.rand()],k=0)
    return a,b,c
def cleb(N):
    a = np.random.rand(N)
    c22 = np.random.rand()
    c33 = np.random.rand()
    c11 = a[1]*a[2]/(a[1]-a[2])*(c22/a[2]-c33/a[1]-(c22-c33)/a[0])
    b11 = np.random.rand()
    b22 = b11
    b33 = b11
    print(a[0])
    b = np.diag([b11,b22,b33],k=0)
    c = np.diag([c11,c22,c33],k=0)
    return a,b,c
def LSK(N):
    a = np.random.rand(N)
    b22 = np.random.rand()
    b33 = np.random.rand()
    b11 = a[1]*a[2]/(a[1]-a[2])*(b22/a[2]-b33/a[1]-(b22-b33)/a[0])
    c11 = np.random.rand()
    c22 = c11-(b22-b33)**2/a[0]+(b33-b11)**2/a[2]
    c33 = -(b22-b33)**2/a[0]+(b11-b22)**2/a[2]
    b = np.diag([b11,b22,b33],k=0)
    c = np.diag([c11,c22,c33],k=0)
    return a,b,c
def hamiltonian(a,b,c,W,theta):
    N = 3;
    H1 = 0;
    H2 = 0;
    H3 = 0;
    m=[W[2][1],W[0][2],W[1][0]];
    p=[theta[2][1],theta[0][2],theta[1][0]];
    for i in range(0,N):
        H1=H1+a[i]*m[i]**2
    for k in range(0,N):
        for j in range(0,N):
            H2=H2+b[k][j]*(p[k]*m[j]+m[k]*p[j])
            H3=H3+c[k][j]*p[k]*p[j]
    H=1/2*(H1+H2+H3)
    return H
def LieBracket(A,B):
    return A@B-B@A
def inertia(W,theta,a,B,C):
    m=[W[2][1],W[0][2],W[1][0]]
    m = np.array(m)
    p=[theta[2][1],theta[0][2],theta[1][0]]
    p = np.array(p)
    w=np.zeros(3)
    Bp=B@p
    for i in range(0,3):
        w[i]=a[i]*m[i]+Bp[i]
    u=B@m+C@p
    M1 = np.array([[0,-w[2],w[1]],[w[2],0,-w[0]],[-w[1],w[0],0]])
    M2 = np.array([[0,-u[2],u[1]],[u[2],0,-u[0]],[-u[1],u[0],0]])
    return M1,M2
def fixed_point_iteration(W,theta,NMAT,h,a,B,C):
    W0 = np.zeros((NMAT,NMAT))
    theta0 = np.zeros((NMAT,NMAT))
    for i in range(1,15):
        M1,M2=inertia(W0,theta0,a,B,C)
        W1 = W0-(-W+W0-h/2*LieBracket(W0,M1)-h/2*LieBracket(theta0,M2)-h**2/4*(M1@W0@M1+M2@theta0@M1+M1@theta0@M2))
        theta1 = theta0-(-theta+theta0-h/2*LieBracket(theta0,M1)-h**2/4*M1@theta0@M1)
        W0 = W1
        theta0 = theta1
    tildeW = (W1-W1.T)/2
    tildeTheta = (theta1-theta1.T)/2
    return tildeW,tildeTheta
#############
h = 1e-1
NMAT = 3
itermax = 10000
#a,b,c=LSK(3)
a = np.array([0.2955,0.3329,0.4671])
b = np.diag([0.8422,0.8422,0.8422],k=0)
c = np.diag([0.9227,0.6482,0.0252],k=0)
#a = np.array([0.5495,0.4852,0.8905])
#b = np.diag([0.7823,0.7990,0.7343],k=0)
#c = np.diag([0.0513,0.0485,-0.0073],k=0)
B = (b+b.T)/2;
C = (c+c.T)/2;
#W0=initial_values_so(NMAT)
#theta0 = initial_values_so(NMAT)
W0 = np.array([[0,-0.7586,0.5855],[0.7586,0,0.2858],[-0.5855,-0.2858,0]])
theta0 = np.array([[0,-0.7921,-0.3530],[0.7921,0,0.4980],[0.3530,-0.4980,0]])
#W0 = np.array([[0,0.9958,-0.0881],[-0.9958,0,-0.0253],[0.0881,0.0253,0]])
#theta0 = np.array([[0,0.7082,-0.7060],[-0.7082,0,-0.0086],[0.7060,0.0086,0]])
lambda0 = np.sort(np.imag(np.linalg.eig(theta0)[0]))[0]
cas0 = np.real(np.trace(W0@theta0))
H0 = hamiltonian(a,b,c,W0,theta0)
dlambda = np.zeros(itermax)
tphys = np.zeros(itermax)
cas = np.zeros(itermax)
VAR_HAM = np.zeros(itermax)
W = W0
theta = theta0
COMPW = np.zeros((itermax,NMAT,NMAT))
COMPTh = np.zeros((itermax,NMAT,NMAT))
for i in range(0,NMAT):
        for j in range(0,NMAT):
            COMPW[0][i][j] = W[i][j]
            COMPTh[0][i][j] = theta[i][j]
for i in range(0,itermax-1):
    tphys[i]=h*i
    VAR_HAM[i]=hamiltonian(a,b,c,W,theta)-H0
    dlambda[i] = np.sort(np.imag(np.linalg.eig(theta)[0]))[0]-lambda0
    cas[i] = np.real(np.trace(W@theta)-cas0)
    tildeW,tildeTheta = fixed_point_iteration(W,theta,NMAT,h,a,B,C)
    M1,M2=inertia(tildeW,tildeTheta,a,B,C)
    W = W+h*LieBracket(tildeW,M1)+h*LieBracket(tildeTheta,M2)
    theta = theta+h*LieBracket(tildeTheta,M1)
    for k in range(0,NMAT):
        for j in range(0,NMAT):
            COMPW[i+1][k][j] = W[k][j]
            COMPTh[i+1][k][j] = theta[k][j]
tphys[itermax-1]=h*(itermax-1)
VAR_HAM[itermax-1]=hamiltonian(a,b,c,W,theta)-H0;
dlambda[itermax-1] = np.sort(np.imag(np.linalg.eig(theta)[0]))[0]-lambda0
cas[itermax-1] = np.real(np.trace(W@theta)-cas0);
plt.rcParams['font.family'] = 'Times New Roman'
plt.rc('font', size=18) 
csfont = {'fontname':'Times New Roman'}
fig1 = plt.figure(figsize=(8,5))
plt.plot(tphys,dlambda)
plt.xlabel("t",**csfont)
plt.title(r"Spectrum of $\Theta$ variation",**csfont)
plt.savefig("spec_Kir_cleb.eps")
fig2 = plt.figure(figsize=(8,5))
plt.plot(tphys,cas)
plt.xlabel("t")
plt.title("Cross-helicity variation")
plt.savefig("cas_Kir_cleb.eps")
fig3 = plt.figure(figsize=(8,5))
plt.plot(tphys,VAR_HAM)
plt.xlabel("t")
plt.title("Hamiltonian variation")
plt.savefig("ham_Kir_LSK.eps")
W12=np.zeros(itermax)
W23=np.zeros(itermax)
for i in range(0,itermax):
    W12[i]=COMPW[i][0][2]
    W23[i]=COMPW[i][2][1]
fig4 = plt.figure(figsize=(9,6))
plt.plot(W12,W23,linewidth=0.5)
plt.xlabel(r"$W_{13}$")
plt.ylabel(r"$W_{32}$",rotation=0)
plt.title(r"Phase portrait for $W$")
plt.savefig("phaseW_Kir_LSK.eps")
W12=np.zeros(itermax)
W23=np.zeros(itermax)
for i in range(0,itermax):
    W12[i]=COMPTh[i][0][2]
    W23[i]=COMPTh[i][2][1]
fig5 = plt.figure(figsize=(9,6))
plt.plot(W12,W23,linewidth=0.3)
plt.xlabel(r"$\Theta_{13}$")
plt.ylabel(r"$\Theta_{32}$",rotation=0)
plt.title(r"Phase portrait for $\Theta$")
plt.savefig("phaseTh_Kir_LSK.eps")
print(W+W.T)
    
    
    