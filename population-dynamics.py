import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import random as rd
import time
time.time()

wks_s = np.array([1.0, 1.96, 2.425, 3.63])           #% da população de monócitos ao longo da gravidez em relação a não grávida
wks_z = np.array([1.0, 0.81, 0.64, 0.57])            #% da população de linfócitos ao longo da gravidez em relação a não grávida
wks_IFN = np.array([1.0, 0.5963, 0.5611, 1.0704])    #% da produção de linfócitos-T com o mesmo estímulo
a = np.array([0.001, 0.002, 0.003])                  #taxa de invasão bem-sucedida em um monócito susceptível
wks_mu = 80 * wks_s                                  #monócitos produzidos/dia.uL
wks_eta = 0.265 * wks_z                              #linfo-T produzidos/dia.uL para equilíbrio de 2000 linfo-T na ausência de infecção
alpha = 1/3                                          #1/período de vida de um monócito em dias
beta = 1/0.5                                         #1/período de infecção de um monócito
gamma = 0.8                                          #taxa de liberação de vírus
k = 20                                               #taxa de multiplicação de vírus
nu = 0.001                                           #taxa de eliminação de monócito infectado
delta = 1/365                                        #1/período de vida de linfo-T
wks_c = 0.01 * wks_IFN                               #estímulo de produção de linfócitos-T pela densidade de monócitos infectados
wks_d = 0.03 * wks_IFN                               #estímulo de produção de linfócitos-T pelos contatos com monócitos infectados
wks_beta1 = beta + wks_eta * nu / delta
wks_c1 = wks_c + wks_d * wks_eta / delta

h = 0.01                                             #parâmetro de precisão para método numérico de Runge-Kutta
h1 = h/2                                             #parâmetro intermediário de precisão para método numérico de Runge-Kutta
t = np.arange(0, 10+h, h)                            #vetor de tempo, até 10 dias
time = len(t)
pregnancy = len(wks_s)
kutta = 4

s = np.zeros((pregnancy,time,kutta))
v = np.zeros((pregnancy,time,kutta))
i = np.zeros((pregnancy,time,kutta))
z = np.zeros((pregnancy,time,kutta))

s_1 = np.zeros((pregnancy+9,time,kutta))
v_1 = np.zeros((pregnancy+9,time,kutta))
i_1 = np.zeros((pregnancy+9,time,kutta))
z_1 = np.zeros((pregnancy+9,time,kutta))

ds = np.zeros((pregnancy,time - 1,kutta))
dv = np.zeros((pregnancy,time - 1,kutta))
di = np.zeros((pregnancy,time - 1,kutta))
dz = np.zeros((pregnancy,time - 1,kutta))

incremento_s = np.zeros((pregnancy,time))
incremento_i = np.zeros((pregnancy,time))
incremento_v = np.zeros((pregnancy,time))
incremento_z = np.zeros((pregnancy,time))

dados_2 = np.zeros((pregnancy+9,time,kutta))
dados_3 = np.zeros((pregnancy+9,time,kutta))
dados_4 = np.zeros((pregnancy+9,time,kutta))
dados_5 = np.zeros((pregnancy+9,time,kutta))

s[:,0,0] = 250 * wks_s
i[:,0,0] = 10
v[:,0,0] = 165
z[:,0,0] = 2000 * wks_z

dados_x = np.zeros((pregnancy+55,time,kutta))

def graph_generator(eixo_x,eixo_y,a_x,t,dados,nome):
    fig, ax = plt.subplots()
    plt.subplots_adjust(right=0.97, top=0.97)
    plt.xlabel(eixo_x, fontsize = 12, labelpad = 4)
    plt.ylabel(eixo_y,  fontsize = 12, labelpad = 6)
    plt.plot(t,dados[0,:,0],'#0000ff',label = 'Pre-pregnancy',linewidth = '1.8')
    plt.plot(t,dados[1,:,0],'#00df00',label='8 weeks',linewidth = '1.8')
    plt.plot(t,dados[2,:,0],'#ff00ff',label='16 weeks',linewidth = '1.8')
    plt.plot(t,dados[3,:,0],'#ff4500',label='24 weeks',linewidth = '1.8')
    leg = ax.legend()
    plt.grid(color = 'gray', linestyle = '--', linewidth = 0.5)
    name = nome + '[' + str(a_x) + '].png'
    plt.savefig(name,format='png')
    return 0

def graph_generator_2(eixo_x,eixo_y,t,dados,nome):
    fig, ax = plt.subplots()
    plt.subplots_adjust(right=0.97, top=0.97)
    plt.xlabel(eixo_x, fontsize = 12, labelpad = 4)
    plt.ylabel(eixo_y,  fontsize = 12, labelpad = 6)
    plt.plot(t,dados[0,:,0],'#0000ff',label = 'Stage 1',linewidth = '1.8')
    plt.plot(t,dados[1,:,0],'#30C636',label='Stage 2',linewidth = '1.8')
    plt.plot(t,dados[2,:,0],'#cc0000',label='Stage 3',linewidth = '1.8')
    leg = ax.legend()
    plt.grid(color = 'gray', linestyle = '--', linewidth = 0.5)
    name = nome + '.png'
    plt.savefig(name,format='png')
    return 0

def graph_generator_3(eixo_x,eixo_y,a_x,t,dados,nome):
    fig, ax = plt.subplots()
    plt.subplots_adjust(right=0.97, top=0.97)
    plt.xlabel(eixo_x, fontsize = 12, labelpad = 4)
    plt.ylabel(eixo_y,  fontsize = 12, labelpad = 6)
    plt.plot(t,dados[0,:,0],'#0000ff',label = 'Pre-pregnancy',linewidth = '1.8')
    leg = ax.legend()
    plt.grid(color = 'gray', linestyle = '--', linewidth = 0.5)
    name = nome + '[' + str(a_x) + '].png'
    plt.savefig(name,format='png')
    return 0

def stages_graphs(a_x,dados,flag):
    if flag == 'i_m':
        a_x+= 12
    elif flag == 'z_p':
        a_x+= 24
    elif flag == 'l_p':
        a_x+= 36
    dados_x[a_x,:,0] = dados[0,:,0]
    dados_x[a_x+3,:,0] = dados[1,:,0] 
    dados_x[a_x+6,:,0] = dados[2,:,0] 
    dados_x[a_x+9,:,0] = dados[3,:,0] 

    return

for w in range(1):
    print(w)
    for n in range(time - 1):
        for r in range(kutta):
            ds[w,n,r] = wks_mu[w] - alpha * s[w,n,r] - a[2] * s[w,n,r] * v[w,n,r]
            di[w,n,r] = a[2] * s[w,n,r] * v[w,n,r] - beta * i[w,n,r] - nu * i[w,n,r] * z[w,n,r]
            dv[w,n,r] = k * i[w,n,r] - gamma * v[w,n,r] - a[2] * s[w,n,r] * v[w,n,r]
            dz[w,n,r] = wks_eta[w] + wks_c[w] * i[w,n,r] + wks_d[w] * i[w,n,r] * z[w,n,r] - delta * z[w,n,r]

            if r < (kutta - 1):
                s[w,n,r+1] = s[w,n,r] + h1 * ds[w,n,r]
                i[w,n,r+1] = i[w,n,r] + h1 * di[w,n,r]
                v[w,n,r+1] = v[w,n,r] + h1 * dv[w,n,r]
                z[w,n,r+1] = z[w,n,r] + h1 * dz[w,n,r]

        incremento_s[w,n] = (h/6) * (ds[w,n,0] + 2 * ds[w,n,1] + 2 * ds[w,n,2] + ds[w,n,3])
        incremento_i[w,n] = (h/6) * (di[w,n,0] + 2 * di[w,n,1] + 2 * di[w,n,2] + di[w,n,3])
        incremento_v[w,n] = (h/6) * (dv[w,n,0] + 2 * dv[w,n,1] + 2 * dv[w,n,2] + dv[w,n,3])
        incremento_z[w,n] = (h/6) * (dz[w,n,0] + 2 * dz[w,n,1] + 2 * dz[w,n,2] + dz[w,n,3])

        s[w,n+1,0] = s[w,n,0] + incremento_s[w,n]
        i[w,n+1,0] = i[w,n,0] + incremento_i[w,n]
        v[w,n+1,0] = v[w,n,0] + incremento_v[w,n]
        z[w,n+1,0] = z[w,n,0] + incremento_z[w,n]

graph_generator_3('Time (days)',u'Monocytes/\u03bcL',0,t,s,'Susceptible_monocyte_population')
graph_generator_3('Time (days)',u'Monocytes/\u03bcL',0,t,i,'Population_of_infected_monocytes')
graph_generator_3('Time (days)',u'Viral particles/\u03bcL',0,t,v,'Zika_virus_population')
graph_generator_3('Time (days)',u'T-Lymphocytes/\u03bcL',0,t,z,'T-lymphocyte_population')