# -*- coding: utf-8 -*-
import pylab as p

ploti=0
epsilon_0=8.85e-12
k=1/(4*p.pi*epsilon_0)
g=0,-9.81

angulo_equilibrio=[]

#criando uma barra
massa_barra=5e-3
comprimento_barra_cm=5e-2
comprimento_barra=100
densidade_barra=comprimento_barra_cm/comprimento_barra
comprimento_suporte=10
barra1=p.empty((comprimento_barra,2))
barra2=p.empty((comprimento_barra,2))

cargas=[]

#carga de 1e-2 a 6e-1
for carga_p in range(1,60):
    carga=carga_p*1e-2
    cargas.append(carga)
    sigma=carga/comprimento_barra
    sigma2=sigma*sigma
    ksigmaden2=k*sigma2*densidade_barra*densidade_barra

    angulos=range(0,90,5)
    torques=[]

    for angulo_normal in angulos:
        #angulo_normal=30
        angulo_barra=(90-angulo_normal)*p.pi/180
        barra2_ini=comprimento_barra*p.cos(angulo_barra)+comprimento_suporte,comprimento_barra*p.sin(angulo_barra)
        barra1_ini=(comprimento_barra*p.cos(angulo_barra),comprimento_barra*p.sin(angulo_barra))
        centro_massa_1=0,0
        centro_massa_2=0,0

        for i in range(comprimento_barra):
            barra1[i][0]=i*p.cos(angulo_barra)
            barra1[i][1]=i*p.sin(angulo_barra)
            barra2[i][0]=barra2_ini[0] + i*p.cos(angulo_barra)
            barra2[i][1]=barra2_ini[1] - i*p.sin(angulo_barra)
            centro_massa_1+=barra1[i]
            centro_massa_2+=barra2[i]

        torque1=0

        for bar1 in barra1:
            forca=0,0
            r=bar1 - barra1_ini
            for bar2 in barra2:
                delta_bar=bar2-bar1
                distancia=p.linalg.norm(delta_bar)
                forca+=ksigmaden2*delta_bar/(distancia**3)
            torque1+=p.cross(r,forca)


        #centro de massa e torque
        barra_transposta1=barra1.transpose()
        barra_transposta2=barra2.transpose()
        centro_massa_1/=comprimento_barra
        centro_massa_2/=comprimento_barra
        #print barra1_ini-centro_massa_1,g,p.cross(barra1_ini-centro_massa_1,g),angulo_barra
        torque1_g=p.cross(barra1_ini-centro_massa_1,g)*massa_barra #sai negativo
        #p.plot(barra_transposta1[0],barra_transposta1[1],lw=0.5)
        #p.plot(centro_massa_1[0],centro_massa_1[1],'r.')
        #print torque1,torque1_g,torque1+torque1_g
        torques.append(torque1+torque1_g)

    torques=abs(p.array(torques))
    angulo_equilibrio.append(angulos[p.where(abs(torques) == min(abs(torques)))[0]])
    print angulo_equilibrio[-1]
    #p.plot(angulos,torques)
p.plot(cargas,angulo_equilibrio)
p.show()

if ploti==1:
    p.plot(barra_transposta1[0],barra_transposta1[1],lw=0.5)
    p.plot(barra_transposta2[0],barra_transposta2[1],lw=0.5)
    p.plot(centro_massa_x1,centro_massa_y1,'r.')
    p.plot(centro_massa_x2,centro_massa_y2,'r.')
    p.show()
