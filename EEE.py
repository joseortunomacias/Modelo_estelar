# -*- coding: utf-8 -*-

import numpy as np
from math import *


#*******************
#Ritmo de generación de energía
#*******************
def calc_epsilon(Z, T):
	# ciclo PP
	X1PP = X #fracción de masa del hidrógeno
	X2PP = X #fracción de masa de elementos pesados

	if (T*10.0 <= 4.0):
		epsilon1PP = 0.0
		nuPP = 0.0
	elif(T*10.0 <= 6.0 and T*10.0 > 4.0):
		epsilon1PP = pow(10.0, -6.84)
		nuPP = 6.0
	elif(T*10.0 <= 9.0 and T*10.0 > 6.0):
		epsilon1PP = pow(10.0, -6.04)
		nuPP = 5.0
	elif(T*10.0 <= 11.0 and T*10.0 > 9.0):
		epsilon1PP = pow(10.0, -5.56)
		nuPP = 4.5
	elif(T*10.0 <= 16.0 and T*10.0 > 11.0):
		epsilon1PP = pow(10.0, -5.02)
		nuPP = 4.0
	elif(T*10.0 <= 24.0 and T*10.0 > 16.0):
		epsilon1PP = pow(10.0, -4.40)
		nuPP = 3.5
	elif( T*10.0 >= 24.0):
		epsilon1PP =  pow(10.0, -4.40)
		nuPP = 3.5
		print 'Error: T supera el límite de (T/10) < 24 -->', T*10.0
	epsilonPP = epsilon1PP*X1PP*X2PP*pow(T*10, nuPP)

	#ciclo CNO
	X1CNO = X #fracción de masa del hidrógeno
	X2CNO = 1.0/3.0*Z #fracción de masa de elementos pesados
	if (T*10.0 <= 12.0):
		epsilon1CNO = 0.0
		nuCNO = 0.0
	elif(T*10.0 <= 16.0 and T*10.0 > 12.0):
		epsilon1CNO = pow(10.0, -22.2)
		nuCNO = 20.0
	elif(T*10.0 <= 24.0 and T*10.0 > 16.0):
		epsilon1CNO = pow(10.0, -19.8)
		nuCNO = 18.0
	elif(T*10.0 <= 31.0 and T*10.0 > 24.0):
		epsilon1CNO = pow(10.0, -17.1)
		nuCNO = 16.0
	elif(T*10.0 <= 36.0 and T*10.0 > 31.0):
		epsilon1CNO = pow(10.0, -15.6)
		nuCNO = 15.0
	elif(T*10.0 <= 50.0 and T*10.0 > 36.0):
		epsilon1CNO = pow(10.0, -12.5)
		nuCNO = 13.0
	epsilonCNO = epsilon1CNO*X1CNO*X2CNO*pow(T*10, nuCNO)

	#Nos quedamos con el ciclo que proporcione mayor ritmo de generación energética (epsilon)
	if (epsilonPP >= epsilonCNO):
		epsilon = epsilonPP
		epsilon1 = epsilon1PP
		nu = nuPP
		X1 = X1PP
		X2 = X2PP
		Tipo_ciclo = 'PP'
	else:
		epsilon = epsilonCNO
		epsilon1 = epsilon1CNO
		nu = nuCNO
		X1 = X1CNO
		X2 = X2CNO
		Tipo_ciclo = 'CNO'
	return epsilon, epsilon1, nu, X1, X2, Tipo_ciclo


#*****************
#Declaración de funciones
#*****************
#función1
def func1(h, Z, mu, r, P, T, M, L): #OJO: las variables son vectores de tres componentes
	epsilon, epsilon1, nu, X1, X2, Tipo_ciclo = calc_epsilon(Z, T[2])

	dM0 = 0.01523*mu*pow(r[0],2.0)*P[0]/T[0]
	dM1 = 0.01523*mu*pow(r[1],2.0)*P[1]/T[1]
	dM2 = 0.01523*mu*pow(r[2],2.0)*P[2]/T[2]
	dM =  np.array([dM0,dM1,dM2])

	dP0 = -8.084*mu*M[0]*P[0]/T[0]/pow(r[0],2.0)
	dP1 = -8.084*mu*M[1]*P[1]/T[1]/pow(r[1],2.0)
	dP2 = -8.084*mu*M[2]*P[2]/T[2]/pow(r[2],2.0)
	dP =  np.array([dP0,dP1,dP2])

	dL0 = (0.01845*epsilon1*X1*X2*pow(10.0,nu)*mu**2.0)*((P[0])**2)*pow(T[0],nu-2.0)*((r[0])**2)
	dL1 = (0.01845*epsilon1*X1*X2*pow(10.0,nu)*mu**2.0)*((P[1])**2)*pow(T[1],nu-2.0)*((r[1])**2)
	dL2 = (0.01845*epsilon1*X1*X2*pow(10.0,nu)*mu**2.0)*((P[2])**2)*pow(T[2],nu-2.0)*((r[2])**2)
	dL =  np.array([dL0,dL1,dL2])

	dT0 = -(0.01679*Z*(1.0+X)*(mu**2))*pow(P[0],2.0)*pow(T[0],-8.5)*pow(r[0],-2.0)*L[0]
	dT1 = -(0.01679*Z*(1.0+X)*(mu**2))*pow(P[1],2.0)*pow(T[1],-8.5)*pow(r[1],-2.0)*L[1]
	dT2 = -(0.01679*Z*(1.0+X)*(mu**2))*pow(P[2],2.0)*pow(T[2],-8.5)*pow(r[2],-2.0)*L[2]
	dT =  np.array([dT0,dT1,dT2])

	return dM, dP, dL, dT

#función2
def func2 (h, r, P, T, M, L, dM, dP, dL, dT):

	delta_i_1 = h*dP[2] - h*dP[1]
	delta_i_2 = h*dP[2] - 2*h*dP[1] + h*dP[0]
	P1_est = P[2] + h*dP[2] + 0.5*delta_i_1 + 5.0/12.0*delta_i_2  #P (i+1)

	delta_i_1 = h*dT[2] - h*dT[1]
	delta_i_2 = h*dT[2] - 2*h*dT[1] + h*dT[0]
	T1_est = T[2] + h*dT[2] + 0.5*delta_i_1 + 5.0/12.0*delta_i_2  #T (i+1)
	return P1_est, T1_est

#funcion3
def func3 (h, mu, r, P, T, M, L, dM, dP, dL, dT, P1_est, T1_est):

	dM1 = 0.01523*mu*P1_est/T1_est*pow(r[2]+h, 2.0) # dM (i+1)
	delta_1i_1 = h*dM1 -h*dM[2]
	M1_cal = M[2] + h*dM1 - 0.5*delta_1i_1 
	return dM1, M1_cal

#función 4
def func4(h, mu, r, P, T, M, L,dM, dP, dL, dT,  P1_est, T1_est, dM1, M1_cal):

	dP1 = -8.084*mu*P1_est/T1_est*pow(r[2]+h, -2.0)*M1_cal #dP (i+1)
	delta_1i_1 = h*dP1 -h*dP[2]
	P1_cal = P[2] + h*dP1 - 0.5*delta_1i_1

	return dP1, P1_cal

#funcion 6
def func6(h, Z, mu, r, P, T, M, L, dM, dP, dL, dT, P1_est, T1_est, dP1, P1_cal):

	epsilon, epsilon1, nu, X1, X2, Tipo_ciclo = calc_epsilon(Z, T[2])

	dL1 = 0.01845*epsilon1*X1*X2*pow(10.0, nu)*((mu)**2.0)*((P1_cal)**2.0)*pow(T1_est,nu-2.0)*pow(r[2]+h,2.0)
	delta_1i_1 = h*dL1 - h*dL[2]
	delta_1i_2 = h*dL1 - 2*h*dL[2] + h*dL[1]
	L1_cal = L[2] + h*dL1 - 0.5*delta_1i_1 -1.0/12.0*delta_1i_2

	return dL1, L1_cal, Tipo_ciclo

#funcion 7 
def func7(h, Z, mu, r, P, T, M, L, dM, dP, dL, dT, P1_est, T1_est, dP1, P1_cal, dL1, L1_cal):

	dT1 = -0.01679*Z*(1+X)*((mu)**2.0)*pow(P1_cal,2.0)*pow(T1_est,-8.5)*L1_cal*pow(r[2]+h,-2.0)
	delta_1i_1 = h*dT1 - h*dT[2]
	T1_cal = T[2] + h*dT1 - 0.5*delta_1i_1

	return dT1, T1_cal




def Modelo_interior(X, Y,R_total, L_total, M_total):

	#*******************
	# Valores iniciales
	#*******************
	Z = 1.0-X-Y
	mu = 1.0/(2.0*X+3.0/4.0*Y+0.5*Z)

	R_ini = 0.9*R_total
	h = -R_ini/100.0 # el menos es porque integramos desde la superficie, sería positivo si empezáramos desde el centro

	#tres primeras capas
	r0 = R_ini
	r1 = r0 + h
	r2 = r1 + h


	A1 = 1.9022*mu*M_total
	T0 = A1*(1.0/r0 - 1/R_total)
	T1 = A1*(1.0/r1 - 1/R_total)
	T2 = A1*(1.0/r2 - 1/R_total)



	A2 = 10.645*sqrt(1.0/(mu*Z*(1.0+X))*M_total/L_total)
	P0 = A2*pow(T0,4.25)
	P1 = A2*pow(T1,4.25)
	P2 = A2*pow(T2,4.25)


	File.write('#**************************************************************************************\n')
	File.write('Valores iniciales---->Radio: ')
	File.write(str(R_total))
	File.write('       Luminosidad: ')
	File.write(str(L_total ))
	File.write('\n')
	File.write('#**************************************************************************************\n \n')

	File.write('-')
	File.write( '	Inicio')
	File.write('	')
	File.write( str(0))
	File.write('	')
	File.write( str(r0))
	File.write('	')
	File.write( str(P0)) 
	File.write('	')	
	File.write(str(T0))
	File.write('	')
	File.write( str(L_total))
	File.write('	')
	File.write( str(M_total))
	File.write('	')
	File.write(str( 0.0))
	File.write('\n') 
	File.write('-')
	File.write('	')
	File.write( 'Inicio')
	File.write('	')
	File.write( str(1))
	File.write('	')
	File.write( str(r1))
	File.write('	')
	File.write( str(P1))
	File.write('	')
	File.write(str(T1))
	File.write('	')
	File.write( str(L_total))
	File.write('	')
	File.write( str(M_total))
	File.write('	')
	File.write(str( 0.0))
	File.write('\n') 
	File.write('-')
	File.write('	')
	File.write( 'Inicio')
	File.write('	')
	File.write( str(2))
	File.write('	')
	File.write( str(r2))
	File.write('	')
	File.write( str(P2))
	File.write('	') 
	File.write(str(T2))
	File.write('	')
	File.write( str(L_total))
	File.write('	')
	File.write(str( M_total))
	File.write('	')
	File.write(str( 0.0))
	File.write('\n') 


	#************************
	# FASE radiativa A.1.1
	#************************

	r = np.array([r0, r1, r2])
	P = np.array([P0, P1, P2])
	T = np.array([T0, T1, T2])
	M = np.array([M_total, M_total, M_total])
	L = np.array([L_total, L_total, L_total])

	i = 2 #número de la siguiente fase
	while True:
		dM, dP, dL, dT = func1(h, Z, mu, r, P, T, M, L)
		P1_est, T1_est = func2(h, r, P, T, M, L, dM, dP, dL, dT)

		while True: 
			while True:
				dP1, P1_cal = func4(h, mu, r, P, T, M, L,dM, dP, dL, dT,  P1_est, T1_est, 0.0, M_total)
				if ( abs((P1_cal - P1_est)/P1_cal) < 0.0001 ):
					break
				P1_est = P1_cal

			dT1, T1_cal = func7(h, Z, mu, r, P, T, M, L, dM, dP, dL, dT, P1_est, T1_est, dP1, P1_cal, 0.0, L_total)
			if (abs((T1_cal - T1_est)/T1_cal) < 0.0001):
				break
			T1_est = T1_cal

		dM1, M1_cal = func3 (h, mu, r, P, T, M, L, dM, dP, dL, dT, P1_cal, T1_cal)

		if ( abs((M1_cal - M_total)/M_total) > 0.0001 ): 
			break
		#vuelves a meter valores de r, P, T, M, L para la capa (i+1)
		r[0] = r[1]
		r[1] = r[2]
		r[2] = r[2] + h

		P[0] = P[1]
		P[1] = P[2]
		P[2] = P1_est

		T[0] = T[1]
		T[1] = T[2]
		T[2] = T1_est
		i = i + 1 

		File.write('-')
		File.write('	')
		File.write( 'A.1.1.')
		File.write('	')
		File.write( str(i))
		File.write('	')
		File.write( str(r[2]))
		File.write('	')
		File.write( str(P[2]))
		File.write('	') 
		File.write(str(T[2]))
		File.write('	')
		File.write( str(L_total))
		File.write('	')
		File.write(str( M_total))
		File.write('	')
		File.write(str( 0.0))
		File.write('\n') 


	#************************
	# FASE radiativa A.1.2
	#************************

	while True:
		dM, dP, dL, dT = func1(h, Z, mu, r, P, T, M, L)
		P1_est, T1_est = func2(h, r, P, T, M, L, dM, dP, dL, dT)

		while True: 
			while True:
				dM1, M1_cal = func3(h, mu, r, P, T, M, L, dM, dP, dL, dT, P1_est, T1_est)
				dP1, P1_cal = func4(h, mu, r, P, T, M, L, dM, dP, dL, dT,  P1_est, T1_est, dM1, M1_cal)
				if ( abs((P1_cal - P1_est)/P1_cal) < 0.0001 ):
					break
				P1_est = P1_cal

			dT1, T1_cal = func7(h, Z, mu, r, P, T, M, L, dM, dP, dL, dT, P1_est, T1_est, dP1, P1_cal, 0.0, L_total)
			if (abs((T1_cal - T1_est)/T1_cal) < 0.0001):
				break
			T1_est = T1_cal

		dL1, L1_cal, Tipo_ciclo = func6(h, Z, mu, r, P, T, M, L, dM, dP, dL, dT, P1_cal, T1_cal, dP1, P1_cal)

		if ( abs((L1_cal - L_total)/L1_cal) > 0.0001 ):
			break
		#vuelves a meter valores de r, P, T, M, L para la capa (i+1)
		r[0] = r[1]
		r[1] = r[2]
		r[2] = r[2] + h

		P[0] = P[1]
		P[1] = P[2]
		P[2] = P1_est

		T[0] = T[1]
		T[1] = T[2]
		T[2] = T1_est

		M[0] = M[1]
		M[1] = M[2]
		M[2] = M1_cal
		i = i + 1 

		File.write(Tipo_ciclo)
		File.write('	')
		File.write( 'A.1.2.')
		File.write('	')
		File.write( str(i))
		File.write('	')
		File.write( str(r[2]))
		File.write('	')
		File.write( str(P[2]))
		File.write('	') 
		File.write(str(T[2]))
		File.write('	')
		File.write( str(L_total))
		File.write('	')
		File.write(str( M[2]))
		File.write('	')
		File.write(str( 0.0))
		File.write('\n') 


	#************************
	# FASE radiativa A.1.3
	#************************

	while True:
		dM, dP, dL, dT = func1(h, Z, mu, r, P, T, M, L)
		P1_est, T1_est = func2(h, r, P, T, M, L, dM, dP, dL, dT)

		while True: 
			while True:
				dM1, M1_cal = func3(h, mu, r, P, T, M, L, dM, dP, dL, dT, P1_est, T1_est)
				dP1, P1_cal = func4(h, mu, r, P, T, M, L, dM, dP, dL, dT,  P1_est, T1_est, dM1, M1_cal)
				if ( abs((P1_cal - P1_est)/P1_cal) < 0.0001 ):
					break
				P1_est = P1_cal
			dL1, L1_cal, Tipo_ciclo = func6(h, Z, mu, r, P, T, M, L, dM, dP, dL, dT, P1_est, T1_est, dP1, P1_cal)
			dT1, T1_cal = func7(h, Z, mu, r, P, T, M, L, dM, dP, dL, dT, P1_est, T1_est, dP1, P1_cal, dL1, L1_cal)
			if (abs((T1_cal - T1_est)/T1_cal) < 0.0001):
				break
			T1_est = T1_cal

		n1 = (T1_cal/P1_cal)*(dP1/dT1)# n1 es el coeficiente n+1

		if ( n1 <= 2.5 ):
			r_rad_convec_1 = r[1] + h     #Definimos el radio anterior a la convección
			r_rad_convec_2 = r[2] + h     #Definimos el radio primero de convección
			n1_r_rad_convec_2 = n1        #guardamos el valor de n+1 para la capa ya convectiva
			break
		elif( n1 >= 2.5 ):
			n1_r_rad_convec_1 = n1        #El último valor que se guarde, será el n+1 de la última capa radiativa


		#vuelves a meter valores de r, P, T, M, L para la capa (i+1)
		r[0] = r[1]
		r[1] = r[2]
		r[2] = r[2] + h

		P[0] = P[1]
		P[1] = P[2]
		P[2] = P1_est

		T[0] = T[1]
		T[1] = T[2]
		T[2] = T1_est

		M[0] = M[1]
		M[1] = M[2]
		M[2] = M1_cal

		L[0] = L[1]
		L[1] = L[2]
		L[2] = L1_cal
		i = i + 1 

		File.write(Tipo_ciclo)
		File.write('	')
		File.write( 'A.1.3.')
		File.write('	')
		File.write( str(i))
		File.write('	')
		File.write( str(r[2]))
		File.write('	')
		File.write( str(P[2]))
		File.write('	') 
		File.write(str(T[2]))
		File.write('	')
		File.write( str(L[2]))
		File.write('	')
		File.write(str( M[2]))
		File.write('	')
		File.write(str( n1))
		File.write('\n') 

	#Hallamos el radio en donde se produce el paso de radiativo a convección
	#donde n+1 = 2.50 (interpolamos entre las dos capas del entorno)
	#usamos para ello el valor de n+1 para cada capa
	#Como   (n1_r_rad_convec_1)*alpha + (1-alpha)*(n1_r_rad_convec_2) = 2.5 tenemos:
	alpha = (2.5 - n1_r_rad_convec_2)/(n1_r_rad_convec_1 - n1_r_rad_convec_2)
	r_rad_convec = alpha*r_rad_convec_1 + (1.0-alpha)*r_rad_convec_2 #radio exacto donde se produce el paso de rad a convec





	#************************
	# FASE Convectivo 4.2
	#************************
	K = P1_cal/(pow(T1_cal,2.5))
	contador = 0 # para que el siguiente bucle if sólo se ejecute una vez
	while True:

		epsilon, epsilon1, nu, X1, X2, Tipo_ciclo = calc_epsilon(Z, T[2])                

		dM0 = 0.01523*mu*K*pow(T[0],1.5)*(r[0]**2.0)
		dM1 = 0.01523*mu*K*pow(T[1],1.5)*(r[1]**2.0)
		dM2 = 0.01523*mu*K*pow(T[2],1.5)*(r[2]**2.0)
		dM =  np.array([dM0,dM1,dM2])

		dP0 = -8.084*mu*M[0]*K*pow(T[0],1.5)/pow(r[0],2.0)
		dP1 = -8.084*mu*M[1]*K*pow(T[1],1.5)/pow(r[1],2.0)
		dP2 = -8.084*mu*M[2]*K*pow(T[2],1.5)/pow(r[2],2.0)
		dP =  np.array([dP0,dP1,dP2])

		dL0 = (0.01845*epsilon1*X1*X2*pow(10.0,nu)*(mu**2.0))*(K**2)*pow(T[0],3.0+nu)*((r[0])**2.0)
		dL1 = (0.01845*epsilon1*X1*X2*pow(10.0,nu)*(mu**2.0))*(K**2)*pow(T[1],3.0+nu)*((r[1])**2.0)
		dL2 = (0.01845*epsilon1*X1*X2*pow(10.0,nu)*(mu**2.0))*(K**2)*pow(T[2],3.0+nu)*((r[2])**2.0)
		dL =  np.array([dL0,dL1,dL2])

		dT0 = (-3.234*mu)*M[0]*pow(r[0],-2.0)
		dT1 = (-3.234*mu)*M[1]*pow(r[1],-2.0)
		dT2 = (-3.234*mu)*M[2]*pow(r[2],-2.0)
		dT =  np.array([dT0,dT1,dT2])
		P1_est, T1_est = func2(h, r, P, T, M, L, dM, dP, dL, dT) #P1_est no nos sirve porque es para radiativo

		while True: 

			P1_est = K*pow((abs(T1_est)),2.5)	
			dM1, M1_cal = func3(h, mu, r, P, T, M, L, dM, dP, dL, dT, P1_est, T1_est)
	
			dT1 = -3.234*mu*M1_cal/((r[2]+h)**2.0)
			delta_1i_1 = h*dT1 - h*dT[2]
			T1_cal = T[2] + h*dT1 - 0.5*delta_1i_1

			if (T1_est < 0.0):     # para evitar el bucle el infinito cuando la temperatura sea negativa
				T1_est = 0.0001
				T1_cal = 0.0001

			if ((abs(((T1_cal) - (T1_est))/T1_cal) < 0.0001) and T1_est > 0.0):
				break
			T1_est = (T1_cal)

		P1_cal = K*pow(T1_est,2.5)
		dL1, L1_cal, Tipo_ciclo = func6(h, Z, mu, r, P, T, M, L, dM, dP, dL, dT, P1_est, T1_est, dP1, P1_cal)

		#vuelves a meter valores de r, P, T, M, L para la capa (i+1)
		r[0] = r[1]
		r[1] = r[2]
		r[2] = r[2] + h

		P[0] = P[1]
		P[1] = P[2]
		P[2] = P1_est

		T[0] = T[1]
		T[1] = T[2]
		T[2] = T1_est

		M[0] = M[1]
		M[1] = M[2]
		M[2] = M1_cal

		L[0] = L[1]
		L[1] = L[2]
		L[2] = L1_cal
		i = i + 1 

		n1 = (T1_cal/P1_cal)*(dP1/dT1)# n1 es el coeficiente n+1

		File.write(Tipo_ciclo)
		File.write('	')
		File.write( 'A.2.')
		File.write('	')
		File.write( str(i))
		File.write('	')
		File.write( str(r[2]))
		File.write('	')
		File.write( str(P[2]))
		File.write('	') 
		File.write(str(T[2]))
		File.write('	')
		File.write( str(L[2]))
		File.write('	')
		File.write(str( M[2]))
		File.write('	')
		File.write(str( 0.0))
		File.write('\n') 



		if (contador == 0):
			#interpolamos para hallar las propiedades en r_rad_convec
			r_final1 = alpha*r[1] + (1.0-alpha)*r[2]
			P_final1 = alpha*P[1] + (1.0-alpha)*P[2]
			T_final1 = alpha*T[1] + (1.0-alpha)*T[2]
			M_final1 = alpha*M[1] + (1.0-alpha)*M[2]
			L_final1 = alpha*L[1] + (1.0-alpha)*L[2]
			contador = 1 #para que no vuelva a entrar en el if

		if ( (r[2]) <= 0.15 ): # evita r = 0 (encontraría problemas porque hay cálculos donde r está en el denominador)
			break

	#*************************************************************************************************************************
	#**********************INTEGRACION DESDE EL CENTRO************************************************************************
	#*************************************************************************************************************************

	for j in range (0,40):

		#*******************
		# Valores iniciales
		#*******************
		Z = 1.0-X-Y
		mu = 1.0/(2.0*X+3.0/4.0*Y+0.5*Z)

		h = R_ini/100.0 # el menos es porque integramos desde la superficie, sería positivo si empezáramos desde el centro

		#tres primeras capas
		epsilon, epsilon1, nu, X1, X2, Tipo_ciclo = calc_epsilon(Z, Tc[j])                        # =?)??=)=?=?= CICLO CN??

		r0 = 0.00001 # para evitar problemas en los denominadores
		r1 = 0.0 + h
		r2 = r1 + h

		T0 = Tc[j]
		T1 = Tc[j] -0.008207*(mu**2.0)*K*pow(Tc[j],1.5)*(r1**2.0)
		T2 = Tc[j] -0.008207*(mu**2.0)*K*pow(Tc[j],1.5)*(r2**2.0)

		P0 = K*pow(T0,2.5)
		P1 = K*pow(T1,2.5)
		P2 = K*pow(T2,2.5)

		M0 = 0.005077*mu*K*pow(Tc[j],1.5)*(r0**3.0)
		M1 = 0.005077*mu*K*pow(Tc[j],1.5)*(r1**3.0)
		M2 = 0.005077*mu*K*pow(Tc[j],1.5)*(r2**3.0)

		L0 = 0.00615*epsilon1*X1*X2*pow(10.0,nu)*((mu*K)**2.0)*pow(Tc[j],3.0+nu)*(r0**3.0)
		L1 = 0.00615*epsilon1*X1*X2*pow(10.0,nu)*((mu*K)**2.0)*pow(Tc[j],3.0+nu)*(r1**3.0)
		L2 = 0.00615*epsilon1*X1*X2*pow(10.0,nu)*((mu*K)**2.0)*pow(Tc[j],3.0+nu)*(r2**3.0)

		File.write('#Integración desde centro:\n')
		File.write('#Temperatura:  ')
		File.write(str(Tc[j]))
		File.write('\n')
		File.write(Tipo_ciclo)
		File.write('	')
		File.write( 'A.2.')
		File.write('	')
		File.write( str(0))
		File.write('	')
		File.write( str(r0))
		File.write('	')
		File.write( str(P0))
		File.write('	') 
		File.write(str(T0))
		File.write('	')
		File.write( str(L0))
		File.write('	')
		File.write(str( M0))
		File.write('	')
		File.write(str( 0.0))
		File.write('\n') 

		File.write(Tipo_ciclo)
		File.write('	')
		File.write( 'A.2.')
		File.write('	')
		File.write( str(1))
		File.write('	')
		File.write( str(r1))
		File.write('	')
		File.write( str(P1))
		File.write('	') 
		File.write(str(T1))
		File.write('	')
		File.write( str(L1))
		File.write('	')
		File.write(str( M1))
		File.write('	')
		File.write(str( 0.0))
		File.write('\n')

		File.write(Tipo_ciclo)
		File.write('	')
		File.write( 'A.2.')
		File.write('	')
		File.write( str(2))
		File.write('	')
		File.write( str(r2))
		File.write('	')
		File.write( str(P2))
		File.write('	') 
		File.write(str(T2))
		File.write('	')
		File.write( str(L2))
		File.write('	')
		File.write(str( M2))
		File.write('	')
		File.write(str( 0.0))
		File.write('\n')

		r = np.array([r0, r1, r2])
		P = np.array([P0, P1, P2])
		T = np.array([T0, T1, T2])
		M = np.array([M0, M1, M2])
		L = np.array([L0, L1, L2])
		i = 2
		#************************
		# FASE Convectiva 4.2****
		#************************

		while True:

			epsilon, epsilon1, nu, X1, X2, Tipo_ciclo = calc_epsilon(Z, T[2])                   

			dM0 = 0.01523*mu*K*pow(T[0],1.5)*(r[0]**2.0)
			dM1 = 0.01523*mu*K*pow(T[1],1.5)*(r[1]**2.0)
			dM2 = 0.01523*mu*K*pow(T[2],1.5)*(r[2]**2.0)
			dM =  np.array([dM0,dM1,dM2])

			dP0 = -8.084*mu*M[0]*K*pow(T[0],1.5)/pow(r[0],2.0)
			dP1 = -8.084*mu*M[1]*K*pow(T[1],1.5)/pow(r[1],2.0)
			dP2 = -8.084*mu*M[2]*K*pow(T[2],1.5)/pow(r[2],2.0)
			dP =  np.array([dP0,dP1,dP2])

			dL0 = (0.01845*epsilon1*X1*X2*pow(10.0,nu)*(mu**2.0))*(K**2)*pow(T[0],3.0+nu)*((r[0])**2.0)
			dL1 = (0.01845*epsilon1*X1*X2*pow(10.0,nu)*(mu**2.0))*(K**2)*pow(T[1],3.0+nu)*((r[1])**2.0)
			dL2 = (0.01845*epsilon1*X1*X2*pow(10.0,nu)*(mu**2.0))*(K**2)*pow(T[2],3.0+nu)*((r[2])**2.0)
			dL =  np.array([dL0,dL1,dL2])

			dT0 = (-3.234*mu)*M[0]*pow(r[0],-2.0)
			dT1 = (-3.234*mu)*M[1]*pow(r[1],-2.0)
			dT2 = (-3.234*mu)*M[2]*pow(r[2],-2.0)
			dT =  np.array([dT0,dT1,dT2])
			P1_est, T1_est = func2(h, r, P, T, M, L, dM, dP, dL, dT) #P1_est no nos sirve porque es para radiativo

			while True: 
				#print 'TEAT', T1_est, T[2], dT[2], r[2]	
				P1_est = K*pow(T1_est,2.5)	
				dM1, M1_cal = func3(h, mu, r, P, T, M, L, dM, dP, dL, dT, P1_est, T1_est)

				dT1 = -3.234*mu*M1_cal/((r[2]+h)**2.0)
				delta_1i_1 = h*dT1 - h*dT[2]
				T1_cal = T[2] + h*dT1 - 0.5*delta_1i_1


				if (abs((T1_cal - T1_est)/T1_cal) < 0.0001):
					break
				T1_est = T1_cal

			P1_cal = K*pow(T1_cal,2.5)
			dL1, L1_cal, Tipo_ciclo = func6(h, Z, mu, r, P, T, M, L, dM, dP, dL, dT, P1_est, T1_est, dP1, P1_cal)

			#vuelves a meter valores de r, P, T, M, L para la capa (i+1)
			r[0] = r[1]
			r[1] = r[2]
			r[2] = r[2] + h

			P[0] = P[1]
			P[1] = P[2]
			P[2] = P1_est

			T[0] = T[1]
			T[1] = T[2]
			T[2] = T1_est

			M[0] = M[1]
			M[1] = M[2]
			M[2] = M1_cal

			L[0] = L[1]
			L[1] = L[2]
			L[2] = L1_cal
			i = i + 1 

			File.write(Tipo_ciclo)
			File.write('	')
			File.write( 'A.2.')
			File.write('	')
			File.write( str(i))
			File.write('	')
			File.write( str(r[2]))
			File.write('	')
			File.write( str(P[2]))
			File.write('	') 
			File.write(str(T[2]))
			File.write('	')
			File.write( str(L[2]))
			File.write('	')
			File.write(str( M[2]))
			File.write('	')
			File.write(str( 0.0))
			File.write('\n')

			if ( (r[2]) >= r_rad_convec ):
				r_final2 = alpha*r[2] + (1.0-alpha)*r[1] #cogemos los valores interpolados en la capa que separa radiativo y convectivo
				P_final2 = alpha*P[2] + (1.0-alpha)*P[1]
				T_final2 = alpha*T[2] + (1.0-alpha)*T[1]
				M_final2 = alpha*M[2] + (1.0-alpha)*M[1]
				L_final2 = alpha*L[2] + (1.0-alpha)*L[1] 
				break


		error_r = (r_final2 - r_final1)/r_final1
		error_P = (P_final2 - P_final1)/P_final1
		error_T = (T_final2 - T_final1)/T_final1
		error_M = (M_final2 - M_final1)/M_final1
		error_L = (L_final2 - L_final1)/L_final1
		error_total = sqrt(error_r**2.0 + error_P**2.0 + error_T**2.0 + error_M**2.0 + error_L**2.0)
		Errores[j] = error_total
		minimum_errores = min(Errores)
		Tc_minimum_error = Tc[np.argmin(Errores)]
	return minimum_errores, Tc_minimum_error




#***************************************************************
#*************************MAIN**********************************
#***************************************************************
print '********************************************************'
print '**************Modelo de Interior Estelar****************'
print '********************************************************\n'
print('Introduzca nº de columnas (= nº líneas) de la tabla ')
cl =input ('de Radios y Luminosidades de la estrella \n')#lo llamaremos 'cl'
cl = int(cl)

Radio_ini = input ('Introduzca el valor inicial del radio:\n')*1.0
Luminosidad_ini = input ('Introduzca el valor inicial de la Luminosidad:\n')*1.0
M_total = input ('Introduzca el valor de la Masa total:\n')*1.0
Rintervalo = input ('Introduzca el valor del intervalo en R a estudiar:\n')*1.0
Lintervalo = input ('Introduzca el valor del intervalo en L a estudiar:\n')*1.0
X = input ('Introduzca el valor de la abundancia de H (X):\n')*1.0
Y = input ('Introduzca el valor de la abundancia de He (Y):\n')*1.0


L_total = np.zeros(cl)
R_total = np.zeros(cl)
for i in range(0,cl):
	R_total[i] = (Radio_ini - Rintervalo/2.0) + (Rintervalo/(cl*1.0))*i #Estudiamos radios entre (R_ini - Rintervalo/2.0) y (R_ini + Rintervalo/2.0)
	L_total[i] = (Luminosidad_ini - Lintervalo/2.0) + (Lintervalo/(cl*1.0))*i  #Estudiamos Lumin entre (L_ini - Lintervalo/2.0) y (L_ini + Lintervalo/2.0)

ErroresTotal = np.zeros((cl,cl)) #Donde almacenaré los errores para cada radio y luminosidad
TcTotal = np.zeros((cl,cl)) 	 #Donde almacena las temperaturas con mínimo error para cada combinación de R y L

Tc = np.zeros(40)	#creamos un vector con 40 temperaturas a probar, paara hallar la que proporcione menor error
Errores = np.zeros(40)  #errores para cada temperatura a un R y L dados
for i in range (0,40):
	Tc[i] = 1.8 + (0.4/40.0)*i #Estudiamos temperaturas entre 1.8 y 2.2


File = open('SalidaEEE.txt','w')
for m in range (0,cl): 		#Ejecutamos la función 'Modelo_interior' para cada Luminosidad (m) y para cada Radio (n)
	for n in range(0,cl):
			ErroresTotal[m,n] , TcTotal[m,n] = Modelo_interior(X,Y,R_total[n],L_total[m],M_total)

File.close() # para cerrar el archivo usado por la función 'Modelo_interior' y volverlo a reescribir sólo con los datos de R y L que nos dan error mínimo 


#Mostramos los resultados en la terminal
print 'Las luminosidades y radios estudiados son:'
print L_total
print R_total
print 'La tabla de errores es'
print ErroresTotal

print  'El error relativo es ', np.amin(ErroresTotal)

fila = np.argmin(ErroresTotal)/cl #división entera para que me de el primer dígito
columna = np.argmin(ErroresTotal)%cl #división entera para que me de el segundo dígito que es el resto de dicha división

print 'Tc es ', TcTotal[fila,columna]
print 'R_total es ',  R_total[columna]
print 'L_total es' , L_total[fila]


#Para escribir en un archivo los datos de las capas para la simulación que produce menor error:
File = open('SalidaEEE.txt','w')

File.write('#**************************************************************************************\n')
File.write('#**************************************************************************************\n \n')
File.write('La luminosidades estudiadas son: ')
File.write(str(L_total))
File.write('\nLos radios estudiados son: ')
File.write(str(R_total ))
File.write('\nLos errores para la tabla cuyas columnas son los radios y las líneas las luminosidades:\n')
File.write(str(ErroresTotal))
File.write('\n\nEl error relativo es: ')
File.write(str(np.amin(ErroresTotal)))
File.write ('\n\nTc es ')
File.write(str(TcTotal[fila,columna]))
File.write('\nEl mínimo de error se encuentra para:')
File.write ('\nR_total es \n')
File.write(str(  R_total[columna]))
File.write ('\nL_total es \n' )
File.write(str( L_total[fila]))
File.write('#**************************************************************************************\n\n\n\n')

File.write('#E	fase	i	r	P	T	L	M	n+1\n') 
Modelo_interior(X,Y,R_total[columna],L_total[fila],M_total)


File.close()


