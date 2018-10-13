import parser  
import matplotlib.pyplot as plt
import sympy as sym
from sympy.parsing.sympy_parser import parse_expr
sym.init_printing()
y, t = sym.symbols('y t')

def euler(param):	
	y0 = float(param[1])
	t0 = float(param[2])	
	h0 = float(param[3])
	n = int(param[4])
	f = parse_expr(param[5])
	lista_t = []
	lista_y = []
	print('y(',t0,') =',y0)
	print('h =', h0)
	for j in range(n+1):		
		print("%2.d %.7f" % (j, y0))
		lista_y.append(y0)
		lista_t.append(t0)		
		k = float(f.subs(y, y0).subs(t, t0))
		y0 = y0 + h0 * k
		t0 = t0 + h0
	return lista_t, lista_y

def euler_inverso(param):
	y0 = float(param[1])
	t0 = float(param[2])	
	h0 = float(param[3])
	n = int(param[4])
	f = parse_expr(param[5])
	lista_t = []
	lista_y = []
	print('y(',t0,') =',y0)
	print('h =', h0)
	for j in range(n+1):
		print("%2.d %.7f" % (j, y0))
		lista_y.append(y0)
		lista_t.append(t0)
		k = float(f.subs(y, (y0 + float(f.subs(y, y0).subs(t, t0)) * h0)).subs(t, t0+h0))
		y0 = y0 + h0 * k
		t0 = t0 + h0
	return lista_t, lista_y

def euler_aprimorado(param):
	print('Metodo de Euler Aprimorado')
	yk0 = float(param[1])
	tk0 = float(param[2])	
	h0 = float(param[3])
	n = int(param[4])
	f = parse_expr(param[5])
	lista_t = []
	lista_y = []
	print('y(',tk0,') =',yk0)
	print('h =', h0)
	for j in range(n+1):
		
		print("%2.d %.7f" % (j, yk0))
		lista_y.append(yk0)
		lista_t.append(tk0)
		
		t0 = tk0
		y0 = yk0
		k1 = float(f.subs(y, y0).subs(t, t0))
		
		t0 = tk0 + h0 
		y0 = yk0 + h0 * k1
		k2 = float(f.subs(y, y0).subs(t, t0))
		
		yk0 = yk0 + (h0/2) * (k1 + k2)
		tk0 = tk0 + h0
	return lista_t, lista_y

def runge_kutta(param):
	print('Metodo de Runge Kutta')
	yk = float(param[1])
	tk = float(param[2])	
	h = float(param[3])
	n = int(param[4])
	f = parse_expr(param[5])
	lista_t = []
	lista_y = []
	print('y(',tk,') =',yk)
	print('h =', h)
	for j in range(n+1):
		
		print("%2.d %.7f" % (j, yk))
		lista_y.append(yk)
		lista_t.append(tk)
		
		t0 = tk
		y0 = yk
		k1 = float(f.subs(y, y0).subs(t, t0))
			
		t0 = tk + h * 0.5
		y0 = yk + h * k1 * 0.5
		k2 = float(f.subs(y, y0).subs(t, t0))
	
		t0 = tk + 0.5 * h
		y0 = yk + 0.5 * h * k2
		k3 = float(f.subs(y, y0).subs(t, t0))
	
		t0 = tk + h
		y0 = yk + h * k3
		k4 = float(f.subs(y, y0).subs(t, t0))
	
		yk = yk + ((h/6) * (k1 + (2 * k2) + (2 * k3) + k4))
		tk = tk + h
	return lista_t, lista_y

def adam_bashforth(param, t_inic, y_inic, ordem):
	lista_t = t_inic
	lista_y = y_inic
	f = parse_expr(param[ordem + 4])
	h =  float(param[ordem + 2])
	n =  int(param[ordem + 3])
	
	if(ordem == 2):
		b0 = -1/2
		b1 = 3/2	

		for j in range(ordem , n + 1):
			k0 = float(f.subs(y, float(lista_y[j-2])).subs(t, float(lista_t[j-2])))
			k1 = float(f.subs(y, float(lista_y[j-1])).subs(t, float(lista_t[j-1])))
			k2 = float(lista_y[j-1]) + b1*h*k1 + b0*h*k0			
			t2 = float(lista_t[j-1] + h)
			lista_y.append(k2)
			lista_t.append(t2)
			print("%2.d %.7f" % (j, float(lista_y[j])))

	elif (ordem == 3):		
		b0 =  5/12
		b1 = -4/3
		b2 = 23/12	
		 		
		for j in range(ordem , n + 1):
			k0 = float(f.subs(y, float(lista_y[j-3])).subs(t, float(lista_t[j-3])))
			k1 = float(f.subs(y, float(lista_y[j-2])).subs(t, float(lista_t[j-2])))
			k2 = float(f.subs(y, float(lista_y[j-1])).subs(t, float(lista_t[j-1])))
			k3 = float(lista_y[j-1]) + b2*h*k2 + b1*h*k1 + b0*h*k0
			 			
			t2 = float(lista_t[j-1] + h)
			lista_y.append(k3)
			lista_t.append(t2)
			print("%2.d %.7f" % (j, float(lista_y[j])))
	elif (ordem == 4):
		b0 =  -3/8
		b1 = 37/24
		b2 = -59/24	
		b3 = 55/24 		
		for j in range(ordem , n + 1):
			k0 = float(f.subs(y, float(lista_y[j-4])).subs(t, float(lista_t[j-4])))
			k1 = float(f.subs(y, float(lista_y[j-3])).subs(t, float(lista_t[j-3])))
			k2 = float(f.subs(y, float(lista_y[j-2])).subs(t, float(lista_t[j-2])))
			k3 = float(f.subs(y, float(lista_y[j-1])).subs(t, float(lista_t[j-1])))
			k4 = float(lista_y[j-1]) + b3*h*k3 + b2*h*k2 + b1*h*k1 + b0*h*k0

			t2 = float(lista_t[j-1] + h)
			lista_y.append(k4)
			lista_t.append(t2)
			print("%2.d %.7f" % (j, float(lista_y[j])))
	elif (ordem == 5):		
		b0 =  float(251/720)
		b1 = float(-637/360)
		b2 = float(109/30)	
		b3 = float(-1387/360) 		
		b4 = float(1901/720)
		for j in range(ordem , n + 1):
			k0 = float(f.subs(y, float(lista_y[j-5])).subs(t, float(lista_t[j-5])))
			k1 = float(f.subs(y, float(lista_y[j-4])).subs(t, float(lista_t[j-4])))
			k2 = float(f.subs(y, float(lista_y[j-3])).subs(t, float(lista_t[j-3])))
			k3 = float(f.subs(y, float(lista_y[j-2])).subs(t, float(lista_t[j-2])))
			k4 = float(f.subs(y, float(lista_y[j-1])).subs(t, float(lista_t[j-1])))
			k5 = float(lista_y[j-1]) + b4*h*k4 + b3*h*k3 + b2*h*k2 + b1*h*k1 + b0*h*k0

			t2 = float(lista_t[j-1] + h)
			lista_y.append(k5)
			lista_t.append(t2)
			print("%2.d %.7f" % (j, float(lista_y[j])))
	elif (ordem == 6):
		b0 =  -95/288
		b1 = 959/480
		b2 = -3649/720
		b3 = 4991/720
		b4 = -2641/480
		b5 = 4277/1440
		for j in range(ordem , n + 1):
			k0 = float(f.subs(y, float(lista_y[j-6])).subs(t, float(lista_t[j-6])))
			k1 = float(f.subs(y, float(lista_y[j-5])).subs(t, float(lista_t[j-5])))
			k2 = float(f.subs(y, float(lista_y[j-4])).subs(t, float(lista_t[j-4])))
			k3 = float(f.subs(y, float(lista_y[j-3])).subs(t, float(lista_t[j-3])))
			k4 = float(f.subs(y, float(lista_y[j-2])).subs(t, float(lista_t[j-2])))
			k5 = float(f.subs(y, float(lista_y[j-1])).subs(t, float(lista_t[j-1])))
			k6 = float(lista_y[j-1]) + b5*h*k5 + b4*h*k4 + b3*h*k3 + b2*h*k2 + b1*h*k1 + b0*h*k0

			t2 = float(lista_t[j-1] + h)
			lista_y.append(k6)
			lista_t.append(t2)
			print("%2.d %.7f" % (j, float(lista_y[j])))
	elif (ordem == 7):
		b0 =  19087/60480
		b1 = -5603/2520
		b2 = 135713/20160
		b3 = -10754/945
		b4 = 235183/20160
		b5 = -18637/2520
		b6 = 198721/60480
		for j in range(ordem , n + 1):
			k0 = float(f.subs(y, float(lista_y[j-7])).subs(t, float(lista_t[j-7])))
			k1 = float(f.subs(y, float(lista_y[j-6])).subs(t, float(lista_t[j-6])))
			k2 = float(f.subs(y, float(lista_y[j-5])).subs(t, float(lista_t[j-5])))
			k3 = float(f.subs(y, float(lista_y[j-4])).subs(t, float(lista_t[j-4])))
			k4 = float(f.subs(y, float(lista_y[j-3])).subs(t, float(lista_t[j-3])))
			k5 = float(f.subs(y, float(lista_y[j-2])).subs(t, float(lista_t[j-2])))
			k6 = float(f.subs(y, float(lista_y[j-1])).subs(t, float(lista_t[j-1])))
			k7 = float(lista_y[j-1]) + b6*h*k6 + b5*h*k5 + b4*h*k4 + b3*h*k3 + b2*h*k2 + b1*h*k1 + b0*h*k0

			t2 = float(lista_t[j-1] + h)
			lista_y.append(k7)
			lista_t.append(t2)
			print("%2.d %.7f" % (j, float(lista_y[j])))
	elif (ordem == 8):
		b0 =  -5257/17280
		b1 = 32863/13440
		b2 = -115747/13440
		b3 = 2102243/120960
		b4 = -296053/13440
		b5 = 242653/13440
		b6 = -1152169/120960
		b7 = 16083/4480
		for j in range(ordem , n + 1):
			k0 = float(f.subs(y, float(lista_y[j-8])).subs(t, float(lista_t[j-8])))
			k1 = float(f.subs(y, float(lista_y[j-7])).subs(t, float(lista_t[j-7])))
			k2 = float(f.subs(y, float(lista_y[j-6])).subs(t, float(lista_t[j-6])))
			k3 = float(f.subs(y, float(lista_y[j-5])).subs(t, float(lista_t[j-5])))
			k4 = float(f.subs(y, float(lista_y[j-4])).subs(t, float(lista_t[j-4])))
			k5 = float(f.subs(y, float(lista_y[j-3])).subs(t, float(lista_t[j-3])))
			k6 = float(f.subs(y, float(lista_y[j-2])).subs(t, float(lista_t[j-2])))
			k7 = float(f.subs(y, float(lista_y[j-1])).subs(t, float(lista_t[j-1])))
			k8 = float(lista_y[j-1]) + b7*h*k7 + b6*h*k6 + b5*h*k5 + b4*h*k4 + b3*h*k3 + b2*h*k2 + b1*h*k1 + b0*h*k0

			t2 = float(lista_t[j-1] + h)
			lista_y.append(k8)
			lista_t.append(t2)
			print("%2.d %.7f" % (j, float(lista_y[j])))
	return lista_t, lista_y

def adam_multon(param, t_inic, y_inic, ordem):
	lista_t = t_inic
	lista_y = y_inic
	f = parse_expr(param[ordem + 4])
	h =  float(param[ordem + 2])
	n =  int(param[ordem + 3])
	
	if(ordem == 2):
		b0 = -1/2
		b1 = 3/2	

		for j in range(ordem , n + 1):
			k0 = float(f.subs(y, float(lista_y[j-2])).subs(t, float(lista_t[j-2])))
			k1 = float(f.subs(y, float(lista_y[j-1])).subs(t, float(lista_t[j-1])))
			k2 = float(lista_y[j-1]) + b1*h*k1 + b0*h*k0			
			t2 = float(lista_t[j-1] + h)
			lista_y.append(k2)
			lista_t.append(t2)
			print("%2.d %.7f" % (j, float(lista_y[j])))

	elif (ordem == 3):		
		b0 =  5/12
		b1 = -4/3
		b2 = 23/12	
		 		
		for j in range(ordem , n + 1):
			k0 = float(f.subs(y, float(lista_y[j-3])).subs(t, float(lista_t[j-3])))
			k1 = float(f.subs(y, float(lista_y[j-2])).subs(t, float(lista_t[j-2])))
			k2 = float(f.subs(y, float(lista_y[j-1])).subs(t, float(lista_t[j-1])))
			k3 = float(lista_y[j-1]) + b2*h*k2 + b1*h*k1 + b0*h*k0
			 			
			t2 = float(lista_t[j-1] + h)
			lista_y.append(k3)
			lista_t.append(t2)
			print("%2.d %.7f" % (j, float(lista_y[j])))
	elif (ordem == 4):
		b0 =  -3/8
		b1 = 37/24
		b2 = -59/24	
		b3 = 55/24 		
		for j in range(ordem , n + 1):
			k0 = float(f.subs(y, float(lista_y[j-4])).subs(t, float(lista_t[j-4])))
			k1 = float(f.subs(y, float(lista_y[j-3])).subs(t, float(lista_t[j-3])))
			k2 = float(f.subs(y, float(lista_y[j-2])).subs(t, float(lista_t[j-2])))
			k3 = float(f.subs(y, float(lista_y[j-1])).subs(t, float(lista_t[j-1])))
			k4 = float(lista_y[j-1]) + b3*h*k3 + b2*h*k2 + b1*h*k1 + b0*h*k0

			t2 = float(lista_t[j-1] + h)
			lista_y.append(k4)
			lista_t.append(t2)
			print("%2.d %.7f" % (j, float(lista_y[j])))
	elif (ordem == 5):		
		b0 =  float(251/720)
		b1 = float(-637/360)
		b2 = float(109/30)	
		b3 = float(-1387/360) 		
		b4 = float(1901/720)
		for j in range(ordem , n + 1):
			k0 = float(f.subs(y, float(lista_y[j-5])).subs(t, float(lista_t[j-5])))
			k1 = float(f.subs(y, float(lista_y[j-4])).subs(t, float(lista_t[j-4])))
			k2 = float(f.subs(y, float(lista_y[j-3])).subs(t, float(lista_t[j-3])))
			k3 = float(f.subs(y, float(lista_y[j-2])).subs(t, float(lista_t[j-2])))
			k4 = float(f.subs(y, float(lista_y[j-1])).subs(t, float(lista_t[j-1])))
			k5 = float(lista_y[j-1]) + b4*h*k4 + b3*h*k3 + b2*h*k2 + b1*h*k1 + b0*h*k0

			t2 = float(lista_t[j-1] + h)
			lista_y.append(k5)
			lista_t.append(t2)
			print("%2.d %.7f" % (j, float(lista_y[j])))
	elif (ordem == 6):
		b0 =  -95/288
		b1 = 959/480
		b2 = -3649/720
		b3 = 4991/720
		b4 = -2641/480
		b5 = 4277/1440
		for j in range(ordem , n + 1):
			k0 = float(f.subs(y, float(lista_y[j-6])).subs(t, float(lista_t[j-6])))
			k1 = float(f.subs(y, float(lista_y[j-5])).subs(t, float(lista_t[j-5])))
			k2 = float(f.subs(y, float(lista_y[j-4])).subs(t, float(lista_t[j-4])))
			k3 = float(f.subs(y, float(lista_y[j-3])).subs(t, float(lista_t[j-3])))
			k4 = float(f.subs(y, float(lista_y[j-2])).subs(t, float(lista_t[j-2])))
			k5 = float(f.subs(y, float(lista_y[j-1])).subs(t, float(lista_t[j-1])))
			k6 = float(lista_y[j-1]) + b5*h*k5 + b4*h*k4 + b3*h*k3 + b2*h*k2 + b1*h*k1 + b0*h*k0

			t2 = float(lista_t[j-1] + h)
			lista_y.append(k6)
			lista_t.append(t2)
			print("%2.d %.7f" % (j, float(lista_y[j])))
	elif (ordem == 7):
		b0 =  19087/60480
		b1 = -5603/2520
		b2 = 135713/20160
		b3 = -10754/945
		b4 = 235183/20160
		b5 = -18637/2520
		b6 = 198721/60480
		for j in range(ordem , n + 1):
			k0 = float(f.subs(y, float(lista_y[j-7])).subs(t, float(lista_t[j-7])))
			k1 = float(f.subs(y, float(lista_y[j-6])).subs(t, float(lista_t[j-6])))
			k2 = float(f.subs(y, float(lista_y[j-5])).subs(t, float(lista_t[j-5])))
			k3 = float(f.subs(y, float(lista_y[j-4])).subs(t, float(lista_t[j-4])))
			k4 = float(f.subs(y, float(lista_y[j-3])).subs(t, float(lista_t[j-3])))
			k5 = float(f.subs(y, float(lista_y[j-2])).subs(t, float(lista_t[j-2])))
			k6 = float(f.subs(y, float(lista_y[j-1])).subs(t, float(lista_t[j-1])))
			k7 = float(lista_y[j-1]) + b6*h*k6 + b5*h*k5 + b4*h*k4 + b3*h*k3 + b2*h*k2 + b1*h*k1 + b0*h*k0

			t2 = float(lista_t[j-1] + h)
			lista_y.append(k7)
			lista_t.append(t2)
			print("%2.d %.7f" % (j, float(lista_y[j])))
	elif (ordem == 8):
		b0 =  -5257/17280
		b1 = 32863/13440
		b2 = -115747/13440
		b3 = 2102243/120960
		b4 = -296053/13440
		b5 = 242653/13440
		b6 = -1152169/120960
		b7 = 16083/4480
		for j in range(ordem , n + 1):
			k0 = float(f.subs(y, float(lista_y[j-8])).subs(t, float(lista_t[j-8])))
			k1 = float(f.subs(y, float(lista_y[j-7])).subs(t, float(lista_t[j-7])))
			k2 = float(f.subs(y, float(lista_y[j-6])).subs(t, float(lista_t[j-6])))
			k3 = float(f.subs(y, float(lista_y[j-5])).subs(t, float(lista_t[j-5])))
			k4 = float(f.subs(y, float(lista_y[j-4])).subs(t, float(lista_t[j-4])))
			k5 = float(f.subs(y, float(lista_y[j-3])).subs(t, float(lista_t[j-3])))
			k6 = float(f.subs(y, float(lista_y[j-2])).subs(t, float(lista_t[j-2])))
			k7 = float(f.subs(y, float(lista_y[j-1])).subs(t, float(lista_t[j-1])))
			k8 = float(lista_y[j-1]) + b7*h*k7 + b6*h*k6 + b5*h*k5 + b4*h*k4 + b3*h*k3 + b2*h*k2 + b1*h*k1 + b0*h*k0

			t2 = float(lista_t[j-1] + h)
			lista_y.append(k8)
			lista_t.append(t2)
			print("%2.d %.7f" % (j, float(lista_y[j])))
	return lista_t, lista_y


arq = open('input.txt', 'r')
for linha in arq:
	param = linha.split()
	if(param[0]=="euler"):
		print('Metodo de Euler')
		graphic_t, graphic_y  = euler(param)
		plt.plot(graphic_t, graphic_y)
		plt.show()
	elif(param[0]=="euler_aprimorado"):
		print('Metodo de Euler Aprimorado')
		graphic_t, graphic_y  = euler_aprimorado(param)
		plt.plot(graphic_t, graphic_y)
		plt.show()
	elif(param[0]=="euler_inverso"):
		print('Metodo de Euler Inverso')
		graphic_t, graphic_y  = euler_inverso(param)
		plt.plot(graphic_t, graphic_y)
		plt.show()
	elif(param[0]=="runge_kutta"):
		print('Metodo de Runge Kutta')
		graphic_t, graphic_y  = runge_kutta(param)
		plt.plot(graphic_t, graphic_y)
		plt.show()
	elif("adam_bashforth" in str(param[0])):
		ordem = int(param[len(param)-1])
		if("euler_aprimorado" in str(param[0])):
			print("Adam_Bashforth by Euler_Aprimorado")
			param_euler_aprimorado = param[0:6]			
			param_euler_aprimorado[4]= int(param[6])-1
			t_inic,y_inic = euler_aprimorado(param_euler_aprimorado)
			param_aux = []
			param_aux[0:1] = param[0:1]
			param_aux[1:ordem] = y_inic
			param_aux[ordem+1:ordem+4]=param[2:7]
			param = param_aux
		elif("euler_inverso" in str(param[0])):
			print("Adam_Bashforth by Euler_Inverso")
			param_euler_inverso = param[0:6]			
			param_euler_inverso[4]= int(param[6])-1
			t_inic,y_inic = euler_inverso(param_euler_inverso)
			param_aux = []
			param_aux[0:1] = param[0:1]
			param_aux[1:ordem] = y_inic
			param_aux[ordem+1:ordem+4]=param[2:7]
			param = param_aux
		elif("runge_kutta" in str(param[0])):
			print("Adam_Bashforth by Runge_Kutta")
			param_rk = param[0:6]			
			param_rk[4]= int(param[6])-1
			t_inic,y_inic = runge_kutta(param_rk)
			param_aux = []
			param_aux[0:1] = param[0:1]
			param_aux[1:ordem] = y_inic
			param_aux[ordem+1:ordem+4]=param[2:7]
			param = param_aux
		elif("euler" in str(param[0])):
			print("Adam_Bashforth by Euler")
			param_euler = param[0:6]			
			param_euler[4]= int(param[6])-1
			t_inic,y_inic = euler(param_euler)
			param_aux = []
			param_aux[0:1] = param[0:1]
			param_aux[1:ordem] = y_inic
			param_aux[ordem+1:ordem+4]=param[2:7]
			param = param_aux
		else:#adam_bashforth passando os parametros como parametros :)				
			print("Metodo Adan-Bashforth ordem", ordem)
			print('y(',param[ordem+1],') =',param[1])
			print('h =', param[ordem+2])
			t_inic = []
			y_inic = []
			h_aux = 0			
			t0 = float(param[ordem+1])
			i = 0
			#y_inic = param[ 1 :len(param)-5]
			for i in range(ordem): #inicializacao dos parametros				
				y_inic.append(float(param[i+1]))
				t_inic.append(float(t0+h_aux))
				h_aux = float(param[len(param)-4])+h_aux
				print("%2.d %.7f" % (i, y_inic[i]))						
		graphic_t, graphic_y  = adam_bashforth(param, t_inic, y_inic, ordem )
		plt.plot(graphic_t, graphic_y)	
		plt.show()
	elif("adam_multon" in str(param[0])):
		ordem = int(param[len(param)-1])
		if("euler_aprimorado" in str(param[0])):
			print("Adam_Multon by Euler_Aprimorado")
			param_euler_aprimorado = param[0:6]			
			param_euler_aprimorado[4]= int(param[6])-1
			t_inic,y_inic = euler_aprimorado(param_euler_aprimorado)
			param_aux = []
			param_aux[0:1] = param[0:1]
			param_aux[1:ordem] = y_inic
			param_aux[ordem+1:ordem+4]=param[2:7]
			param = param_aux
		elif("euler_inverso" in str(param[0])):
			print("Adam_Multon by Euler_Inverso")
			param_euler_inverso = param[0:6]			
			param_euler_inverso[4]= int(param[6])-1
			t_inic,y_inic = euler_inverso(param_euler_inverso)
			param_aux = []
			param_aux[0:1] = param[0:1]
			param_aux[1:ordem] = y_inic
			param_aux[ordem+1:ordem+4]=param[2:7]
			param = param_aux
		elif("runge_kutta" in str(param[0])):
			print("Adam_Multon by Runge_Kutta")
			param_rk = param[0:6]			
			param_rk[4]= int(param[6])-1
			t_inic,y_inic = runge_kutta(param_rk)
			param_aux = []
			param_aux[0:1] = param[0:1]
			param_aux[1:ordem] = y_inic
			param_aux[ordem+1:ordem+4]=param[2:7]
			param = param_aux
		elif("euler" in str(param[0])):
			print("Adam_Multon by Euler")
			param_euler = param[0:6]			
			param_euler[4]= int(param[6])-1
			t_inic,y_inic = euler(param_euler)
			param_aux = []
			param_aux[0:1] = param[0:1]
			param_aux[1:ordem] = y_inic
			param_aux[ordem+1:ordem+4]=param[2:7]
			param = param_aux
		else:#adam_bashforth passando os parametros como parametros :)				
			print("Metodo Adan-Multon ordem", ordem)
			print('y(',param[ordem+1],') =',param[1])
			print('h =', param[ordem+2])
			t_inic = []
			y_inic = []
			h_aux = 0			
			t0 = float(param[ordem+1])
			i = 0
			#y_inic = param[ 1 :len(param)-5]
			for i in range(ordem): #inicializacao dos parametros				
				y_inic.append(float(param[i+1]))
				t_inic.append(float(t0+h_aux))
				h_aux = float(param[len(param)-4])+h_aux
				print("%2.d %.7f" % (i, y_inic[i]))						
		graphic_t, graphic_y  = adam_multon(param, t_inic, y_inic, ordem )
		plt.plot(graphic_t, graphic_y)	
		plt.show()

	print()
arq.close()
