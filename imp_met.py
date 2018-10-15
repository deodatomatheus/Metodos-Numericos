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
	
	wrt.write('y('+ str(t0) + ') = '+ str(y0)+'\n')
	wrt.write('h =' + str(h0) + '\n')
	print('y(',t0,') =',y0)
	print('h =', h0)
	
	for j in range(n+1):		
		wrt.write('{:2d}'.format(j) + ' {:.7f}'.format(y0) + '\n')
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
	wrt.write('y('+ str(t0) + ') = '+ str(y0)+'\n')
	wrt.write('h =' + str(h0) + '\n')
	print('y(',t0,') =',y0)
	print('h =', h0)
	for j in range(n+1):
		print("%2.d %.7f" % (j, y0))
		wrt.write('{:2d}'.format(j) + ' {:.7f}'.format(y0) + '\n')
		lista_y.append(y0)
		lista_t.append(t0)
		k = float(f.subs(y, (y0 + float(f.subs(y, y0).subs(t, t0)) * h0)).subs(t, t0+h0))
		y0 = y0 + h0 * k
		t0 = t0 + h0
	return lista_t, lista_y

def euler_aprimorado(param):
	yk0 = float(param[1])
	tk0 = float(param[2])	
	h0 = float(param[3])
	n = int(param[4])
	f = parse_expr(param[5])
	lista_t = []
	lista_y = []
	wrt.write('y('+ str(tk0) + ') = '+ str(yk0)+'\n')
	wrt.write('h =' + str(h0) + '\n')
	print('y(',tk0,') =',yk0)
	print('h =', h0)
	for j in range(n+1):
		
		print("%2.d %.7f" % (j, yk0))
		wrt.write('{:2d}'.format(j) + ' {:.7f}'.format(yk0) + '\n')
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
	yk = float(param[1])
	tk = float(param[2])	
	h = float(param[3])
	n = int(param[4])
	f = parse_expr(param[5])
	lista_t = []
	lista_y = []
	wrt.write('y('+ str(tk) + ') = '+ str(yk)+'\n')
	wrt.write('h =' + str(h) + '\n')
	print('y(',tk,') =',yk)
	print('h =', h)
	for j in range(n+1):
		wrt.write('{:2d}'.format(j) + ' {:.7f}'.format(yk) + '\n')
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
			wrt.write('{:2d}'.format(j) + ' {:.7f}'.format(float(lista_y[j])) + '\n')
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
			wrt.write('{:2d}'.format(j) + ' {:.7f}'.format(float(lista_y[j])) + '\n')
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
			wrt.write('{:2d}'.format(j) + ' {:.7f}'.format(float(lista_y[j])) + '\n')
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
			wrt.write('{:2d}'.format(j) + ' {:.7f}'.format(float(lista_y[j])) + '\n')
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
			wrt.write('{:2d}'.format(j) + ' {:.7f}'.format(float(lista_y[j])) + '\n')
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
			wrt.write('{:2d}'.format(j) + ' {:.7f}'.format(float(lista_y[j])) + '\n')
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
			wrt.write('{:2d}'.format(j) + ' {:.7f}'.format(float(lista_y[j])) + '\n')
			print("%2.d %.7f" % (j, float(lista_y[j])))
	return lista_t, lista_y

def adam_moulton(param, t_inic, y_inic, ordem):
	lista_t = t_inic
	lista_y = y_inic
	f = parse_expr(param[ordem + 4])
	h =  float(param[ordem + 2])
	n =  int(param[ordem + 3])
	
	if(ordem == 2):
		b0 = 1/2
		b1 = 1/2	

		for j in range(ordem , n + 1):
			k0 = float(f.subs(y, float(lista_y[j-1])).subs(t, float(lista_t[j-1])))
			k1 = float(f.subs(y, (lista_y[j-1] + float(f.subs(y, lista_y[j-1]).subs(t, lista_t[j-1])) * h)).subs(t, lista_t[j-1]+h))
			k2 = float(lista_y[j-1]) + b1*h*k1 + b0*h*k0			
			t2 = float(lista_t[j-1] + h)
			lista_y.append(k2)
			lista_t.append(t2)
			wrt.write('{:2d}'.format(j) + ' {:.7f}'.format(float(lista_y[j])) + '\n')
			print("%2.d %.7f" % (j, float(lista_y[j])))

	elif (ordem == 3):		
		b0 =  -1/12
		b1 = 2/3
		b2 = 5/12	
		 		
		for j in range(ordem , n + 1):
			k0 = float(f.subs(y, float(lista_y[j-2])).subs(t, float(lista_t[j-2])))
			k1 = float(f.subs(y, float(lista_y[j-1])).subs(t, float(lista_t[j-1])))
			k2 = float(f.subs(y, (lista_y[j-1] + float(f.subs(y, lista_y[j-1]).subs(t, lista_t[j-1])) * h)).subs(t, lista_t[j-1]+h))
			k3 = float(lista_y[j-1]) + b2*h*k2 + b1*h*k1 + b0*h*k0
			 			
			t2 = float(lista_t[j-1] + h)
			lista_y.append(k3)
			lista_t.append(t2)
			print("%2.d %.7f" % (j, float(lista_y[j])))
			wrt.write('{:2d}'.format(j) + ' {:.7f}'.format(float(lista_y[j])) + '\n')
	elif (ordem == 4):
		b0 = 1/24
		b1 = -5/24
		b2 = 19/24	
		b3 = 3/8 		
		for j in range(ordem , n + 1):
			k0 = float(f.subs(y, float(lista_y[j-3])).subs(t, float(lista_t[j-3])))
			k1 = float(f.subs(y, float(lista_y[j-2])).subs(t, float(lista_t[j-2])))
			k2 = float(f.subs(y, float(lista_y[j-1])).subs(t, float(lista_t[j-1])))
			k3 = float(f.subs(y, (lista_y[j-1] + float(f.subs(y, lista_y[j-1]).subs(t, lista_t[j-1])) * h)).subs(t, lista_t[j-1]+h))
			k4 = float(lista_y[j-1]) + b3*h*k3 + b2*h*k2 + b1*h*k1 + b0*h*k0

			t2 = float(lista_t[j-1] + h)
			lista_y.append(k4)
			lista_t.append(t2)
			print("%2.d %.7f" % (j, float(lista_y[j])))
			wrt.write('{:2d}'.format(j) + ' {:.7f}'.format(float(lista_y[j])) + '\n')
	elif (ordem == 5):		
		b0 = -19/720
		b1 = 53/360
		b2 = -11/30
		b3 = 323/360 		
		b4 = 251/720
		for j in range(ordem , n + 1):
			k0 = float(f.subs(y, float(lista_y[j-4])).subs(t, float(lista_t[j-4])))
			k1 = float(f.subs(y, float(lista_y[j-3])).subs(t, float(lista_t[j-3])))
			k2 = float(f.subs(y, float(lista_y[j-2])).subs(t, float(lista_t[j-2])))
			k3 = float(f.subs(y, float(lista_y[j-1])).subs(t, float(lista_t[j-1])))
			k4 = float(f.subs(y, (lista_y[j-1] + float(f.subs(y, lista_y[j-1]).subs(t, lista_t[j-1])) * h)).subs(t, lista_t[j-1]+h))
			k5 = float(lista_y[j-1]) + b4*h*k4 + b3*h*k3 + b2*h*k2 + b1*h*k1 + b0*h*k0
			t2 = float(lista_t[j-1] + h)
			lista_y.append(k5)
			lista_t.append(t2)
			print("%2.d %.7f" % (j, float(lista_y[j])))
			wrt.write('{:2d}'.format(j) + ' {:.7f}'.format(float(lista_y[j])) + '\n')
	elif (ordem == 6):
		b0 = 3/160
		b1 = -173/1440
		b2 = 241/720
		b3 = -133/240
		b4 = 1427/1440
		b5 = 95/288
		for j in range(ordem , n + 1):
			k0 = float(f.subs(y, float(lista_y[j-5])).subs(t, float(lista_t[j-5])))
			k1 = float(f.subs(y, float(lista_y[j-4])).subs(t, float(lista_t[j-4])))
			k2 = float(f.subs(y, float(lista_y[j-3])).subs(t, float(lista_t[j-3])))
			k3 = float(f.subs(y, float(lista_y[j-2])).subs(t, float(lista_t[j-2])))
			k4 = float(f.subs(y, float(lista_y[j-1])).subs(t, float(lista_t[j-1])))
			k5 = float(f.subs(y, (lista_y[j-1] + float(f.subs(y, lista_y[j-1]).subs(t, lista_t[j-1])) * h)).subs(t, lista_t[j-1]+h))
			k6 = float(lista_y[j-1]) + b5*h*k5 + b4*h*k4 + b3*h*k3 + b2*h*k2 + b1*h*k1 + b0*h*k0

			t2 = float(lista_t[j-1] + h)
			lista_y.append(k6)
			lista_t.append(t2)
			print("%2.d %.7f" % (j, float(lista_y[j])))
			wrt.write('{:2d}'.format(j) + ' {:.7f}'.format(float(lista_y[j])) + '\n')
	elif (ordem == 7):
		b0 = -863/60480
		b1 = 263/2520
		b2 = -6737/20160
		b3 = 586/945
		b4 = -15487/20160
		b5 = 2713/2520
		b6 = 19087/60480
		for j in range(ordem , n + 1):
			k0 = float(f.subs(y, float(lista_y[j-6])).subs(t, float(lista_t[j-6])))
			k1 = float(f.subs(y, float(lista_y[j-5])).subs(t, float(lista_t[j-5])))
			k2 = float(f.subs(y, float(lista_y[j-4])).subs(t, float(lista_t[j-4])))
			k3 = float(f.subs(y, float(lista_y[j-3])).subs(t, float(lista_t[j-3])))
			k4 = float(f.subs(y, float(lista_y[j-2])).subs(t, float(lista_t[j-2])))
			k5 = float(f.subs(y, float(lista_y[j-1])).subs(t, float(lista_t[j-1])))
			k6 = float(f.subs(y, (lista_y[j-1] + float(f.subs(y, lista_y[j-1]).subs(t, lista_t[j-1])) * h)).subs(t, lista_t[j-1]+h))
			k7 = float(lista_y[j-1]) + b6*h*k6 + b5*h*k5 + b4*h*k4 + b3*h*k3 + b2*h*k2 + b1*h*k1 + b0*h*k0

			t2 = float(lista_t[j-1] + h)
			lista_y.append(k7)
			lista_t.append(t2)
			print("%2.d %.7f" % (j, float(lista_y[j])))
			wrt.write('{:2d}'.format(j) + ' {:.7f}'.format(float(lista_y[j])) + '\n')
	elif (ordem == 8):
		b0 = 275/24192
		b1 = -11351/120960
		b2 = 1537/4480
		b3 = -88547/120960
		b4 = 123133/120960
		b5 = -4511/4480
		b6 = 139849/120960
		b7 = 5257/17280
		for j in range(ordem , n + 1):
			k0 = float(f.subs(y, float(lista_y[j-7])).subs(t, float(lista_t[j-7])))
			k1 = float(f.subs(y, float(lista_y[j-6])).subs(t, float(lista_t[j-6])))
			k2 = float(f.subs(y, float(lista_y[j-5])).subs(t, float(lista_t[j-5])))
			k3 = float(f.subs(y, float(lista_y[j-4])).subs(t, float(lista_t[j-4])))
			k4 = float(f.subs(y, float(lista_y[j-3])).subs(t, float(lista_t[j-3])))
			k5 = float(f.subs(y, float(lista_y[j-2])).subs(t, float(lista_t[j-2])))
			k6 = float(f.subs(y, float(lista_y[j-1])).subs(t, float(lista_t[j-1])))
			k7 = float(f.subs(y, (lista_y[j-1] + float(f.subs(y, lista_y[j-1]).subs(t, lista_t[j-1])) * h)).subs(t, lista_t[j-1]+h))
			k8 = float(lista_y[j-1]) + b7*h*k7 + b6*h*k6 + b5*h*k5 + b4*h*k4 + b3*h*k3 + b2*h*k2 + b1*h*k1 + b0*h*k0

			t2 = float(lista_t[j-1] + h)
			lista_y.append(k8)
			lista_t.append(t2)
			print("%2.d %.7f" % (j, float(lista_y[j])))
			wrt.write('{:2d}'.format(j) + ' {:.7f}'.format(float(lista_y[j])) + '\n')
	return lista_t, lista_y

def formula_inversa(param, t_inic, y_inic, ordem):
	lista_t = t_inic
	lista_y = y_inic
	f = parse_expr(param[ordem + 4])
	h =  float(param[ordem + 2])
	n =  int(param[ordem + 3])
	
	if(ordem == 2):
		b0 = 4/3
		b1 = -1/3
		b2 = 2/3
		for j in range(ordem , n + 1):
			k0 = lista_y[j-1]
			k1 = lista_y[j-2]
			k2 = float(f.subs(y, (lista_y[j-1] + float(f.subs(y, lista_y[j-1]).subs(t, lista_t[j-1])) * h)).subs(t, lista_t[j-1]+h))
			k3 = b0*k0 + b1*k1 + b2*k2			
			t2 = float(lista_t[j-1] + h)
			lista_y.append(k3)
			lista_t.append(t2)
			print("%2.d %.7f" % (j, float(lista_y[j])))
			wrt.write('{:2d}'.format(j) + ' {:.7f}'.format(float(lista_y[j])) + '\n')
	elif (ordem == 3):		
		b0 = 18/11
		b1 = -9/11
		b2 = 2/11	
		b3 = 6/11
		for j in range(ordem , n + 1):
			k0 = lista_y[j-1]
			k1 = lista_y[j-2]
			k2 = lista_y[j-3]
			k3 = float(f.subs(y, (lista_y[j-1] + float(f.subs(y, lista_y[j-1]).subs(t, lista_t[j-1])) * h)).subs(t, lista_t[j-1]+h))
			k4 = b0*k0 + b1*k1 + b2*k2 + b3*k3
			t2 = float(lista_t[j-1] + h)
			lista_y.append(k4)
			lista_t.append(t2)
			print("%2.d %.7f" % (j, float(lista_y[j])))
			wrt.write('{:2d}'.format(j) + ' {:.7f}'.format(float(lista_y[j])) + '\n')
	elif (ordem == 4):
		b0 = 48/25
		b1 = -36/25
		b2 = 16/25	
		b3 = -3/25
		b4 = 12/25	
		for j in range(ordem , n + 1):
			k0 = lista_y[j-1]
			k1 = lista_y[j-2]
			k2 = lista_y[j-3]
			k3 = lista_y[j-4]
			k4 = float(f.subs(y, (lista_y[j-1] + float(f.subs(y, lista_y[j-1]).subs(t, lista_t[j-1])) * h)).subs(t, lista_t[j-1]+h))
			k5 = b0*k0 + b1*k1 + b2*k2 + b3*k3 + b4*k4*h
			t2 = float(lista_t[j-1] + h)
			lista_y.append(k5)
			lista_t.append(t2)
			print("%2.d %.7f" % (j, float(lista_y[j])))
			wrt.write('{:2d}'.format(j) + ' {:.7f}'.format(float(lista_y[j])) + '\n')
	elif (ordem == 5):		
		b0 = 300/137
		b1 = -300/137
		b2 = 200/137	
		b3 = -75/137
		b4 = 12/137
		b5 = 60/137	
		for j in range(ordem , n + 1):
			k0 = lista_y[j-1]
			k1 = lista_y[j-2]
			k2 = lista_y[j-3]
			k3 = lista_y[j-4]
			k4 = lista_y[j-5]
			k5 = float(f.subs(y, (lista_y[j-1] + float(f.subs(y, lista_y[j-1]).subs(t, lista_t[j-1])) * h)).subs(t, lista_t[j-1]+h))
			k6 = b0*k0 + b1*k1 + b2*k2 + b3*k3 + b4*k4 + b5*k5*h
			t2 = float(lista_t[j-1] + h)
			lista_y.append(k6)
			lista_t.append(t2)
			print("%2.d %.7f" % (j, float(lista_y[j])))
			wrt.write('{:2d}'.format(j) + ' {:.7f}'.format(float(lista_y[j])) + '\n')
	elif (ordem == 6):
		b0 = 360/147
		b1 = -450/147
		b2 = 400/147	
		b3 = -225/147
		b4 = 72/147
		b5 = 10/147	
		b6 = 60/147
		for j in range(ordem , n + 1):
			k0 = lista_y[j-1]
			k1 = lista_y[j-2]
			k2 = lista_y[j-3]
			k3 = lista_y[j-4]
			k4 = lista_y[j-5]
			k5 = lista_y[j-6]
			k6 = float(f.subs(y, (lista_y[j-1] + float(f.subs(y, lista_y[j-1]).subs(t, lista_t[j-1])) * h)).subs(t, lista_t[j-1]+h))
			k7 = b0*k0 + b1*k1 + b2*k2 + b3*k3 + b4*k4 + b5*k5+ b6*k6*h
			t2 = float(lista_t[j-1] + h)
			lista_y.append(k7)
			lista_t.append(t2)
			print("%2.d %.7f" % (j, float(lista_y[j])))
			wrt.write('{:2d}'.format(j) + ' {:.7f}'.format(float(lista_y[j])) + '\n')
	return lista_t, lista_y

arq = open('input.txt', 'r')
wrt = open('output.txt', 'w')
for linha in arq:
	param = linha.split()
	if(param[0]=="euler"):
		wrt.write('Metodo de Euler\n')
		print('Metodo de Euler')
		graphic_t, graphic_y  = euler(param)

		plt.xlabel('t')
		plt.ylabel('y')
		plt.title('Metodo de Euler')
		plt.plot(graphic_t, graphic_y)
		plt.show()
	elif(param[0]=="euler_aprimorado"):
		wrt.write('Metodo de Euler Aprimorado\n')
		print('Metodo de Euler Aprimorado')
		graphic_t, graphic_y  = euler_aprimorado(param)
		plt.xlabel('t')
		plt.ylabel('y')
		plt.title('Metodo de Euler Aprimorado')
		plt.plot(graphic_t, graphic_y)
		plt.show()
	elif(param[0]=="euler_inverso"):
		wrt.write('Metodo de Euler Inverso\n')
		print('Metodo de Euler Inverso')
		graphic_t, graphic_y  = euler_inverso(param)
		plt.xlabel('t')
		plt.ylabel('y')
		plt.title('Metodo de Euler Inverso')
		plt.plot(graphic_t, graphic_y)
		plt.show()
	elif(param[0]=="runge_kutta"):
		wrt.write('Metodo de Runge Kutta\n')
		print('Metodo de Runge Kutta')
		graphic_t, graphic_y  = runge_kutta(param)
		plt.xlabel('t')
		plt.ylabel('y')
		plt.title('Metodo de Runge Kutta')
		plt.plot(graphic_t, graphic_y)
		plt.show()
	elif("adam_bashforth" in str(param[0])):
		ordem = int(param[len(param)-1])
		if("euler_aprimorado" in str(param[0])):
			plt.title('Metodo de Adam Bashforth por Euler Aprimorado')
			wrt.write('Adam_Bashforth by Euler_Aprimorado\n')
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
			plt.title('Metodo de Adam Bashforth por Euler Inverso')
			wrt.write('Adam_Bashforth by Euler_Inverso\n')
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
			plt.title('Metodo de Adam Bashforth por Runge Kutta')
			wrt.write('Adam_Bashforth by Runge_Kutta\n')
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
			plt.title('Metodo de Adam Bashforth por Euler')
			wrt.write('Adam_Bashforth by Euler\n')
			print("Adam_Bashforth by Euler")
			param_euler = param[0:6]			
			param_euler[4]= int(param[6])-1
			t_inic,y_inic = euler(param_euler)
			param_aux = []
			param_aux[0:1] = param[0:1]
			param_aux[1:ordem] = y_inic
			param_aux[ordem+1:ordem+4]=param[2:7]
			param = param_aux
		else:		
			plt.title('Metodo de Adam Bashforth')	
			wrt.write('Metodo Adan-Bashforth ordem'+ str(ordem) + '\n')			
			wrt.write('y('+ str(param[ordem+1]) + ') = '+ str(param[1])+'\n')
			wrt.write('h =' + str(param[ordem+2]) + '\n')
			print("Metodo Adan-Bashforth ordem", ordem)
			print('y(',param[ordem+1],') =',param[1])
			print('h =', param[ordem+2])
			t_inic = []
			y_inic = []
			h_aux = 0			
			t0 = float(param[ordem+1])
			i = 0			
			for i in range(ordem): 			
				y_inic.append(float(param[i+1]))
				t_inic.append(float(t0+h_aux))
				h_aux = float(param[len(param)-4])+h_aux
				print("%2.d %.7f" % (i, y_inic[i]))
				wrt.write('{:2d}'.format(i) + ' {:.7f}'.format(y_inic[i]) + '\n')			
		graphic_t, graphic_y  = adam_bashforth(param, t_inic, y_inic, ordem )
		plt.xlabel('t')
		plt.ylabel('y')
		plt.plot(graphic_t, graphic_y)	
		plt.show()
	elif("adam_moulton" in str(param[0])):
		ordem = int(param[len(param)-1])
		if("euler_aprimorado" in str(param[0])):
			plt.title('Metodo de Adam Moulton por Euler Aprimorado')
			wrt.write("Adam_Moulton by Euler_Aprimorado\n")
			print("Adam_Moulton by Euler_Aprimorado")
			param_euler_aprimorado = param[0:6]			
			param_euler_aprimorado[4]= int(param[6])-1
			t_inic,y_inic = euler_aprimorado(param_euler_aprimorado)
			param_aux = []
			param_aux[0:1] = param[0:1]
			param_aux[1:ordem] = y_inic
			param_aux[ordem+1:ordem+4]=param[2:7]
			param = param_aux
		elif("euler_inverso" in str(param[0])):
			plt.title('Metodo de Adam Moulton por Euler Inverso')
			wrt.write("Adam_Moulton by Euler_Inverso\n")
			print("Adam_Moulton by Euler_Inverso")
			param_euler_inverso = param[0:6]			
			param_euler_inverso[4]= int(param[6])-1
			t_inic,y_inic = euler_inverso(param_euler_inverso)
			param_aux = []
			param_aux[0:1] = param[0:1]
			param_aux[1:ordem] = y_inic
			param_aux[ordem+1:ordem+4]=param[2:7]
			param = param_aux
		elif("runge_kutta" in str(param[0])):
			plt.title('Metodo de Adam Moulton por Runge Kutta')
			wrt.write("Adam_Moulton by Runge_Kutta\n")
			print("Adam_Moulton by Runge_Kutta")
			param_rk = param[0:6]			
			param_rk[4]= int(param[6])-1
			t_inic,y_inic = runge_kutta(param_rk)
			param_aux = []
			param_aux[0:1] = param[0:1]
			param_aux[1:ordem] = y_inic
			param_aux[ordem+1:ordem+4]=param[2:7]
			param = param_aux
		elif("euler" in str(param[0])):
			plt.title('Metodo de Adam Moulton por Euler')
			wrt.write("Adam_Moulton by Euler\n")
			print("Adam_Moulton by Euler")
			param_euler = param[0:6]			
			param_euler[4]= int(param[6])-1
			t_inic,y_inic = euler(param_euler)
			param_aux = []
			param_aux[0:1] = param[0:1]
			param_aux[1:ordem] = y_inic
			param_aux[ordem+1:ordem+4]=param[2:7]
			param = param_aux
		else:				
			plt.title('Metodo de Adam Moulton')
			wrt.write('Metodo Adan-Moulton ordem'+ str(ordem) + '\n')			
			wrt.write('y('+ str(param[ordem+1]) + ') = '+ str(param[1])+'\n')
			wrt.write('h =' + str(param[ordem+2]) + '\n')
			print("Metodo Adan-Moulton ordem", ordem)
			print('y(',param[ordem+1],') =',param[1])
			print('h =', param[ordem+2])
			t_inic = []
			y_inic = []
			h_aux = 0			
			t0 = float(param[ordem+1])
			i = 0			
			for i in range(ordem): #inicializacao dos parametros				
				y_inic.append(float(param[i+1]))
				t_inic.append(float(t0+h_aux))
				h_aux = float(param[len(param)-4])+h_aux
				print("%2.d %.7f" % (i, y_inic[i]))
				wrt.write('{:2d}'.format(i) + ' {:.7f}'.format(y_inic[i]) + '\n')						
		graphic_t, graphic_y  = adam_moulton(param, t_inic, y_inic, ordem )
		plt.xlabel('t')
		plt.ylabel('y')

		plt.plot(graphic_t, graphic_y)	
		plt.show()
	elif("formula_inversa" in str(param[0])):
		ordem = int(param[len(param)-1])
		if("euler_aprimorado" in str(param[0])):
			plt.title('Metodo da Formula Inversa por Euler Aprimorado')
			wrt.write("formula_inversa by Euler_Aprimorado\n")
			print("formula_inversa by Euler_Aprimorado")
			param_euler_aprimorado = param[0:6]			
			param_euler_aprimorado[4]= int(param[6])-1
			t_inic,y_inic = euler_aprimorado(param_euler_aprimorado)
			param_aux = []
			param_aux[0:1] = param[0:1]
			param_aux[1:ordem] = y_inic
			param_aux[ordem+1:ordem+4]=param[2:7]
			param = param_aux
		elif("euler_inverso" in str(param[0])):
			plt.title('Metodo da Formula Inversa por Euler Inverso')
			wrt.write("formula_inversa by Euler_Inverso\n")
			print("formula_inversa by Euler_Inverso")
			param_euler_inverso = param[0:6]			
			param_euler_inverso[4]= int(param[6])-1
			t_inic,y_inic = euler_inverso(param_euler_inverso)
			param_aux = []
			param_aux[0:1] = param[0:1]
			param_aux[1:ordem] = y_inic
			param_aux[ordem+1:ordem+4]=param[2:7]
			param = param_aux
		elif("runge_kutta" in str(param[0])):
			plt.title('Metodo da Formula Inversa por Runge Kutta')
			wrt.write("formula_inversa by Runge_Kutta\n")
			print("formula_inversa by Runge_Kutta")
			param_rk = param[0:6]			
			param_rk[4]= int(param[6])-1
			t_inic,y_inic = runge_kutta(param_rk)
			param_aux = []
			param_aux[0:1] = param[0:1]
			param_aux[1:ordem] = y_inic
			param_aux[ordem+1:ordem+4]=param[2:7]
			param = param_aux
		elif("euler" in str(param[0])):
			plt.title('Metodo da Formula Inversa por Euler')
			wrt.write("formula_inversa by Euler\n")
			print("formula_inversa by Euler")
			param_euler = param[0:6]			
			param_euler[4]= int(param[6])-1
			t_inic,y_inic = euler(param_euler)
			param_aux = []
			param_aux[0:1] = param[0:1]
			param_aux[1:ordem] = y_inic
			param_aux[ordem+1:ordem+4]=param[2:7]
			param = param_aux
		else:	
			plt.title('Metodo da Formula Inversa')
			wrt.write('Metodo formula_inversa ordem'+ str(ordem) + '\n')			
			wrt.write('y('+ str(param[ordem+1]) + ') = '+ str(param[1])+'\n')
			wrt.write('h =' + str(param[ordem+2]) + '\n')
			print("Metodo formula_inversa ordem", ordem)
			print('y(',param[ordem+1],') =',param[1])
			print('h =', param[ordem+2])
			t_inic = []
			y_inic = []
			h_aux = 0			
			t0 = float(param[ordem+1])
			i = 0
			for i in range(ordem): #inicializacao dos parametros				
				y_inic.append(float(param[i+1]))
				t_inic.append(float(t0+h_aux))
				h_aux = float(param[len(param)-4])+h_aux
				print("%2.d %.7f" % (i, y_inic[i]))
				wrt.write('{:2d}'.format(i) + ' {:.7f}'.format(y_inic[i]) + '\n')						
		graphic_t, graphic_y  = formula_inversa(param, t_inic, y_inic, ordem )
		plt.xlabel('t')
		plt.ylabel('y')
		plt.plot(graphic_t, graphic_y)	
		plt.show()
	wrt.write('\n')	
wrt.close()	
arq.close()
