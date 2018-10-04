import parser  
import matplotlib.pyplot as plt
import sympy as sym
from sympy.parsing.sympy_parser import parse_expr
sym.init_printing()
y, t = sym.symbols('y t')

def euler(param):
	print('Metodo de Euler')
	y = float(param[1])
	t = float(param[2])	
	h = float(param[3])
	n = int(param[4])
	f = parser.expr(param[5]).compile()
	lista_t = []
	lista_y = []
	print('y(',t,') =',y)
	print('h =', h)
	for j in range(n+1):
		print(j, " " , y)
		lista_y.append(y)
		lista_t.append(t)
		k = eval(f)
		y = y + h * k
		t = t + h
	return lista_t, lista_y

def euler_inverso(param):
	print('Metodo de Euler Inverso')
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
		print(j, " " , y0)
		lista_y.append(y0)
		lista_t.append(t0)

		k = float(f.subs(y, (y0 + float(f.subs(y, y0).subs(t, t0)) * h0)).subs(t, t0+h0))
		y0 = y0 + h0 * k
		t0 = t0 + h0
	return lista_t, lista_y

def euler_aprimorado(param):
	print('Metodo de Euler Aprimorado')
	yk = float(param[1])
	tk = float(param[2])	
	h = float(param[3])
	n = int(param[4])
	f = parser.expr(param[5]).compile()
	lista_t = []
	lista_y = []
	print('y(',tk,') =',yk)
	print('h =', h)
	for j in range(n+1):
		print(j, " " , yk)
		lista_y.append(yk)
		lista_t.append(tk)
		
		t = tk
		y = yk
		k1 = eval(f)
		
		t = tk + h 
		y = yk + h * k1
		k2 = eval(f)

		yk = yk + (h/2) * (k1 + k2)
		tk = tk + h
	return lista_t, lista_y

def runge_kutta(param):
	print('Metodo de Runge Kutta')
	yk = float(param[1])
	tk = float(param[2])	
	h = float(param[3])
	n = int(param[4])
	f = parser.expr(param[5]).compile()
	lista_t = []
	lista_y = []
	print('y(',tk,') =',yk)
	print('h =', h)
	for j in range(n+1):
		print(j, " " , yk)
		lista_y.append(yk)
		lista_t.append(tk)
		
		t = tk
		y = yk
		k1 = eval(f)
		
		t = tk + h * 0.5
		y = yk + h * k1 * 0.5
		k2 = eval(f)

		t = tk + 0.5 * h
		y = yk + 0.5 * h * k2
		k3 = eval(f)

		t = tk + h
		y = yk + h * k3
		k4 = eval(f)

		yk = yk + ((h/6) * (k1 + (2 * k2) + (2 * k3) + k4))
		tk = tk + h
	return lista_t, lista_y

arq = open('input.txt', 'r')
for linha in arq:
	param = linha.split()
	if(param[0]=="euler"):
		graphic_t, graphic_y  = euler(param)
		plt.plot(graphic_t, graphic_y)
		plt.show()
	elif(param[0]=="euler_aprimorado"):
		graphic_t, graphic_y  = euler_aprimorado(param)
		plt.plot(graphic_t, graphic_y)
		plt.show()
	elif(param[0]=="euler_inverso"):
		graphic_t, graphic_y  = euler_inverso(param)
		plt.plot(graphic_t, graphic_y)
		plt.show()
	elif(param[0]=="runge_kutta"):
		graphic_t, graphic_y  = runge_kutta(param)
		plt.plot(graphic_t, graphic_y)
		plt.show()
	print()
arq.close()
