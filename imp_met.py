import parser  
import matplotlib.pyplot as plt

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
	print()
arq.close()
