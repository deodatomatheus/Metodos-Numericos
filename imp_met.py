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

		
arq = open('input.txt', 'r')
for linha in arq:
	param = linha.split()
	if(param[0]=="euler"):
		graphic_t, graphic_y  = euler(param)
		plt.plot(graphic_t, graphic_y)
		plt.show()

arq.close()
