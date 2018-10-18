# Projeto de Métodos Numéricos Computacionais 
Este é um projeto pertencente a Matheus Deodato (e-mail mdo2@cin.ufpe.br)
foi desenvolvido como projeto da disciplina de Métodos Numéricos Computacionais.<br/>
O programa executa o método de Euler, Euler Inverso, Euler Aprimorado, Runge-Kutta (4 ordem), Adam-Bashforth(da 2º a 8º ordem),
Adam-Moulton(da 1º a 7º ordem) e a Fórmula Inversa (da 2º a 8º ordem).
## Requisitos de Execução
1. Ambiente Linux
2. Python 3
3. Biblioteca Python sympy
4. Biblioteca Python numpy
5. Biblioteca Python matplotlib
#### Para executar:</br>
###### ` $ python3 imp_met.py `

## Entrada do Programa
As entradas do programa devem estar em um arquivo que deve ficar no diretório do programa. O arquivo de entrada 
deve se chamar **input.txt**, e deve ser escrito seguindo a sintaxe especificada adiante.</br>

###### `nome_do_metodo t(0) y(0) [y(1) .. y(ordem)] h n f [ordem] /n ` </br>

Pode-se ter várias entradas no arquivo.

 - nome_do_metodo : Os possíveis nomes para os métodos, são:</br>
  eulers_runge_kutta<br/>
  euler</br>
  euler_inverso</br>
  euler_aprimorado</br>
  runge_kutta</br>
  adam_bashforth</br>
  adam_bashforth_by_euler</br>
  adam_bashforth_by_euler_inverso</br>
  adam_bashforth_by_euler_aprimorado</br>
  adam_bashforth_by_runge_kutta</br>
  adam_multon</br>
  adam_multon_by_euler</br>
  adam_multon_by_euler_inverso</br>
  adam_multon_by_euler_aprimorado</br>
  adam_multon_by_runge_kutta</br>
  formula_inversa </br>
  formula_inversa_by_euler</br>
  formula_inversa_by_euler_inverso</br>
  formula_inversa_by_euler_aprimorado</br>
  formula_inversa_by_runge_kutta</br>
- t(0): Valor inicial do t</br>
- y(0): Valor inicial de y</br>
- [y(1) .. y(ordem)] : Valores Iniciais usados em Adam-Bashforth, Adam-Moulton e Formula Inversa</br>
- h : Tamanho do passo</br>
- n : Quantidade de passos</br>
- f : função</br>
- [ordem] : ordem do método (usado em Adam-Bashforth, Adam-Moulton e Fórmula Inversa, cada um com suas variações)</br>

#### Exemplo de Entrada:
```
eulers_runge_kutta 0 0 0.1 20 1-t+4*y
euler 0 0 0.1 20 1-t+4*y
runge_kutta 0 0 0.1 20 1-t+4*y
adam_bashforth 0.0 0.1 0.23 0.402 0.6328 0 0.1 20 1-t+4*y 5
adam_bashforth_by_euler 0 0 0.1 20 1-t+4*y 6
adam_multon 0.0 0.1 0.23 0.402 0.6328 0 0.1 20 1-t+4*y 6
formula_inversa 0.0 0.1 0.23 0.402 0.6328 0 0.1 20 1-t+4*y 5


```

## Saída do Programa
O programa exibe um gráfico (o grafico de eulers_runge_kutta é plotado de uma vez) para cada função dada no arquivo, e guarda o resultado no diretório do programa
como um arquivo chamado **`output.txt`**.




