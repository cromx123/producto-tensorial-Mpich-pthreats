Calculadora de producto tensorial con algoritmo paralelo en memoria privada y en memoria compartida
-------------------------
Compilacion: mpicc -o t2.exe t2.c
-------------------------
Ejecucion: mpirun -np k ./t2.exe -N -P -O data.txt
donde k = Numero de procesos
      N = Numero de hilos
      P = {V (Particion vertical), H (Particion horizontal)}
      O = {S (Silent), V (Verbose)}
Se incluyen data.txt y data2.txt como datos de prueba.
------------------
Maquina de prueba:

CPU: Ryzen 5 2600x
RAM: 16 gb ram
OS: Ubuntu 22.04.3
COMPILADORES: gcc 11.4.0, MPICH 4.0

