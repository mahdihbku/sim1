# sim1
## Compilation
```bash
gcc -O3 simulation.c -lcrypto -pthread -lm

gcc -O3 simulation3.c -lcrypto -pthread -lm
```
## Usage
```bash
./a.out <m> <n> <nbr_hash_fncts> <p_size> <nbr_threads>
./a.out 8000 1000 5 13 4

./a.out <m> <n> <k> <p_size> <car_threads> <server_threads>
./a.out 16000 2000 7 13 1 4
```
