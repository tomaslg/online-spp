# online-spp
Implements Thompson Sampling Algorithm using a Gaussian Field distribution assumption over arcs velocity


# Run real life instance to replicate the results in the paper:

Windows: 
```
FOR /L %i IN (0,1,99) DO python main_from_paths.py 150 %i
```
Linux: 
```
for number in {0..99}; do python main_from_paths.py 150 $number; done
```
## Data:
https://github.com/felipelagos/beijing-instance.git

# Run artificial instances to replicate the results in the paper:

```
python main.py 50 51 1 0 1 500
python main.py 50 51 1 0 3 500
python main.py 50 51 1 0 4 500
python main.py 50 51 2 0 1 500
python main.py 50 51 2 0 3 500
python main.py 50 51 2 0 4 500
python main.py 50 51 3 0 1 500
python main.py 50 51 3 0 3 500
python main.py 50 51 3 0 4 500
python main.py 50 51 4 0 1 500
python main.py 50 51 4 0 3 500
python main.py 50 51 4 0 4 500
```
