# online-spp
Implements Thompson Sampling Algorithm using a Gaussian Field distribution assumption over arcs velocity

Execute the codes inside the project folder.

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

Folder containing the codes should look like:

/path_to_folder/online-spp -> project folder

/path_to_folder/data/new_large/nodes.csv: format -> id;latitude;longitude

/path_to_folder/data/new_large/arcs.csv: format -> id;node1;node2;valores


# Run artificial instances to replicate the results in the paper:


Windows: 
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
python main.py 150 151 1 0 1 20
python main.py 150 151 1 0 3 20
python main.py 150 151 1 0 4 20
python main.py 150 151 2 0 1 20
python main.py 150 151 2 0 3 20
python main.py 150 151 2 0 4 20
python main.py 150 151 3 0 1 20
python main.py 150 151 3 0 3 20
python main.py 150 151 3 0 4 20
python main.py 150 151 4 0 1 20
python main.py 150 151 4 0 3 20
python main.py 150 151 4 0 4 20

```
Linux: 
```
for number1 in {1..4}; do for number2 in {1,3,4}; do python main.py 50 51 $number1 0 $number2 500; done; done
for number1 in {1..4}; do for number2 in {1,3,4}; do python main.py 150 151 $number1 0 $number2 20; done; done
```
