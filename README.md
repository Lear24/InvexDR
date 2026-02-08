# InvexDR

This code of simulation part for paper "Distributed prediction under heterogeneity with unidentifiable parameter" 

## Requirements

`addpath(genpath(('./')))`

## simulation setting
- coef1 - sig1 - m10 - n300; p=10/16;
- ex1: `m = sin(2*xb1)*exp(xb2);`  set1 - p16  m10 
- ex2: `m = 3*xb1./(1+(1+xb2).^2);` set2 - p10  m5
- `optionsKer = 'mave'
- `lamALL = 1/1.2/1.6;` for heterogeneity parameter  theta = pi/3, pi/4, pi/8.

## comparison method
- Our method: InvexDR
- 1. MAVE
- 2. NR-B: The separate NR-iteration with the identity block constraints. Initialization with MAVE
- 3. NR-O: The separate NR-iteration with the orthogonal constraints.Initialization with MAVE (with `mNR-O`)
- 4. InvexDR-B:  Our method under identity block constraints.
- 5. InvexDR-O: Our method under orthogonal constraints.
- 6. pNR-B: the pooling NR-iteration with  the identity block constraints. Initialization with pMAVE.
- 7. pNR-O: the pooling NR-iteration with  the orthogonal constraints.Initialization with pMAVE.
 

## Simulation 

- 1. Vary n: [100, 200, 300, 400, 500];  sig = 1; coef1;
	- ex1-p16-m10; 
	- ex2-p10-m5;  
- 2. Vary m:  [2, 4, 6, 10, 16];   with n = 400; 
	- ex1-p16; show 4;
	- ex2-p10; show 5;
- 3. Vary sig:  [0.8, 1, 1.2] with n 400
	- m =10; n=300; p=10;
- 4. Vary sd: [0.3, 0.5,  0.7]; 
	- m =10; n=300; p = 10; sig = 1;

## Code 

- Lib
- opt
- simulaiton
  - var_n
  - var_m
  - var_sig
  - var_sd
