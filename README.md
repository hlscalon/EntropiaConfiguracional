# EntropiaConfiguracional
Implementação do cálculo da entropia configuracional

### Instalação

- Rodar script instalar.sh ou

- Compilar nauty:
	- Entrar em src/cpp/libs/nauty
	- ./configure
	- alterar Makefile:
		- adicionar opção -fPIC em CFLAGS
	- make
	- alterar nauty.h:
		- substituir bloco `#define _FILE_OFFSET_BITS` por:
		
		```
		#ifndef _FILE_OFFSET_BITS
		#define _FILE_OFFSET_BITS 0
		#if _FILE_OFFSET_BITS == 64
		#define _LARGEFILE_SOURCE
		#else
		#undef _FILE_OFFSET_BITS
		#endif
		#endif
		```

- Mudanças Aboria:
	- Entrar em src/cpp/libs/Aboria/third-party/nanoflann/nanoflann.hpp
	- Alterar WORDSIZE para WORDSIZE_NANNOFLANN

### Rodar o programa

- Argumentos:

```
<$1>: from measurement (Y or N)
```

- Se sim (Y)

```
<$2>: measurement filepath
```

- Se não (N)

```
<$2>: xyz filepath
<$3>: file type (V = Vink file or N = Normal)
<$4>: covalent radii cut off
<$5>: c
<$6>: initial n
<$7>: final n
<$8>: calculate (Y or N)
```

- Possíveis mudanças necessárias para rodar o programa:

    - Alterar número de posições aleatórias geradas, (padrão `m = n^2N`)
        - Alterar função `getNumberRandomPositions` no arquivo `src/python/configurational_entropy.py`
    - Alterar eixos do gráfico gerado, (padrão `-10, 10`):
        - Alterar função `calculateConfigurationalEntropy` no arquivo `src/python/configurational_entropy.py`
	- `plt.axis([n1, n2, -10, 10])`

#### Exemplo

- python2 main.py N arquivos_xyz/SILICA_SAMPLES/3k/3k-total-2-ase.xyz V 1.12 5 20 22 N
- python2 main.py Y medicoes/med_fcc2048.xyz_01_1_2019_00_00_00.ce
