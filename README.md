# EntropiaConfiguracional
Implementação do cálculo da entropia configuracional

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
