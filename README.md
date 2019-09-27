MGranul is an implementation of granulometry based on cross-correlation coded with OpenCV 2.4+. In this repo, there are three folder: `cpp` (the implementation itself), `python` (a bind for Python), `docker` (DockerFile to build a container with OpenCV2 and MGranul).

## Dependencies
To build fully MGranul are required these dependencies:
```
 - git;
 - cmake;
 - make;
 - g++;
 - Python3.6;
 - python3-dev;
```

Please, check online on how to get all these dependencies installed.

## Building Opencv2

MGranul relies on OpenCV2 to provide several cv algorithms. There is a script that will install and save the built dependencies in a folder in the same folder of MGranul. To run this script, open your terminal and type:

```
$ cd cpp/
$ ./build_opencv2.sh
```

If you prefer it is possible to use a version previously installed of OpenCV2. But, there is no guarantee of working.

## Building MGranul

After the OpenCV2 build, you are ready to build the MGranul, just type:

```
$ cd cpp/
$ ./build.sh
```

It will compile and will create a `bin/mgranul` file. This file is a shell script that maps the local OpenCV2 dependency and the `exe_granul`, the actual executable. To run msgranul:

```
$ bin/mgranul

< MGranul.exe: Programas para granulometria multi-formas v1.2.1>
_______________________________________________________________________________
Programas:
  Kernel   - Gera imagem das mascaras a partir de kernel.cfg
  Correla  - Maximos locais (.ho1) das correlacoes com mascaras multiformas
  Mostra   - Mostra maximos locais (.ho1) sem filtrar
  Filtra   - Filtra maximos locais (.ho1) e mostra
  Relat    - Le .ho2 e gera .rel
  Msgranul - Filtra maximos locais (.ho1) e mostra
  Msgranul_kmeans - Filtra maximos locais (.ho1) e mostra
...............................................................................
```

## Building Python bind

In order to build Python bind, we must have to have MGranul built in the `cpp` and other dependencies. To install all dependencies just run:

```
$ pip3 install numpy
$ cd python/
$ ./build_python.sh
```

The python module will be put in `module`.

## Platform

MSGranul was successfully built on macOS, Ubuntu, Debian, and Fedora.
