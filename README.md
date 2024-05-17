## Установка membranesolver

Программа с которой вам предстоит работать представленна целью `model` в проекте https://boogie.inm.ras.ru/liogky/membranemodel . Проект имеет очень много зависимостей (более 10 библиотек):
- [INMOST](https://github.com/INMOST-DEV/INMOST) базовая версия, для чтения xml-файлов и для решения линейных систем
- [Ani3d](https://sourceforge.net/projects/ani3d/) модули для построения сеток
- [CasAdi](https://github.com/casadi/casadi) используется для работы с символьными выражениями
- [Rpoly](https://github.com/sweeneychris/RpolyPlusPlus) предоставляет эффективный решатель полиномиальных уравнений
- [Eigen3](https://eigen.tuxfamily.org/) используется для работы с плотными матрицами
- [GSL](https://www.gnu.org/software/gsl/) используется для вычисления эллиптических интегралов и решения полиномиальных уравнений
- [CGAL](https://www.cgal.org/download.html) входит в ядро модели и предоставлет структуру для треугольной сетки
- [nlohmann_json](https://github.com/nlohmann/json) используется для работы с json форматом некоторыми частями проекта
- [Boost](https://www.boost.org/) требуется для CGAL
- [SUNDIALS](https://computing.llnl.gov/projects/sundials) предоставляет нелинейную решалку kinsol
- [blas](https://ru.wikipedia.org/wiki/Basic_Linear_Algebra_Subprograms) и [lapack](https://ru.wikipedia.org/wiki/LAPACK) - требуются для Ani3d
- стандартная библиотека языка Фортан - требуются для Ani3d
- [GL](https://packages.ubuntu.com/search?keywords=libgl-dev) - графическая библиотека для Ani3d
- [gmp](https://gmplib.org/) и [mpfr](https://www.mpfr.org/) - требуются для  CGAL, библиотеки для работы с числами произвольной точности

К счастью, вам не придётся возится с установкой всего этого зоопарка.
Я подготовил и оттестировал скрипты, которые автоматически выкачивают и устанавливают все необходимые зависимости.

Итак, 
1. Устанавливаем поддерживаемые пакетным менеджером зависимости:
```bash
$ sudo apt -y install cmake pkg-config git gfortran\
	libblas-dev liblapack-dev libgl-dev libboost-all-dev\
	libgmp-dev libmpfr-dev libeigen3-dev
```
(Для запуска на образе `ubuntu` в `docker` ещё требуется установить компилятор, что можно сделать добавив в команду выше либу `build-essential`)

2. Скачиваем код программного пакета и распаковываем его. Далее заходим в папку с проектом. Если вы все сделали верно, то команда `ls` должна выдать что-то такое
```bash
$ ls
AVSim  CMakeLists.txt  README.txt  benchmarks  cmake  config.cmake.in  examples
```
3. Далее выполняем сборку посредством CMake (она может занять порядка 15 минут!):
```bash
$ mkdir build; cd build
$ cmake -DDOWNLOAD_Ani3d=ON -DDOWNLOAD_CGAL=ON -DDOWNLOAD_GSL=ON -DDOWNLOAD_casadi=ON -DDOWNLOAD_inmost=ON -DDOWNLOAD_nlohmann_json=ON -DDOWNLOAD_rpoly=ON -DDOWNLOAD_sundials=ON -DCMAKE_BUILD_TYPE=Release ..
$ cmake --build . --target model -j8
# model - название программы, ключ -j8 ускоряет установку используя 8 потоков
```
4. Если всё завершилось нормально, то по пути `benchmarks/general` должен был появится исполняемый файл `model`. Можно убедится в его работоспособности следующим скриптом
```bash
$ ./benchmarks/general/model --help
This program is designed to simulate the behavior of a thin shell in a membrane or shell approximation.
Copyright (C) 2023, the Marchuk Institute of Numerical Mathematics of the Russian Academy of Sciences.

 Command line options: 
  -cf, --config     FILE    <Configuration file>
 ...  
```

## Установка пакета для моделирования индентирования материала

Установить [anifem++](https://github.com/INMOST-DEV/INMOST-FEM) можно выполнив следующие действия в терминале:
1. Клонируем репозиторий
```bash
git clone https://github.com/INMOST-DEV/INMOST-FEM.git
```
2. Устанавливаем поддерживаемые пакетным менеджером зависимости
```bash
cd INMOST-FEM; mkdir build; cd build

cmake -DWITH_INMOST=ON -DDOWNLOAD_inmost=ON -DWITH_KINSOL=ON 
-DDOWNLOAD_sundials=ON -DWITH_EIGEN=ON -DDOWNLOAD_eigen3=ON 
-DWITH_CASADI=ON -DDOWNLOAD_casadi=ON -DCOMPILE_EXAMPLES=ON 
-DCOMPILE_TESTS=ON ../
```
3. Собираем проект при помощи CMake
```bash
cmake --build .

cmake -DCMAKE_INSTALL_PREFIX=/home/student/libs/anifem++_install ..

cmake --install .
```
