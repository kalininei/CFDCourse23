1. Сборка
=========

1.1 Windows, VisualStudio
-------------------------
Проверено на Windows 10, Visual Studio 2019

- установить cmake версии выше 3.0
https://github.com/Kitware/CMake/releases/download/v3.24.2/cmake-3.24.2-windows-x86_64.msi
Нужно, чтобы cmake.exe был в системных путях. Если ставить через msi, он об этом спросит при установке.

- cоздать папку build в корне проекта

- скопировать скрипт winbuild64.bat в папку build. Далее вносить изменения
  только в скопированном файле.

- скрипт написан для версии Visual Studio 2019. Если используется другая версия,
  изменить в скрипте значение переменной CMGenerator на соответствующие вашей версии.
  Значения для разных версий Visual Studio написаны ниже

SET CMGenerator="Visual Studio 17 2022"
SET CMGenerator="Visual Studio 16 2019"
SET CMGenerator="Visual Studio 15 2017"
SET CMGenerator="Visual Studio 14 2015"

- запустить скрипт winbuild64.bat из папки build. Нуже доступ к интернету.
  В процессе будет скачано около 200Мб пакетов, поэтому первый запуск может занять время

- После сборки в папке build появится проект VisualStudio cfdapp.sln
  Сборка бинарных файлов будет осуществлятся в папку build/bin.
  Выходные файлы при использовании функции from_output_path будут писаться в build/output

- Для проверки можно собрать и запустить программу poisson.exe

1.2 Linux, GCC
--------------
- в консоли
mkdir build
cd build
cmake ..
make
- программы соберутся в папку build/bin
