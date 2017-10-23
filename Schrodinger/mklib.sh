#!/bin/bash

# My cpp lib make scrip

# Colors
echo_black() { echo "$(tput setaf 0)$*$(tput setaf 7)"; }
echo_darkred() { echo "$(tput setaf 1)$*$(tput setaf 7)"; }
echo_green() { echo "$(tput setaf 2)$*$(tput setaf 7)"; }
echo_yellow() { echo "$(tput setaf 3)$*$(tput setaf 7)"; }
echo_purple() { echo "$(tput setaf 4)$*$(tput setaf 7)"; }
echo_pink() { echo "$(tput setaf 5)$*$(tput setaf 7)"; }
echo_cyan() { echo "$(tput setaf 6)$*$(tput setaf 7)"; }
echo_white() { echo "$(tput setaf 7)$*$(tput setaf 7)"; }
echo_grey() { echo "$(tput setaf 8)$*$(tput setaf 7)"; }
echo_red() { echo "$(tput setaf 9)$*$(tput setaf 7)"; }

clear

# Select parameters
program="$1"  #
library="$2"  # MetComp
compiler="$3" # clang, gcc, intel
opt="$4"      # Optimization: 0, 1, 2, 3
debug="-debug full"

echo_green "Compiling lib$library and $program with $compiler $debug -O$opt"

printf '\n'
echo_purple "Creating directories..."
mkdir -p obj lib bin

# Compile library
printf '\n'
echo_purple "Compiling library to object file..."
if [ $compiler = "gcc" ]; then
  g++ -Wall -c $debug -O$opt ./src/$library.cpp -std=c++17 -I ./inc -o ./obj/$library.o
elif [ $compiler = "clang" ]; then
  clang -v -c -O$opt ./src/$library.cpp -std=c++17 -I ./inc -o ./obj/$library.o
elif [ $compiler = "intel" ]; then
  time icpc -Wall -c -O$opt ./src/$library.cpp -o ./obj/$library.o
fi

# Generate library file
printf '\n'
echo_purple "Generating library archive file..."
if [ $compiler = "gcc" ] ||  [ $compiler = "intel" ]; then
  time ar rcs ./lib/lib$library.a ./obj/$library.o
elif [ $compiler = "clang" ]; then
  libtool -static ./obj/$library.o -o ./lib/lib$library.a
fi

# Compile program and link
printf '\n'
if [ $compiler = "gcc" ]; then
  echo_purple "Compiling $program and linking to $library..."
  g++ -Wall -O$opt $program.cpp -std=c++17 -I ./inc -L ./lib -static -l$library -o ./bin/$program
elif [ $compiler = "intel" ]; then
  echo_purple "Compiling $program and linking to $library..."
  time icpc -Wall -O$opt -L ./lib -I ./inc $program.cpp ./lib/lib$library.a -o ./bin/$program
elif [ $compiler = "clang" ]; then
  echo_purple "Compiling $program to object file..."
  clang -v -c $debug -O$opt ./src/$program.cpp -std=c++17 -I ./inc -o ./obj/$program.o
  echo_purple "Linking $program and $library..."
  ld -v ./obj/$program.o -framework CoreFoundation -lSystem -L. -l$library -o ./bin/$program
fi

# Run test program
printf '\n'
echo_purple "Running $program..."
printf '\n'
time ./bin/$program

printf '\n'
