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

# Select parameters
program="$1"
options="$2"

clear

echo_green "Compiling $program with $options"

printf '\n'
echo_purple "Creating directories..."
mkdir -p obj lib bin

# Compile headers
printf '\n'
echo_purple "Precompiling headers to object files and creating library archive..."
headers=./inc/*.h
for header in $headers
do
 file=${header:6:-2}
  printf '\n'

  echo_cyan "$file"

  time icpc -Wall -c $options ./src/$file.cpp -I ./inc -o ./obj/$file.o
  ar rcs ./lib/libMetComp.a ./obj/$file.o
done

# Compile program and link
printf '\n'
echo_purple "Compiling $program and linking to library..."
icpc -Wall $options -L ./lib -I ./inc ./src/$program.cpp ./lib/libMetComp.a -o ./bin/$program

# Run test program
printf '\n'
echo_purple "Running $program..."
printf '\n'
#time ./bin/$program
#./bin/$program | gnuplot -p

printf '\n'
