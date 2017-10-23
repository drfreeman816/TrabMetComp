#!/bin/sh
#/usr/bin/env bash

# My cpp lib make script

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

echo_green "Cleaning folder..."
rm -Rf bin lib obj
mkdir bin lib obj

# Compile library to object code
echo_green "Compiling library to object file..."
g++ ./src/MetComp.cpp -c -std=c++17 -I ./inc -o ./obj/MetComp.o

# Generate library file
echo_green "Generating library archive file..."
ar rcs ./lib/libMetComp.a ./obj/MetComp.o

# Use library in test program
echo_green "Compiling main program to object file..."
g++ main.cpp -std=c++17 -I ./inc -L ./lib -static -lMetComp -o ./bin/main

# Run test program
echo_green "Running test program..."
./bin/main
