package main

import (
	"fmt"
)

func mult() {
	var one float64
	var two float64

	fmt.Println("What is number one for multiplying?")
	fmt.Scan(&one)

	clearScreen()

	fmt.Println("Number two now please.")
	fmt.Scan(&two)

	clearScreen()

	var total float64 = one * two

	println("Your number is:", total)
}
