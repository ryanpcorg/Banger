package main

import (
	"fmt"
)

func subtract() {
	var one float64
	var two float64

	fmt.Println("What is number one?")
	fmt.Scan(&one)
	clearScreen()

	fmt.Println("What is number two?")
	fmt.Scan(&two)
	clearScreen()

	var total float64 = one - two
	fmt.Print(total)
}
