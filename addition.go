package main

import (
	"fmt"
	"time"
)

func addition() {
	var one float64
	var two float64

	fmt.Println("What is number 1 for adding?")
	fmt.Scan(&one)
	clearScreen()

	fmt.Println("What is number 2?")
	fmt.Scan(&two)
	clearScreen()

	var total float64 = one + two
	fmt.Print("Your number is:", total)
	time.Sleep(3000)
	clearScreen()
}
