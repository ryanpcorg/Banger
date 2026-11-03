package main

import (
	"fmt"
	"os"
	"os/exec"
	"runtime"
)

func clearScreen() {
	var cmd *exec.Cmd
	if runtime.GOOS == "windows" {
		cmd = exec.Command("cmd", "/c", "cls")
	} else {
		cmd = exec.Command("clear")
	}
	cmd.Stdout = os.Stdout
	cmd.Run()
}

func main() {
	var one float64
	var two float64

	fmt.Println("What is number one for multiplying?")
	fmt.Scan(&one)

	clearScreen()

	fmt.Println("Number two now please.")
	fmt.Scan(&two)

	clearScreen()

	var total float64 = one * two

	println(total)
}
