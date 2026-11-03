package main

import (
	"fmt"
	"os"
	"os/exec"
)

var choice string

func runthego() {
	var target string
	switch choice {
	case "1":
		mult()
	case "2":
		addition()
	default:
		return
	}
	cmd := exec.Command("go", "run", target)
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	cmd.Stdin = os.Stdin
	err := cmd.Run()
	if err != nil {
		fmt.Println("Error running file:", err)
	}
}

func main() {

	fmt.Println("What option do you want to do today?")
	fmt.Println("1: Multiply")
	fmt.Println("2: Addition")
	fmt.Println("3: Subtraction")
	fmt.Println("4: Division")
	fmt.Scan(&choice)
	clearScreen()
	runthego()
}
