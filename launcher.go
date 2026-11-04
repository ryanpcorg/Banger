package main

import (
	"bufio"
	"fmt"
	"os"
	"strings"
)

func main() {
	scanner := bufio.NewScanner(os.Stdin)

	println("██████╗  █████╗ ███╗   ██╗ ██████╗ ███████╗██████╗ ")
	println("██╔══██╗██╔══██╗████╗  ██║██╔════╝ ██╔════╝██╔══██╗")
	println("██████╔╝███████║██╔██╗ ██║██║  ███╗█████╗  ██████╔╝")
	println("██╔══██╗██╔══██║██║╚██╗██║██║   ██║██╔══╝  ██╔══██╗")
	println("██████╔╝██║  ██║██║ ╚████║╚██████╔╝███████╗██║  ██║")
	println("╚═════╝ ╚═╝  ╚═╝╚═╝  ╚═══╝ ╚═════╝ ╚══════╝╚═╝  ╚═╝")
	println("---------------------------------------------------")
	println("The calculator your mum hates and you dad thinks is BS")

	for {
		fmt.Print("Calculate > ")
		if !scanner.Scan() {
			break
		}

		input := strings.TrimSpace(scanner.Text())

		if input == "exit" {
			break
		}

		if input == "" {
			continue
		}

		result := EvaluateExpression(input)

		fmt.Printf("Result    = %v\n\n", result)
	}
}
