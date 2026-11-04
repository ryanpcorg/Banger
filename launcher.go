package main

import (
	"bufio"
	"fmt"
	"os"
	"os/exec"
	"runtime"
	"strconv"
	"strings"
	"time"
)

func printBanner() {
	println("██████╗  █████╗ ███╗   ██╗ ██████╗ ███████╗██████╗ ")
	println("██╔══██╗██╔══██╗████╗  ██║██╔════╝ ██╔════╝██╔══██╗")
	println("██████╔╝███████║██╔██╗ ██║██║  ███╗█████╗  ██████╔╝")
	println("██╔══██╗██╔══██║██║╚██╗██║██║   ██║██╔══╝  ██╔══██╗")
	println("██████╔╝██║  ██║██║ ╚████║╚██████╔╝███████╗██║  ██║")
	println("╚═════╝ ╚═╝  ╚═╝╚═╝  ╚═══╝ ╚═════╝ ╚══════╝╚═╝  ╚═╝")
	println("---------------------------------------------------")
	println("The calculator your mum hates and you dad thinks is BS (Ballistically Stupid)")
	println("Commands: 'clear' or 'cls' to clear screen | 'plot <expr> [min] [max]' | 'exit' to quit\n")
}

func clearTerminal() {
	if runtime.GOOS == "windows" {
		cmd := exec.Command("cmd", "/c", "cls")
		cmd.Stdout = os.Stdout
		_ = cmd.Run()
	} else {
		fmt.Print("\033[H\033[2J")
	}
}

func main() {
	scanner := bufio.NewScanner(os.Stdin)
	time.Sleep(100 * time.Millisecond)
	printBanner()

	for {
		fmt.Print("Calculate > ")
		if !scanner.Scan() {
			break
		}

		input := strings.TrimSpace(scanner.Text())

		if input == "exit" || input == "quit" {
			break
		}

		if input == "clear" || input == "cls" {
			clearTerminal()
			printBanner()
			continue
		}

		if input == "" {
			continue
		}

		if strings.HasPrefix(strings.ToLower(input), "plot ") || strings.HasPrefix(strings.ToLower(input), "graph ") {
			parts := strings.Fields(input)
			if len(parts) >= 2 {
				expr := parts[1]
				xMin := -10.0
				xMax := 10.0
				if len(parts) >= 3 {
					if v, err := strconv.ParseFloat(parts[2], 64); err == nil {
						xMin = v
					}
				}
				if len(parts) >= 4 {
					if v, err := strconv.ParseFloat(parts[3], 64); err == nil {
						xMax = v
					}
				}
				graphAscii := RenderGraph(expr, xMin, xMax, 60, 20)
				fmt.Print(graphAscii)
				fmt.Println()
				continue
			}
		}

		result := EvaluateExpression(input)
		fmt.Printf("Result    = %v\n\n", result)
	}
}
