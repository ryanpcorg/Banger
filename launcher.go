package main

import (
	"fmt"
	"os"
	"os/exec"
	"runtime"
)

var filename string

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

func runthego() {
	cmd := exec.Cmd("go", "run", filename)
}

func main() {

	fmt.Println("What option do you want to do today?")
	fmt.Println("1: Multiply")
	fmt.Scan(&filename)

}
