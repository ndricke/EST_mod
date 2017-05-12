# Basic binary test printing
from termcolor import colored

def print_result(conditional):
    print colored("PASS", "green") if conditional else colored("FAIL", "red")
