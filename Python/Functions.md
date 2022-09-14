---
title: Functions
author: Federico Cesari
tags: [Doc]
date: 15.09.2022
---
# Functions
Functions *name* a peace of code like a variable does with strings, numbers... and they *take arguments* in a similar way of [[Inputs#sys argv|argv]]. With these two abilities function are able to define **mini-scripts**.
A function is declared using `def` keyword and it is concluded with dedentation.

In this case our `function` has no declared arguments. `def function("Arg")`.
```python
def function():
	print("This is my function!")

function()
```

Instead in this case I specify an argument.
```python
def print_my_name(my_name):
	print(f"My name is {my_name}.")

print_my_name("otto")
```
