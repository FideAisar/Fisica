---
title: Inputs
author: Federico Cesari
tags: [Doc]
date: 04.09.2022
---

# Different input methods
## Input
Asks the user for an input. It accepts strings only.
```python
x = input()
print("...")
print(x)
```

*terminal:*
```bash
> this is my input
> ...
> this is my input
```
## sys.argv
`argv` stands for *argument variable*, it consists in a list of strings which contain the command-line arguments passed to the script. To use `sys.argv`, you will first have to import the sys module.
```python
from sys import argv



```
