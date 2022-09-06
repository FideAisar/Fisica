---
title: Managing text files
author: Federico Cesari
tags: [Doc]
date: 03.09.2022
---
# Text files
## Open function
	Syntax: open(filename, 'mode')
```python
open(filename, 'r') # Opens a file for reading
open(filename, 'a') # Opens a file for appending
open(filename, 'w') # Opens a file for writing
open(filename, 'x') # Creates the specified file

open(filename) # Opens the file in read mode by default
```

## Methods
	Syntax: open(filename, 'mode').method(param)

| Name     | Function                                                   | Param                                                         |
| -------- | ---------------------------------------------------------- | ------------------------------------------------------------- |
| close    | close the file                                             | ()                                                            |
| read     | read the file                                              | ('bite to read of the file [-1 by default = the whole file]') |
| readline | read a specific line                                       | ()                                                            |
| truncate | empties the file                                           | ()                                                             |
| write    | write "stuff" in the file                                  | ('stuff')                                                              |
| seek     | Moves the read/write location to the beginning of the file | (0)                                                              |
- - - 
*ex:*
- method.py
- test.txt

```python
open("test.txt", 'w').write("Line 1")
```
terminal:
```bash
> python method.py
```
*test.txt*
```txt
Line 1
```
- - -