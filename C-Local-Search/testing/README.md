# How to test your problem

- Edit the file mk/problem2.mk to indicate which API functions you
have implemented:
    - Check that all variables are correctly set (the defaults should just
    work):
        - `PROBPATH`: path to the source code of your problem
        - `PROBOBJ`: name of the `.o` file
        - `PROBINIT`: path to the problem init file (e.g.
          `problem_init/default_init.c`)
        - Set all the other variables to `1` or `0` to indicate whether the
          function with the same name should be tested or not, respectively.

- Run the Makefile indicating the problem name (must be the name given in the
`.mk` file), e.g.: 
```
make PROBNAME=problem2
```
This will create an executable file named `test-problem2`.

- Run the tests (and provide any additional arguments required, as defined in
`initProblem()`):
```
`./test-fakeproblem
```

