# Debugging IMPACT with CSSH

## Introduction

A CSSH- and gdb-based debugging tool for IMPACT.

## How to use it

Put 

```
+  INTEGER :: gdb
+
+  gdb = 0
+
+  do while (gdb == 0)
+    call sleep(2)
+  end do
```

into impact.f90 under "INCLUDE 'mpif.h'" (a seperate impact_debugger.f90 should be created for this). All kalman module subroutine calls are currently commented out as they are causing memory errors. In targets.mk make sure you are using the flags "-O0 -g".

Run

```
./debug_impact.sh
```

in a bash, which will launch impact and start CSSH. Within each CSSH xterm change to this directory and execute

``` 
source xterm_init.sh
```

after that run 

```
./xterm_sort_pids.sh
```
 
only in a single shell and after that again 

``` 
source xterm_attach.sh
```

in all CSSH xterms. Then a gdb instance will be attached to a different MPI process in each xterm. You can then 

```
up 2
set variable gdb = 1
```
which will unblock the sleeping loop and then before running `continue` define any breakpoints or others as usual.
