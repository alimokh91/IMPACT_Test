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

in a bash, which will launch IMPACT using MPI and start CSSH. Using CSSH's master terminal perform the following steps within each CSSH xterm

``` 
cd /path/to/your/IMPACT/debug
source xterm_attach.sh
```
 
which for each xterm will attach a gdb instance to a different MPI IMPACT process. Within gdb you then 

```
source gdb_init
```

which will unblock the sleeping loop and as of now stop execution before in the MPI_INIT statement of IMPACT's main routine. 
