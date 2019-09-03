# Debugging IMPACT with CSSH

## Introduction

A CSSH- and gdb-based debugging tool for IMPACT.

## What's required to make it work within IMPACT

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

into impact.f90 under "INCLUDE 'mpif.h'" - this was realized in impact_debugger.f90. All kalman module subroutine calls are currently commented out as they are causing memory errors. In targets.mk make sure you are using the flags "-O0 -g".

## How to use it

### Preliminary setup

Make sure you're running an SSH server on the target machine, e.g. by using 

```
systemctl status ssh 
```

on Ubuntu 18.04. If it is not running, start it using

```
systemctl start ssh
``` 

### Running the debugger

Run

```
./debug_impact.sh
```

in a bash, which will launch the debugger version of IMPACT using MPI and start CSSH. Using CSSH's master terminal perform the following steps within each CSSH xterm

``` 
cd /path/to/your/IMPACT/debug
source xterm_attach.sh
```
 
which for each xterm will attach a gdb instance to a different MPI IMPACT process. Within gdb you then 

```
source gdb_init
```

which will unblock the sleeping loop and as of now stop execution before in the MPI_INIT statement of IMPACT's main routine. 

### Using gdb

A good summary is available from [here](https://cs.brown.edu/courses/cs033//docs/guides/gdb.pdf) from the [Brown CS course](https://cs.brown.edu/courses/cs033//).
