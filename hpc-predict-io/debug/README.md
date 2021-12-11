# Debugging hpc-predict-io with CSSH

## Introduction

A CSSH- and gdb-based debugging tool for hpc-predict-io. This launches hpc-predict-io over MPI and then creates a parallel set of gdb sessions each attached to a different MPI rank that are centrally controlled through the CSSH master terminal.

## How to use it

### Preliminary setup

Make sure you're running an SSH server on the target machine. This requires installing the package first, e.g.

```
sudo apt-get install ssh-server
```

on Ubuntu 18.04. Check if it is running, e.g. by using 

```
systemctl status ssh 
```

on Ubuntu 18.04. If it is not running, start it using

```
systemctl start ssh
``` 

### Running the debugger

Run the shell script

```
./debug_hpc_predict_io.sh
```

which will launch the debugger version of hpc-predict-io (see below) using MPI and start CSSH. Using CSSH's master terminal (on which the cursor is focused by default) perform the following steps (automatically executed within each CSSH xterm)

``` 
cd /path/to/your/hpc_predict_io/debug
source xterm_attach.sh
```
 
which for each xterm will attach a gdb instance to a different MPI hpc-predict-io process. Within gdb you then 

```
source gdb_init
```

which will unblock the MPI processes (from the sleeping loop) and stop execution at the MPI_INIT statement of fortran/test/mr_io_test_impact_mri's main routine. 

### Using gdb

A good summary is available from [here](https://cs.brown.edu/courses/cs033//docs/guides/gdb.pdf) from the [Brown CS course](https://cs.brown.edu/courses/cs033//).

### Terminating

Quit by just killing the process that runs `debug_hpc_predict_io.sh` - this will automatically terminate CSSH and clean up hpc-predict-io's MPI processes.

## Implementation details: Debugger version of hpc-predict-io

In fortran/test/mr_io_test_impact_mri.f90, the lines

```
  INTEGER :: gdb

  gdb = 0

  do while (gdb == 0)
    call sleep(2)
  end do
```

were added to the executable to let the processes spin until a gdb processes from an xterm attaches (in xterm_attach.sh) and sets the value of gdb /= 0 (in gdb_init). The program then continues until the breakpoint defined in gdb_init. All kalman module subroutine calls are currently commented out as they are causing memory errors. Whenever debugging with gdb, make sure that in targets.mk you are using the flags "-O0 -g".

## Known issues

The ptrace system call can be deactivated by default on Ubuntu, disabling gdb to attach to a running process. This can be undone by setting

```
kernel.yama.ptrace_scope = 0
```

in the file /etc/sysctl.d/10-ptrace.conf and running `sysctl --system` to reload configuration files. Use `sysctl -a` to display the currently active configuration. 
