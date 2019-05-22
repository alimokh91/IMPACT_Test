#! /usr/bin/python
#**************************************************************************************************
#* PyMPACT                                                                                        *
#* by Barna Becsek, ARTORG Center, University of Bern, (barna.becsek@artorg.unibe.ch)             *
#* October 2015                                                                                   *
#**************************************************************************************************


#from __future__ import division
from mpi4py import MPI #MPI_init() is called upon import
#import h5py as HDF5 # This is still done in Fortran
from datetime import datetime
import impact 
import impactTimeloop

#==================================================================================================
#=== Initialization ===============================================================================
#==================================================================================================

#--- Initialize MPI -------------------------------------------------------------------------------
impact.mod_vars.rank = MPI.COMM_WORLD.Get_rank()
# DEBUGGING:
#MPI.COMM_WORLD.Set_name('Test1')
#impact.usr_func.print_fcomm_size() # This routine ouputs the name set above. Thus, it's the same commi.

#--- Initialize HDF5 ------------------------------------------------------------------------------
impact.usr_func.init_hdf5()

#--- Set alarm if queue.txt is being used ---------------------------------------------------------
impact.mod_lib.init_alarm()

#--- Read and set configuration / topology --------------------------------------------------------
impact.configuration()

#--- Test of input parameters ---------------------------------------------------------------------
impact.mod_test.test_parameter()

#--- Setup ----------------------------------------------------------------------------------------
 #--- General ----------------------------------------------------------------------------------
impact.mod_setup.init_general()

 #--- MPI --------------------------------------------------------------------------------------
impact.mod_setup.init_parallel()

 #--- Type of boundary conditions --------------------------------------------------------------
impact.mod_setup.init_boundaries()

 #--- Limits of indices ------------------------------------------------------------------------
impact.mod_setup.init_limits()

#--- Physical coordinates -------------------------------------------------------------------------
impact.mod_geometry.coordinates()

#--- Determine differential coefficients ----------------------------------------------------------
impact.mod_coeffs.fd_coeffs()

#--- Get stencil of operator div(grad( )), order of convergence is 2 ------------------------------
impact.mod_coeffs.get_stencil()
impact.mod_coeffs.get_stencil_helm()

#--- Get interpolation coefficients (multigrid, order of convergence is 2) ------------------------
impact.mod_coeffs.interp_coeffs()
impact.mod_coeffs.interp_coeffs_helm()
impact.mod_coeffs.restr_coeffs()
impact.mod_coeffs.restr_coeffs_helm()

#--- Test differential coefficients ---------------------------------------------------------------
if impact.mod_vars.write_test_yes:
	impact.mod_coeffs.test_coeffs()

#--- Weights for stop criterion -------------------------------------------------------------------
impact.mod_coeffs.get_weights()

#--- Beta (for analysis) --------------------------------------------------------------------------
impact.mod_lib.get_beta()

#==================================================================================================


#==================================================================================================
#=== Main task ====================================================================================
#==================================================================================================
#--- DNS / time integration -----------------------------------------------------------------------
# Todo: disassemble this routine, so that the loop is done inside Python. This will then yield the
#       desired interface to couple with other functions. OR: provide arguments to timeintegration()
#       which e.g. tell it how far to advance in time (and from where).
if impact.mod_vars.task == 1:
	impactTimeloop.timeintegration()
#--- read and average fields ----------------------------------------------------------------------
if impact.mod_vars.task == 2:
	impact.postprocess()
#--- solve eigenvalue problem --- (task == 3, now deprecated / removed) ---------------------------
if impact.mod_vars.task == 3:
	if impact.mod_vars.rank == 0:
		print 'Task 3 is no longer available. Please read information concerning update on 06 May 2013.'

#--- analyze matrices -----------------------------------------------------------------------------
if impact.mod_vars.task == 4:
	analyze_matrix(0)

#==================================================================================================

#==================================================================================================
#=== Finalize =====================================================================================
#==================================================================================================
if impact.mod_vars.rank == 0:
	print ''
	print 'DONE ...'
	print ''

#--- close HDF5 -----------------------------------------------------------------------------------
impact.usr_func.finl_hdf5()

#--- close MPI ------------------------------------------------------------------------------------
#impact.usr_func.finl_mpi()
#==================================================================================================





