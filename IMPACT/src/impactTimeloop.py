# This is the main timeintegration loop in IMPACT. It is translated into python and calls wrapped Fortran routines.

#from __future__ import division
from mpi4py import MPI
from datetime import datetime
import impact

def timeintegration():
#--- start vector ---------------------------------------------------------------------------------
# already necessary here in order to be able to add more concentrations (BCs, however, are initialized below).
	impact.initial_conditions_vel()

	if impact.mod_vars.restart == 0:
		impact.mod_vars.time          = impact.mod_vars.time_start
		impact.mod_vars.time_out_vect = impact.mod_vars.time_start
		impact.mod_vars.time_out_scal = impact.mod_vars.time_start

		impact.mod_vars.dtime         = 0.0
		impact.mod_vars.timestep      = 0

		impact.mod_vars.new_dtime      = 1 # or true?
		impact.mod_vars.write_out_vect = 1 # or true?
		impact.mod_vars.write_out_scal = 1 # or true?

		impact.mod_vars.write_count = 0

		if impact.mod_vars.dtime_out_vect == 0.0:
			impact.mod_vars.write_out_vect = 0 # or false?
		if impact.mod_vars.dtime_out_scal == 0.0:
			impact.mod_vars.write_out_scal = 0 # or false?
	else:
		impact.mod_inout.read_restart()
		if impact.mod_vars.dtime_out_scal != 0.0:
			impact.read_restart_stats()
		if (impact.mod_vars.rank == 0) and (impact.mod_vars.write_stout_yes):
			print("             time = %13.5E" % (impact.mod_vars.time))
			print("         timestep = %5i" % (impact.mod_vars.timestep))
			
			print("            dtime = %13.5E" % (impact.mod_vars.dtime))
			print("    time_out_vect = %13.5E" % (impact.mod_vars.time_out_vect))
			print("    time_out_scal = %13.5E" % (impact.mod_vars.time_out_scal))

			print("        new_dtime = %s" % bool(impact.mod_vars.new_dtime))
			print("   write_out_vect = %s" % bool(impact.mod_vars.write_out_vect))
			print("   write_out_scal = %s" % bool(impact.mod_vars.write_out_scal))

			print("      write_count = %5i" % (impact.mod_vars.write_count))

	impact.mod_vars.timestep_old  = impact.mod_vars.timestep
	impact.mod_vars.dtime_old     = impact.mod_vars.dtime
	impact.mod_vars.dtime_average = 0.0
	impact.mod_vars.finish_yes    = False

	#--- Determine null-spaces -----------------------------------------------------------------------
	# Is here because the correction vector "th" needs allocation and definition only after "configuration"
	# (whereafter "initial_conditions_" are called subsequently, s.o.)
	# Alternative: create a subroutine for "th" ...
	if impact.mod_vars.nullspace_yes:
		impact.mod_coeffs.get_stencil_transp()
		impact.mod_solvers.get_nullspace()
		impact.mod_coeffs.get_stencil() # TEST!!! not nice! better pack into impact.f90 and change order ...
	
	#--- Initialize BCs ------------------------------------------------------------------------------
	impact.mod_lib.init_bc()

	#--- Test divergence-freeness --------------------------------------------------------------------
	impact.mod_test.test_divergence()

	#--- Open several files --------------------------------------------------------------------------
	impact.open_stats()

	#--- Recreate file for intermediate termination of timeintegration -------------------------------
	if impact.mod_vars.rank == 0:
		f = open('send_signal.txt','w')
		f.write("%i\n" % 0)
		f.write("%i\n" % impact.mod_vars.n_timesteps)
		f.write("%i"   % impact.mod_vars.time_end)
		f.flush()
		f.close()
	
	if (impact.mod_vars.rank == 0) and impact.mod_vars.write_stout_yes:
		print
		print("================================================================================")
		print("================================================================================")
		print("================================================================================")
		print
		print("--------------------------- START TIME-INTEGRATION -----------------------------")
		print
		print("                    Re = %12.5E" % (impact.mod_vars.re))
		print("         box resolution: %4i x %4i x %4i" % (impact.mod_dims.m1,impact.mod_dims.m2,impact.mod_dims.m3))
		print("         box dimension : %12.5E x %12.5E x %12.5E" % (impact.mod_vars.l1,impact.mod_vars.l2,impact.mod_vars.l3))
		print("                  epsU = %12.5E" % (impact.mod_vars.epsu))
		print("================================================================================")
		print("================================================================================")
		print("================================================================================")
		print
	
	#--- Start time measurement ----------------------------------------------------------------------
	if impact.mod_vars.rank == 0:
		ctime = datetime.now()
		impact.mod_vars.day  = ctime.day
		impact.mod_vars.hour = ctime.hour
		impact.mod_vars.minu = ctime.minute
		impact.mod_vars.sec  = ctime.second
		impact.mod_vars.msec = ctime.microsecond/1000.0 # there is no milisecond member

		impact.mod_vars.elatime = impact.mod_vars.msec + 1000*(impact.mod_vars.sec + 60*(impact.mod_vars.minu + 60*impact.mod_vars.hour))
		f99 = open("test_wallclocktime_restart%s.txt" % (''.join(impact.mod_vars.restart_char)),'w')
		f99.write("Begin time integration at %2i.%2i.%4i, %2i:%2i:%2i.%3i\n" % (ctime.day,ctime.month,ctime.year,ctime.hour,ctime.minute,ctime.second,ctime.microsecond/1000.0))
		f99.flush()
	
	#--- Write-out -----------------------------------------------------------------------------------
	if impact.mod_vars.write_out_scal:
		impact.compute_stats()
	if impact.mod_vars.write_out_vect:
		impact.mod_inout.write_fields()
	#=================================================================================================
	

	#=================================================================================================
	#=== Time integration ============================================================================
	#=================================================================================================
	while True:

		impact.mod_lib.get_dtime()
		
		if (impact.mod_vars.rank == 0) and impact.mod_vars.write_stout_yes:
			if impact.mod_vars.timeint_mode == 0:
				print
				print("================================================================================")
				print("================================================================================")
			print("time step = %8i ; time = %25.17E ; dtime = %25.17E" % (impact.mod_vars.timestep,impact.mod_vars.time,impact.mod_vars.dtime))
			if impact.mod_vars.timeint_mode == 0:
				print("================================================================================")
				print("================================================================================")
		
		if (impact.mod_vars.rank == 0) and impact.mod_vars.log_iteration_yes:
			f10 = open("log_iterations.txt",'w')
			f10.write("\n")
			f10.write("================================================================================\n")
			f10.write("================================================================================\n")
			f10.write("time step = %8i ; time = %25.17E ; dtime = %25.17E\n" % (impact.mod_vars.timestep,impact.mod_vars.time,impact.mod_vars.dtime))
			f10.write("================================================================================\n")
			f10.write("================================================================================\n")

		#=========================================================================================

		for impact.mod_vars.substep in range(1,impact.mod_vars.rk_steps+1): # +1 because Python excludes upper bound

			if (impact.mod_vars.rank == 0) and impact.mod_vars.write_stout_yes and (impact.mod_vars.timeint_mode == 0):
				print
				print("================================================================================")
				print("Runge-Kutta sub-step: %2i" % (impact.mod_vars.substep))
				print("================================================================================")
			if (impact.mod_vars.rank == 0) and impact.mod_vars.log_iteration_yes:
				f10.write("\n")
				f10.write("================================================================================\n")
				f10.write("Runge-Kutta sub-step: %2i\n" % (impact.mod_vars.substep))
				f10.write("================================================================================\n")
			
			#--- Time ------------------------------------------------------------------------
			if impact.mod_vars.substep == 1:
				impact.mod_vars.subtime = impact.mod_vars.time + impact.mod_vars.dtime* impact.mod_vars.ark[0] # (1) in fortran
			if impact.mod_vars.substep == 2:
				impact.mod_vars.subtime = impact.mod_vars.time + impact.mod_vars.dtime*(impact.mod_vars.ark[0]+impact.mod_vars.ark[1]+impact.mod_vars.brk[1])
			if impact.mod_vars.substep == 3:
				impact.mod_vars.subtime = impact.mod_vars.time + impact.mod_vars.dtime
					
			
			#--- Ghost cell update (for RHS) -------------------------------------------------
			impact.mod_exchange.exchange_all_all(True,impact.mod_vars.vel,impact.mod_vars.n1,impact.mod_vars.n2,impact.mod_vars.n3, \
							          impact.mod_vars.b1u,impact.mod_vars.b2u,impact.mod_vars.b3u, \
								  impact.mod_vars.b1l,impact.mod_vars.b2l,impact.mod_vars.b3l) # slightly different syntax


			
			#--- Interpolate advection velocity + update ghost cells -------------------------
			# vel(:,:,:,i) --> worki(:,:,:)
			impact.mod_diff.interpolate_vel(False) # TEST!!! Has partially already been done at time-step definition!

			#--- IBM (bbecsek 2015) ----------------------------------------------------------
			if impact.mod_vars.ib_on:
				if impact.mod_vars.timestep != 0:
					impact.mod_ibm.interpolate_vel_to_ib()
					if impact.mod_vars.rank == 0:
						impact.mod_ibm.update_boundary()
						impact.mod_ibm.calculate_displacements()
						impact.mod_ibm.fe_setup_triangular_3node()
				if impact.mod_vars.rank == 0:
					impact.mod_ibm.compute_force()

			#--- rhs (optionally override Neumann-BC) ----------------------------------------
			impact.mod_rhs.rhs_vel()


			#--- Helmholtz-multiplier --------------------------------------------------------
			impact.mod_vars.multl = impact.mod_vars.thetal \
					*(impact.mod_vars.ark[impact.mod_vars.substep-1]+impact.mod_vars.brk[impact.mod_vars.substep-1]) \
					*impact.mod_vars.dtime / impact.mod_vars.re

			#--- Rescale (efficiency, initial guess) -----------------------------------------
			if (not impact.mod_vars.init_pre[impact.mod_vars.substep-1]):
				impact.mod_vars.pre[(impact.mod_vars.s1p-1-impact.mod_vars.b1l+1):(impact.mod_vars.n1p-impact.mod_vars.b1l+1), \
						    (impact.mod_vars.s2p-1-impact.mod_vars.b2l+1):(impact.mod_vars.n2p-impact.mod_vars.b2l+1), \
						    (impact.mod_vars.s3p-1-impact.mod_vars.b3l+1):(impact.mod_vars.n3p-impact.mod_vars.b3l+1)] \
				= impact.mod_vars.pre[(impact.mod_vars.s1p-1-impact.mod_vars.b1l+1):(impact.mod_vars.n1p-impact.mod_vars.b1l+1), \
						      (impact.mod_vars.s2p-1-impact.mod_vars.b2l+1):(impact.mod_vars.n2p-impact.mod_vars.b2l+1), \
						      (impact.mod_vars.s3p-1-impact.mod_vars.b3l+1):(impact.mod_vars.n3p-impact.mod_vars.b3l+1)] \
				* (impact.mod_vars.ark[impact.mod_vars.substep-1] + impact.mod_vars.brk[impact.mod_vars.substep-1]) \
				* impact.mod_vars.dtime # upper bounds are not substracted one because python excludes upper bounds; account for lower negative index "bil" plus one for fortran 0 index

			#--- Solver ----------------------------------------------------------------------
			if (impact.mod_vars.timeint_mode == 1) or (impact.mod_vars.thetal == 1.0):
				impact.mod_solvers.explicit()
			else:
				if impact.mod_vars.twostep_yes:
					impact.mod_solvers.twostep()
				else:
					impact.mod_solvers.outer_iteration()

			#--- Physical pressure -----------------------------------------------------------
			impact.mod_vars.pre[(impact.mod_vars.s1p-1-impact.mod_vars.b1l+1):(impact.mod_vars.n1p-impact.mod_vars.b1l+1), \
					    (impact.mod_vars.s2p-1-impact.mod_vars.b2l+1):(impact.mod_vars.n2p-impact.mod_vars.b2l+1), \
					    (impact.mod_vars.s3p-1-impact.mod_vars.b3l+1):(impact.mod_vars.n3p-impact.mod_vars.b3l+1)] \
			= impact.mod_vars.pre[(impact.mod_vars.s1p-1-impact.mod_vars.b1l+1):(impact.mod_vars.n1p-impact.mod_vars.b1l+1), \
					      (impact.mod_vars.s2p-1-impact.mod_vars.b2l+1):(impact.mod_vars.n2p-impact.mod_vars.b2l+1), \
					      (impact.mod_vars.s2p-1-impact.mod_vars.b3l+1):(impact.mod_vars.n3p-impact.mod_vars.b3l+1)] \
			/ (impact.mod_vars.ark[impact.mod_vars.substep-1] + impact.mod_vars.brk[impact.mod_vars.substep-1]) \
			/ impact.mod_vars.dtime # upper bounds are not substracted one because python excludes upper bounds; account for lower negative index "bil" plus one for fortran 0 index

			#--- Fill undefined corners / edges ----------------------------------------------
			impact.mod_lib.fill_corners(impact.mod_vars.pre,impact.mod_vars.n1,impact.mod_vars.n2,impact.mod_vars.n3, \
							          impact.mod_vars.b1u,impact.mod_vars.b2u,impact.mod_vars.b3u, \
								  impact.mod_vars.b1l,impact.mod_vars.b2l,impact.mod_vars.b3l) # slightly different syntax

		#--- for loop finish
		#=========================================================================================
		impact.mod_vars.timestep = impact.mod_vars.timestep + 1
		impact.mod_vars.time     = impact.mod_vars.time + impact.mod_vars.dtime

		#--- Read send_signal.txt ----------------------------------------------------------------
		impact.mod_lib.check_signal()


		#--- Fix pressure level ------------------------------------------------------------------
		impact.mod_lib.level_pressure()


		#--- Write out ---------------------------------------------------------------------------
		if impact.mod_vars.write_out_scal:
			impact.compute_stats()
		if impact.mod_vars.write_out_vect:
			impact.mod_inout.write_fields()


		#-----------------------------------------------------------------------------------------
		if (impact.mod_vars.rank == 0) and impact.mod_vars.log_iteration_yes:
			f10.close()

		if (impact.mod_vars.rank == 0) and impact.mod_vars.write_stout_yes and (impact.mod_vars.timeint_mode == 0):
			print
		if (impact.mod_vars.rank == 0) and impact.mod_vars.write_stout_yes and (impact.mod_vars.timeint_mode == 0): # why twice?
			print
		
		impact.usr_func.mpi_bcast_fort() # This is in fortran because all the MPI communication is done there
						 # but it could be moved to Python in the future.

		if impact.mod_vars.finish_yes:
			break
	#--- for loop finish
	#=================================================================================================


	#--- Terminate time measurement ------------------------------------------------------------------
	if impact.mod_vars.rank == 0:
		ctime = datetime.now()
		impact.mod_vars.hour = ctime.hour
		impact.mod_vars.minu = ctime.minute
		impact.mod_vars.sec  = ctime.second
		impact.mod_vars.msec = ctime.microsecond/1000.0

		if ctime.day != impact.mod_vars.day:
			# Remark: only valid for jobs <= 24h
			impact.mod_vars.elatime = impact.mod_vars.msec + 1000*(impact.mod_vars.sec + 60*(impact.mod_vars.minu + 60*impact.mod_vars.hour)) - impact.mod_vars.elatime + 24*60*60*1000
		else:
			impact.mod_vars.elatime = impact.mod_vars.msec + 1000*(impact.mod_vars.sec + 60*(impact.mod_vars.minu + 60*impact.mod_vars.hour)) - impact.mod_vars.elatime

		f99.write("Finish time integration at %2i.%2i.%4i, %2i:%2i:%2i:%3i\n" % (ctime.day,ctime.month,ctime.year,ctime.hour,ctime.minute,ctime.second,ctime.microsecond/1000.0))
		f99.write("elapsed time [sec] %13.5E" % float(impact.mod_vars.elatime/1000.0))
		f99.close()
	
	#--- Write restart -------------------------------------------------------------------------------
	impact.mod_vars.restart = impact.mod_vars.restart + 1
	impact.mod_inout.write_restart()
	impact.write_restart_stats()


	#--- Evaluate iteration stats --------------------------------------------------------------------
	impact.mod_lib.iteration_stats()

	#--- Close several files -------------------------------------------------------------------------
	impact.close_stats()


	return
