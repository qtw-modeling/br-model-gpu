all: compile_cpu compile_gpu

compile_cpu:
	pgcc main_br.c CurrentNa.c CurrentK.c CurrentX.c CurrentS.c GateVarCommonForm.c extra.c -o exec_br_cpu

compile_gpu:
	pgcc -acc -fast -Minfo=accel -ta=tesla main_br.c CurrentNa.c CurrentK.c CurrentX.c CurrentS.c GateVarCommonForm.c extra.c -o exec_br_gpu

phase:
	pgcc  WritePhaseDataIntoFile.c CurrentNa.c CurrentK.c CurrentX.c CurrentS.c GateVarCommonForm.c -o exec_br_phase

clear_output:
	rm output/*.vtk
