all: compile_cpu compile_gpu

compile_cpu:
	pgcc  main_br.c CurrentNa.c CurrentK.c CurrentX.c CurrentS.c GateVarCommonForm.c -o exec_br_cpu

compile_gpu:
	pgcc -acc -fast -Minfo=accel -ta=nvidia main_br.c CurrentNa.c CurrentK.c CurrentX.c CurrentS.c GateVarCommonForm.c -o exec_br_gpu

#run:
#	./exec_br_gpu

clear_output:
	rm output/*.vtk
