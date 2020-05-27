#all: compile_cpu run

#compile_cpu:
#	g++ -std=c++11 -g main_br.cpp CurrentNa.cpp CurrentK.cpp CurrentX.cpp CurrentS.cpp GateVarCommonForm.cpp -o exe_br

compile_gpu:
	pgcc -acc -Minfo=accel -ta=nvidia main_br.c CurrentNa.c CurrentK.c CurrentX.c CurrentS.c GateVarCommonForm.c -o exec_br_gpu

run:
	./exe_br_gpu

clear_output:
	rm output/*.vtk
