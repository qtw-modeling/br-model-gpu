all: compile_cpu run

compile_cpu:
	g++ -std=c++11 -g main_br.cpp CurrentNa.cpp CurrentK.cpp CurrentX.cpp CurrentS.cpp GateVarCommonForm.cpp -o exe_br

#compile_gpu:
#	pgc++ main_yni.cpp -acc -Minfo=accel -ta=nvidia -o yni_exec

run:
	./exe_br

clear_output:
	rm output/*.vtk