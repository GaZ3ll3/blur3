all: setenv bash lib lib2d lib2d_sphere lib3d

setenv:
	python util/setenv.py 
	
bash: setenv
	bash setenv.bash
lib:
	cd lib && make 
	
lib2d:
	cd lib/2d && make
	
lib2d_sphere:
	cd lib/2d_sphere && make
	
lib3d:
	cd lib/3d && make
	
clean: cleanlib cleanlib2d cleanlib3d cleanlib2d_sphere cleanenv

cleanlib:
	cd lib && make clean
cleanlib2d:
	cd lib/2d && make clean
	
cleanlib2d_sphere:
	cd lib/2d_sphere && make clean
	
cleanlib3d:
	cd lib/3d && make clean

cleanenv:
	rm -f setenv.*
