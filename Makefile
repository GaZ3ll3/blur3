all: 
	cd libs/qhull && make && cd ../meshgencpp && make && cd ../../src && make

clean: cleanlibs cleansrc

cleanlibs:
	cd libs && make clean
	
cleansrc:
	cd src && make clean
