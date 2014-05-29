/*
 * Media.cpp
 *
 *  Created on: May 28, 2014
 *      Author: lurker
 */

#include "Media.h"

Media::Media(SpatialMesh& Space, AngularMesh& Angle) {
	// TODO Auto-generated constructor stub
	ua.resize(Space.nt);
	us.resize(Space.nt);
	for (register int i = 0; i < Space.nt; i++){
		// linear  DG
		ua[i].resize(3);
		us[i].resize(3);
	}
	std::fstream uafile("./ua.txt", std::ios_base::in);
	std::fstream usfile("./us.txt", std::ios_base::in);
	for (register int i = 0; i < Space.nt; i++){
		for (register int j = 0 ; j < 3; j++){
			uafile >> ua[i][j];
			usfile >> us[i][j];
		}

	}
	uafile.close();
	usfile.close();

}

Media::~Media() {
	// TODO Auto-generated destructor stub
}

void Media::Show(){
	for (std::size_t i = 0 ; i < ua.size() ; i++){
		for (int j = 0; j < 3 ; j++){
			std::cout << "ua[" << i << "][" << j << "] = " << ua[i][j] << std::endl;
			std::cout << "us[" << i << "][" << j << "] = " << us[i][j] << std::endl;
		}
	}
}
