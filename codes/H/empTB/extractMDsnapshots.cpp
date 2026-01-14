#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>

using namespace std;

int main(int argc, char* argv[]){                         // extracts posistions from MD and writes it into outfile   


	// defining variables

	int atoms_per_unitcell = 12;
		
	ofstream fout;
	string folder1 = argv[1];
	string folder2 = argv[2];
	int num_unitcells = stoi(argv[3]);
	int first_snapshot = stoi(argv[4]);
	int last_snapshot = stoi(argv[5]);
	int N_snapshots = stoi(argv[6]);
	int num_atoms = atoms_per_unitcell * (num_unitcells * num_unitcells * num_unitcells);
	string infile = folder1 + "/position.lammpstrj", outfile = folder2 + "/snapshots/traj", temp;
	//infile.at(17) = outfile.at(14) = argv[1][1];            // name of the dum

	ifstream fin (infile.c_str());

	cout << infile << "\n";
	cout << outfile << "\n";

	double data[3];

	//while (frame < timesteps + 1){
	//for (int frame=first_snapshot; frame<=last_snapshot; frame++){
	int frame = first_snapshot;
	for (int i = 0; i < N_snapshots; ++i) {
		if (N_snapshots > 1) {
    		frame = first_snapshot + round(i * double(last_snapshot - first_snapshot)/(N_snapshots - 1));
		}


		for (int i=0; i<9; ++i){
			getline(fin,temp);                           // getline reads a line from inputfile fin and stores it into temp
			//if (frame == 0){
			//	cout << temp << '\n';
			//}
		}
			
		string outfile1 = outfile + to_string(frame) + ".xyz";
		//outfile.at(6) = frame/1000 + 48;
		//outfile.at(7) = (frame%1000)/100 + 48;
		//outfile.at(8) = (frame%100)/10 + 48;
		//outfile.at(9) = frame%10 + 48;



		fout.open(outfile1.c_str());
		for (int i=0; i<num_atoms; ++i){     // 49152 = 12*16^3    12 atomes per unit cell and 16^3 unit cells
			fin >> temp >> temp >> temp >> data[0] >> data[1] >> data[2];
			fout << data[0] << " " << data[1] << " " << data[2] << '\n';
			}
		getline(fin,temp);
		fout.close();

		}
	
	return 0;
	}
