#include "main.h"
#include <math.h>
#include "measures.h"
#include "robust_filesystem.h"
#include "initialization.h"
#include "montecarlo.h"
#include <fstream>

namespace fs = std::filesystem;
size_t L, n_steps;
double T;
bool restart;


// Function to generate random double between min and max
double randomDouble(double min, double max) {
    return min + (max - min) * (double)rand() / RAND_MAX;
}


int main(int argc, char *argv[]) {

    srand(static_cast<unsigned int>(time(NULL)));
    std::vector <Node> Lattice;
    struct H_parameters Hp{};
    struct MC_parameters MC{};
    int NSTART=0;
    int N;



    std::string directory_read;
    std::string directory_write;


    if(argc == 14) {
        L = static_cast<size_t>(std::atoi(argv[1]));
        MC.n_steps = static_cast<size_t>(std::atoi(argv[2]));
        MC.transient = static_cast<size_t>(std::atoi(argv[3]));
        MC.tau = static_cast<size_t>(std::atoi(argv[4]));
        T = std::atof(argv[5]);
        restart = std::atof(argv[6]);
        Hp.K = std::atof (argv[7]);
        Hp.J1 = std::atof (argv[8]);
        Hp.J2 = std::atof (argv[9]);
        Hp.e = std::atof(argv[10]);
        MC.theta_box= std::atof (argv[11]);
        paths_dir::DIR_IN = directory_read = argv[12];
        paths_dir::DIR_OUT = directory_write = argv[13];
    }
    else{
        myhelp(argc, argv);
    }

    N=L*L;

    std::cout<< L << std::endl;
    std::cout<<"Numero di step "<< n_steps << std::endl;
    std::cout<< T << std::endl;
    std::cout<< Hp.K << std::endl;
    std::cout<< Hp.J1 << std::endl;
    std::cout<< Hp.J2 << std::endl;
    std::cout<< Hp.e <<std::endl;
    std::cout<< restart << std::endl;
    std::cout<<MC.theta_box<<std::endl;
    std::cout<< directory_read << std::endl;
    std::cout<< directory_write << std::endl;

    //Declaration of lattice structure
    Lattice.resize(N);
    //Initialization of the Lattice
    initialize_Lattice (Lattice, directory_read, directory_write, restart, Hp);

    // Run simulation --->>
    mainloop (Lattice, T, MC, N, Hp, directory_write, NSTART);

    save_lattice(Lattice, directory_write, std::string("Lattice_dir_fin.txt"));

    return 0;
}

void mainloop(std::vector <Node> &Site, double T, struct MC_parameters &MC, size_t N, struct H_parameters &Hp, std::string directory_write, int NSTART) {

    std::vector<double> magnetization_values;
    Measures mis;


    std::string Filename_magnetization1=(directory_write+"/Single_Magnetization.txt");
    std::string Filename_energy=(directory_write+"/Energy.txt");
    std::string Filename_helicity1=(directory_write+"/Helicity_modulus1.txt");
    std::string Filename_helicity2=(directory_write+"/Helicity_modulus2.txt");
    std::string Filename_trsb_magn=(directory_write+"/trsb_magnetization.txt");
    //std::string Filename_vortices1=(directory_write+"/Vortices.txt");
    //std::string Filename_vortices2=(directory_write+"/Vortices.txt");

    std::ofstream File_Magetization1 (Filename_magnetization1);
    std::ofstream File_Energy (Filename_energy);
    std::ofstream File_helicity1 (Filename_helicity1);
    std::ofstream File_helicity2 (Filename_helicity2);
    std::ofstream File_trsb_magn (Filename_trsb_magn);
    //std::ofstream File_vortices1 (Filename_vortices1);
    //std::ofstream File_vortices2 (Filename_vortices2);

    if(NSTART==0){
        /**Thermalization**/
        std::cout<< "Thermalization" << std::endl;
        for (int t = 0; t < MC.transient; t++) {
            metropolis(Site, MC, Hp, T);
        }
    }

    for (int nM = NSTART; nM<MC.n_steps; nM++) {

        for (int t = 0; t < MC.tau; t++) {
            metropolis(Site, MC, Hp,  T);
        }

        //Measures
        mis.reset();

        energy(mis, Hp, Site);
        single_magnetization(Site, mis, N);   //In this function we are calculating the magnetization of both layers separately
        trsb_magnetization(mis, Site, N);
        helicity_modulus(Hp, Site, mis, N);   //Both lattices
        //vortex(Site, mis, N) ;

        File_Energy        << mis.E << std::endl;
        File_Magetization1 << mis.m_phase[0] << " " << mis.m_phase[1] << std::endl;
        File_helicity1     << mis.Jd[0] << " " << mis.Ic[0] << std::endl;
        File_helicity2     << mis.Jd[1] << " " << mis.Ic[1] << std::endl;
        File_trsb_magn     << mis.trsb_m << std::endl;
        //File_vortices1     <<mis.vortices[0] << " " << mis.antivortices[0] <<std::endl;
        //File_vortices2     <<mis.vortices[1] <<" " << mis.antivortices[1] <<std::endl;

    }

    File_Energy.close();
    File_Magetization1.close();
    File_helicity1.close();
    File_helicity2.close();
    File_trsb_magn.close();
    //File_vortices1.close();
    //File_vortices2.close();
}


void myhelp(int argd, char **argu) {
    int i;
    fprintf(stderr,"Errore nei parametri su linea di comando; hai scritto:\n");
    for (i=0;i<argd;i++) fprintf(stderr," %s",argu[i]);
    fprintf(stderr,"\n");
    fprintf(stderr,"%s <L> <n_steps> <T> <restart> <K> <J1> <J2> <DIR_IN> <DIR_OUT> \n",argu[0]);
    exit (EXIT_FAILURE);
}

