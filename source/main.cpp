#include "main.h"
#include <cmath>
#include "measures.h"
#include "initialization.h"
#include "montecarlo.h"
#include "parallel_temp.h"
#include <fstream>
#include <iostream>

namespace fs = std::filesystem;
size_t L;
double T;
bool restart;

int main(int argc, char *argv[]) {

    //srand(static_cast<unsigned int>(time(NULL)));
    std::vector <Node> Lattice;
    struct H_parameters Hp{};
    struct MC_parameters MC{};
    struct PT_parameters PTp;
    struct PTroot_parameters PTroot;
    double my_beta=0.244;
    int my_ind=0;
    int NSTART=0;
    size_t N;

    std::string directory_read;
    std::string directory_write;


    if(argc == 17) {
        L = static_cast<size_t>(std::strtol(argv[1], nullptr,  10));
        MC.n_steps = static_cast<int>(std::strtol(argv[2], nullptr, 10));
        MC.transient = static_cast<int>(std::strtol(argv[3], nullptr, 10));
        MC.tau = static_cast<int>(std::strtol(argv[4], nullptr, 10));
        T = std::strtod(argv[5], nullptr);
        restart = std::strtol(argv[6], nullptr, 10);
        Hp.K = std::strtod(argv[7], nullptr);
        Hp.J1 = std::strtod(argv[8], nullptr);
        Hp.J2 = std::strtod(argv[9], nullptr);
        Hp.e = std::strtod(argv[10], nullptr);
        Hp.beta_high = std::strtod(argv[11], nullptr);
        Hp.beta_low = std::strtod(argv[12], nullptr);
        MC.theta_box= std::strtod(argv[13], nullptr);
        MC.theta_box_A = std::strtod(argv[14], nullptr);
        paths_dir::DIR_IN = directory_read = argv[15];
        paths_dir::DIR_OUT = directory_write = argv[16];
    }
    else{
        myhelp(argc, argv);
    }

    N=L*L;

    std::cout<< L << std::endl;
    std::cout<<"Numero di step "<< MC.n_steps << std::endl;
    std::cout<< T << std::endl;
    std::cout<< Hp.K << std::endl;
    std::cout<< Hp.J1 << std::endl;
    std::cout<< Hp.J2 << std::endl;
    std::cout<< Hp.e <<std::endl;
    std::cout << Hp.beta_high <<std::endl;
    std::cout << Hp.beta_low <<std::endl;
    std::cout<< restart << std::endl;
    std::cout<<MC.theta_box<<std::endl;
    std::cout<<MC.theta_box_A<<std::endl;
    std::cout<< directory_read << std::endl;
    std::cout<< directory_write << std::endl;

    //Declaration of lattice structure
    Lattice.resize(N);

    MPI_Init(nullptr, nullptr); /* START MPI */
/*DETERMINE RANK OF THIS PROCESSOR*/
    MPI_Comm_rank(MPI_COMM_WORLD, &PTp.rank);
/*DETERMINE TOTAL NUMBER OF PROCESSORS*/
    MPI_Comm_size(MPI_COMM_WORLD, &PTp.np);

    if(PTp.rank == PTp.root) {
        //Initialization ranks arrays
        initialize_PTarrays( PTp, PTroot, Hp);
    }
    MPI_Scatter(PTroot.beta.data(), 1, MPI_DOUBLE, &my_beta, 1, MPI_DOUBLE, PTp.root, MPI_COMM_WORLD);
    MPI_Scatter(PTroot.rank_to_ind.data(), 1, MPI_INT, &my_ind, 1, MPI_INT, PTp.root, MPI_COMM_WORLD);

    printf("I'm rank %d and this is my beta %lf\n", PTp.rank, my_beta);

    //Initialization of the Lattice
    initialize_Lattice (Lattice, directory_read, directory_write, restart, Hp );



    // Run simulation --->>
    mainloop (Lattice, T, MC, my_ind, my_beta, PTp, PTroot, N, Hp, directory_write, NSTART);

    save_lattice(Lattice, directory_write, std::string("Lattice_dir_fin.txt"), std::string("A_dir_fin.txt"), Hp);

    return 0;
}

void mainloop(std::vector <Node> &Site, double T_temp, struct MC_parameters &MC, int &my_ind, double &my_beta, struct PT_parameters PTp, struct PTroot_parameters PTroot, size_t N, struct H_parameters &Hp, const std::string directory_write, int NSTART) {

    std::vector<double> magnetization_values;
    Measures mis;


    std::string Filename_magnetization1=(directory_write+"/Single_Magnetization.txt");
    std::string Filename_energy=(directory_write+"/Energy.txt");
    std::string Filename_kin_energy=(directory_write+"/Kin_Energy.txt");
    std::string Filename_joseph_energy=(directory_write+"/Joseph_Energy.txt");
    std::string Filename_B_energy=(directory_write+"/B_Energy.txt");
    std::string Filename_helicity1=(directory_write+"/Helicity_modulus1.txt");
    std::string Filename_helicity2=(directory_write+"/Helicity_modulus2.txt");
    std::string Filename_trsb_magn=(directory_write+"/trsb_magnetization.txt");
    //std::string Filename_vortices1=(directory_write+"/Vortices.txt");
    //std::string Filename_vortices2=(directory_write+"/Vortices.txt");
    std::string Filename_dual_stiff=(directory_write+"/Dual_Stiffness.txt");


    std::ofstream File_Magetization1 (Filename_magnetization1);
    std::ofstream File_Energy (Filename_energy);
    std::ofstream File_Kin_Energy (Filename_kin_energy);
    std::ofstream File_Joseph_Energy (Filename_joseph_energy);
    std::ofstream File_B_Energy (Filename_B_energy);
    std::ofstream File_helicity1 (Filename_helicity1);
    std::ofstream File_helicity2 (Filename_helicity2);
    std::ofstream File_trsb_magn (Filename_trsb_magn);
    //std::ofstream File_vortices1 (Filename_vortices1);
    //std::ofstream File_vortices2 (Filename_vortices2);
    std::ofstream File_dual_stiff (Filename_dual_stiff);

    if(NSTART==0){
        /**Thermalization**/
        std::cout<< "Thermalization" << std::endl;
        for (int t = 0; t < MC.transient; t++) {
            metropolis(Site, MC, Hp, T_temp);
        }
    }

    for (int nM = NSTART; nM<MC.n_steps; nM++) {

        for (int t = 0; t < MC.tau; t++) {
            metropolis(Site, MC, Hp, T_temp);
        }

        //Measures
        mis.reset();

        energy(mis, Hp, Site);
        single_magnetization(Site, mis, N);   //In this function we are calculating the magnetization of both layers separately
        trsb_magnetization(mis, Site);
        helicity_modulus(Hp, Site, mis, N);   //Both lattices
        if(Hp.e == 1) {
            dual_stiffness(mis, Site);
        }
        mis.my_rank=PTp.rank;
        MPI_Barrier(MPI_COMM_WORLD);

        //Parallel Tempering swap
        parallel_temp(mis.E, my_beta, my_ind, PTp, PTroot);


        File_Energy        << mis.E << std::endl;
        File_Kin_Energy    << mis.E_kinetic << std::endl;
        File_Joseph_Energy << mis.E_josephson << std::endl;
        File_B_Energy      << mis.E_B << std::endl;
        File_Magetization1 << mis.m_phase[0] << " " << mis.m_phase[1] << std::endl;
        File_helicity1     << mis.Jd[0] << " " << mis.Ic[0] << std::endl;
        File_helicity2     << mis.Jd[1] << " " << mis.Ic[1] << std::endl;
        File_trsb_magn     << mis.trsb_m << std::endl;
        //File_vortices1     <<mis.vortices[0] << " " << mis.antivortices[0] <<std::endl;
        //File_vortices2     <<mis.vortices[1] <<" " << mis.antivortices[1] <<std::endl;
        File_dual_stiff    << mis.dual_stiff_Z << std::endl;

    }

    File_Energy.close();
    File_Kin_Energy.close();
    File_Joseph_Energy.close();
    File_B_Energy.close();
    File_Magetization1.close();
    File_helicity1.close();
    File_helicity2.close();
    File_trsb_magn.close();
    //File_vortices1.close();
    //File_vortices2.close();
    File_dual_stiff.close();
}


void myhelp(int argd, char **argu) {
    int i;
    fprintf(stderr,"Errore nei parametri su linea di comando; hai scritto:\n");
    for (i=0;i<argd;i++) fprintf(stderr," %s",argu[i]);
    fprintf(stderr,"\n");
    fprintf(stderr,"%s <L> <n_steps> <transient> <tau> <T> <restart> <K> <J1> <J2> <e> <theta_box> <DIR_IN> <DIR_OUT> \n",argu[0]);
    exit (EXIT_FAILURE);
}

