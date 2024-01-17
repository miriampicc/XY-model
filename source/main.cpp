#include "main.h"
#include <math.h>
#include "measures.h"
#include "robust_filesystem.h"
#include "initialization.h"

// Initialize random spin directions
std::vector<double> spins1;
//double* spin_1 = double[N];
std::vector<double> spins2;
std::vector<double> cc;
size_t L, n_steps;
double T;
bool restart;
size_t N;

// Function to generate random double between min and max
double randomDouble(double min, double max) {
    return min + (max - min) * (double)rand() / RAND_MAX;
}


int main(int argc, char *argv[]) {

    srand(static_cast<unsigned int>(time(NULL)));
    //srand(2);
    std::string directory_read;
    std::string directory_write;
    struct H_parameters Hp{};

    if(argc == 10) {
        L = static_cast<size_t>(std::atoi(argv[1]));
        n_steps = static_cast<size_t>(std::atoi(argv[2]));
        T = std::atof(argv[3]);
        restart = std::atof(argv[4]);
        Hp.K = std::atof (argv[5]);
        Hp.J1 = std::atof (argv[6]);
        Hp.J2 = std::atof (argv[7]);
        paths_dir::DIR_IN = directory_read = argv[8];
        paths_dir::DIR_OUT = directory_write = argv[9];
    }
    else{
        myhelp(argc, argv);
    }

    N=L*L;
    spins1.resize(N);
    spins2.resize(N);
    cc.resize(n_steps);
    std::cout<< L << std::endl;
    std::cout<< n_steps << std::endl;
    std::cout<< T << std::endl;
    std::cout<< Hp.K << std::endl;
    std::cout<< Hp.J1 << std::endl;
    std::cout<< Hp.J2 << std::endl;
    std::cout<< restart << std::endl;
    std::cout<< directory_read << std::endl;
    std::cout<< directory_write << std::endl;

    // Initialize random spin directions

    std::string Filename_spin_1_init=(directory_write+"/Spin_1_dir_init.txt");
    std::string Filename_spin_1_fin=(directory_write+"/Spin_1_dir_fin.txt");
    std::string Filename_spin_2_init=(directory_write+"/Spin_2_dir_init.txt");
    std::string Filename_spin_2_fin=(directory_write+"/Spin_2_dir_fin.txt");

    if (restart == 1) {
        std::ifstream inputFile1 (directory_read + "/Spin_1_dir_fin.txt");
        if (!inputFile1.is_open()) {
            std::cerr << "Error opening the file for Input." << std::endl;
            return 1;
        }
        int index = 0;
        while (inputFile1 >> spins1[index]) {
            index++;
        }
        inputFile1.close();
        std::ofstream File_Spin_1(Filename_spin_1_init);
        if(File_Spin_1.is_open()) {
            for (int i = 0; i < N; i++) {
                File_Spin_1 << spins1[i] << std::endl;
            }
        }
        File_Spin_1.close();

        //second lattice
        std::ifstream inputFile2 (directory_read + "/Spin_2_dir_fin.txt");
        if (!inputFile2.is_open()) {
            std::cerr << "Error opening the file for Input." << std::endl;
            return 1;
        }
        int index2 = 0;
        while (inputFile2 >> spins1[index2]) {
            index2++;
        }
        inputFile2.close();
        std::ofstream File_Spin_2(Filename_spin_2_init);
        if(File_Spin_2.is_open()) {
            for (int i = 0; i < N; i++) {
                File_Spin_2 << spins2[i] << std::endl;
            }
        }
        File_Spin_2.close();

    } else {
        for (int i = 0; i < N; i++) {
                spins1[i] = rn::uniform_real_box(0, 2 * M_PI);
                std::ofstream File_Spin_1(Filename_spin_1_init);
                File_Spin_1 << spins1[i] << std::endl;

                spins2[i] = rn::uniform_real_box(0, 2 * M_PI);
                std::ofstream File_Spin_2(Filename_spin_2_init);
                File_Spin_2 << spins2[i] << std::endl;
            }
        }


    // Run simulation --->>
    mainloop (spins1, spins2, T, n_steps, N, Hp, directory_write);

    std::ofstream File_Spin_1_fin(Filename_spin_1_fin);
    if(File_Spin_1_fin.is_open()) {
        for (int i = 0; i < N; i++) {
            File_Spin_1_fin << spins1[i] << std::endl;
        }
    }
    File_Spin_1_fin.close();

    std::ofstream File_Spin_2_fin(Filename_spin_2_fin);
    if(File_Spin_2_fin.is_open()) {
        for (int i = 0; i < N; i++) {
            File_Spin_2_fin << spins2[i] << std::endl;
        }
    }
    File_Spin_2_fin.close();

    return 0;
}

void mainloop(std::vector<double> & spins1, std::vector<double> & spins2, double T, int n_steps, size_t N, struct H_parameters &Hp, std::string directory_write) {

    double thetabox1 = M_PI/4.;
    double thetabox2 = M_PI/4.;
    std::vector<double> magnetization_values;
    double acc_ideal = 0.5;
    double acc1 = 0;
    double acc2 = 0;
    Measures mis;
    int acc_rate1 = 0;
    int acc_rate2 = 0;

    std::cout<< "aperto!"<< std::endl;
    std::cout<< "The temperature is: " << T << std::endl;

    //Aprire i file output

    std::string Filename_magnetization1=(directory_write+"/Magnetization1.txt");

    std::string Filename_magnetization2=(directory_write+"/Magnetization2.txt");
    std::string Filename_energy=(directory_write+"/Energy.txt");
    std::string Filename_restart=(directory_write+"/Restart.txt");
    std::string Filename_helicity1=(directory_write+"/Helicity_modulus1.txt");
    std::string Filename_helicity2=(directory_write+"/Helicity_modulus2.txt");
    std::string Filename_trsb_magn=(directory_write+"/trsb_magnetization.txt");
    std::string Filename_vortices1=(directory_write+"/Vortices1.txt");
    std::string Filename_antivortices1=(directory_write+"/Antivortices1.txt");
    std::string Filename_vortices2=(directory_write+"/Vortices2.txt");
    std::string Filename_antivortices2=(directory_write+"/Antivortices2.txt");


    std::ofstream File_Magetization1 (Filename_magnetization1);

    std::ofstream File_Magetization2 (Filename_magnetization2);
    std::ofstream File_Energy (Filename_energy);
    std::ofstream File_restart (Filename_restart);
    std::ofstream File_helicity1 (Filename_helicity1);
    std::ofstream File_helicity2 (Filename_helicity2);
    std::ofstream File_trsb_magn (Filename_trsb_magn);
    std::ofstream File_vortices1 (Filename_vortices1);
    std::ofstream File_antivortices1 (Filename_antivortices1);
    std::ofstream File_vortices2 (Filename_vortices2);
    std::ofstream File_antivortices2 (Filename_antivortices2);

    if (File_Energy.is_open()){

        std::cout<< "aperto!"<< std::endl;

    }

    //std::string magn1_log = "";

    for (size_t step = 0; step < n_steps; step++) { //so, we know that we have n_steps montecarlo
        acc_rate1 = 0;
        acc_rate2 = 0;
        int k =0;

        for (size_t n_up = 0; n_up < N; n_up++) {
            mc_step(spins1, spins2, T, acc_rate1, thetabox1, Hp);   //and, here, for each monte carlo step, we have to update each spin
            mc_step(spins2, spins1, T, acc_rate2, thetabox2, Hp);
        }

        mis.reset();
        energy(spins1, spins2,  mis);
        magnetization (spins1, mis, N, k);
        vortex(spins1, mis, N, k) ;
        helicity_modulus (spins1, mis, N, k);
        k=1;
        magnetization (spins2, mis, N, k);
        vortex(spins2, mis, N, k) ;
        helicity_modulus (spins2, mis, N, k);
        trsb_magnetization(spins1, spins2,  mis, N);


        //std::cout<<"Step   "<< step << std::endl;

        //std::cout<<mis.Ic<<std::endl;    CORRETTO fino a qua

        //std::cout<<"Step:   "<<step  << " Fictitious magnetisation:" <<mis.Mp<<std::endl;

        //magn1_log += mis.M1;
        //magn1_log += "\n";

        File_Energy << mis.E << std::endl;
        File_Magetization1 << mis.M1 << std::endl;
        File_Magetization2 << mis.M2 << std::endl;
        File_helicity1 << mis.Jd1 << " " << mis.Ic1 << std::endl;
        File_helicity2 << mis.Jd2 << " " << mis.Ic2 << std::endl;
        File_trsb_magn << mis.Mp << std::endl;
        File_vortices1 <<mis.n_vort1 <<std::endl;
        File_antivortices1 << mis.n_antivort1 <<std::endl;
        File_vortices2 <<mis.n_vort2 <<std::endl;
        File_antivortices2 << mis.n_antivort2 <<std::endl;


        //CORRETTO fino a qua

        acc1 = (acc_rate1) / ((double)N);
        acc2 = (acc_rate2) / ((double)N);

        double new_thetabox1 = thetabox1 * (0.5 + 0.5 * (acc1 / acc_ideal));
        //std::cout<< "thetabox: "<< thetabox <<  " new thetabox: "<< new_thetabox << " acc: " << acc<< std::endl;
        thetabox1 = new_thetabox1;

        if (thetabox1 < 0.0) {
            thetabox1 += 2.0 * M_PI;
        } else if (thetabox1 > 2.0 * M_PI) {
            thetabox1 -= 2.0 * M_PI;
        }

        double new_thetabox2 = thetabox2 * (0.5 + 0.5 * (acc2 / acc_ideal));
        //std::cout<< "thetabox: "<< thetabox <<  " new thetabox: "<< new_thetabox << " acc: " << acc<< std::endl;
        thetabox2 = new_thetabox2;

        if (thetabox2 < 0.0) {
            thetabox2 += 2.0 * M_PI;
        } else if (thetabox2 > 2.0 * M_PI) {
            thetabox2 -= 2.0 * M_PI;
        }

    }
    //std::string Filename_magnetization1=();

    /*std::ofstream File_Magetization1 (directory_write+"/Magnetization1.txt");
    File_Magetization1 << magn1_log << std::endl;
    File_Magetization1.close(); */

    File_Energy.close();
    File_Magetization1.close();
    File_Magetization2.close();
    File_helicity1.close();
    File_helicity2.close();
    File_trsb_magn.close();
    File_restart.close();
    File_vortices1.close();
    File_antivortices1.close();
    File_vortices2.close();
    File_antivortices2.close();
}

// Define Monte Carlo step
void mc_step(std::vector<double> & spins, std::vector<double> & spins_not, double T, int &acc, double thetabox, struct H_parameters &Hp) {

    int i = rn::uniform_integer_box(0,L-1); //prendi solo un numero random tra 0 e N (posso ricares sia i che j)
    int j =rn::uniform_integer_box(0,L-1);

    double d_spin = randomDouble(-thetabox, thetabox);

    std::vector<double> spins_new = spins;

    spins_new[i+j*L] += d_spin;

    //Check if the spin angles is between 0 and 2pi, because if it is not, the number will just increase.
    if (spins_new[i + j * L] < 0.0) {
        spins_new[i + j * L] += 2.0 * M_PI;
    } else if (spins_new[i + j * L] > 2.0 * M_PI) {
        spins_new[i + j * L] -= 2.0 * M_PI;
    }

    //std::cout<< i << " " << j << " " << " " << d_spin << " " << spins[i+j*L] << " " << spins_new[i+j*L] << std::endl;

    double delta_e = local_energy(spins_new, spins_not, Hp) - local_energy(spins, spins_not, Hp); //make a local energy function

    if (delta_e < 0) {
        spins[i+j*L] = spins_new[i+j*L] ;
        acc += 1;
    } else {
        double k = rn::uniform_real_box(0,1);
        if (exp(-delta_e / T) >  k) {
            spins[i+j*L] = spins_new[i+j*L] ;
            acc += 1;
        }
    }
}

double local_energy(const std::vector<double>& spins, const std::vector<double>& spins_not, struct H_parameters &Hp) {
    double sum_cosines_s = 0.0;
    double sum_cosines_n = 0.0;
    double interaction = 0.0;

    for (size_t i = 0; i < L; i++) {
        for (size_t j = 0; j < L; j++) {
            sum_cosines_s += (cos(spins[i+j*L] - spins[((i + 1) % L)+j*L]) +
                           cos(spins[i+j*L] - spins[i+((j + 1)% L)*L]));
            sum_cosines_n += cos(spins_not[i+j*L] - spins_not[((i + 1) % L)+j*L]) +
                            cos(spins_not[i+j*L] - spins_not[i+((j + 1)% L)*L]);
            interaction +=  cos(2*(spins[i+j*L]-spins_not[i+j*L]));
        }
    }
    return - Hp.J1 * sum_cosines_s - Hp.J2 * sum_cosines_n + Hp.K * interaction;
}

void myhelp(int argd, char** argu) {
    int i;
    fprintf(stderr,"Errore nei parametri su linea di comando; hai scritto:\n");
    for (i=0;i<argd;i++) fprintf(stderr," %s",argu[i]);
    fprintf(stderr,"\n");
    fprintf(stderr,"%s <L> <n_steps> <T> <restart> <DIR_IN> <DIR_OUT> \n",argu[0]);
    exit (EXIT_FAILURE);
}

