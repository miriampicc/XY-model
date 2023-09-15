#include "main.h"
#include "measures.h"
#include "robust_filesystem.h"

// Initialize random spin directions
std::vector<double> spins;
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
    //srand(time(NULL));
    srand(static_cast<unsigned int>(time(NULL)));
    //srand(2);
    std::string directory_read;
    std::string directory_write;

    if(argc == 7) {
        L = static_cast<size_t>(std::atoi(argv[1]));
        n_steps = static_cast<size_t>(std::atoi(argv[2]));
        T = std::atof(argv[3]);
        restart = std::atof(argv[4]);
        paths_dir::DIR_IN = directory_read = argv[5];
        paths_dir::DIR_OUT = directory_write = argv[6];
    }
    else{
        myhelp(argc, argv);
    }

    N=L*L;
    spins.resize(N);
    cc.resize(n_steps);
    std::cout<< L << std::endl;
    std::cout<< n_steps << std::endl;
    std::cout<< T << std::endl;
    std::cout<< restart << std::endl;
    std::cout<< directory_read << std::endl;
    std::cout<< directory_write << std::endl;

    // Initialize random spin directions

    std::string Filename_spin_init=(directory_write+"/Spin_dir_init.txt");
    std::string Filename_spin_fin=(directory_write+"/Spin_dir_fin.txt");

    if (restart == 1) {
        std::ifstream inputFile (directory_read + "/Spin_dir_fin.txt");
        if (!inputFile.is_open()) {
            std::cerr << "Error opening the file for Input." << std::endl;
            return 1;
        }
        int index = 0;
        while (inputFile >> spins[index]) {
            index++;
        }
        inputFile.close();
        std::ofstream File_Spin(Filename_spin_init);
        if(File_Spin.is_open()) {
            for (int i = 0; i < N; i++) {
                File_Spin << i << " " << spins[i] << std::endl;
            }
        }
        File_Spin.close();

    } else {
        for (int i = 0; i < N; i++) {
                spins[i] = rn::uniform_real_box(0, 2 * M_PI);
                std::ofstream File_Spin(Filename_spin_init);
                File_Spin << i << " " << spins[i ] << std::endl;
            }
        }


    // Run simulation --->>
    mainloop (spins, T, n_steps, N, directory_write);

    std::ofstream File_Spin_fin(Filename_spin_fin);
    if(File_Spin_fin.is_open()) {
        for (int i = 0; i < N; i++) {
            File_Spin_fin << i << " " << spins[i] << std::endl;
        }
    }
    File_Spin_fin.close();

    return 0;
}

void mainloop(std::vector<double> & spins, double T, int n_steps, size_t N, std::string directory_write) {

    double thetabox = M_PI/4.;
    std::vector<double> magnetization_values;
    double acc_ideal = 0.5;
    double acc = 0;
    Measures mis;
    int acc_rate = 0;

    std::cout<< "aperto!"<< std::endl;

    //Aprire i file output

    std::string Filename_magnetization=(directory_write+"/Magnetization.txt");
    std::string Filename_energy=(directory_write+"/Energy.txt");
    std::string Filename_restart=(directory_write+"/Restart.txt");
    std::string Filename_helicity=(directory_write+"Helicity modulus.txt");
    std::string Filename_vortices=(directory_write+"/Vortices.txt");


    std::ofstream File_Magetization (Filename_magnetization);
    std::ofstream File_Energy (Filename_energy);
    std::ofstream File_restart (Filename_restart);
    std::ofstream File_helicity (Filename_helicity);
    std::ofstream File_vortices (Filename_vortices);

    if (File_Energy.is_open()){

        std::cout<< "aperto!"<< std::endl;

    }

    for (size_t step = 0; step < n_steps; step++) {
        acc_rate = 0;

        for (size_t n_up = 0; n_up < N; n_up++) {
            mc_step(spins, T, acc_rate, thetabox);
        }

        mis.reset();
        energy(spins, mis);
        magnetization (spins, mis, N);
        helicity_modulus (spins, mis, N);
        vortex(spins, mis, N) ;

        //std::cout<<mis.Ic<<std::endl;    CORRETTO fino a qua

        File_Energy << step << " " << mis.E << std::endl;
        File_Magetization << step << " " << mis.M << std::endl;
        File_helicity << mis.Jd << " " << mis.Ic << std::endl;
        File_vortices <<mis.n_vort << "  " << mis.n_antivort <<std::endl;

        //CORRETTO fino a qua

        acc = (acc_rate) / ((double)N);

        double new_dthetabox = 2 * thetabox * (0.5 * (acc / acc_ideal));

        thetabox += new_dthetabox;
        thetabox/=2.;

        if (thetabox < 0.0) {
            thetabox += 2.0 * M_PI;
        } else if (thetabox > 2.0 * M_PI) {
            thetabox -= 2.0 * M_PI;
        }
    }

    File_Energy.close();
    File_Magetization.close();
    File_helicity.close();
    File_restart.close();
    File_vortices.close();
}

// Define Monte Carlo step
void mc_step(std::vector<double> & spins, double T, int &acc, double thetabox) {

    int i = rand() % L;
    int j = rand() % L;

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

    double delta_e = local_energy(spins_new) - local_energy(spins); //make a local energy function

    if (delta_e < 0) {
        spins[i+j*L] = spins_new[i+j*L] ;
        acc += 1;
    } else {
        double k = randomDouble(0, 1);
        if (exp(-delta_e / T) >  k) {
            spins[i+j*L] = spins_new[i+j*L] ;
            acc += 1;
        }
    }
}

double local_energy(const std::vector<double>& spins) {
    double sum_cosines = 0.0;

    for (size_t i = 0; i < L; i++) {
        for (size_t j = 0; j < L; j++) {
            sum_cosines += cos(spins[i+j*L] - spins[((i + 1) % L)+j*L]) +
                           cos(spins[i+j*L] - spins[i+((j + 1)% L)*L]);
        }
    }
    return -sum_cosines;
}

void myhelp(int argd, char** argu) {
    int i;
    fprintf(stderr,"Errore nei parametri su linea di comando; hai scritto:\n");
    for (i=0;i<argd;i++) fprintf(stderr," %s",argu[i]);
    fprintf(stderr,"\n");
    fprintf(stderr,"%s <L> <n_steps> <T> <restart> <DIR_IN> <DIR_OUT> \n",argu[0]);
    exit (EXIT_FAILURE);
}

