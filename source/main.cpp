#include "main.h"
#include "measures.h"

// Initialize random spin directions
std::vector<double> spins;
std::vector<double> cc;
size_t L, n_steps;
double T;
bool restart;
int N, n_vort_start, n_vort_end;

// Function to generate random double between min and max
double randomDouble(double min, double max) {
    return min + (max - min) * (double)rand() / RAND_MAX;
}

double wrapToPi(double angle) {

    if (angle > M_PI) {
        angle -= 2*M_PI;
    }   else if (angle < -M_PI) {
        angle += 2*M_PI;
        }
    return angle;
}

int main(int argc, char *argv[]) {
    //srand(time(NULL));
    srand(static_cast<unsigned int>(time(NULL)));
    //srand(2);

    L = static_cast<size_t>(std::atoi(argv[1]));
    n_steps = static_cast<size_t>(std::atoi(argv[2]));
    T = std::atof(argv[3]);
    restart = std::atof(argv[4]);
    N=L*L;
    spins.resize(N);
    cc.resize(n_steps);
    std::cout<< L << std::endl;
    std::cout<< n_steps << std::endl;
    std::cout<< T << std::endl;
    std::cout<< restart << std::endl;


    // Initialize random spin directions


    std::ofstream outputFile1 ("/Users/mirimi/CLionProjects/XY_model/results.txt");
    std::ofstream outputFile2 ("/Users/mirimi/CLionProjects/XY_model/results_m.txt");
    std::ofstream outputFile3 ("/Users/mirimi/CLionProjects/XY_model/spin_dir.txt");
    std::ofstream outputFile5 ("/Users/mirimi/CLionProjects/XY_model/helicity_modulus.txt");


    // Check if the files are opened successfully
    if (!outputFile1.is_open()) {
        std::cerr << "Error opening the file 1." << std::endl;
        return 1;
    }

    if (!outputFile2.is_open()) {
        std::cerr << "Error opening the file 2." << std::endl;
        return 1;
    }

    if (!outputFile3.is_open()) {
        std::cerr << "Error opening the file 3." << std::endl;
        return 1;
    }

    if (!outputFile5.is_open()) {
        std::cerr << "Error opening the file 5." << std::endl;
        return 1;
    }



    if (restart == 1) {


        std::ifstream inputFile ("/Users/mirimi/CLionProjects/XY_model/spin_dir_f.txt");
        if (!inputFile.is_open()) {
            std::cerr << "Error opening the file 4 for input." << std::endl;
            return 1;
        }

        int index = 0;
        while (inputFile >> spins[index]) {
            index++;
        }
        inputFile.close();

        for (int i = 0; i < spins.size(); i++){

            //std::cout<<" "<< spins[i] <<std::endl;
            //std::cerr << "Non sta stampando niente" << std::endl;
        }

    } else {
        for (int i = 0; i < L; i++) {
            for (int j = 0; j < L; j++) {

                spins[i + j * L] = randomDouble(0, 2 * M_PI);
                outputFile3 << i << " " << spins[i + j * L] << std::endl;
            }
        }
    }

    //Let us obtain the number of vortices at the beginning of the simulation

    n_vort_start = vortex(spins) ;
    std::cout<<"Number of vortices start :  "<< n_vort_start <<std::endl;


    // Run simulation --->> void mainloop()

    std::ofstream outputFile4 ("/Users/mirimi/CLionProjects/XY_model/spin_dir_f.txt");
    if (!outputFile4.is_open()) {
        std::cerr << "Error opening the file 4." << std::endl;
        return 1;
    }

    double thetabox = M_PI/4.;
    std::vector<double> magnetization_values;
    double acc_ideal = 0.5;
    double acc = 0;
    double Jp;
    Measures mis;
    int acc_rate = 0;

    for (size_t step = 0; step < n_steps; step++) {
        acc_rate = 0;

        for (size_t n_up = 0; n_up < N; n_up++) {
            mc_step(spins, T, acc_rate, thetabox);
        }

        mis.reset();
        energy(spins, mis);
        magnetization (spins, mis, N);
        helicity_modulus (spins, mis, N);

        //std::cout<<mis.Ic<<std::endl;    CORRETTO fino a qua

        outputFile1 << step << " " << mis.E << std::endl;
        outputFile2 << step << " " << mis.M << std::endl;
        outputFile5 << mis.Jd << " " << mis.Ic << std::endl;

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

    outputFile5.close();

    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {

            outputFile4 << spins[i+j*L] << std::endl;
        }
    }

    n_vort_end = vortex(spins) ;
    std::cout<<"Number of vortices end :  "<< n_vort_end <<std::endl;


    std::ifstream inputFile5 ("/Users/mirimi/CLionProjects/XY_model/helicity_modulus.txt");
    if (!inputFile5.is_open()) {
        std::cerr << "Error opening the file 5 for input." << std::endl;
        return 1;
    }

    int index = 0;
    while (inputFile5 >> cc[index]) {
        index++;
    }
    inputFile5.close();

    for (int i=0; i< cc.size(); i++) {
        std::cout<< cc[i] <<std::endl;
        }

    std::cout<<"cc vector size is:   "<< cc.size() << std::endl;

    //Jp = J_param ( cc, T);

    std::cout<<"J_p:  "<< Jp <<std::endl;

    outputFile1.close();
    outputFile2.close();
    outputFile3.close();
    outputFile4.close();

    return 0;
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

int vortex (std::vector<double> spins) {

    double sq; //this should represent the sum of the angle of every square
    int i, j;
    int n_plus, n_minus, n_tot;
    double phi_1, phi_1_fin,  phi_2, phi_2_fin, phi_3, phi_3_fin, phi_4, phi_4_fin;


    n_plus = 0;
    n_minus = 0;
    n_tot = 0;

    for ( i=0; i<L; i++ ) {

        for ( j=0; j<L; j++ ) {

            sq = 0 ;
            phi_1 = 0;
            phi_2 = 0;
            phi_3 = 0;
            phi_4 = 0;

            if ( j == (L-1) & i != (L-1)){

                phi_1 = spins[(i+1)+j*L] - spins[i+j*L] ;
                phi_1_fin = wrapToPi(phi_1);

                phi_2 = spins[(i+1)+(j+1-L)*L] - spins[(i+1)+j*L] ;
                phi_2_fin = wrapToPi(phi_2);

                phi_3 = spins[i+(j+1-L)*L] - spins[(i+1)+(j+1-L)*L] ;
                phi_3_fin = wrapToPi(phi_3);

                phi_4 = spins[i+j*L] - spins[i+(j+1-L)*L] ;
                phi_4_fin = wrapToPi(phi_4);

            } else if ( j != (L-1) & i == (L-1) ) {

                phi_1 = spins[(i+1-L)+j*L] - spins[i+j*L] ;
                phi_1_fin = wrapToPi(phi_1);

                phi_2 = spins[(i+1-L)+(j+1)*L] - spins[(i+1-L)+j*L] ;
                phi_2_fin = wrapToPi(phi_2);

                phi_3 = spins[i+(j+1)*L] - spins[(i+1-L)+(j+1)*L] ;
                phi_3_fin = wrapToPi(phi_3);

                phi_4 = spins[i+j*L] - spins[i+(j+1)*L] ;
                phi_4_fin = wrapToPi(phi_4);


            } else if ( j == (L-1) & i == (L-1) ) {

                phi_1 = spins[(i+1-L)+j*L] - spins[i+j*L] ;
                phi_1_fin = wrapToPi(phi_1);

                phi_2 = spins[(i+1-L)+(j+1-L)*L] - spins[(i+1-L)+j*L] ;
                phi_2_fin = wrapToPi(phi_2);

                phi_3 = spins[i+(j+1-L)*L] - spins[(i+1-L)+(j+1-L)*L] ;
                phi_3_fin = wrapToPi(phi_3);

                phi_4 = spins[i+j*L] - spins[i+(j+1-L)*L] ;
                phi_4_fin = wrapToPi(phi_4);


            } else {

                phi_1 = spins[(i+1)+j*L] - spins[i+j*L] ;

                // Wrap phi_1 to the range [-π, π]
                phi_1_fin = wrapToPi(phi_1);


                phi_2 = spins[(i+1)+(j+1)*L] - spins[(i+1)+j*L] ;
                phi_2_fin = wrapToPi(phi_2);


                phi_3 = spins[i+(j+1)*L] - spins[(i+1)+(j+1)*L] ;
                phi_3_fin = wrapToPi(phi_3);

                phi_4 = spins[i+j*L] - spins[i+(j+1)*L] ;
                phi_4_fin = wrapToPi(phi_4);

            }

            sq = phi_1_fin + phi_2_fin + phi_3_fin + phi_4_fin ;

            sq /= 2*M_PI;

            if ( sq == 1 ) {
                n_plus ++;
            } else if ( sq == -1) {
                n_minus++;
            }

        }
    }

    n_tot = n_plus + n_minus ;

    return n_tot;
}

/*double J_param (std::vector<double> cc, double T) {

    int i;
    double J;
    double sum_sq = 0, sum_cc = 0;
    double mean_sq, mean_cc;

    for (i=0; i < cc.size(); i++){

        sum_sq += cc[i] * cc[i];
        sum_cc += cc[i];

    }

    mean_sq = sum_sq / n_steps ;
    mean_cc = sum_cc / n_steps ;

    J = N * (mean_sq - mean_cc * mean_cc) / T ;

    return J;
}*/

