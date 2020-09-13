#include<iostream>
#include<math.h>
#include<stdlib.h>
#include<fstream>
#include<random>
#include<cassert>
#include<array>

std::ofstream DATA("data.dat",std::ios::out);

const int L = 6;
const int D = 3;
const int n = std::pow(L,D);

double T = 12.0;
const double T_min = 0.5;

const double dT = 0.1;
const int sweeps = 1000;
const int transient = 1000;
//normalization for averaging
double norm = (1.0/double(sweeps*n));

//samples uniform distribution to get random value between 0 and 1
std::random_device rd;
std::mt19937_64 gen(rd());
double gen_rand(){
    std::uniform_real_distribution<double> dis(0.0, 1.0);
    return dis(gen);
}

//returns random integer in specified range
int gen_int(int min, int max){
    std::uniform_int_distribution<int> dis(min,max);
    return dis(gen);
}

//class for an individual spin
struct Site{
    int index;
    int state;
    int energy;
    
    //each pair has the following data: (index of site this site is coupled to, value of coupling)
    //Coupling goes one way and there are no redundant couplings
    std::vector<std::pair<int,int>> coupled_to;

    //holds info on which are coupled to a particular site
    std::vector<std::pair<int,int>> coupled_from;
};

//class for entire lattice
struct Lattice{

    //holds all the sites
    std::array<std::array<std::array<Site, L>, L>,L> lat;

    //returns index given an ordered pair
    int index(int i, int j, int k) {
        //return ((m+L)%L)*L + ((i+L)%L);
        return (i % L)*L*L + ((j)%L)*L + ((k)%L);

    }

    //returns (x,y,z) ordered pair when given an index
    std::vector<int> to_pos(int i){
        std::vector<int> a;

        a.push_back((int)floor(i/(L))%L);
        a.push_back(i%L);
        a.push_back(floor(i/(L*L)));
    
    return a;
    }

    //used in initializtion
    int gen_couple(){
        //return gen_rand() < 0.5 ? -1 : 1; //use this one for spin glass
        return 1; //use this one for ferromagnet
    }

    void initialize(){
        //sets spins to +- 1
        for(int i=0; i<L; i++){
            for(int j=0; j<L; j++){
                for(int k=0; k<L; k++){
                    if(gen_rand()>0.5)
                        lat[i][j][k].state = 1;

                    else    
                        lat[i][j][k].state = -1;
                }
            }
        }
        
        //assigns each spin an index
        int counter = 0;
        for(int i=0; i<L; i++){
            for(int j=0; j<L; j++){
                for(int k=0; k<L; k++){
                    lat[i][j][k].index = counter;
                    counter++;
                }
            }
        }
        
        //fix this
        //sets couplings
        int a = 0;
        int b = 0;
        int c = 0;
        int d = 0;
        int e = 0;
        int f = 0;

        for(int i=0; i<L; i++){
            for(int j=0; j<L; j++){
                for(int k=0; k<L; k++){
                    a = index(i,j,k+1);    
                    b = index(i,j+1,k);
                    c = index(i+1,j,k);
                    d = gen_couple();
                    e = gen_couple();
                    f = gen_couple();

                    //sites our chosen site is coupled to
                    lat[i][j][k].coupled_to.emplace_back(a,d);
                    lat[i][j][k].coupled_to.emplace_back(b,e);
                    lat[i][j][k].coupled_to.emplace_back(c,f);

                    //sites coupled to our chosen site
                    lat[to_pos(a)[0]][to_pos(a)[1]][to_pos(a)[2]].coupled_from.emplace_back(index(i,j,k),d);
                    lat[to_pos(b)[0]][to_pos(b)[1]][to_pos(b)[2]].coupled_from.emplace_back(index(i,j,k),e);
                    lat[to_pos(c)[0]][to_pos(c)[1]][to_pos(c)[2]].coupled_from.emplace_back(index(i,j,k),f);
                }
            }
        }
    }

    //picks an integer suitible for use as an index at random
    int choose_site(){
        return gen_int(0,n-1);
    }

    //flips spin
    void flip(Site& s){
        s.state = -s.state;
    }

    //returns energy of arbitrary site
    int e_site(Site& s){
        double e = 0;
        for(int i=0; i<D; i++){
            //adds coupled_to and coupled_from elements
            e += -1.0*s.state*(s.coupled_to[i].second*lat[to_pos(s.coupled_to[i].first)[0]][to_pos(s.coupled_to[i].first)[1]][to_pos(s.coupled_to[i].first)[2]].state + s.coupled_from[i].second*lat[to_pos(s.coupled_from[i].first)[0]][to_pos(s.coupled_from[i].first)[1]][to_pos(s.coupled_from[i].first)[2]].state);
        }   
        //s.energy = e;
        return e;
    }

    bool test_flip(Site s, int &de){
        //change in energy for specific site
        de=-2.0*e_site(s);
        
        //flips due to lower energy
        if (de<0)
            return true;
    
        //flips due to heat bath
        else if (gen_rand()<exp(-de/T))
            return true;
        
        else
            return false;
    }

    void transient_results(){
        int de=0;

        for(int i=0; i<transient; i++){
            for(int j=0; j<n; j++){
                for(int k=0; k<L; k++){
                    int location = choose_site();
                    int x = to_pos(location)[0];
                    int y = to_pos(location)[1];
                    int z = to_pos(location)[2];
                    
                    if (test_flip(lat[x][y][z],de)){
                        flip(lat[x][y][z]);
                    }
                }
            }
        }
    }

    int total_magnetization(){
        int m=0;
        for(int i=0; i<L; i++){
            for(int j=0; j<L; j++){
                for(int k=0; k<L; k++){
                    m += lat[i][j][k].state;  
                }
            }
        }
        return m;
    }

    int total_energy(){
        int e=0;
        for(int i=0; i<L; i++){
            for(int j=0; j<L; j++){
                for(int k=0; k<L; k++){
                    e += e_site(lat[i][j][k]);
                }
            }
        }
        return e;
    }
    
    //this is only for the 2D case
    void output(){
        std::cout << std::string( 100, '\n' );
        for(int i=0; i<L; i++){
            for(int j=0; j<L; j++){
                if (lat[i][j][0].state<0)
                    std::cout<<"- ";
                else
                    std::cout<<"+ ";
            }
            std::cout << std::endl;
        }
    }

};

int main(){

    Lattice spins;
    spins.initialize();

    double E=0,E_sq=0,E_sq_avg=0,E_avg=0,E_tot=0,E_tot_sq=0;
    double M=0,M_sq=0,M_sq_avg=0,M_avg=0,M_tot=0,M_tot_sq=0;
    double M_abs=0,M_abs_avg=0,Mq_avg=0,M_abs_tot=0,mqtot=0;
    int de = 0;
    int x = 0;
    int y = 0;
    int z = 0;

    for( ; T>=T_min; T-=dT){
        std::cout << T << std::endl;
        spins.transient_results();

        //observables adopt equilibrated lattice configurations values 
        M = spins.total_magnetization();
        M_abs = abs(spins.total_magnetization());
        E = spins.total_energy();
        
        //reset at each temp
        E_tot = 0;
        E_tot_sq = 0;
        M_tot = 0;
        M_tot_sq = 0;
        M_abs_tot = 0;
        mqtot = 0;
        
        //Metropolis
        for(int i=0; i<sweeps; i++){
            for(int j=0; j<n; j++){

                int location = spins.choose_site();
                x = spins.to_pos(location)[0];
                y = spins.to_pos(location)[1];
                z = spins.to_pos(location)[2];

                if (spins.test_flip(spins.lat[x][y][z],de)){
                    spins.flip(spins.lat[x][y][z]);
               
                    //adjust observables
                    E += 2.0*de;
                    M += 2.0*spins.lat[x][y][z].state;
                    M_abs += abs(spins.lat[x][y][z].state);
                }
            }

            //keep summation of observables
            E_tot += E/2.0;
            //E_tot += E;
            
            //so as not to count the energy for each spin twice //Are the factors needed with the rewrite?
            E_tot_sq+=E/2.0*E/2.0;
            //E_tot_sq += E*E;
            M_tot += M;
            M_tot_sq += M*M;
            mqtot += M*M*M*M;
            M_abs_tot += (sqrt(M*M));
        }

        //average observables
        E_avg = E_tot*norm;
        E_sq_avg = E_tot_sq*norm;
        M_avg = M_tot*norm;
        M_sq_avg = M_tot_sq*norm;
        M_abs_avg = M_abs_tot*norm; 
        Mq_avg = mqtot*norm;

        DATA<<T<<                                       //temperature 
       "\t"<<M_avg<<"\t"<<M_abs_avg<<"\t"<<M_sq_avg<<     //<M>;<|M|>;<M^2> per spin
       "\t"<<(M_sq_avg-(M_avg*M_avg*n))/(T)<<            //susceptibility per spin (X)
       "\t"<<(M_sq_avg-(M_abs_avg*M_abs_avg*n))/(T)<<      //susceptibility per spin (X’)
       "\t"<<E_avg<<"\t"<<E_sq_avg<<                     //<E>;<E^2> per spin
       "\t"<<(E_sq_avg-(E_avg*E_avg*n))/(T*T)<<          //heat capacity (C) per spin
       "\t"<<1.0-(Mq_avg/(3.0*M_sq_avg*M_sq_avg))<<std::endl;            //Binder cumulant (U_L)

       spins.output();
    }
}