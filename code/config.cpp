#include "header.h"

//Create container for main parameters and read them from the parser 
Config::Config(int pid, int nproc, int argc, char *argv[], bool &exit, string print_file = ""):
    pid(pid),
    master(pid == 0),
    nproc(nproc),
    L_scale(1.),
    T_scale(1.)
{
    
    //Read parameters
    OptionsParser args(argc, argv);
    args.AddOption(&Lx, "-Lx", "--Lx",
                   "Length of domain in x direction.");
    args.AddOption(&Ly, "-Ly", "--Ly",
                   "Length of domain in y direction.");
    
    args.AddOption(&dt_init, "-dt", "--time_step",
                   "Initial time step.");
    args.AddOption(&t_final, "-t_f", "--t_final",
                   "Final time.");
    args.AddOption(&vis_steps_max, "-v_s", "--visualization_steps",
                   "Visualize every n-th timestep.");
    
    args.AddOption(&refinements, "-ref", "--refinements",
                   "Number of total uniform refinements.");
    args.AddOption(&order, "-o", "--order",
                   "Finite element order (polynomial degree) or -1 for isoparametric space.");
    args.AddOption(&abstol_solver, "-abstol_s", "--abstol_s",
                   "Absolute tolerance of solver.");
    args.AddOption(&reltol_solver, "-reltol_s", "--reltol_s",
                   "Relative tolerance of solver.");
    args.AddOption(&iter_solver, "-iter_s", "--iter_s",
                   "Iterations of solver.");
    args.AddOption(&abstol_sundials, "-abstol_sun", "--abstol_sun",
                   "Absolute tolerance of SUNDIALS.");
    args.AddOption(&reltol_sundials, "-reltol_sun", "--reltol_sun",
                   "Relative tolerance of SUNDIALS.");
    args.AddOption(&Epsilon, "-eps", "--epsilon",
                   "Small constant.");
    
    args.AddOption(&Rx, "-Rx", "--Rx",
                   "Coordinate x of vortex.");
    args.AddOption(&Ry, "-Ry", "--Ry",
                   "Coordinate y of vortex.");
    args.AddOption(&Sigma, "-sigma", "--sigma",
                   "Standard deviation of distribution.");
    args.AddOption(&Gamma, "-gamma", "--gamma",
                   "Intensity of distribution.");
    args.AddOption(&Viscosity, "-visc", "--visc",
                   "Kinematic Viscosity.");
    
    //Check if parameters were read correctly and print results
    args.Parse();
    if (!args.Good()){
        if (master) args.PrintUsage(cout);
        exit = true;
    } else {
        if (master){
            args.PrintOptions(cout);
            cout << "\n";

            ofstream out;
            out.open(print_file, std::ios::trunc);
            args.PrintOptions(out);
            out << "\n";
            out.close();
        }
    }
};     

//Change dimension scales in the parameters
void Config::Adimentionalize(double new_T_scale, double new_L_scale, string print_file = ""){

    //Set change factors
    double ref_T = T_scale/new_T_scale;
    double ref_L = L_scale/new_L_scale;

    //Time variables
    dt_init *= ref_T;
    t_final *= ref_T;

    //Space variables
    Lx *= ref_L;
    Ly *= ref_L;
    Rx *= ref_L;
    Ry *= ref_L;
    Sigma *= ref_L;

    //Other variables
    Gamma *= ref_L*ref_L/ref_T;
    Viscosity *= ref_L*ref_L/ref_T;

    //Update scale
    T_scale = new_T_scale;
    L_scale = new_L_scale;

    //Print characteristic variables
    if (master){
        cout << "Time scale: " << T_scale << " s\n"
             << "Lenght scale: " << L_scale << " cm\n"
             << "Reynolds number: " << 1/Viscosity << "\n"
             << "Vortex ratio: " << Sigma/Rx << "\n\n";

        ofstream out;
        out.open(print_file, std::ios::app);
        out << "Time scale: " << T_scale << " s\n"
            << "Lenght scale: " << L_scale << " cm\n"
            << "Reynolds number: " << 1/Viscosity << "\n"
            << "Vortex ratio: " << Sigma/Rx << "\n\n";
        out.close();
    }
}
