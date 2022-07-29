#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "mfem.hpp"

using namespace std;
using namespace mfem;

//Main variables for the program
struct Config{
    
    //Initialization
    Config(int pid, int nproc, int argc, char *argv[], bool &exit, string print_file = ""):
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

        
        //Check if parameters were read correctly
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

    //MPI parameters
    bool master;
    int nproc;
    int pid;

    //Time and visualization parameters
    double dt_init;
    double t_final;
    int vis_steps_max;

    //FEM parameters
    int refinements;
    int order;
    double reltol_solver;
    double abstol_solver;
    int iter_solver;
    double reltol_sundials;
    double abstol_sundials;
    double Epsilon;

    //Vortex parameters
    double Lx, Ly;    
    double Rx, Ry;    
    double Sigma;     
    double Gamma;     
    double Viscosity; 

    //Dimention scales
    double L_scale;   
    double T_scale;   

    void Adimentionalize(double new_T_scale, double new_L_scale, string print_file = ""){

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
};

//Solver for the temperature and salinity field
class EvolutionOperator : public TimeDependentOperator{
    public:
        //Initialization of the solver
        EvolutionOperator(Config config, ParMesh *pmesh):
            config(config),
            pmesh(pmesh),
            M_solver(MPI_COMM_WORLD), C_solver(MPI_COMM_WORLD)
        {
            //Create finite element spaces
            fecH1 = new H1_FECollection(config.order, pmesh->Dimension());
            fespaceH1 = new ParFiniteElementSpace(pmesh, fecH1);
            
            //Initialize operator
            height = width = fespaceH1->GetTrueVSize();
            t = 0.;
            type = EXPLICIT;
            eval_mode = NORMAL;

            //Initialize ODE solver
            ode_solver = new ForwardEulerSolver;
            ode_solver->Init(*this);

            //Initialize variables
            vorticity.SetSpace(fespaceH1); vorticity = 0.;
            velocity_0.SetSpace(fespaceH1); velocity_0 = 0.;
            velocity_1.SetSpace(fespaceH1); velocity_1 = 0.;

            Vorticity.SetSize(fespaceH1->GetTrueVSize()); Vorticity = 0.;
            Velocity_0.SetSize(fespaceH1->GetTrueVSize()); Velocity_0 = 0.;
            Velocity_1.SetSize(fespaceH1->GetTrueVSize()); Velocity_1 = 0.;
            Z.SetSize(fespaceH1->GetTrueVSize()); Z = 0.;
            B_0.SetSize(fespaceH1->GetTrueVSize()); B_0 = 0.;
            B_1.SetSize(fespaceH1->GetTrueVSize()); B_1 = 0.;

            //Set boundary dofs
            /****
             * Define essential boundary conditions
             *                 2
             *         /---------------\
             *         |               |
             *         |               |
             *       3 |               | 1
             *         |               | 
             *         |               |
             *         \---------------/
             *                 0         
             ****/
            Array<int> ess_tdof;
            Array<int> ess_bdr(pmesh->bdr_attributes.Max()); 
            ess_bdr = 1; ess_bdr[3] = 0;
            fespaceH1->GetEssentialTrueDofs(ess_bdr, ess_tdof);

            //Set initial and boundary conditions
            FunctionCoefficient initial_vorticity([=](const Vector &x){
                    return (config.Gamma/(M_PI*pow(config.Sigma, 2)))*exp(-pow(hypot(x(0)-config.Rx, x(1)-config.Ry)/config.Sigma, 2));
                    }); 
            ConstantCoefficient Zero(0.);
           
            vorticity.ProjectCoefficient(initial_vorticity);
            vorticity.ProjectBdrCoefficient(Zero, ess_bdr);
            vorticity.GetTrueDofs(Vorticity);

            //Coefficients
            FunctionCoefficient coeff_r([](const Vector &x){return x(0);});
            FunctionCoefficient coeff_nu_r([=](const Vector &x){return config.Viscosity*x(0);});
            FunctionCoefficient coeff_neg_r([](const Vector &x){return -x(0);});

            FunctionCoefficient coeff_inv_r([=](const Vector &x){return pow(x(0)+config.Epsilon, -1);});
            FunctionCoefficient coeff_nu_inv_r([=](const Vector &x){return config.Viscosity*pow(x(0)+config.Epsilon, -1);});

            //Create bilinear forms
            ParBilinearForm m(fespaceH1);
            m.AddDomainIntegrator(new MassIntegrator(coeff_r));
            m.Assemble();
            m.EliminateEssentialBC(ess_bdr);
            m.Finalize();
            M = m.ParallelAssemble();

            ParBilinearForm k(fespaceH1);
            k.AddDomainIntegrator(new DiffusionIntegrator(coeff_nu_r));
            k.AddDomainIntegrator(new MassIntegrator(coeff_nu_inv_r));
            k.Assemble();
            k.EliminateEssentialBC(ess_bdr, Operator::DIAG_ZERO);
            k.Finalize();
            K = k.ParallelAssemble();

            ParBilinearForm c0(fespaceH1);
            c0.AddDomainIntegrator(new DiffusionIntegrator(coeff_r));
            c0.AddDomainIntegrator(new MassIntegrator(coeff_inv_r));
            c0.Assemble();
            c0.EliminateEssentialBC(ess_bdr);
            c0.Finalize();
            C0 = c0.ParallelAssemble();
            C0_e = c0.ParallelAssembleElim();

            ParBilinearForm c1(fespaceH1);
            c1.AddDomainIntegrator(new DiffusionIntegrator(coeff_r));
            c1.Assemble();
            c1.EliminateEssentialBC(ess_bdr);
            c1.Finalize();
            C1 = c1.ParallelAssemble();
            C1_e = c1.ParallelAssembleElim();

            ParBilinearForm d0(fespaceH1);
            d0.AddDomainIntegrator(new DerivativeIntegrator(coeff_neg_r, 1));
            d0.AddDomainIntegrator(new MassIntegrator);
            d0.Assemble();
            d0.Finalize();
            D0 = d0.ParallelAssemble(); 

            ParBilinearForm d1(fespaceH1);
            d1.AddDomainIntegrator(new DerivativeIntegrator(coeff_r, 0));
            d1.Assemble();
            d1.Finalize();
            D1 = d1.ParallelAssemble(); 

            M_prec.SetPrintLevel(0);
            M_solver.SetPrintLevel(0);
            M_solver.SetPreconditioner(M_prec);
            M_solver.SetTol(config.reltol_solver);
            M_solver.SetAbsTol(config.abstol_solver);
            M_solver.SetMaxIter(config.iter_solver);

            M_solver.SetOperator(*M); 
            M_prec.SetOperator(*M);

            C_prec.SetPrintLevel(0);
            C_solver.SetPrintLevel(0);
            C_solver.SetPreconditioner(C_prec);
            C_solver.SetTol(config.reltol_solver);
            C_solver.SetAbsTol(config.abstol_solver);
            C_solver.SetMaxIter(config.iter_solver);

            //Get Initial Velocity
            D0->Mult(Vorticity, B_0);
            C0->EliminateBC(*C0_e, ess_tdof, Velocity_0, B_0);
            C_solver.SetOperator(*C0);
            C_prec.SetOperator(*C0);
            C_solver.Mult(B_0, Velocity_0);

            D1->Mult(Vorticity, B_1);
            C1->EliminateBC(*C1_e, ess_tdof, Velocity_1, B_1);
            C_solver.SetOperator(*C1);
            C_prec.SetOperator(*C1);
            C_solver.Mult(B_1, Velocity_1);
        };

        //Update of the solver on each iteration
        void Step(double &t, double &dt){
            ode_solver->Step(Vorticity, t, dt);
            //SolveVelocity();
        }
        //void SolveVelocity(){
        //    Vorticity_v.GetBlock(0) = Vorticity;
        //    D->Mult(Vorticity_v, B);
        //    C_solver.Mult(B, Velocity);
        //}
        //void SetParameters();

        //Time-evolving functions
        virtual void Mult(const Vector &X, Vector &dX_dt) const {
            Z = 0.; dX_dt = 0.;
            K->Mult(-1., X, 1., Z);
            M_solver.Mult(Z, dX_dt);
        }
        //virtual int SUNImplicitSetup(const Vector &X, const Vector &RHS, int j_update, int *j_status, double scaled_dt);
	    //virtual int SUNImplicitSolve(const Vector &X, Vector &X_new, double tol);
        
        ParGridFunction *GetVorticity() {return &vorticity;}
        ParGridFunction *GetVelocity_0() {return &velocity_0;}
        ParGridFunction *GetVelocity_1() {return &velocity_1;}

        void UpdateVorticity() {
            vorticity.SetFromTrueDofs(Vorticity);
        }
        void UpdateVelocity() {
            velocity_0.SetFromTrueDofs(Velocity_0);
            velocity_1.SetFromTrueDofs(Velocity_1);
        }

        virtual ~EvolutionOperator(){
            delete fecH1;
            delete fespaceH1;

            delete M;
            delete K;
            delete C0;
            delete C1;
            delete C0_e;
            delete C1_e;
            delete D0;
            delete D1;
        };
    protected:

        //Global parameters
        Config config;

        //FEM parameters
        ParMesh *pmesh = NULL;
        FiniteElementCollection *fecH1 = NULL;
        ParFiniteElementSpace *fespaceH1 = NULL;

        //ODE parameters
        ODESolver *ode_solver = NULL;

        //Auxiliar grid functions
        ParGridFunction vorticity, velocity_0, velocity_1;

        Vector Vorticity, Velocity_0, Velocity_1;
        Vector B_0, B_1;
        mutable Vector Z;

        HypreParMatrix *M = NULL, *K = NULL;

        HypreParMatrix *C0 = NULL, *C0_e = NULL;
        HypreParMatrix *C1 = NULL, *C1_e = NULL;
        HypreParMatrix *D0 = NULL, *D1   = NULL;

        HyprePCG M_solver, C_solver;
        HypreBoomerAMG M_prec, C_prec;

};

ParMesh* CreateMesh(Config config, string print_file = ""){
    int nx = ceil(config.Lx/config.Sigma);
    int ny = ceil(config.Ly/config.Sigma);
    Mesh mesh = Mesh::MakeCartesian2D(nx, ny, Element::TRIANGLE, false, config.Lx, config.Ly);
    int dim = mesh.Dimension();
    
    //Calculate how many serial refinements are needed
    //More than 1000 cells per processor
    int elements = mesh.GetNE();
    int min_elements = 1000.*config.nproc;
    double serial_refinements, parallel_refinements;
    if (min_elements > elements)
        serial_refinements = min(config.refinements, (int)floor(log(min_elements/elements)/(dim*log(2.))));
    else
        serial_refinements = 0;
    parallel_refinements = config.refinements - serial_refinements;
    
    //Refine mesh (serial)
    for (int ii = 0; ii < serial_refinements; ii++)
        mesh.UniformRefinement();
    
    //Make mesh (parallel), delete the serial
    ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, mesh);
    mesh.Clear();
    
    //Refine mesh (parallel)
    for (int ii = 0; ii < parallel_refinements; ii++)
        pmesh->UniformRefinement();
    
    //Calculate minimum size of elements
    double null, h_min;
    pmesh->GetCharacteristics(h_min, null, null, null);
    elements = pmesh->GetNE();

    if (config.master){
    cout << "\nMesh characteristics:\n"
         << "Serial refinements: " << serial_refinements << "\n"
         << "Parallel refinements: " << parallel_refinements << "\n"
         << "Total refinements: " << config.refinements << "\n"
         << "Number of Elements: " << elements << "\n"
         << "Mesh Size: " << h_min << " (" << h_min*config.L_scale << " mm)\n\n";

    ofstream out;
    out.open(print_file, std::ios::app);
    out << "\nMesh characteristics:\n"
        << "Serial refinements: " << serial_refinements << "\n"
        << "Parallel refinements: " << parallel_refinements << "\n"
        << "Total refinements: " << config.refinements << "\n"
        << "Number of Elements: " << elements << "\n"
        << "Mesh Size: " << h_min << " (" << h_min*config.L_scale << " mm)\n\n";
    out.close();
    }

    return pmesh;
}

int main(int argc, char *argv[]){

    //Create MPI session
    Mpi::Init(argc, argv);
    //Hypre::Init();

    //Initialize program parameters
    bool exit = false;
    Config config(Mpi::WorldRank(), Mpi::WorldSize(), argc, argv, exit, "results/graph/settings.txt");
    if (exit) return 1.;
    config.Adimentionalize(config.Sigma*config.Sigma/config.Gamma, config.Sigma, "results/graph/settings.txt");

    //Create mesh
    ParMesh *pmesh = CreateMesh(config, "results/graph/settings.txt");

    //Initialize operator
    EvolutionOperator evo_oper(config, pmesh);

    ParGridFunction *vorticity = evo_oper.GetVorticity();
    ParGridFunction *velocity_0 = evo_oper.GetVelocity_0();
    ParGridFunction *velocity_1 = evo_oper.GetVelocity_1();

    //Open the paraview output and print initial state
    string folder = "results/graph"; 
    ParaViewDataCollection paraview(folder, pmesh);
    paraview.SetDataFormat(VTKFormat::BINARY);
    paraview.SetLevelsOfDetail(config.order);
    paraview.RegisterField("Vorticity", vorticity);
    paraview.RegisterField("Velocity_0", velocity_0);
    paraview.RegisterField("Velocity_1", velocity_1);
    paraview.SetCycle(0);
    paraview.SetTime(0);
    paraview.Save();

    //Time evolution parameters
    double t = 0.;
    double dt = config.dt_init;
    int iteration = 0;
    int vis_iteration = 0;
    int vis_print = 0.;
    int vis_steps = config.vis_steps_max;
    bool last = false;

    if (config.master){
        cout.precision(4);
        cout << left << setw(12)
             << "--------------------------------------------------------------\n"
             << left << setw(12)
             << "Step" << setw(12)
             << "Printed" << setw(12)
             << "Dt(s)" << setw(12)
             << "Time(s)" << setw(12)
             << "Progress"
             << left << setw(12)
             << "\n--------------------------------------------------------------\n"
             << left << setw(12)
             << 0 << setw(12)
             << vis_print << setw(12)
             << dt*config.T_scale << setw(12)
             << t*config.T_scale  << setw(12)
             << "0%" << "\r";
    }

    for (iteration = 1, vis_iteration = 1; !last; iteration++, vis_iteration++){

        //Update iteration parameters
        last = (t >= config.t_final - 1e-8*config.dt_init);
        dt = min(dt, config.t_final - t);
        
        //Perform the time_step
        evo_oper.Step(t, dt);
        
        //Update visualization steps
        vis_steps = (dt == config.dt_init) ? config.vis_steps_max : int((config.dt_init/dt)*config.vis_steps_max);
        
        //Print visualization on certain steps
        if (last || vis_steps <= vis_iteration){
            //Update parameters
            vis_iteration = 0;
            vis_print++;

            //Update information
            evo_oper.UpdateVorticity();
            evo_oper.UpdateVelocity();
        
            //Print fields
            paraview.SetCycle(vis_print);
            paraview.SetTime(t*config.T_scale);
            paraview.Save();
        }
        
        //Print the system state
        double percentage = 100*t/config.t_final;
        string progress = to_string((int)percentage)+"%";
        if (config.master){
            cout.flush();
            cout.precision(4);
            cout << left << setw(12)
                 << iteration << setw(12)
                 << vis_print << setw(12)
                 << dt*config.T_scale << setw(12)
                 << t*config.T_scale  << setw(12)
                 << progress << "\r";
        }
    }

    //Free memory
    delete pmesh;

    return 0;
}
