#pragma once

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
    Config(int pid, int nproc, int argc, char *argv[], bool &exit, string print_file);
    void Adimentionalize(double new_T_scale, double new_L_scale, string print_file);

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
};

//Solver for the time evolution
class EvolutionOperator : public TimeDependentOperator{
    public:
        //Initialization of the solver
        EvolutionOperator(Config config, ParMesh *pmesh);
        void Setup();

        //Update of the solver on each iteration
        void Step(double &t, double &dt);
        void SolveVelocity();

        //Time-evolving functions
        virtual void Mult(const Vector &X, Vector &dX_dt) const;
        virtual int SUNImplicitSetup(const Vector &X, const Vector &RHS, int j_update, int *j_status, double scaled_dt);
        virtual int SUNImplicitSolve(const Vector &X, Vector &X_new, double tol);

        //Transfer information functions
        ParGridFunction* GetVorticity() {return &vorticity;}
        ParGridFunction* GetVelocity() {return &velocity;}
        void UpdateVorticity() {vorticity.SetFromTrueDofs(Vorticity);}
        void UpdateVelocity() {velocity.SetFromTrueDofs(Velocity);}

        virtual ~EvolutionOperator();
    protected:

        //Global parameters
        Config config;

        //FEM parameters
        ParMesh *pmesh = NULL;
        FiniteElementCollection *fecH1 = NULL, *fecH1_v = NULL;
        ParFiniteElementSpace *fespaceH1 = NULL, *fespaceH1_v = NULL;

        Array<int> ess_bdr;
        Array<int> ess_tdof;
        Array<int> ess_tdof_v;

        //ODE parameters
        CVODESolver *ode_solver = NULL;

        //Linear variables
        ParGridFunction vorticity, velocity, vorticity_v, velocity_v;

        Vector Vorticity, Velocity, Vorticity_v, Velocity_v;
        mutable Vector Z, B;

        //Bilinear variables
        HypreParMatrix *M = NULL, *K = NULL, *T = NULL;
        HypreParMatrix *C = NULL, *D0 = NULL, *D1 = NULL;

        //Solvers
        HyprePCG M_solver, C_solver;
        GMRESSolver T_solver;
        HypreBoomerAMG M_prec, T_prec, C_prec;
};

ParMesh* CreateMesh(Config config, string print_file);
