#include "header.h"

//Create time evolution operator and configure initial variables
EvolutionOperator::EvolutionOperator(Config config, ParMesh *pmesh):
    config(config),
    pmesh(pmesh),
    M_solver(MPI_COMM_WORLD), T_solver(MPI_COMM_WORLD), C_solver(MPI_COMM_WORLD)
{
    //Create finite element spaces
    fecH1 = new H1_FECollection(config.order, pmesh->Dimension());
    fespaceH1 = new ParFiniteElementSpace(pmesh, fecH1);
    
    fecH1_v = new H1_FECollection(config.order, pmesh->Dimension());
    fespaceH1_v = new ParFiniteElementSpace(pmesh, fecH1_v, pmesh->Dimension());

    //Initialize linear variables
    vorticity.SetSpace(fespaceH1); vorticity = 0.;
    velocity.SetSpace(fespaceH1_v); velocity = 0.;
    vorticity_v.SetSpace(fespaceH1_v); vorticity_v = 0.;
    velocity_v.SetSpace(fespaceH1_v); velocity_v = 0.;

    Vorticity.SetSize(fespaceH1->GetTrueVSize()); Vorticity = 0.;
    Velocity.SetSize(fespaceH1_v->GetTrueVSize()); Velocity = 0.;
    Vorticity_v.SetSize(fespaceH1_v->GetTrueVSize()); Vorticity_v = 0.;
    Velocity_v.SetSize(fespaceH1_v->GetTrueVSize()); Velocity_v = 0.;
    Z.SetSize(fespaceH1->GetTrueVSize()); Z = 0.;
    B.SetSize(fespaceH1_v->GetTrueVSize()); B = 0.;
    
    //Initialize operator
    height = width = fespaceH1->GetTrueVSize();
    t = 0.;
    type = IMPLICIT;
    eval_mode = NORMAL;

    //Initialize ODE solver
    ode_solver = new CVODESolver(MPI_COMM_WORLD, CV_BDF);
    ode_solver->Init(*this);
    ode_solver->SetSStolerances(config.reltol_sundials, config.abstol_sundials);
    ode_solver->SetMaxStep(config.dt_init);
    ode_solver->SetStepMode(CV_ONE_STEP);

    //Initialize linear solvers
    C_prec.SetPrintLevel(0);
    C_solver.SetPrintLevel(0);
    C_solver.SetPreconditioner(C_prec);
    C_solver.SetTol(config.reltol_velocity);
    C_solver.SetAbsTol(config.abstol_velocity);
    C_solver.SetMaxIter(config.iter_velocity);

    M_prec.SetPrintLevel(0);
    M_solver.SetPrintLevel(0);
    M_solver.SetPreconditioner(M_prec);
    M_solver.SetTol(config.reltol_explicit);
    M_solver.SetAbsTol(config.abstol_explicit);
    M_solver.SetMaxIter(config.iter_explicit);

    T_prec.SetPrintLevel(0);
    T_solver.SetPrintLevel(0);
    T_solver.SetPreconditioner(T_prec);
    T_solver.SetTol(config.reltol_implicit);
    T_solver.SetAbsTol(config.abstol_implicit);
    T_solver.SetMaxIter(config.iter_implicit);
};

void EvolutionOperator::Setup(){
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
    ess_bdr.SetSize(pmesh->bdr_attributes.Max()); 
    ess_bdr = 1; ess_bdr[3] = 0;
    fespaceH1->GetEssentialTrueDofs(ess_bdr, ess_tdof);
    fespaceH1_v->GetEssentialTrueDofs(ess_bdr, ess_tdof_v);

    //Set initial and boundary conditions
    FunctionCoefficient initial_vorticity([=](const Vector &x){
            return (config.Gamma/(M_PI*pow(config.Sigma, 2)))*exp(-pow(hypot(x(0)-config.Rx, x(1)-config.Ry)/config.Sigma, 2));
            }); 
    ConstantCoefficient Zero(0.); 
    vorticity.ProjectCoefficient(initial_vorticity);
    vorticity.ProjectBdrCoefficient(Zero, ess_bdr);
    vorticity.GetTrueDofs(Vorticity);

    Vector zero_v(pmesh->Dimension()); zero_v = 0.;
    VectorConstantCoefficient Zero_v(zero_v);
    velocity.ProjectCoefficient(Zero_v);
    velocity.GetTrueDofs(Velocity);

    //Coefficients for further integrators
    FunctionCoefficient coeff_r([](const Vector &x){return x(0);});
    VectorFunctionCoefficient coeff_inv_r_hat(pmesh->Dimension(), [=](const Vector &x, Vector &f){
            f = 0.;
            f(0) = pow(x(0)+config.Epsilon, -1);
            });
    
    //Create bilinear forms for velocity solver
    ParBilinearForm c(fespaceH1_v);
    c.AddDomainIntegrator(new VectorDiffusionIntegrator(coeff_r));
    c.AddDomainIntegrator(new VectorMassIntegrator(coeff_inv_r_hat));
    c.Assemble();
    c.EliminateEssentialBC(ess_bdr);
    c.Finalize();
    C = c.ParallelAssemble();

    ParMixedBilinearForm d0(fespaceH1, fespaceH1_v);
    d0.AddDomainIntegrator(new GradientIntegrator(coeff_r));
    d0.Assemble();
    d0.Finalize();
    D0 = d0.ParallelAssemble();

    ParBilinearForm d1(fespaceH1_v);
    d1.AddDomainIntegrator(new VectorMassIntegrator);
    d1.Assemble();
    d1.Finalize();
    D1 = d1.ParallelAssemble();

    C_prec.SetOperator(*C);
    C_solver.SetOperator(*C); 

    //Calculate initial velocity
    SolveVelocity();

    //Create bilinear forms for vorticity evolution
    ParBilinearForm m(fespaceH1);
    m.AddDomainIntegrator(new MassIntegrator(coeff_r));
    m.Assemble();
    m.EliminateEssentialBC(ess_bdr);
    m.Finalize();
    M = m.ParallelAssemble();

    M_prec.SetOperator(*M);
    M_solver.SetOperator(*M); 
}

//Perform a time step in the system evolution
void EvolutionOperator::Step(double &t, double &dt){
    ode_solver->Step(Vorticity, t, dt);
    SolveVelocity();
}

//Calculate velocity corresponding to the current vorticity
void EvolutionOperator::SolveVelocity(){
    //Setup RHS vorticity-dependent
    UpdateVorticity();
    VectorArrayCoefficient coeff_vorticity_v(pmesh->Dimension());
    coeff_vorticity_v.Set(0, new GridFunctionCoefficient(&vorticity), true);
    vorticity_v.ProjectCoefficient(coeff_vorticity_v);
    vorticity_v.GetTrueDofs(Vorticity_v);

    B = 0.;
    D0->Mult(1., Vorticity, 1., B);
    D1->Mult(1., Vorticity_v, 1., B);
    B.SetSubVector(ess_tdof_v, 0.);

    //Solve for the rotated velocity and get the true velocity
    C_solver.Mult(B, Velocity_v);
    velocity_v.SetFromTrueDofs(Velocity_v);
    VectorGridFunctionCoefficient coeff_velocity_v(&velocity_v);
    MatrixFunctionCoefficient coeff_rot(pmesh->Dimension(), [](const Vector &x, DenseMatrix &f){ f = 0.; f(0,1) = -1.; f(1,0) = 1.;});
    MatrixVectorProductCoefficient coeff_rot_velocity_v(coeff_rot, coeff_velocity_v);
    velocity.ProjectCoefficient(coeff_rot_velocity_v);
    velocity.GetTrueDofs(Velocity);

    //Update K matrix
    FunctionCoefficient coeff_r([](const Vector &x){return x(0);});
    FunctionCoefficient coeff_nu_r([=](const Vector &x){return config.Viscosity*x(0);});
    FunctionCoefficient coeff_nu_inv_r([=](const Vector &x){return config.Viscosity*pow(x(0)+config.Epsilon, -1);});
    VectorFunctionCoefficient coeff_neg_r_hat(pmesh->Dimension(), [](const Vector &x, Vector &f){
            f = 0.; f(0) = -1.;
            });
    VectorGridFunctionCoefficient coeff_velocity(&velocity);
    ScalarVectorProductCoefficient coeff_r_velocity(coeff_r, coeff_velocity);
    InnerProductCoefficient coeff_neg_vr(coeff_neg_r_hat, coeff_velocity);

    if (K) delete K;
    ParBilinearForm k(fespaceH1);
    k.AddDomainIntegrator(new DiffusionIntegrator(coeff_nu_r));
    k.AddDomainIntegrator(new MassIntegrator(coeff_nu_inv_r));
    k.AddDomainIntegrator(new ConvectionIntegrator(coeff_r_velocity));
    k.AddDomainIntegrator(new MassIntegrator(coeff_neg_vr));
    k.Assemble();
    k.EliminateEssentialBC(ess_bdr, Operator::DIAG_ZERO);
    k.Finalize();
    K = k.ParallelAssemble();
}

//From  M(dX_dt) + K(X) = 0
//Solve M(dX_dt) + K(X) = 0 for dX_dt
void EvolutionOperator::Mult(const Vector &X, Vector &dX_dt) const {
    Z = 0.; dX_dt = 0.;
    K->Mult(-1., X, 1., Z);
    M_solver.Mult(Z, dX_dt);
}

//Setup the ODE Jacobian T = M + dt*K
int EvolutionOperator::SUNImplicitSetup(const Vector &X, const Vector &RHS, int j_update, int *j_status, double scaled_dt){
    if (T) delete T;
    T = Add(1., *M, scaled_dt, *K);
    T_prec.SetOperator(*T);
    T_solver.SetOperator(*T);
    *j_status = 1;
    return 0;
}

//From  M(dX_dt) + K(X) = 0
//Solve M(X_new - X) + dt*K(X_new) = 0 for X_new
int EvolutionOperator::SUNImplicitSolve(const Vector &X, Vector &X_new, double tol){
    Z = 0.; X_new = X;
    M->Mult(X, Z);          
    T_solver.Mult(Z, X_new);
    return 0; 
}

//Free memory
EvolutionOperator::~EvolutionOperator(){
    delete fecH1;
    delete fecH1_v;
    delete fespaceH1;
    delete fespaceH1_v;

    delete ode_solver;

    delete M;
    delete K;
    delete T;
    delete C;
    delete D0;
    delete D1;
};
