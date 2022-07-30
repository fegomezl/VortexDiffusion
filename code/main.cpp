#include "header.h"

int main(int argc, char *argv[]){

    //Create MPI session
    Mpi::Init(argc, argv);

    //Initialize program parameters
    bool exit = false;
    Config config(Mpi::WorldRank(), Mpi::WorldSize(), argc, argv, exit, "results/graph/settings.txt");
    if (exit) return 1.;
    config.Adimentionalize(config.Sigma*config.Sigma/config.Gamma, config.Sigma, "results/graph/settings.txt");

    //Create mesh
    ParMesh *pmesh = CreateMesh(config, "results/graph/settings.txt");

    //Initialize operator and fields
    EvolutionOperator evo_oper(config, pmesh);

    ParGridFunction *vorticity = evo_oper.GetVorticity(); evo_oper.UpdateVorticity();
    ParGridFunction *velocity = evo_oper.GetVelocity(); evo_oper.UpdateVelocity();

    //Open the paraview output and print initial state
    string folder = "results/graph"; 
    ParaViewDataCollection paraview(folder, pmesh);
    paraview.SetDataFormat(VTKFormat::BINARY);
    paraview.SetLevelsOfDetail(config.order);
    paraview.RegisterField("Vorticity", vorticity);
    paraview.RegisterField("Velocity", velocity);
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

    //Start program check
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
