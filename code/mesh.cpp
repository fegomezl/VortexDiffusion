#include "header.h"

//Create or read mesh
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

    //Print important variables
    if (config.master){
    cout << "Mesh characteristics:\n"
         << "Serial refinements: " << serial_refinements << "\n"
         << "Parallel refinements: " << parallel_refinements << "\n"
         << "Total refinements: " << config.refinements << "\n"
         << "Number of Elements: " << elements << "\n"
         << "Mesh Size: " << h_min << " (" << h_min*config.L_scale << " mm)\n\n";

    ofstream out;
    out.open(print_file, std::ios::app);
    out << "Mesh characteristics:\n"
        << "Serial refinements: " << serial_refinements << "\n"
        << "Parallel refinements: " << parallel_refinements << "\n"
        << "Total refinements: " << config.refinements << "\n"
        << "Number of Elements: " << elements << "\n"
        << "Mesh Size: " << h_min << " (" << h_min*config.L_scale << " mm)\n\n";
    out.close();
    }

    return pmesh;
}
