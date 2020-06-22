#include "../include/linear_convection_2d.h"

std::string int_to_string (unsigned int value, const unsigned int digits) {
    std::string lc_string = std::to_string(value);  
    
    if (lc_string.size() < digits) {
        // We have to add the padding zeroes in front of the number
        const unsigned int padding_position = (lc_string[0] == '-')
                                                ?
                                                1
                                                :
                                                0;

        const std::string padding(digits - lc_string.size(), '0');
        lc_string.insert(padding_position, padding);
    }
         
    return lc_string;
}


// Plotting commands 

void LinearConvection2D ::plot_matplotlib() const {
    
    std::ofstream mpl_x;
    std::ofstream mpl_y;
    std::ofstream mpl_u;
    
    mpl_x.open ("plots/x.dat");
    mpl_y.open ("plots/y.dat");
    mpl_u.open ("plots/u.dat");

    if ( !(mpl_x.is_open()) || !(mpl_y.is_open()) || !(mpl_u.is_open()) ) {
        
        ThrowFileOpeningError();
    }
    
    else {
        

        for (unsigned int i=0; i<Params.N_x; i++)  {
            mpl_x << x[i] << std::endl; 
        }

        for (unsigned int j=0; j<Params.N_y; j++)  {
            mpl_y << y[j] << std::endl; 
        }

        for (unsigned int i=0; i<Params.N_x; i++) {
		    for (unsigned int j=0; j<Params.N_y; j++) {
                mpl_u << U(i, j) << " " ;
		    }
		    
            mpl_u << std::endl; 
	    }
    }

    mpl_u.close();
}

// ======================================================================================

// GNU Plot 

void LinearConvection2D::plot_gnuplot() const {
    
    std::ofstream gpl_u;

    gpl_u.open ("plots/gpl_u.dat");

    if ( !(gpl_u.is_open()) ) {
        
        ThrowFileOpeningError(); 
    }

    else {
    
        
        for (unsigned int i=0; i < Params.N_x; i++) {
		    for (unsigned int j=0; j < Params.N_y; j++) {
                gpl_u << x[i] << " " << y[j] << " " << U(i, j) << std::endl;
	    	}
		    gpl_u << std::endl;
        }
    }

    gpl_u.close();
}

// ======================================================================================

void LinearConvection2D::plot_tecplot(unsigned int i, const unsigned int digits) const {

    std::ofstream tpl;
    const std::string filename = "../plots/plot_" + int_to_string (i, digits) + ".dat";
    tpl.open (filename);
    tpl.flags( std::ios::dec | std::ios::scientific );
	tpl.precision(6);

    if ( !(tpl.is_open()) ) {
        ThrowFileOpeningError();
    }

    else {
        
        tpl << "TITLE = \"Linear Convection 2D\" " << std::endl
		    << "VARIABLES = \"x\", \"y\", \"U\" " << std::endl;
        tpl << "Zone I = " << Params.N_y << " J = " << Params.N_x << std::endl; 

        for (unsigned int i=0; i < Params.N_x; i++) {
		    for (unsigned int j=0; j < Params.N_y; j++) {
                tpl << x[i] << "\t" << y[j] << "\t" << U(i, j) << std::endl;
		    }
	    }
    }

    tpl.close(); 
}

void LinearConvection2D::plot_vtk(double time, unsigned int i, const unsigned int digits) const {

    std::ofstream vtk;
    const std::string filename = "../plots/plot_" + int_to_string (i, digits) + ".vtk";
    vtk.open (filename);
    vtk.flags( std::ios::dec | std::ios::scientific );
	vtk.precision(6);

    if ( !(vtk.is_open()) ) {
        ThrowFileOpeningError();
    }
    
    else {
	
		vtk << "# vtk DataFile Version 2.0" << "\n"; 
        vtk << "2D Linear Convection" << "\n";
        vtk << "ASCII" << "\n";
        vtk << "\nDATASET STRUCTURED_GRID" << "\n"; 
        vtk << "\nFIELD FieldData 1" << "\n"; 
        vtk << "TIME 1 1 double" << "\n"; 
        vtk << time << "\n"; 
        vtk << "\nDIMENSIONS " <<  Params.N_y << " " << Params.N_x << " " << 1 << "\n"; 
        vtk << "POINTS " << Params.N_x*Params.N_y <<  " double" << "\n";
        vtk << "\n"; 
        
        for (unsigned int i = 0; i < Params.N_x; i++) {
            for (unsigned int j = 0; j < Params.N_y; j++) {
                vtk << x[i] << " " << y[j] << " " << 0.0 << "\n"; 
            }
        }
        
        vtk << "\nPOINT_DATA " << Params.N_x*Params.N_y << "\n"; 
        
        vtk << "\nSCALARS U double 1" << "\n"; 
        vtk << "LOOKUP_TABLE default" << "\n"; 
        vtk << "\n";
        
        for (unsigned int i = 0; i < Params.N_x; i++) {
            for (unsigned int j = 0; j < Params.N_y; j++) {
                vtk << U(i,j) << "\n"; 
            }
        }
	}

    vtk.close(); 
}


void LinearConvection2D::write_restart_file(double time, unsigned int counter) const {

    std::ofstream restart;
    restart.open ("../plots/restart.dat");
    restart.flags( std::ios::dec | std::ios::scientific );
	restart.precision(8);



    if ( !(restart.is_open())) {
        ThrowFileOpeningError();
    }

    else {

        restart << time << std::endl; 
        restart << counter << std::endl; 

        for (unsigned int i=0; i < Params.N_x; ++i) {
		    for (unsigned int j=0; j < Params.N_y; ++j) {
                restart << U(i, j) << std::endl;
		    }
	    }
    }

    restart.close();
}
