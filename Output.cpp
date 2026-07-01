#include "Output.h"
#include <list>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <fstream>

using namespace std;
using namespace Const;

void Output::write2file_Object()
{


    MPI_Status status;
    MPI_File fh;

    // -------------------------------------------------
    // Open MPI file
    // -------------------------------------------------
    MPI_File_open(
        MPI_COMM_WORLD,
        "results/Object.dat",
        MPI_MODE_CREATE | MPI_MODE_WRONLY,
        MPI_INFO_NULL,
        &fh
    );

       // -------------------------------------------------
    // Fixed line size for MPI offsets
    // -------------------------------------------------
    const int lineSize = 64;
 // -------------------------------------------------
    // Each MPI rank writes its own contiguous block
    // -------------------------------------------------
    MPI_Offset displace =
        (MPI_Offset)domain.mpi_rank *
        Nr * Nz_phys * lineSize;

    MPI_File_set_view(
        fh,
        displace,
        MPI_CHAR,
        MPI_CHAR,
        "native",
        MPI_INFO_NULL
    );

// -------------------------------------------------
    // Output buffer
    // -------------------------------------------------
    char buf[lineSize];

    // -------------------------------------------------
    // Loop over local domain
    // i -> radial
    // k -> axial(z)
    // -------------------------------------------------
    for (int k = 0; k < Nz_phys; k++)
    {
        for (int i = 0; i < Nr; i++)
        {
        // clear buffer
            char buf[lineSize] = {0};
            memset(buf, 0, lineSize-1);
            vectorCyl Vcyl;
            vector3d  V;
            
             // -----------------------------------------
            // Physical coordinates
            // -----------------------------------------
            Vcyl.r     = i * domain.dh[0];
            Vcyl.theta = 0.0;

            // global z coordinate
            Vcyl.zee =
                (k + domain.mpi_rank * (Nz_phys-1))
                * domain.dh[2];

            // convert to Cartesian for visualization
            V = Cyl2Cart(Vcyl);

            // -----------------------------------------
            // Write formatted data
            // -----------------------------------------
            int nchar = snprintf(
                buf,
                lineSize-2,
                "%.4e %.4e %.4e",
                double(Vcyl.zee), 
                double(Vcyl.r),
                double(domain.object[i + k*Nr]) 
                );
            // safety
			std:: fill(std::begin(buf)+nchar, std::end(buf)-1,' ');
            // newline at end
			buf[lineSize-1] = '\n';
			// -----------------------------------------
            // Fixed-size MPI write
            // IMPORTANT:
            // always write lineSize bytes
            // -----------------------------------------
            MPI_File_write(
                fh,
                buf,
                lineSize,
                MPI_CHAR,
                &status
            );
        	}}
        	// -------------------------------------------------
    // Close file
    // -------------------------------------------------
    MPI_File_close(&fh);

}

void Output::write2file_SField2d(Field* F, char* fileName )
{
    MPI_Status status;
    MPI_File fh;

    // -------------------------------------------------
    // Open MPI file
    // -------------------------------------------------
    MPI_File_open(
        MPI_COMM_WORLD,
        fileName,
        MPI_MODE_CREATE | MPI_MODE_WRONLY,
        MPI_INFO_NULL,
        &fh
    );

    // -------------------------------------------------
    // Fixed line size for MPI offsets
    // -------------------------------------------------
    const int lineSize = 64;
 // -------------------------------------------------
    // Each MPI rank writes its own contiguous block
    // -------------------------------------------------
    MPI_Offset displace =
        (MPI_Offset)domain.mpi_rank *
        Nr * Nz_phys * lineSize;

    MPI_File_set_view(
        fh,
        displace,
        MPI_CHAR,
        MPI_CHAR,
        "native",
        MPI_INFO_NULL
    );

    // -------------------------------------------------
    // Output buffer
    // -------------------------------------------------
    char buf[lineSize];

    // -------------------------------------------------
    // Loop over local domain
    // i -> radial
    // k -> axial(z)
    // -------------------------------------------------
    for (int k = 0; k < Nz_phys; k++)
    {
        for (int i = 0; i < Nr; i++)
        {
            // clear buffer
            char buf[lineSize] = {0};
            memset(buf, 0, lineSize-1);
            vectorCyl Vcyl;

            // -----------------------------------------
            // Physical coordinates
            // -----------------------------------------
            Vcyl.r     = i * domain.dh[0];
            Vcyl.theta = 0.0;
            // global z coordinate
            Vcyl.zee =
                 (k + domain.mpi_rank * (Nz_phys-1))
                * domain.dh[2];

            // -----------------------------------------
            // Write formatted data
            // -----------------------------------------
            int nchar = snprintf(
                buf,
                lineSize-2,
                "%.4e %.4e %.4e", 
                double(Vcyl.zee),
                double(Vcyl.r),
                double(F->data[i +(k+1)*Nr]) 
                );
            // safety
			std:: fill(std::begin(buf)+nchar, std::end(buf)-1,' ');
            // newline at end
			buf[lineSize-1] = '\n';
			// -----------------------------------------
            // Fixed-size MPI write
            // IMPORTANT:
            // always write lineSize bytes
            // -----------------------------------------
            MPI_File_write(
                fh,
                buf,
                lineSize,
                MPI_CHAR,
                &status
            );
        	}}
        	// -------------------------------------------------
    // Close file
    // -------------------------------------------------
    MPI_File_close(&fh);
}


void Output::write2file_VField2d(FieldCyl* F, const char* fileName)
{
    MPI_Status status;
    MPI_File fh;

    // -------------------------------------------------
    // Open MPI file
    // -------------------------------------------------
    MPI_File_open(
        MPI_COMM_WORLD,
        fileName,
        MPI_MODE_CREATE | MPI_MODE_WRONLY,
        MPI_INFO_NULL,
        &fh
    );

    // -------------------------------------------------
    // Fixed line size for MPI offsets
    // -------------------------------------------------
    const int lineSize = 64;

    // -------------------------------------------------
    // Each MPI rank writes its own contiguous block
    // -------------------------------------------------
   MPI_Offset displace =
        (MPI_Offset)domain.mpi_rank *
        Nr * Nz_phys * lineSize;

    MPI_File_set_view(
        fh,
        displace,
        MPI_CHAR,
        MPI_CHAR,
        "native",
        MPI_INFO_NULL
    );

    // -------------------------------------------------
    // Output buffer
    // -------------------------------------------------
    char buf[lineSize];

    // -------------------------------------------------
    // Loop over local domain
    // i -> radial
    // k -> axial(z)
    // -------------------------------------------------
    for (int k = 0; k < Nz_phys; k++)
    {
        for (int i = 0; i < Nr; i++)
        {
            // clear buffer
            char buf[lineSize] = {0};
            memset(buf, 0, lineSize-1);

            vectorCyl Vcyl;

           // -----------------------------------------
            // Physical coordinates
            // -----------------------------------------
            Vcyl.r     = i * domain.dh[0];
            Vcyl.theta = 0.0;
            // global z coordinate
            Vcyl.zee =
                 (k + domain.mpi_rank * (Nz_phys-1))
                * domain.dh[2];

            // -----------------------------------------
            // Write formatted data
            // -----------------------------------------
            int nchar = snprintf(
                buf,
                lineSize-2,
                "%.4e %.4e %.4e %.4e %.4e",
                double(Vcyl.zee),
                double(Vcyl.r),
                double(F->data[i+(k+1)*Nr].r),
                double(F->data[i+(k+1)*Nr].theta),
                double(F->data[i+(k+1)*Nr].zee)
            );

            // safety
			std:: fill(std::begin(buf)+nchar, std::end(buf)-1,' ');
            // newline at end
			buf[lineSize-1] = '\n';
            
            // -----------------------------------------
            // Fixed-size MPI write
            // IMPORTANT:
            // always write lineSize bytes
            // -----------------------------------------
            MPI_File_write(
                fh,
                buf,
                lineSize,
                MPI_CHAR,
                &status
            );
        }
    }

    // -------------------------------------------------
    // Close file
    // -------------------------------------------------
    MPI_File_close(&fh);
}
/*
void Output::write2file_VField2d(FieldCyl* F, char* fileName )
{
MPI_Status   status; 
MPI_File fh;
MPI_File_open (MPI_COMM_WORLD, fileName, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh); 

MPI_Offset displace = domain.mpi_rank*Nr*Nzee*64*2  ; // start of the view for each processor 

MPI_File_set_view (fh , displace , MPI_CHAR, MPI_CHAR, "native" ,MPI_INFO_NULL);

char buf[64];
int j ;
	for (int k=0; k<Nzee; k++)
		{
        		for (int i=0; i<Nr; i++)
            		{
			buf[64] = {0};
			j = round(90.*PI/180. / domain.dh[1]);
			buf[50] = {0};
			vectorCyl Vcyl;
			vector3d V;
			Vcyl.r = i* domain.dh[0];	
			Vcyl.theta = j* domain.dh[1];	
			Vcyl.zee = (k + domain.mpi_rank*(Nzee-1))* domain.dh[2];
			V= Cyl2Cart(Vcyl);

			int nchar =  snprintf(buf,61, "%.4e %.4e %.4e %.4e %.4e", 
				double( V.y),
			 		double (V.z), 
						double( F->data[i + Nr*j +(k+1)*Nr*Ntheta].r), 
							double( F->data[i + Nr*j +(k+1)*Nr*Ntheta].theta), 
								double( F->data[i + Nr*j +(k+1)*Nr*Ntheta].zee)
						);

			std::fill(std::begin(buf)+nchar, std::end(buf)-1, ' ');
			buf[63]= '\n';

			MPI_File_write(fh,buf,strlen(buf), MPI_CHAR, &status);
		
			j = round(270.*PI/180. / domain.dh[1]);
			buf[50] = {0};
			
			Vcyl.r = i* domain.dh[0];	
			Vcyl.theta = j* domain.dh[1];	
			Vcyl.zee = (k + domain.mpi_rank*(Nzee-1))* domain.dh[2];
			V= Cyl2Cart(Vcyl);

			nchar =  snprintf(buf,61, "%.4e %.4e %.4e %.4e %.4e", 
				double( V.y),
			 		double (V.z), 
						double( F->data[i + Nr*j +(k+1)*Nr*Ntheta].r), 
							double( F->data[i + Nr*j +(k+1)*Nr*Ntheta].theta), 
								double( F->data[i + Nr*j +(k+1)*Nr*Ntheta].zee)
						);

			std::fill(std::begin(buf)+nchar, std::end(buf)-1, ' ');
			buf[63]= '\n';

			MPI_File_write(fh,buf,strlen(buf), MPI_CHAR, &status);
        	}}

MPI_File_close(&fh); 
//}

*/

void Output:: outputPart(
	    std::vector<std::unique_ptr<Species>> &speciesList,
	    const char* fileName
		)
{		
 	MPI_Status status;
    MPI_File fh;

    // -------------------------------------------------
    // Open MPI file
    // -------------------------------------------------
    MPI_File_open(
        MPI_COMM_WORLD,
        fileName,
        MPI_MODE_CREATE | MPI_MODE_WRONLY,
        MPI_INFO_NULL,
        &fh
    );

    // -------------------------------------------------
    // Fixed line size for MPI offsets
    // -------------------------------------------------
    const int lineSize = 64;
    
	double sizePart=0;

	for (auto &speciesPtr : speciesList)
		{
    	KineticSpecies* ks =
        	dynamic_cast<KineticSpecies*>(
        	    speciesPtr.get()
        	);

	    if (ks == nullptr)
	        continue;

	    sizePart += ks->particles.size();
		}
		
	std :: cout<<sizePart <<" = size part \t" <<domain.mpi_rank<<" \n";

	int *tmpOff;
	tmpOff = new int[domain.mpi_size];
	for (int i =0; i <domain.mpi_size; i ++)
		{
		tmpOff [i] = sizePart;
		}

	for (int i =0; i <domain.mpi_size; i ++)
		{
		std :: cout<< tmpOff [i] << " = tmp \t";
		}
	 std ::cout <<domain.mpi_rank<<"\n";

	MPI_Alltoall(&tmpOff[0], 1, MPI_INT, tmpOff , 1, MPI_INT, MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);

	for (int i =0; i <domain.mpi_size; i ++)
		{
		std :: cout<< tmpOff [i] << " = tmp \t";
		}
 	std ::cout<<"\n";
	int offSet_t=0;
	for (int i =0; i <domain.mpi_rank; i ++)
		{
		offSet_t = offSet_t +tmpOff[i];
		}

	std ::cout <<offSet_t<<"\t"<<domain.mpi_rank<<" = my offset\n";

	 // -------------------------------------------------
    // Each MPI rank writes its own contiguous block
    // -------------------------------------------------
    MPI_Offset displace =offSet_t*lineSize;

    MPI_File_set_view(
        fh,
        displace,
        MPI_CHAR,
        MPI_CHAR,
        "native",
        MPI_INFO_NULL
    );
    // -------------------------------------------------
    // Output buffer
    // -------------------------------------------------
    char buf[lineSize];
	int k=0;
    for (auto &speciesPtr : speciesList)
		{
    	KineticSpecies* ks =
        	dynamic_cast<KineticSpecies*>(
        	    speciesPtr.get()
        	);

	    if (ks == nullptr)
	        continue;
		std :: cout<< "writing "<<ks->particles.size()<< "  particle from "<< domain.mpi_rank<< " =rank \n";
	for (int i=0; i<ks->particles.size(); i++)
		{	
		// clear buffer
            char buf[lineSize] = {0};
            memset(buf, 0, lineSize-1);
		// -----------------------------------------
            // Write formatted data
            // -----------------------------------------
            int nchar = snprintf(
                buf,
                lineSize-2,
                "%.4e %.4e %.4e %.4e", 
                double(ks->particles[i].pos.y),
                	double (ks->particles[i].pos.z), 
                		double (ks->particles[i].pos.x), 
                			double(k) 
            );

            // safety
			std:: fill(std::begin(buf)+nchar, std::end(buf)-1,' ');
            // newline at end
			buf[lineSize-1] = '\n';


			// -----------------------------------------
            // Fixed-size MPI write
            // IMPORTANT:
            // always write lineSize bytes
            // -----------------------------------------
            MPI_File_write(
                fh,
                buf,
                lineSize,
                MPI_CHAR,
                &status
            );
        	}
	k=k+1;
	}
MPI_File_close(&fh); 	    
}

void Output::outputPart2(
        KineticSpecies &species,
        const char* fileName
        )
{
    MPI_Status status;
    MPI_File fh;

    //--------------------------------------------------
    // open MPI file
    //--------------------------------------------------

    MPI_File_open(
        MPI_COMM_WORLD,
        fileName,
        MPI_MODE_CREATE | MPI_MODE_WRONLY,
        MPI_INFO_NULL,
        &fh
    );

    //--------------------------------------------------
    // fixed line size
    //--------------------------------------------------

    const int lineSize = 128;

    //--------------------------------------------------
    // local particle count
    //--------------------------------------------------

    int localCount =
        static_cast<int>(
            species.particles.size()
        );

    //--------------------------------------------------
    // gather counts from all ranks
    //--------------------------------------------------

    std::vector<int> counts(
        domain.mpi_size,
        0
    );

    MPI_Allgather(
        &localCount,
        1,
        MPI_INT,
        counts.data(),
        1,
        MPI_INT,
        MPI_COMM_WORLD
    );

    //--------------------------------------------------
    // compute offset
    //--------------------------------------------------

    int offsetParticles = 0;

    for (int r=0; r<domain.mpi_rank; r++)
    {
        offsetParticles += counts[r];
    }

    MPI_Offset displacement =
        static_cast<MPI_Offset>(
            offsetParticles
        ) * lineSize;

    //--------------------------------------------------
    // set MPI file view
    //--------------------------------------------------

    MPI_File_set_view(
        fh,
        displacement,
        MPI_CHAR,
        MPI_CHAR,
        "native",
        MPI_INFO_NULL
    );

    //--------------------------------------------------
    // write particles
    //--------------------------------------------------

    char buf[lineSize];

    for (size_t p=0;
         p<species.particles.size();
         p++)
    {
        const Particle& part =
            species.particles[p];

        //--------------------------------------------------
        // clear buffer
        //--------------------------------------------------

        char buf[lineSize] = {0};
        memset(buf, 0, lineSize-1);

        //--------------------------------------------------
        // format line
        //--------------------------------------------------

        int nchar = snprintf(
            buf,
            lineSize-2,
            "%.4e %.4e %.4e %.4e %.4e %.4e %.d",
            part.pos.x,
            part.pos.y,
            part.pos.z,
            part.vel.x,
            part.vel.y,
            part.vel.z,
            part.STATUS
        );

        //--------------------------------------------------
        // newline
        //--------------------------------------------------

       // safety
			std:: fill(std::begin(buf)+nchar, std::end(buf)-1,' ');
            // newline at end
			buf[lineSize-1] = '\n';

        //--------------------------------------------------
        // MPI write
        //--------------------------------------------------

        MPI_File_write(
            fh,
            buf,
            lineSize,
            MPI_CHAR,
            &status
        );
    }

    //--------------------------------------------------
    // close file
    //--------------------------------------------------

    MPI_File_close(&fh);
}

