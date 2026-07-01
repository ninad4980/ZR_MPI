/*definitions for species functions*/
#include "Species.h"
#include "InputPar.h"
//#include "PotentialSolver.h"	/*for vec*/
#include <mpi.h>

/*adds a new particle, rewinding velocity by half dt*/
void KineticSpecies::addParticle(
    const vector3d& position,
    vector3d velocity
)
{
    double lc[3];

    domain.Posv3d2lc(lc, position);

    //--------------------------------------------------
    // inside physical domain only
    //--------------------------------------------------

    if (lc[0] < 0.0 || lc[0] >= domain.Nr-1)
        return;

    if (lc[2] < 0.0 || lc[2] >= domain.Nz_phys-1)
        return;

    //--------------------------------------------------
    // gather electric field
    //--------------------------------------------------

    vectorCyl ef_cyl =
        domain.Estatic.gather(
            lc,
            0.0,
            domain.dh[0]
        );

    //--------------------------------------------------
    // convert E to Cartesian
    //--------------------------------------------------

    vectorCyl pos_c = Cart2Cyl(position);

    vector3d ef =
        cylV2cart(ef_cyl, pos_c);

    //--------------------------------------------------
    // velocity rewind
    //--------------------------------------------------

    velocity =
        velocity
        -
        ef * (charge/mass)
        * (0.5*dt);

    //--------------------------------------------------
    // create particle
    //--------------------------------------------------

    Particle part;

    part.pos    = position;
    part.vel    = velocity;
    part.STATUS = PARTICLE_ACTIVE;

    particles.emplace_back(part);
}

void KineticSpecies::applyDomainConditions()
{
    sendL = 0;
    sendR = 0;

    for (Particle &part : particles)
    {
        if (part.STATUS != PARTICLE_ACTIVE)
            continue;

        double z = part.pos.z;

        //--------------------------------------------------
        // crossing right MPI boundary
        //--------------------------------------------------

        if (z > domain.rd[2] &&
            domain.mpi_rank < domain.mpi_size-1)
        {
            part.STATUS = PARTICLE_SEND_RIGHT;
            sendR++;
        }

        //--------------------------------------------------
        // crossing left MPI boundary
        //--------------------------------------------------

        else if (z < domain.r0[2] &&
                 domain.mpi_rank > 0)
        {
            part.STATUS = PARTICLE_SEND_LEFT;
            sendL++;
        }
    }
}
static vector3d computeWallNormal(
    const vector3d& pos,
    Domain& domain
)
{
    double lc[3];

    domain.Posv3d2lc(
        lc,
        pos
    );

    int i =
        static_cast<int>(
            floor(lc[0])
        );

    int k =
        static_cast<int>(
            floor(lc[2])
        );

    vector3d n;

    n.x = 0.0;
    n.y = 0.0;
    n.z = 0.0;

    vectorCyl pc =
        Cart2Cyl(pos);

    //--------------------------------------------------
    // inward radial direction
    //--------------------------------------------------

    double erx = 0.0;
    double ery = 0.0;

    if (pc.r > 1e-12)
    {
        erx = pos.x / pc.r;
        ery = pos.y / pc.r;
    }

    //--------------------------------------------------
    // neighbor toward smaller r
    //--------------------------------------------------

    if (i > 0)
    {
        int id =
            (i-1)
            +
            k*domain.Nr;

        int obj =
            domain.object[id];

        //----------------------------------------------
        // plasma is inward
        //----------------------------------------------

        if (
            obj == OBJ_PLASMA ||
            obj == OBJ_EXTRACTION
        )
        {
            n.x -= erx;
            n.y -= ery;
        }
    }

    //--------------------------------------------------
    // neighbor toward larger r
    //--------------------------------------------------

    if (i < domain.Nr-1)
    {
        int id =
            (i+1)
            +
            k*domain.Nr;

        int obj =
            domain.object[id];

        //----------------------------------------------
        // plasma is outward
        //----------------------------------------------

        if (
            obj == OBJ_PLASMA ||
            obj == OBJ_EXTRACTION
        )
        {
            n.x += erx;
            n.y += ery;
        }
    }

    //--------------------------------------------------
    // neighbor toward -z
    //--------------------------------------------------

    if (k > 0)
    {
        int id =
            i
            +
            (k-1)*domain.Nr;

        int obj =
            domain.object[id];

        if (
            obj == OBJ_PLASMA ||
            obj == OBJ_EXTRACTION
        )
        {
            n.z -= 1.0;
        }
    }

    //--------------------------------------------------
    // neighbor toward +z
    //--------------------------------------------------

    if (k < domain.Nz_phys-1)
    {
        int id =
            i
            +
            (k+1)*domain.Nr;

        int obj =
            domain.object[id];

        if (
            obj == OBJ_PLASMA ||
            obj == OBJ_EXTRACTION
        )
        {
            n.z += 1.0;
        }
    }

    //--------------------------------------------------
    // normalize
    //--------------------------------------------------

    double mag =
        sqrt(
            n.x*n.x +
            n.y*n.y +
            n.z*n.z
        );

    if (mag > 1e-12)
    {
        n.x /= mag;
        n.y /= mag;
        n.z /= mag;
    }
    else
    {
        //--------------------------------------------------
        // fallback
        //--------------------------------------------------

        n.x = 0.0;
        n.y = 0.0;
        n.z = 1.0;
    }

    return n;
}
//--------------------------------------------------
// specular reflection
//--------------------------------------------------

static vector3d specularReflect(
    const vector3d& vel,
    const vector3d& normal
)
{
    double vn =
          vel.x * normal.x
        + vel.y * normal.y
        + vel.z * normal.z;

    vector3d vnew;

    vnew.x =
        vel.x - 2.0 * vn * normal.x;

    vnew.y =
        vel.y - 2.0 * vn * normal.y;

    vnew.z =
        vel.z - 2.0 * vn * normal.z;

    return vnew;
}
//--------------------------------------------------
// diffuse wall reflection
//--------------------------------------------------

static vector3d diffuseReflect(
    const vector3d& vel,
    const vector3d& normal,
    double Twall,
    double mass
)
{
    vector3d vrand;

    double sigma =
        sqrt(Const::kB * Twall / mass);

    //--------------------------------------------------
    // random thermal velocity
    //--------------------------------------------------

    vrand.x =
        sigma * gaussian();

    vrand.y =
        sigma * gaussian();

    vrand.z =
        sigma * gaussian();

    //--------------------------------------------------
    // ensure outward direction
    //--------------------------------------------------

    double vn =
          vrand.x * normal.x
        + vrand.y * normal.y
        + vrand.z * normal.z;

    if (vn < 0.0)
    {
        vrand =
            vrand
            -
            normal * (2.0*vn);
    }

    return vrand;
}

void KineticSpecies::checkBounds()
{
    nLost = 0;
    nFC   = 0;

    double lc[3];

    for (Particle &part : particles)
    {
        if (part.STATUS != PARTICLE_ACTIVE)
            continue;

        //--------------------------------------------------
        // logical coordinates
        //--------------------------------------------------

        domain.Posv3d2lc(lc, part.pos);
		
        int i = static_cast<int>(floor(lc[0]));
		int k = static_cast<int>(floor(lc[2]));

        //--------------------------------------------------
        // outside global mesh
        //--------------------------------------------------
			
		
        if (
            i < 0 ||
            i >= domain.Nr ||
            k < 0 ||
            k >= domain.Nz_phys
        )
        {
            part.STATUS = PARTICLE_LOST;
            nLost++;
            continue;
        }
		 //--------------------------------------------------
        // Faraday Cup
        //--------------------------------------------------

        if (
              part.pos.z >= inputPar.Domain.posFC
        )
        {
            part.STATUS = PARTICLE_COLLECTED;
            nFC++;
            continue;
        }
		
        //--------------------------------------------------
        // geometry lookup
        //--------------------------------------------------

        int obj =
            domain.object[
                i + k*domain.Nr
            ];

        //--------------------------------------------------
        // plasma / extraction region
        //--------------------------------------------------

        if (
            obj == OBJ_PLASMA ||
            obj == OBJ_EXTRACTION
        )
        {
            continue;
        }

        //--------------------------------------------------
        // solid wall / electrode
        //--------------------------------------------------

        if (
            obj == OBJ_WALL ||
            obj == OBJ_ELECTRODE 
        )
        {
            //----------------------------------------------
            // ELECTRONS
            //----------------------------------------------

            if (speciesType == ELECTRON)
            {
                part.STATUS = PARTICLE_LOST;
                nLost++;
            }

            //----------------------------------------------
            // IONS
            //----------------------------------------------

            else if (speciesType == ION)
            {
                //------------------------------------------
                // optional:
                // secondary electron emission
                //------------------------------------------

                /*
                if (rnd() < gammaSEE)
                {
                    createSecondaryElectron(part);
                }
                */

                part.STATUS = PARTICLE_LOST;
                nLost++;
            }

            //----------------------------------------------
            // NEUTRALS
            //----------------------------------------------

            else if (speciesType == NEUTRAL)
			{
			    //--------------------------------------------------
			    // compute wall normal
			    //--------------------------------------------------

			    vector3d normal =
        			computeWallNormal(
			            part.pos,
        			    domain
			        );

			    //--------------------------------------------------
    			// specular reflection
			    //--------------------------------------------------

   				 part.vel =
      			  specularReflect(
            			part.vel,
			            normal
        			);

			    //--------------------------------------------------
    			// move particle outside wall
			    //--------------------------------------------------

			    part.pos =
        			part.pos
        			+
			        normal * (1.1 * domain.dh[0]);
			}

            continue;
        }

        //--------------------------------------------------
        // dielectric
        //--------------------------------------------------

        if (obj == OBJ_DIELECTRIC)
        {
            //----------------------------------------------
            // absorb charged particles
            //----------------------------------------------

            if (
                speciesType == ELECTRON ||
                speciesType == ION
            )
            {
                part.STATUS = PARTICLE_LOST;
                nLost++;
            }

            //----------------------------------------------
            // reflect neutrals
            //----------------------------------------------

            else if (speciesType == NEUTRAL)
            {
                //--------------------------------------------------
			    // compute wall normal
			    //--------------------------------------------------

			    vector3d normal =
        			computeWallNormal(
			            part.pos,
        			    domain
			        );

			    //--------------------------------------------------
    			// specular reflection
			    //--------------------------------------------------

   				 part.vel =
      			  specularReflect(
            			part.vel,
			            normal
        			);

			    //--------------------------------------------------
    			// move particle outside wall
			    //--------------------------------------------------

			    part.pos =
        			part.pos
        			+
			        normal * (1.1 * domain.dh[0]);
            }
        }
    }

    //------------------------------------------------------
    // erase removed particles
    //------------------------------------------------------

    for (
        auto it = particles.begin();
        it != particles.end();
    )
    {
        if (
            it->STATUS == PARTICLE_LOST ||
            it->STATUS ==  PARTICLE_COLLECTED
        )
        {
            it = particles.erase(it);
        }
        else
        {
            ++it;
        }
    }
}
/*compute number density*/
void KineticSpecies::computeGasProperties()
{
    double lc[3];

    den.clear();

    //--------------------------------------------------
    // scatter particle weights
    //--------------------------------------------------

    for (Particle &part : particles)
    {
        if (part.STATUS != PARTICLE_ACTIVE)
            continue;

        domain.Posv3d2lc(lc, part.pos);

        den.scatter(
            lc,
            double(spwt),
            0.0,
            domain.dh[0]
        );
    }

    //--------------------------------------------------
    // accumulate density at boundary
    //--------------------------------------------------
/*	den.accumulateBoundaries(
	domain.mpi_rank,
      domain.mpi_size,
        1
    	);
 */   	
    	
    //--------------------------------------------------
    // convert to number density
    //--------------------------------------------------

	
    den.divideByNodeVolume(
        domain.mpi_rank,
        domain.mpi_size,
        domain.dh[0],
        domain.dh[2]
    );

    //--------------------------------------------------
    // sync again
    //--------------------------------------------------

    den.updateBoundaries(
        domain.mpi_rank,
        domain.mpi_size,
        1
    );
    	
}


void KineticSpecies::transferParticles()
{
    //--------------------------------------------------
    // MPI particle datatype
    //--------------------------------------------------

    MPI_Datatype vec_type;
    MPI_Datatype particle_type;

    // vector3d = 3 doubles
    MPI_Type_contiguous(
        3,
        MPI_DOUBLE,
        &vec_type
    );

    MPI_Type_commit(&vec_type);

    int lengths[3] = {1,1,1};

    MPI_Aint displacements[3];

    displacements[0] =
        offsetof(Particle,pos);

    displacements[1] =
        offsetof(Particle,vel);

    displacements[2] =
        offsetof(Particle,STATUS);

    MPI_Datatype types[3] =
    {
        vec_type,
        vec_type,
        MPI_INT
    };

    MPI_Type_create_struct(
        3,
        lengths,
        displacements,
        types,
        &particle_type
    );

    MPI_Type_commit(&particle_type);

    //--------------------------------------------------
    // Neighbor ranks
    //--------------------------------------------------

    int left  = MPI_PROC_NULL;
    int right = MPI_PROC_NULL;

    if (domain.mpi_rank > 0)
        left = domain.mpi_rank - 1;

    if (domain.mpi_rank < domain.mpi_size-1)
        right = domain.mpi_rank + 1;

    //--------------------------------------------------
    // Reset receive counters
    //--------------------------------------------------

    recvL = 0;
    recvR = 0;

    MPI_Status status;

    //--------------------------------------------------
    // Exchange particle counts
    //--------------------------------------------------

    MPI_Sendrecv(
        &sendL,
        1,
        MPI_INT,
        left,
        0,

        &recvR,
        1,
        MPI_INT,
        right,
        0,

        MPI_COMM_WORLD,
        &status
    );

    MPI_Sendrecv(
        &sendR,
        1,
        MPI_INT,
        right,
        1,

        &recvL,
        1,
        MPI_INT,
        left,
        1,

        MPI_COMM_WORLD,
        &status
    );

    //--------------------------------------------------
    // Allocate buffers
    //--------------------------------------------------

    send_bufL =
        new Particle[std::max(sendL,1)];

    send_bufR =
        new Particle[std::max(sendR,1)];

    recv_bufL =
        new Particle[std::max(recvL,1)];

    recv_bufR =
        new Particle[std::max(recvR,1)];

    //--------------------------------------------------
    // Fill send buffers
    //--------------------------------------------------

    int iL = 0;
    int iR = 0;

    for (
        auto part_it = particles.begin();
        part_it != particles.end();
    )
    {
        Particle &part = *part_it;

        vectorCyl pos_c =
            Cart2Cyl(part.pos);

        //--------------------------------------------------
        // Send LEFT
        //--------------------------------------------------

        if(part.STATUS == PARTICLE_SEND_LEFT)
        	{
            send_bufL[iL++] = part;
    		part_it = particles.erase(part_it);
        	}
        //--------------------------------------------------
        // Send RIGHT
        //--------------------------------------------------

        else if(part.STATUS == PARTICLE_SEND_RIGHT)
			{
    		send_bufR[iR++] = part;
    		part_it = particles.erase(part_it);
			}

        //--------------------------------------------------
        // Keep particle
        //--------------------------------------------------

        else
        {
            ++part_it;
        }
    }

    //--------------------------------------------------
    // Debug info
    //--------------------------------------------------

    /*
    std::cout
        << "Rank "
        << domain.mpi_rank
        << " sendL="
        << sendL
        << " sendR="
        << sendR
        << " recvL="
        << recvL
        << " recvR="
        << recvR
        << "\n";
    */

    //--------------------------------------------------
    // Send LEFT / receive RIGHT
    //--------------------------------------------------

    MPI_Sendrecv(
        (sendL > 0) ? send_bufL : nullptr,
        sendL,
        particle_type,
        left,
        10,

        (recvR > 0) ? recv_bufR : nullptr,
        recvR,
        particle_type,
        right,
        10,

        MPI_COMM_WORLD,
        &status
    );

    //--------------------------------------------------
    // Send RIGHT / receive LEFT
    //--------------------------------------------------

    MPI_Sendrecv(
        (sendR > 0) ? send_bufR : nullptr,
        sendR,
        particle_type,
        right,
        11,

        (recvL > 0) ? recv_bufL : nullptr,
        recvL,
        particle_type,
        left,
        11,

        MPI_COMM_WORLD,
        &status
    );

    //--------------------------------------------------
    // Add received particles
    //--------------------------------------------------
	const double eps =
    1e-10 * domain.dh[2];
    for (int i=0; i<recvL; i++)
	{
    	recv_bufL[i].pos.z += eps;

    	recv_bufL[i].STATUS = PARTICLE_ACTIVE;

    	particles.emplace_back(recv_bufL[i]);
	}

	for (int i=0; i<recvR; i++)
	{
	    recv_bufR[i].pos.z -= eps;

	    recv_bufR[i].STATUS = PARTICLE_ACTIVE;

	    particles.emplace_back(recv_bufR[i]);
	}

    //--------------------------------------------------
    // Cleanup
    //--------------------------------------------------

    delete[] send_bufL;
    delete[] send_bufR;

    delete[] recv_bufL;
    delete[] recv_bufR;

    MPI_Type_free(&particle_type);
    MPI_Type_free(&vec_type);
}
