#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "koenigsberger_macros.h"

using namespace std;

#define local	0
#define remote 1

#define UP1    0
#define UP2		1
#define DOWN1  2
#define DOWN2	3
#define LEFT1	4
#define LEFT2	5
#define RIGHT1 6
#define RIGHT2 7

#define UP		0
#define DOWN 	1
#define LEFT	2
#define RIGHT	3

#define P 		1		//parent
#define L 		2		//Left branch
#define R		3		//Right brach
// helper functions for exponentiation to integer powers
#define P2(x) ((x)*(x))
#define P4(x) ((x)*(x)*(x)*(x))
struct conductance{
	double
		Vm_hm_smc,			///homocellular membrane potential coupling between SMCs
		Vm_hm_ec,			///homocellular membrane potential coupling between ECs
		Ca_hm_smc,			///homocellular Ca coupling between SMCs
		Ca_hm_ec,			///homocellular Ca coupling between ECs
		IP3_hm_smc,			///homocellular IP3 coupling between SMCs
		IP3_hm_ec,			///homocellular IP3 coupling between ECs
		Vm_ht_smc,			///heterocellular membrane potential coupling between SMCs
		Vm_ht_ec,			///heterocellular membrane potential coupling between ECs
		Ca_ht_smc,			///heterocellular Ca coupling between SMCs
		Ca_ht_ec,			///heterocellular Ca coupling between ECs
		IP3_ht_smc,			///heterocellular IP3 coupling between SMCs
		IP3_ht_ec;			///heterocellular IP3 coupling between ECs
	};

typedef struct {
	///General infomation on cell geometry and the geometric primitive constructed.
		double hx_smc, hx_ec,hy_smc, hy_ec,
		requested_length, requested_diameter,corrected_length, corrected_diameter, new_circ;
	int
	///this is a node local information

	///Topology information (fundamental unit or block of cells)
	num_smc_fundblk_circumferentially, num_ec_fundblk_circumferentially,
			num_smc_fundblk_axially, num_ec_fundblk_axially,

			///Total number of ghost cells to be added in the computational in each dimension (circumferentail and axial)
			num_ghost_cells,
			///total grid points axially
			m,
			///total grid points circumferentially
			n,
			/// My coordinates
			coords[2],
			///Coordinates for neighbour tasks
			nbrs[2][4],
			///Node payload information(number of cells laid out on a node)
			num_ec_axially, num_ec_circumferentially, num_smc_axially,
			num_smc_circumferentially, neq_ec_axially, neq_smc_axially,
			///Model related parameters
			///number of equations modelling an EC or SMC
			neq_smc,
			neq_ec,
			///Total number of state variables in the computational domain
			NEQ,
			num_fluxes_smc, num_fluxes_ec, num_coupling_species_smc,
			num_coupling_species_ec,
			///Number of elements added to the Send buffer for sending relevant information on the content of the buffer to receiving task
			added_info_in_send_buf,
			///this is global and local MPI information
			numtasks, universal_rank, sub_universe_numtasks, sub_universe_rank, rank, color,key,
			//Each processor on the edges of each branch contains brach_tag can have one of four values P=parent = 1, L=Left branch = 2, R=Right brach = 3.
			//If branch_tag=0, this implies that the rank is located interior or doesn't  contain a remote neighbour on any other branch.
			branch_tag,
			///variables for remote MPI information (P=parent, L & R = Left & Right branch respectively)
			scheme,offset_P,offset_L,offset_R,flip_array[4],
			///number of elements being sent and received
			num_elements_send_up, num_elements_send_down,
			num_elements_send_left, num_elements_send_right,
			num_elements_recv_up, num_elements_recv_down,
			num_elements_recv_left, num_elements_recv_right;
	///Information for spatial variation in agonist
	double
	min_jplc, max_jplc, gradient, uniform_jplc;
	//Allow two types of communicators to exist, first resulting from comm_split operation on MPI_COMM_WORLD
	//and the other a Cartisian communicator arising from Cart_create operation
	MPI_Comm universe,sub_universe,split_comm,cart_comm;
			}grid_parms;

///Structure to store coupling data received from the neighbouring task.
typedef struct{
            double c,v,I;
}nbrs_data;


typedef struct{
	double	 	*p;				///storage for the state variables corresponding to an SMC.
	double 		*d;				///storage for the first derivative of state variables corresponding to an SMC.
	int		node_row,node_col;	///stores coordinates of the node on which I am located.
	int		my_row,my_col;		///stores my location on the node.
	double*		 A;			    ///stores single cell fluxes
	double*		 B;			    ///stores homogeneous coupling fluxes
	double*		 C;			    ///stores heterogeneous coupling fluxes

	conductance    cpl_cef;
	}celltype1;

typedef struct{
	double	 	*q;				///storage for the state variables corresponding to an SMC.
	double 		*d;				///storage for the first derivative of state variables corresponding to an EC.
	int		node_row,node_col;	///stores coordinates of the node on which I am located.
	int		my_row,my_col;		///stores my location on the node.
	double*		 A;			    ///stores single cell fluxes
	double*		 B;			    ///stores homogeneous coupling fluxes
	double*		 C;			    ///stores heterogeneous coupling fluxes

	double		 JPLC;			///local agonsit concentration  on my GPCR receptor (an ith EC)
	conductance    cpl_cef;
	}celltype2;
#ifdef PARALLEL_IO
typedef struct{
MPI_File logptr, Time, ci, si, vi, wi, Ii, cpCi, cpVi, cpIi, cj,
sj, vj, Ij, cpCj, cpVj, cpIj;
}checkpoint_handle;

#else
typedef struct {
	FILE *logptr, *Time, *ci, *si, *vi, *wi, *Ii, *cpCi, *cpVi, *cpIi, *cj,
	*sj, *vj, *Ij, *cpCj, *cpVj, *cpIj;
}checkpoint_handle;

#endif

typedef struct {
	int key_val;		//am I a bifurcation or a straight segment?
	int m,n;			//row and columns in my MPI_sub_world
	double d,l;			//diameter and length scales
}node;

typedef struct {
	node	internal_info;
	node 	left_child,right_child,parent;
}my_tree;


void check_flag(int, FILE*, const char*);
void* checked_malloc(size_t bytes, FILE*, const char* errmsg);

int couplingParms(int CASE, conductance* cpl_cef);
void Initialize_koeingsberger_smc(grid_parms,double*,celltype1**);
void Initialize_koeingsberger_ec(grid_parms,double*,celltype2**);
void map_GhostCells_to_cells(celltype1**, celltype2**, grid_parms);
void map_solver_to_cells(grid_parms grid,double y[], celltype1** smc, celltype2** ec);

grid_parms communicate_num_recv_elements_to_nbrs(FILE* logptr, grid_parms grid);
void communication_update_sendbuf(FILE* logptr,grid_parms grid,double** sendbuf, celltype1** smc, celltype2** ec);
void communication_update_recvbuf(FILE* logptr,grid_parms grid,double** recvbuf, celltype1** smc, celltype2** ec);
void single_cell(double t,double y[], grid_parms grid, celltype1** smc, celltype2** ec);
void coupling(double t,double y[], grid_parms grid, celltype1** smc, celltype2** ec,conductance cpl_cef);
void determin_source_destination(grid_parms grid, int source[], int dest[]);

void communication_async_send_recv(FILE* logptr, grid_parms grid, double** sendbuf, double** recvbuf,celltype1** smc, celltype2** ec);

///Checkpoint functions.
checkpoint_handle* initialise_checkpoint(int);
void dump_smc(grid_parms, celltype1**, checkpoint_handle*, int);
void dump_ec(grid_parms, celltype2**, checkpoint_handle*, int);
void dump_JPLC(grid_parms, celltype2**, checkpoint_handle*, const char*);
void checkpoint(checkpoint_handle*, grid_parms, double, celltype1**, celltype2**, int);
void final_checkpoint(checkpoint_handle*, grid_parms, double, double);
void dump_rank_info(checkpoint_handle*, conductance,grid_parms);
void dump_smc_with_ghost_cells(grid_parms, celltype1**, checkpoint_handle*, int);
void dump_ec_with_ghost_cells(grid_parms, celltype2**, checkpoint_handle*, int);
void checkpoint_with_ghost_cells(checkpoint_handle*, grid_parms, double, celltype1**, celltype2**, int);

void computeDerivatives(double, double*, double*);
void rksuite_solver_CT(double tnow, double tfinal, double interval, double *y, double* yp,
		int total, double TOL, double* thres, int file_write_per_unit_time,
		checkpoint_handle *check);

void rksuite_solver_UT(double tnow, double tfinal, double interval, double *y, double* yp,
		int total, double TOL, double* thres, int file_write_per_unit_time,
		checkpoint_handle *check);
///These are debugging functions, not used in production runs.
void print_domains(FILE* logptr,grid_parms grid,celltype1** smc, celltype2** ec);
void print_send_buffer(FILE* logptr,grid_parms grid, double** sendbuf);
void print_recv_buffer(FILE* logptr,grid_parms grid, double** recvbuf);
void print_compare(FILE* logptr, double t, double y[],grid_parms grid, celltype1** smc, celltype2** ec);


void communication_update_recvbuf_modified(FILE* logptr,grid_parms grid,double** recvbuf, celltype1** smc, celltype2** ec);




















