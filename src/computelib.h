#ifndef _COMPUTE_LIB_
#define _COMPUTE_LIB_

#include <vector>
#include <mpi.h>
#include <stdio.h>

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#define SRC_LOC __FILE__ ":" TOSTRING(__LINE__)

#define CHECK_MPI_ERROR(fn) \
{ \
int errcode = (fn); \
	if(errcode != MPI_SUCCESS) \
	{ \
		fprintf(stderr, "MPI ERROR: %d; %s.\n", errcode, SRC_LOC); \
		MPI_Abort(MPI_COMM_WORLD, 911); \
	} \
}

/**
 * Helper functions for exponentiation to integer powers.
 * @{ */
#define P2(x) ((x)*(x))
#define P3(x) ((x)*(x)*(x))
#define P4(x) ((x)*(x)*(x)*(x))
/** @} */

/****** marcos for identifying models ******/
#define 		KNBGR			0
#define 		TSK				1

#define UP 0
#define DOWN 1
#define LEFT 2
#define RIGHT 3

#define STRSEG 0 ///< A straight segment.
#define BIF 1 ///< A bifurcation.
#define P 1 ///< Parent.
#define L 2 ///< Left branch.
#define R 3 ///< Right branch.

#define local 0
#define remote 1

#define FILENAME "PhyFiIP3_plotting.csv"

#define PLOTTING 0
#define EXPLICIT_ONLY 1

// Number of elements to plot for Lemon model
#define OUTPUT_PLOTTING_SIZE 14
#define RANK 0
#define EC_COL 3
#define EC_ROW 3
#define SMC_COL 1
#define SMC_ROW 33
extern FILE* var_file;

extern double* plotttingBuffer;
extern int bufferPos;


/**
 * Conductance / coupling coefficients.
 */
struct conductance
{
	double Vm_hm_smc,	///< Homocellular membrane potential coupling between SMCs.
			Vm_hm_ec,	///< Homocellular membrane potential coupling between ECs.
			Ca_hm_smc,	///< Homocellular Ca coupling between SMCs.
			Ca_hm_ec,	///< Homocellular Ca coupling between ECs.
			IP3_hm_smc,	///< Homocellular IP3 coupling between SMCs.
			IP3_hm_ec,	///< Homocellular IP3 coupling between ECs.
			Vm_ht_smc,	///< Heterocellular membrane potential coupling between SMCs.
			Vm_ht_ec,	///< Heterocellular membrane potential coupling between ECs.
			Ca_ht_smc,	///< Heterocellular Ca coupling between SMCs.
			Ca_ht_ec,	///< Heterocellular Ca coupling between ECs.
			IP3_ht_smc,	///< Heterocellular IP3 coupling between SMCs.
			IP3_ht_ec;	///< Heterocellular IP3 coupling between ECs.
	double ec_diffusion[4];
	double smc_diffusion[4];
};

#define DOMAIN_NUM 0
#define DOMAIN_TYPE 1
#define AX_QUADS 2
#define CR_QUADS 3
#define PARENT_DOMAIN_NUM 4
#define LEFT_DOMAIN_NUM 5
#define RIGHT_DOMAIN_NUM 6
#define AX_ECS 7
#define CR_SMCS 8
#define NUM_CONFIG_ELEMENTS 9

// TODO: Initialise and use constants correctly within grid_params.
// const int num_ghost_cells = 2;
typedef struct
{
	///General information on cell geometry and the geometric primitive constructed.
	int num_smc_fundblk_circumferentially;
	int num_ec_fundblk_circumferentially;
	int num_smc_fundblk_axially;
	int num_ec_fundblk_axially;
	int num_ghost_cells;
	int num_fluxes_smc; // Number of SMC ioinic currents to be evaluated for eval of LHS of the d/dt terms of the ODEs.
	int num_fluxes_ec; // Number of EC ioinic currents to be evaluated for eval of LHS of the d/dt terms of the ODEs.
	int num_coupling_species_smc; // Number of SMC coupling species homogenic/heterogenic.
	int num_coupling_species_ec; // Number of SMC coupling species homogenic/heterogenic.
	int neq_smc; // Number of SMC ODEs for a single cell.
	int neq_ec; // Number of EC ODEs for a single cell.

	///Number of elements added to the send buffer for sending relevant information on the content of the buffer to receiving task.
	int added_info_in_send_buf;

	int domain_index;
	int domain_type;
	int boundary_tag;

	int
	/// Global domain information storage with these elements:
	// [DOMAIN_NUM, DOMIN_TYPE, AX_QUADS, CR_QUADS, PARENT_DOMAIN_NUM, LEFT_DOMAIN_NUM, RIGHT_DOMAIN_NUM, AX_ECS, CR_SMCS].
	num_domains, **domain_params,
	///total grid points axially
	m,
	///total grid points circumferentially
	n,
	/// My coordinates
	coords[2],
	///Coordinates for neighbour tasks.
	nbrs[2][4],
	///Node payload information (number of cells laid out on a node).
	num_ec_axially, num_ec_circumferentially, num_smc_axially, num_smc_circumferentially, neq_ec_axially, neq_smc_axially,
	///Total number of state variables in the computational domain
	NEQ,
	///This is global and local MPI information
	num_ranks, universal_rank, /// numtasks = total CPUs in MPI_COMM_WORLD,
	num_ranks_branch, rank_branch, /// tasks = total CPUs in my-subdomain's comm (branch)

	//Each processor on the edges of each branch contains brach_tag can have one of four values P=parent = 1, L=Left branch = 2, R=Right branch = 3.
	//If branch_tag=0, this implies that the rank is located interior or doesn't  contain a remote neighbour on any other branch.
	branch_tag,
	/// Variables for remote MPI information (P=parent, L & R = Left & Right branch respectively).
	offset_P, offset_L, offset_R, flip_array[4],
	/// Number of elements being sent and received.
	num_elements_send_up, num_elements_send_down, num_elements_send_left, num_elements_send_right,
	num_elements_recv_up, num_elements_recv_down, num_elements_recv_left, num_elements_recv_right;

	///Information for spatial variation in agonist.
	double uniform_jplc, stimulus_onset_time;	/// the time when spatially varying agonist kicks in

	// MPI_COMM_WORLD, and branch communicatior.
	MPI_Comm universe, cart_comm;

	int smc_model, ec_model; // These are placeholders for the selection of model to be simulated in each cell.
	int NO_path, cGMP_path;	// Specific for Tsoukias model to signal whether to activate NO and cGMP pathways for vasodilation.

	char solution_dir[1024], time_profiling_dir[1024], config_file[1024];
	bool random;
} grid_parms;

typedef struct
{
	double *vars;		///storage for the state variables corresponding to an SMC.
	double NO, NE, I_stim;		///specific to Tsoukias model
	double* fluxes;			    ///stores single cell fluxes
	double* homo_fluxes;			    ///stores homogeneous coupling fluxes
	double* hetero_fluxes;			    ///stores heterogeneous coupling fluxes
	//conductance cpl_cef;
} SMC_cell;

typedef struct
{
	double* vars;		///storage for the state variables corresponding to an SMC.
	double* fluxes;			    ///stores single cell fluxes
	double* homo_fluxes;			    ///stores homogeneous coupling fluxes
	double* hetero_fluxes;			    ///stores heterogeneous coupling fluxes
	double JPLC;			    ///local agonist concentration  on my GPCR receptor (an ith EC)
	//conductance cpl_cef;
} EC_cell;

class Cell_type {
public:
	Cell_type() {
		_nc = -1;
		_na = -1;
		_nv = -1;
		_nf = -1;
		_ncs = -1;
		this->data.resize(0);
	}
	virtual ~Cell_type() {
		this->destroy();
	}
	void build(int num_circumferentially, 
		       int num_axially, 
		       int num_vars, 
		       int num_fluxes, 
		       int num_coupling_species) {
		_nc = num_circumferentially;
		_na = num_axially;
		_nv = num_vars;
		_nf = num_fluxes;
		_ncs = num_coupling_species;
		this->data.resize(_nc * _na * (_nv + _nf + 2*_ncs));
	}
	void destroy() {
		this->data.resize(0);
	}
	inline double& var(int i, int j, int k) {
		return this->data[this->getFlatIndex(i, j) + k];
	}
	inline double& flux(int i, int j, int k) {
		return this->data[this->getFlatIndex(i, j) + k + _nv];
	}
	inline double& homo_flux(int i, int j, int k) {
		return this->data[this->getFlatIndex(i, j) + k + _nv + _nf];
	}
	inline double& hetero_flux(int i, int j, int k) {
		return this->data[this->getFlatIndex(i, j) + k + _nv + _nf + _ncs];
	}

protected:
	std::vector<double> data;
	int _nc;
	int _na;
	int _nv;
	int _nf;
	int _ncs;
	inline int getFlatIndex(int i, int j) const {
		return (_nv + _nf + 2*_ncs) * (_na*i + j);
	}
};

class EC_type : public Cell_type {
public:
	std::vector<double> jplc;
	void build(int num_circumferentially, 
		       int num_axially, 
		       int num_vars, 
		       int num_fluxes, 
		       int num_coupling_species) {
		Cell_type::build(num_circumferentially, 
		                 num_axially, 
		                 num_vars, 
		                 num_fluxes, 
		                 num_coupling_species);
		this->jplc.resize(_nc * _na);
	}
	void destroy() {
		Cell_type::destroy();
		this->jplc.resize(0);
	}
	inline double& JPLC(int i, int j) {
		return this->jplc[this->getFlatIndex(i, j)];
	}
};

class SMC_type : public Cell_type {
public:
	std::vector<double> no, ne;
	void build(int num_circumferentially, 
		       int num_axially, 
		       int num_vars, 
		       int num_fluxes, 
		       int num_coupling_species) {
		Cell_type::build(num_circumferentially, 
		                 num_axially, 
		                 num_vars, 
		                 num_fluxes, 
		                 num_coupling_species);
		this->no.resize(_nc * _na);
		this->ne.resize(_nc * _na);
	}
	void destroy() {
		Cell_type::destroy();
		this->no.resize(0);		
		this->ne.resize(0);		
	}
	inline double& NO(int i, int j) {
		return this->no[this->getFlatIndex(i, j)];
	}
	inline double& NE(int i, int j) {
		return this->ne[this->getFlatIndex(i, j)];
	}
};


typedef struct
{
	double aggregate_compute, aggregate_comm;
	double aggregate_ec_gather, aggregate_smc_gather;
	double aggregate_ec_write, aggregate_smc_write;
} time_stamps;

//Topology related functions
void set_task_parameters(grid_parms *);
void make_bifucation_cart_grids(grid_parms *);
void make_straight_cart_grid(grid_parms *);
void read_init_ATP(grid_parms *grid, EC_type& ECs);
void set_coupling_parms(int CASE, conductance* cpl_cef);

void determine_source_destination(grid_parms, int*, int*);
void communication_update_recv_size(grid_parms *);
void communication_update_sendbuf(grid_parms, double**, SMC_type& , EC_type& );
void communication_update_recvbuf(grid_parms, double**, SMC_type& , EC_type& );
void communication_async_send_recv(grid_parms, double**, double**, SMC_type& , EC_type& );

//Cell dynamics evaluation handlers. These contain the ODEs for representative models from different sources.
void compute(const grid_parms&, SMC_type& , EC_type& , const conductance& cpl_cef, double, double*, double*);
void compute_implicit(const grid_parms&, SMC_type& , EC_type& , const conductance& cpl_cef, double, double*, double*);
void compute_explicit(const grid_parms&, SMC_type& , EC_type& , const conductance& cpl_cef, double, double*, double*);
void coupling(double, double*, const grid_parms&, SMC_type& , EC_type& , const conductance&);
void coupling_implicit(double, double*, const grid_parms&, SMC_type& , EC_type& , const conductance&);
void coupling_explicit(double, double*, const grid_parms&, SMC_type& , EC_type& , const conductance&);

// Solver wrapper functions.
#ifdef RK_SUITE
void rksuite_solver_CT(double, double, double, double*, double*, int, double, double*, int, char*); //, IO_domain_info*);
#endif
#ifdef ARK_ODE
void arkode_solver(double, double, double, double*, int, double, double, int, char*);
#endif
#ifdef BOOST_ODEINT
void odeint_solver(double, double, double, double*, int, double, double, int, char*);
#endif

int map_solver_output_to_cells(const grid_parms&, double*, SMC_type& , EC_type& );

///These are debugging functions, not used in production runs.
void dump_rank_info(conductance, grid_parms);
void print_domains(FILE*, grid_parms, SMC_type& , EC_type& );
void print_send_buffer(FILE*, grid_parms, double**);
void print_recv_buffer(FILE*, grid_parms, double**);
void print_compare(double, double*, grid_parms, SMC_type& , EC_type& );

void initialize_t_stamp(time_stamps*);
void dump_time_profiling(grid_parms grid, time_stamps* t_stamp);
void dump_time_field(char* file_prefix, grid_parms grid, double field);

/**
 * \brief Catch failed memory allocation.
 */
void* checked_malloc(size_t bytes, const char* errmsg);

#endif
