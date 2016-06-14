#include <malloc.h>
#include <stdlib.h> // malloc

#include "computelib.h"
#include "koenigsberger_model.h"
#include "tsoukias_model.h"

using namespace std;
time_stamps t_stamp;

/**
 * Wrapper around malloc to catch failed memory allocation. If allocation fails MPI_Abort is called.
 *
 * \param bytes Size of requested memory.
 * \param errmsg Message produced in the event of failed memory allocation.
 */
void* checked_malloc(size_t bytes, const char* errmsg)
{
	void *pval = malloc(bytes);
	if (pval == NULL)
	{
		fprintf(stderr, "MEMORY ALLOCATION ERROR: %s\n", errmsg);
		MPI_Abort(MPI_COMM_WORLD, 911);
	}
	return pval;
}

void set_coupling_parms(int CASE, conductance* cpl_cef)
{
	// These values are honestly just trial and error guesses at what makes nice waves...
	// Really need to have a think about what they should be.
	cpl_cef->ec_diffusion[0] = 1;
	cpl_cef->ec_diffusion[1] = 1;
	cpl_cef->ec_diffusion[2] = 1;
	cpl_cef->ec_diffusion[3] = 1;

	cpl_cef->smc_diffusion[0] = 1;
	cpl_cef->smc_diffusion[1] = 1;
	cpl_cef->smc_diffusion[2] = 1;
	cpl_cef->smc_diffusion[3] = 1;


	if(CASE == 1)
	{
		cpl_cef->Vm_hm_smc = 1000.00;
		cpl_cef->Vm_hm_ec = 1000.00;

		cpl_cef->Ca_hm_smc = 0.1;
		cpl_cef->Ca_hm_ec = 0.1;

		cpl_cef->IP3_hm_smc = 0.05;
		cpl_cef->IP3_hm_ec = 0.00;

		cpl_cef->Vm_ht_smc = 50.0;
		cpl_cef->Vm_ht_ec = 50.0;

		cpl_cef->Ca_ht_smc = 0.00;
		cpl_cef->Ca_ht_ec = 0.00;

		cpl_cef->IP3_ht_smc = 0.05;
		cpl_cef->IP3_ht_ec = 0.05;
	}
	else if(CASE == 2)
	{
		cpl_cef->Vm_hm_smc = 1000.00;
		cpl_cef->Vm_hm_ec = 1000.00;

		cpl_cef->Ca_hm_smc = 0.01;
		cpl_cef->Ca_hm_ec = 0.01;

		cpl_cef->IP3_hm_smc = 0.05;
		cpl_cef->IP3_hm_ec = 0.00;

		cpl_cef->Vm_ht_smc = 50.0;
		cpl_cef->Vm_ht_ec = 50.0;

		cpl_cef->Ca_ht_smc = 0.05;
		cpl_cef->Ca_ht_ec = 0.05;

		cpl_cef->IP3_ht_smc = 0.05;
		cpl_cef->IP3_ht_ec = 0.05;
	}
	else if(CASE == 3)
	{
		cpl_cef->Vm_hm_smc = 1000.00;
		cpl_cef->Vm_hm_ec = 1000.00;

		cpl_cef->Ca_hm_smc = 0.05;
		cpl_cef->Ca_hm_ec = 0.05;

		cpl_cef->IP3_hm_smc = 0.05;
		cpl_cef->IP3_hm_ec = 0.05;

		cpl_cef->Vm_ht_smc = 50.0;
		cpl_cef->Vm_ht_ec = 50.0;

		cpl_cef->Ca_ht_smc = 0.05;
		cpl_cef->Ca_ht_ec = 0.05;

		cpl_cef->IP3_ht_smc = 0.05;
		cpl_cef->IP3_ht_ec = 0.05;
	}
	else if(CASE == 4)
	{
		cpl_cef->Vm_hm_smc = 1000.00;
		cpl_cef->Vm_hm_ec = 0.00;

		cpl_cef->Ca_hm_smc = 0.05;
		cpl_cef->Ca_hm_ec = 0.00;

		cpl_cef->IP3_hm_smc = 0.05;
		cpl_cef->IP3_hm_ec = 0.05;

		cpl_cef->Vm_ht_smc = 0.0;
		cpl_cef->Vm_ht_ec = 0.0;

		cpl_cef->Ca_ht_smc = 0.0;
		cpl_cef->Ca_ht_ec = 0.0;

		cpl_cef->IP3_ht_smc = 0.05;
		cpl_cef->IP3_ht_ec = 0.05;
	}
	else if(CASE == 5)
	{
		// Simulating for experiments suggested by Dr. James Kozloski (IBM Watson Centre).
		// The homocellular Ca coupling between SMCs is changed to investigate the effects of
		// strength on coupling on the propagation speed of the spatial waves.
		cpl_cef->Vm_hm_smc = 0.00;
		cpl_cef->Vm_hm_ec = 0.00;

		cpl_cef->Ca_hm_smc = 0.05;
		cpl_cef->Ca_hm_ec = 0.05;

		cpl_cef->IP3_hm_smc = 0.05;
		cpl_cef->IP3_hm_ec = 0.00;

		cpl_cef->Vm_ht_smc = 50.0;
		cpl_cef->Vm_ht_ec = 50.0;

		cpl_cef->Ca_ht_smc = 0.00;
		cpl_cef->Ca_ht_ec = 0.00;

		cpl_cef->IP3_ht_smc = 0.05;
		cpl_cef->IP3_ht_ec = 0.05;
	}
	else if(CASE == 6)
	{
		// Simulating for experiments suggested by Dr. James Kozloski (IBM Watson Centre).
		cpl_cef->Vm_hm_smc = 4000.00;
		cpl_cef->Vm_hm_ec = 1000.00;

		cpl_cef->Ca_hm_smc = 0.05;
		cpl_cef->Ca_hm_ec = 0.05;

		cpl_cef->IP3_hm_smc = 0.05;
		cpl_cef->IP3_hm_ec = 0.00;

		cpl_cef->Vm_ht_smc = 50.0;
		cpl_cef->Vm_ht_ec = 50.0;

		cpl_cef->Ca_ht_smc = 0.00;
		cpl_cef->Ca_ht_ec = 0.00;

		cpl_cef->IP3_ht_smc = 0.05;
		cpl_cef->IP3_ht_ec = 0.05;
	}
	else if(CASE == 7)
	{
		// Simulating for experiments suggested by Dr. James Kazloski (IBM Watson Centre).
		cpl_cef->Vm_hm_smc = 6000.00;
		cpl_cef->Vm_hm_ec = 1000.00;

		cpl_cef->Ca_hm_smc = 0.05;
		cpl_cef->Ca_hm_ec = 0.05;

		cpl_cef->IP3_hm_smc = 0.05;
		cpl_cef->IP3_hm_ec = 0.00;

		cpl_cef->Vm_ht_smc = 50.0;
		cpl_cef->Vm_ht_ec = 50.0;

		cpl_cef->Ca_ht_smc = 0.00;
		cpl_cef->Ca_ht_ec = 0.00;

		cpl_cef->IP3_ht_smc = 0.05;
		cpl_cef->IP3_ht_ec = 0.05;
	}
	else if(CASE == 8)
	{
		// What is this one?
		cpl_cef->Vm_hm_smc = 0.00;
		cpl_cef->Vm_hm_ec = 0.00;

		cpl_cef->Ca_hm_smc = 0.0;
		cpl_cef->Ca_hm_ec = 0.0;

		cpl_cef->IP3_hm_smc = 0.0;
		cpl_cef->IP3_hm_ec = 0.0;

		cpl_cef->Vm_ht_smc = 50.0;
		cpl_cef->Vm_ht_ec = 50.0;

		cpl_cef->Ca_ht_smc = 0.00;
		cpl_cef->Ca_ht_ec = 0.00;

		cpl_cef->IP3_ht_smc = 0.05;
		cpl_cef->IP3_ht_ec = 0.05;
	}
}

// TODO: Move the Tsoukias code to the appropriate location.
// Mapping from state variable vector to cells.
int map_solver_output_to_cells(const grid_parms& grid, 
                               double* y, 
                               SMC_cell* __restrict__ smc, 
                               EC_cell* __restrict__ ec)
{
	int err = 0;
	switch (grid.smc_model)
	{
	case (TSK): {
		int k = 0, offset;
		const int nc = grid.num_smc_circumferentially;
		const int na = grid.num_smc_axially;
		const int ng = grid.num_ghost_cells;
		for (int ij0 = 0; ij0 < nc * na; ij0++) {
			int i = ij0 / (na + ng) + 1;
			int j = ij0 % (na + ng) + 1;
			int ij = i*(na + ng) + j;
			k = (i - 1) * grid.neq_smc_axially;
			smc[ij].vars[smc_Vm] = y[k + ((j - 1) * grid.neq_smc) + smc_Vm];
			smc[ij].vars[smc_d_L] = y[k + ((j - 1) * grid.neq_smc) + smc_d_L];
			smc[ij].vars[smc_f_L] = y[k + ((j - 1) * grid.neq_smc) + smc_f_L];
			smc[ij].vars[smc_p_f] = y[k + ((j - 1) * grid.neq_smc) + smc_p_f];
			smc[ij].vars[smc_p_s] = y[k + ((j - 1) * grid.neq_smc) + smc_p_s];
			smc[ij].vars[smc_q_1] = y[k + ((j - 1) * grid.neq_smc) + smc_q_1];
			smc[ij].vars[smc_q_2] = y[k + ((j - 1) * grid.neq_smc) + smc_q_2];
			smc[ij].vars[smc_p_K] = y[k + ((j - 1) * grid.neq_smc) + smc_p_K];
			smc[ij].vars[smc_Ca_u] = y[k + ((j - 1) * grid.neq_smc) + smc_Ca_u];
			smc[ij].vars[smc_Ca_r] = y[k + ((j - 1) * grid.neq_smc) + smc_Ca_r];
			smc[ij].vars[smc_R_10] = y[k + ((j - 1) * grid.neq_smc) + smc_R_10];
			smc[ij].vars[smc_R_11] = y[k + ((j - 1) * grid.neq_smc) + smc_R_11];
			smc[ij].vars[smc_R_01] = y[k + ((j - 1) * grid.neq_smc) + smc_R_01];
			smc[ij].vars[smc_h_IP3] = y[k + ((j - 1) * grid.neq_smc) + smc_h_IP3];
			smc[ij].vars[smc_R_S_G] = y[k + ((j - 1) * grid.neq_smc) + smc_R_S_G];
			smc[ij].vars[smc_R_S_P_G] = y[k + ((j - 1) * grid.neq_smc) + smc_R_S_P_G];
			smc[ij].vars[smc_G] = y[k + ((j - 1) * grid.neq_smc) + smc_G];
			smc[ij].vars[smc_IP3] = y[k + ((j - 1) * grid.neq_smc) + smc_IP3];
			smc[ij].vars[smc_PIP2] = y[k + ((j - 1) * grid.neq_smc) + smc_PIP2];
			smc[ij].vars[smc_V_cGMP] = y[k + ((j - 1) * grid.neq_smc) + smc_V_cGMP];
			smc[ij].vars[smc_cGMP_i] = y[k + ((j - 1) * grid.neq_smc) + smc_cGMP_i];
			smc[ij].vars[smc_Ca] = y[k + ((j - 1) * grid.neq_smc) + smc_Ca];
			smc[ij].vars[smc_Na_i] = y[k + ((j - 1) * grid.neq_smc) + smc_Na_i];
			smc[ij].vars[smc_K_i] = y[k + ((j - 1) * grid.neq_smc) + smc_K_i];
			smc[ij].vars[smc_Cl_i] = y[k + ((j - 1) * grid.neq_smc) + smc_Cl_i];
			smc[ij].vars[smc_DAG] = y[k + ((j - 1) * grid.neq_smc) + smc_DAG];
		}
		break;
	}
	case (KNBGR): {
		int k = 0, offset;
		const int nc = grid.num_smc_circumferentially;
		const int na = grid.num_smc_axially;
		const int ng = grid.num_ghost_cells;
		for (int ij0 = 0; ij0 < nc * na; ij0++) {
			int i = ij0 / (na + ng) + 1;
            int j = ij0 % (na + ng) + 1;
            int ij = i*(na + ng) + j;
			k = (i - 1) * grid.neq_smc_axially;
			smc[ij].vars[smc_Ca] = y[k + ((j - 1) * grid.neq_smc) + smc_Ca];
			smc[ij].vars[smc_SR] = y[k + ((j - 1) * grid.neq_smc) + smc_SR];
			smc[ij].vars[smc_Vm] = y[k + ((j - 1) * grid.neq_smc) + smc_Vm];
			smc[ij].vars[smc_w] = y[k + ((j - 1) * grid.neq_smc) + smc_w];
			smc[ij].vars[smc_IP3] =y[k + ((j - 1) * grid.neq_smc) + smc_IP3];
		}
		break;
	}
	default: {
		err = 1;
		break;
	}
	}
	switch (grid.ec_model) {
	case (TSK): {
		int k, offset = (grid.neq_smc * grid.num_smc_circumferentially
				* grid.num_smc_axially);

        const int nc = grid.num_ec_circumferentially;
        const int na = grid.num_ec_axially;
        const int ng = grid.num_ghost_cells;
		for (int ij0 = 0; ij0 < nc * na; ij0++) {
			int i = ij0 / (na + ng) + 1;
			int j = ij0 % (na + ng) + 1;
			int ij = i*(na + ng) + j;
			k = offset + ((i - 1) * grid.neq_ec_axially);
			ec[ij].vars[ec_Ca] = y[k + ((j - 1) * grid.neq_ec) + ec_Ca];
			ec[ij].vars[ec_SR] = y[k + ((j - 1) * grid.neq_ec) + ec_SR];
			ec[ij].vars[ec_Vm] = y[k + ((j - 1) * grid.neq_ec) + ec_Vm];
			ec[ij].vars[ec_IP3] = y[k + ((j - 1) * grid.neq_ec) + ec_IP3];
		}
		break;
	}
	case (KNBGR): {
		int k, offset = (grid.neq_smc * grid.num_smc_circumferentially
				* grid.num_smc_axially);

		const int nc = grid.num_ec_circumferentially;
		const int na = grid.num_ec_axially;
		const int ng = grid.num_ghost_cells;
		for (int ij0 = 0; ij0 < nc * na; ij0++) {
			int i = ij0 / (na + ng) + 1;
			int j = ij0 % (na + ng) + 1;
			int ij = i*(na + ng) + j;
			k = offset + ((i - 1) * grid.neq_ec_axially);
			ec[ij].vars[ec_Ca] = y[k + ((j - 1) * grid.neq_ec) + ec_Ca];
			ec[ij].vars[ec_SR] = y[k + ((j - 1) * grid.neq_ec) + ec_SR];
			ec[ij].vars[ec_Vm] = y[k + ((j - 1) * grid.neq_ec) + ec_Vm];
			ec[ij].vars[ec_IP3] = y[k + ((j - 1) * grid.neq_ec) + ec_IP3];
			ec[ij].vars[ec_Gprot] = y[k + ((j - 1) * grid.neq_ec) + ec_Gprot];
		}
		break;
	}
	default: {
		err = 1;
		break;
	}
	}

	return (err);
}

void coupling_implicit(double t, double y[], 
                       const grid_parms& grid, 
                       SMC_cell* __restrict__ smc, 
                       EC_cell* __restrict__ ec, 
                       const conductance& cpl_cef)
{
////******************** HOMOCELLULAR COUPLING *********************/
	int nc = grid.num_smc_circumferentially;
	int na = grid.num_smc_axially;
    const int ng = grid.num_ghost_cells;
#pragma ivdep
#pragma omp parallel for
    for (int ij0 = 0; ij0 < nc * na; ij0++) {
	    int i = ij0 / (na + ng) + 1;
        int j = ij0 % (na + ng) + 1;
        int ij = i*(na + ng) + j;
		const double vSmc_Vm = smc[ij].vars[smc_Vm];
		int up = j - 1, down = j + 1, left = i - 1, right = i + 1;
        const double* __restrict__ vars = smc[ij].vars;
        const double* __restrict__ upVars = smc[i*(na + ng) + up].vars;
        const double* __restrict__ downVars = smc[i*(na + ng) + down].vars;
        const double* __restrict__ leftVars = smc[left*(na + ng) + j].vars;
        const double* __restrict__ rightVars = smc[right*(na + ng) + j].vars;
        const double var = vars[smc_Vm];
        const double upVar = upVars[smc_Vm];
        const double downVar = downVars[smc_Vm];
        const double leftVar =  leftVars[smc_Vm];
        const double rightVar = rightVars[smc_Vm];
        const double* __restrict__ diffusion = cpl_cef.smc_diffusion;
        double* __restrict__ homo_fluxes = smc[ij].homo_fluxes;

		homo_fluxes[cpl_Vm] = -cpl_cef.Vm_hm_smc
			  * (diffusion[0] * (var - upVar)
			  +  diffusion[1] * (var - downVar)
			  +  diffusion[2] * (var - leftVar)
			  +  diffusion[3] * (var - rightVar));

	}	//end ij

    nc = grid.num_ec_circumferentially;
    na = grid.num_ec_axially;
#pragma ivdep
#pragma omp parallel for
    for (int ij0 = 0; ij0 < nc * na; ij0++) {
        int i = ij0 / (na + ng) + 1;
        int j = ij0 % (na + ng) + 1;
        int ij = i*(na + ng) + j;
		int up = j - 1, down = j + 1, left = i - 1, right = i + 1;
        const double* __restrict__ vars =  ec[ij].vars;
        const double* __restrict__ upVars = ec[i*(na + ng) + up].vars;
        const double* __restrict__ downVars = ec[i*(na + ng) + down].vars;
        const double* __restrict__ leftVars = ec[left*(na + ng) + j].vars;
        const double* __restrict__ rightVars = ec[right*(na + ng) + j].vars;
        const double var = vars[ec_Vm];
        const double upVar = upVars[ec_Vm];
        const double downVar = downVars[ec_Vm];
        const double leftVar = leftVars[ec_Vm];
        const double rightVar = rightVars[ec_Vm];
        const double* __restrict__ diffusion = cpl_cef.ec_diffusion;
        double* __restrict__ homo_fluxes = ec[ij].homo_fluxes;

		ec[ij].homo_fluxes[cpl_Vm] = -cpl_cef.Vm_hm_ec
			  * (diffusion[0] * (var - upVar)
			  +  diffusion[1] * (var - downVar)
			  +  diffusion[2] * (var - leftVar)
			  +  diffusion[3] * (var - rightVar));

	}	//end ij

////******************** HETEROCELLULAR COUPLING *********************/
	int offset_smc_circumferentially, offset_ec_axially;

    nc = grid.num_smc_circumferentially;
    na = grid.num_smc_axially;
#pragma omp parallel for
    for (int ij0 = 0; ij0 < nc * na; ij0++) {
	    int i = ij0 / (na + ng) + 1; // x dim
        int j = ij0 % (na + ng) + 1; // y dim
        int ij = i*(na + ng) + j;
		int l = (j - 1) / grid.num_smc_fundblk_axially + 1;
        const double smcIJ = smc[ij].vars[smc_Vm];
		double dummy_smc[3] = { 0.0, 0.0, 0.0 };
#pragma ivdep
		for (int k = 1 + (i - 1) * 5; k <= i * 5; k++) {
			int kl = k*(na + ng) + l;
		    dummy_smc[cpl_Vm] = dummy_smc[cpl_Vm] + (smcIJ - ec[kl].vars[ec_Vm]);
		}
		smc[ij].hetero_fluxes[cpl_Vm] = -cpl_cef.Vm_ht_smc * dummy_smc[cpl_Vm];
	}

    nc = grid.num_ec_circumferentially;
    na = grid.num_ec_axially;
#pragma omp parallel for
        for (int ij0 = 0; ij0 < nc * na; ij0++) {
	    int i = ij0 / (na + ng) + 1;
		int j = ij0 % (na + ng) + 1;
		int ij = i*(na + ng) + j;
        int k = (i - 1) / 5 + 1;
        const double ecIJ = ec[ij].vars[ec_Vm];
		double dummy_ec[3] = { 0.0, 0.0, 0.0 };
#pragma ivdep
		for (int l = 1 + (j - 1) * 13; l <= j * 13; l++) {
			int kl = k*(na + ng) + l;
		    dummy_ec[cpl_Vm] = dummy_ec[cpl_Vm] + (ecIJ - smc[kl].vars[smc_Vm]);
		}
		ec[ij].hetero_fluxes[cpl_Vm] = -cpl_cef.Vm_ht_ec * dummy_ec[cpl_Vm];
	}
}

void coupling_explicit(double t, double y[], 
                       const grid_parms& grid, 
                       SMC_cell* __restrict__ smc, 
                       EC_cell* __restrict__ ec, 
                       const conductance& cpl_cef)
{
////******************** HOMOCELLULAR COUPLING *********************/

    int nc = grid.num_smc_circumferentially;
    int na = grid.num_smc_axially;
    int ng = grid.num_ghost_cells;
#pragma ivdep
#pragma omp parallel for
	for (int ij0 = 0; ij0 < nc * na; ij0++) {
	    int i = ij0 / (na + ng) + 1;
        int j = ij0 % (na + ng) + 1;
        int ij = i*(na + ng) + j;
		int up = j - 1, down = j + 1, left = i - 1, right = i + 1;

        const double* __restrict__ vars = smc[ij].vars;
        const double* __restrict__ upVars = smc[i*(na + ng) + up].vars;
        const double* __restrict__ downVars = smc[i*(na + ng) + down].vars;
        const double* __restrict__ leftVars = smc[left*(na + ng) + j].vars;
        const double* __restrict__ rightVars = smc[right*(na + ng) + j].vars;
        double* __restrict__ homo_fluxes = smc[ij].homo_fluxes;

		homo_fluxes[cpl_Ca] = -cpl_cef.Ca_hm_smc
		                        * (cpl_cef.smc_diffusion[0] * ((vars[smc_Ca] - upVars[smc_Ca]))
					+ (cpl_cef.smc_diffusion[1] * (vars[smc_Ca] - downVars[smc_Ca]))
					+ (cpl_cef.smc_diffusion[2] * (vars[smc_Ca] - leftVars[smc_Ca]))
					+ (cpl_cef.smc_diffusion[3] * (vars[smc_Ca] - rightVars[smc_Ca])));

		homo_fluxes[cpl_IP3] = -cpl_cef.IP3_hm_smc
					* (cpl_cef.smc_diffusion[0] * ((vars[smc_IP3] - upVars[smc_IP3]))
					+ (cpl_cef.smc_diffusion[1] * (vars[smc_IP3] - downVars[smc_IP3]))
					+ (cpl_cef.smc_diffusion[2] * (vars[smc_IP3] - leftVars[smc_IP3]))
					+ (cpl_cef.smc_diffusion[3] * (vars[smc_IP3] - rightVars[smc_IP3])));
	}	//end ij

    nc = grid.num_ec_circumferentially;
    na = grid.num_ec_axially;
#pragma omp parallel for
	for (int ij0 = 0; ij0 < nc * na; ij0++) {
        int i = ij0 / (na + ng) + 1;
        int j = ij0 % (na + ng) + 1;
        int ij = i*(na + ng) + j;
        int up = j - 1, down = j + 1, left = i - 1, right = i + 1;

		const double* __restrict__ vars = ec[ij].vars;
		const double* __restrict__ upVars = ec[i*(na + ng) + up].vars;
        const double* __restrict__ downVars = ec[i*(na + ng) + down].vars;
        const double* __restrict__ leftVars = ec[left*(na + ng) + j].vars;
        const double* __restrict__ rightVars = ec[right*(na + ng) + j].vars;
        double* __restrict__ homo_fluxes = ec[ij].homo_fluxes;

		homo_fluxes[cpl_Ca] = -cpl_cef.Ca_hm_ec
					* cpl_cef.ec_diffusion[0] * ((vars[ec_Ca] - upVars[ec_Ca])
					+ cpl_cef.ec_diffusion[1] * (vars[ec_Ca] - downVars[ec_Ca])
					+ cpl_cef.ec_diffusion[2] * (vars[ec_Ca] - leftVars[ec_Ca])
					+ cpl_cef.ec_diffusion[3] * (vars[ec_Ca] - rightVars[ec_Ca]));

		homo_fluxes[cpl_IP3] = -cpl_cef.IP3_hm_ec
					* cpl_cef.ec_diffusion[0] * ((vars[ec_IP3] - upVars[ec_IP3])
					+ cpl_cef.ec_diffusion[1] * (vars[ec_IP3] - downVars[ec_IP3])
					+ cpl_cef.ec_diffusion[2] * (vars[ec_IP3] - leftVars[ec_IP3])
					+ cpl_cef.ec_diffusion[3] * (vars[ec_IP3] - rightVars[ec_IP3]));

	}	//end ij

////******************** HETEROCELLULAR COUPLING *********************/
	int offset_smc_circumferentially, offset_ec_axially;

    nc = grid.num_smc_circumferentially;
    na = grid.num_smc_axially;
#pragma omp parallel for
    for (int ij0 = 0; ij0 < nc * na; ij0++) {
	    int i = ij0 / (na + ng) + 1; // x dim
	    int j = ij0 % (na + ng) + 1; // y dim
	    int ij = i*(na + ng) + j;
        int l = (j - 1) / grid.num_smc_fundblk_axially + 1;
	    double dummy_smc[3] = { 0.0, 0.0, 0.0 };
        const double* __restrict__ smcIJVars = smc[ij].vars;
        double* __restrict__ hetero_fluxes = smc[ij].hetero_fluxes;
		for (int k = 1 + (i - 1) * 5; k <= i * 5; k++) {
			int kl = k*(na + ng) + l;
            const double* __restrict__ ecKLVars = ec[kl].vars;
		    dummy_smc[cpl_Ca] = dummy_smc[cpl_Ca] + (smcIJVars[smc_Ca] - ecKLVars[ec_Ca]);
			dummy_smc[cpl_IP3] = dummy_smc[cpl_IP3] + (smcIJVars[smc_IP3] - ecKLVars[ec_IP3]);
		}
		hetero_fluxes[cpl_Ca] = -cpl_cef.Ca_ht_smc * dummy_smc[cpl_Ca];
		hetero_fluxes[cpl_IP3] = -cpl_cef.IP3_ht_smc * dummy_smc[cpl_IP3];
	}

    nc = grid.num_ec_circumferentially;
    na = grid.num_ec_axially;
#pragma omp parallel for
    for (int ij0 = 0; ij0 < nc * na; ij0++) {
	    int i = ij0 / (na + ng) + 1;
        int j = ij0 % (na + ng) + 1;
        int ij = i*(na + ng) + j;
		int k = (i - 1) / 5 + 1;
        const double* __restrict__ ecIJVars = ec[ij].vars;
        double* __restrict__ hetero_fluxes = ec[ij].hetero_fluxes;
		double dummy_ec[3] = { 0.0, 0.0, 0.0 };
		for (int l = 1 + (j - 1) * 13; l <= j * 13; l++) {
			int kl = k*(na + ng) + l;
		    const double* __restrict__ smcKLVars = smc[kl].vars;
		    dummy_ec[cpl_Ca] = dummy_ec[cpl_Ca] + (ecIJVars[ec_Ca] - smcKLVars[smc_Ca]);
			dummy_ec[cpl_IP3] = dummy_ec[cpl_IP3] + (ecIJVars[ec_IP3] - smcKLVars[smc_IP3]);
		}
		hetero_fluxes[cpl_Ca] = -cpl_cef.Ca_ht_ec * dummy_ec[cpl_Ca];
	    hetero_fluxes[cpl_IP3] = -cpl_cef.IP3_ht_ec * dummy_ec[cpl_IP3];
	}
}

void coupling(double t, double y[], const grid_parms& grid, 
              SMC_cell* smc, EC_cell* ec, 
              const conductance& cpl_cef)
{
	coupling_implicit(t, y, grid, smc, ec, cpl_cef);
	coupling_explicit(t, y, grid, smc, ec, cpl_cef);
}


void compute(const grid_parms& grid, 
             SMC_cell* smc, 
             EC_cell* ec, 
             const conductance& cpl_cef, 
             double t, double* y, double* f)
{
	int err;

	map_solver_output_to_cells(grid, y, smc, ec);

#if CELL_MODEL == TSK

	tsoukias_smc(grid, smc);
	koenigsberger_ec(grid, ec);

	coupling(t, y, grid, smc, ec, cpl_cef);

	tsoukias_smc_derivatives(f, grid, smc);
	koenigsberger_ec_derivatives(t, f, grid, ec);

#elif CELL_MODEL == KNBGR

	koenigsberger_smc(grid, smc);
	koenigsberger_ec(grid, ec);

	coupling(t, y, grid, smc, ec, cpl_cef);

#if PLOTTING && EXPLICIT_ONLY
	bufferPos = 0;
#endif

	koenigsberger_smc_derivatives(f, grid, smc);
	koenigsberger_ec_derivatives(t, f, grid, ec);

#if PLOTTING && EXPLICIT_ONLY

	if (grid.universal_rank == RANK)
	{
		for (int i = 0; i < OUTPUT_PLOTTING_SIZE; i++)
		{
			fprintf(var_file, "%f,", plotttingBuffer[i]);
		}
		fprintf(var_file, "%f\n", t);
	}
#endif


#endif
}

void compute_implicit(const grid_parms& grid, 
                      SMC_cell* __restrict__ smc, EC_cell* __restrict__ ec, 
                      const conductance& cpl_cef, 
                      double t, double* y, double* f)
{
	map_solver_output_to_cells(grid, y, smc, ec);

#if CELL_MODEL == TSK

	tsoukias_smc(grid, smc);
	koenigsberger_ec(grid, ec);

	coupling(t, y, grid, smc, ec, cpl_cef);

	tsoukias_smc_derivatives(f, grid, smc);
	koenigsberger_ec_derivatives(t, f, grid, ec);

#elif CELL_MODEL == KNBGR

	koenigsberger_smc_implicit(grid, smc);
	koenigsberger_ec_implicit(grid, ec);

	coupling_implicit(t, y, grid, smc, ec, cpl_cef);

	koenigsberger_smc_derivatives_implicit(f, grid, smc);
	koenigsberger_ec_derivatives_implicit(t, f, grid, ec);

#endif
}

void compute_explicit(const grid_parms& grid, 
                      SMC_cell* __restrict__ smc, 
                      EC_cell* __restrict__ ec, 
                      const conductance& cpl_cef, 
                      double t, double* y, double* f)
{
	map_solver_output_to_cells(grid, y, smc, ec);

#if CELL_MODEL == TSK

	tsoukias_smc(grid, smc);
	koenigsberger_ec(grid, ec);

	coupling(t, y, grid, smc, ec, cpl_cef);

	tsoukias_smc_derivatives(f, grid, smc);
	koenigsberger_ec_derivatives(t, f, grid, ec);

#elif CELL_MODEL == KNBGRs

	koenigsberger_smc_explicit(grid, smc);
	koenigsberger_ec_explicit(grid, ec);

	coupling_explicit(t, y, grid, smc, ec, cpl_cef);

	koenigsberger_smc_derivatives_explicit(f, grid, smc);
	koenigsberger_ec_derivatives_explicit(t, f, grid, ec);

#endif
}

#if 0
// TODO: These functions are to be moved to a util_func.cpp or something of that sort.
/************************************************************/
void minimum(double* table, int size, double *value, int *index) {
	///For evaluating minimum of an array.
	*value = table[0];
	for (int i = 0; i < size; i++) {
		if (*value > table[i]) {
			*value = table[i];
			*index = i;
		}
	}
}

/************************************************************/
void maximum(double* table, int size, double *value, int *index) {
	///For evaluating maximum of an array.
	*value = table[0];
	for (int i = 0; i < size; i++) {
		if (*value < table[i]) {
			*value = table[i];
			*index = i;
		}
	}
}

/************************************************************/
void average(double* table, int size, double *value) {
	///For evaluating average of an array.
	*value = 0;
	for (int i = 0; i < size; i++) {
		*value += table[i];
	}
	*value = *value / (double) (size);
}
#endif
