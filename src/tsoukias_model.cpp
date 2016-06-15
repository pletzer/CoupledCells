/*
 * tsoukias_smc_multicell.cpp
 *
 *  Created on: 16/12/2012
 *      Author: mohsinshaikh
 */

#include <math.h>

#include "computelib.h"
#include "tsoukias_model.h"

/**** Standard model parameter values ****/
double z_K = 1.0;
double z_Na = 1.0;
double z_Ca = 2.0;
double z_Cl = -1.0;
//Avogadro's number
double N_Avogadro = 6.022e23;
//Faraday's constant (C/mol)
double F = 96487.0;
//Universal Gas constant (mJ/mol*K)
double R_const = 8.341e3;
//Temperature (K)
double T = 293.0;
//Cell membrane capacitance (pF)
double C_m = 25.0, Am = C_m / 1e6;	//cm^2
//Extracellular concentrations of Ca; Na; K and Cl (mM)
double Ca_o = 2.0, Na_o = 140.0, Cl_o = 129.0, K_o = 5.0;
//Cell volumes (pL)
//cytosolic
double vol_i = 1.0;
//available to free Ca
double vol_Ca = 0.7;
/**** L-type voltage operated Ca2+ channels ****/
//Whole cell permeability (cm/s)
double P_VOCC = 1.88e-5;
/**** Large conductance Ca2+ activated potassium channels ****/
//Single channel permeability (cm^3/s)
double P_BKCa = 3.9e-13;
//Channel density
double N_BKCa = 6.6e6;
//Fast and slow activation time constants (ms)
double tau_pf = 0.84, tau_ps = 35.9;
//Maximum NO-induced V_1/2KCa shift (mV)
double dV_half_KCaNO = 46.3;
//Maximum cGMP-induced V_1/2KCa shift (mV)
double dV_half_KCacGMP = 76.0;
/**** Voltage dependent K channels ****/
//Maximum whole cell conductance (nS)
double g_Kv = 1.35;
//Fast and slow inactivation time constants (ms)
double tau_q1 = 371.00, tau_q2 = 2884.00;
/**** ATP sensetive K channels	****/
//Whole cell conductance (nS)
double g_K_leak = 0.067;
/**** Nonselective cation channels ****/
//DAG concentration for half maximum activation (nM)
double K_NSC = 3000e-6;		//corrected and compared with jsim model
//and changed to mM
//Whole cell Na; K ; Ca permeabilities (cm/s)
double P_NaNSC = 5.11e-7, P_KNSC = 1.06 * P_NaNSC, P_CaNSC = 4.54 * P_NaNSC;
//Constituent activation
double d_NSCmin = 0.0244;
/**** Store operated nonselective cation channels ****/
//Whole cell conductances to Ca and Na (nM)
double g_SOCCa = 0.0083, g_SOCNa = 0.0575;
//Ca concentration in SR for half activation (nM)
double K_SOC = 100e-6;		//changed to mM
/**** Calcium activated chloride channels ****/
//Maximum conductance (nS/pF)
double g_ClCa = 0.23;
//cGMP independent component
double R_ClcGMPmin = 0.0132;
//Hill coefficient
double n_ClCa = 2;
//EC50 for Ca dependent activation (nM)
double K_ClCa = 365.00;
//Hill coefficient;
double n_ClcGMP = 3.3;
//EC50 for cGMP dependent activation (uM)
double K_ClcGMP = 6.4;
/**** Plasma membrane Ca ATPase ****/
//Maximum current (pA)
double I_bar_PMCA = 5.37;
//Michaelis constant (nM)
double K_mPMCA = 170e-6;		//changed to mM
/**** Plasma membrane Na-Ca exchange ****/
//Scaling factor (pS)
double g_NCX = 4.87e-4;		//0.0487;
double gama_NCX = 0.45, d_NCX = 3e-4;
/**** Sodium potasium pump ****/
//Maximum current density (pA/pF)
double I_bar_NaK = 2.3083;
//Binding constants and Hill coefficients (mM) for K_out and Na_intracellular
double K_dKo = 1.6,
		K_dNai = 22.0,
		n_HKo = 1.1,
		n_HNai = 1.7;
/**** Sodium potassium chloride cotransport ****/
//Cotransport coefficient (nmol^2/(J*s*cm^2))
double L_NaKCl = 1.79e-11;		//changed to (mmol^2/(J*s*cm^2))
/**** Sarcoplasmic reticulum ****/
/** (parameters as in Yang(2003) ) **/
//Time constant of the internal diffusion (ms)
double tau_tr = 1000.00;
//Time constant of the diffusion from SR release compartment (ms)
double tau_rel = 0.0333;
//Activation rate constant of RyR (1/ (ms * mM^2))
double Kr1 = 2500.00;
//Inactivation rate constant of RyR (1/ (ms * mM))
double Kr2 = 1.05;
//Unbinding rate constant from activation (1/ms)
double K_r1 = 0.0076;
//Unbinding rate constant from inactivation (1/ms)
double K_r2 = 0.084;
/** (Rest of the parameters as in Tsoukias(2008) ) **/
//Michaelis constant of SR Ca pump (uM)
double K_mUp = 1.0;
//Oscillation and non oscillation parameter
double k_leak = 4.0;
//Maximum SR uptake current (pA)
double I_bar_SERCA = 20.0;//(3.34)*(k_leak+1);               //ranging between 6.68 - 20 pA
//Leakage parameter
double R_leak = 5.35e-5;//1.07e-5*(k_leak);         //ranging between 1.07e-5-5.35e-5
//Binding affinity of calsequestrin (mM)
double K_CSQN = 0.8;
//Concentration of calsequestrin in the release compartment (mM)
double CSQN = 15.0;
//Rate constant of Ca release by IP3R (1/ms)
double I_bar_IP3 = 2880e-6;
//Dissociation constants for Ca activation and inhibitory sites (nM)
double K_actIP3 = 170e-6;		//changed to uM
double K_inhIP3 = 100e-6;		//changed to uM
//Dissociation constant for IP3 binding to IP3R (nM)
double K_IP3 = 120e-6;		//changed to uM
//Rate of Ca binding to the inhibitory site (1/(mM*ms))
double k_on_IP3 = 1.4;
//Volume of uptake compartment  (pl)
double vol_u = 0.07;
//Volume of release compartment (pl)
double vol_r = 0.007;
/**** alpha 1 adrenoreceptor activation and IP3 formation ****/
/** (parameters as in Bennett et al (2005) ) **/
//Total number of receptors
double R_T_G = 2.0e4;
//Unphosphorylated receptor dissociation constanst (mM)
double K_1_G = 0.01;
//Phosphorylated receptor dissociation constant (mM)
double K_2_G = 0.2;
//Receptor recycling rate (1/ms)
double k_r_G = 1.75e-7;
/*kpG declared below */
//Receptor endocytosis rate (1/ms)
double k_e_G = 6e-6;
//Fraction of mobile receptors
double ksi_G = 0.85;
//Total number of G protein molecules
double G_T_G = 1e5;
//IP3 degradation rate (1/ms)
double k_deg_G = 1.25e-3;
//G protein activation rate (1/ms)
double k_a_G = 0.17e-3;
//G protein deactivation rate (1/ms)
double k_d_G = 1.5e-3;
//Total PIP2 molecules
double PIP2_T = 5e7;
//PIP2 replenishment rate (1/ms)
double r_r_G = 0.015e-3;
//Dissociation constant for Ca binding to PLC (mM)
double K_c_G = 0.4e-3;
//Effective signal gain parameter (1/ms)
double alpha_G = 2.781e-8;
//parameter
double gama_G = N_Avogadro * vol_i * 1e-15;
/** (Rest of the parameters as in Tsoukias(2008) ) **/
//Receptor phosphorylation rate (1/s)
double k_p_G = 0.0;	//set it to 0.1(1/ms) while replicating the results in Fig3 of (Tsoukias 2008)
double deltaG_0 = 0;
/**** sGC activation and cGMP formation ****/
//sGC activation and inactivation time constants (s)
double tau_a_sGC = 0.23e3;		//changed to ms
double tau_d_sGC = 10.0e3;		//changed to ms
//Maximum cGMP formation rate (mM/ms)
double V_cGMPmax = 1.26e-7;
//Michaelis Menten Constant (mM)
double K_m_pde = 1e-3, k1_sGC = 2e3,	//Units (1/mM*ms)
		k2_sGC = 0.64e-5,		//Units (1/ms)
		k3_sGC = 4.2,		//Units (1/mM*ms)
		kD_tau_sGC = 0.1e-3,	//Units (1/ms)
		k_1_sGC = 15e-3,		//Units (1/ms)
		k_2_sGC = 0.1e-6,		//Units (1/ms)
		kD_sGC = 0.4e-3,		//Units (1/ms)
		k_4_sGC = 0.1e-3,		//Units (1/ms)
		B5_sGC = k2_sGC / k3_sGC, k_pde_cGMP = 0.0695e-3;
/**** Ca buffering in the cytosol ****/
//Concentration of calmodulin and other buffers (mM) and their dissociation constants (nM)
double S_CM = 0.1,			//mM
		B_F = 0.1,			//mM
		K_d = 260.0e-6,		//nM
		K_dB = 530.0e-6;		//nM

/*******************************************************************************************/
void Initialize_tsoukias_smc(grid_parms grid, double y[], SMC_type&  smc)
/*******************************************************************************************/
{

	int k = 0, offset;

		for (int i = 1; i <= grid.num_smc_circumferentially; i++) {
			for (int j = 1; j <= grid.num_smc_axially; j++) {
				if (i > 1)
					k = ((i - 1) * grid.neq_smc_axially);
				else if (i == 1)
					k = 0;

				y[k + ((j - 1) * grid.neq_smc) + smc_Vm] 	= -59.4;		//mV
				y[k + ((j - 1) * grid.neq_smc) + smc_Ca]	= 68-6; 		//mM
				y[k + ((j - 1) * grid.neq_smc) + smc_Ca_u]	= 0.66;			//mM
				y[k + ((j - 1) * grid.neq_smc) + smc_Ca_r]	= 0.57;			//mM
				y[k + ((j - 1) * grid.neq_smc) + smc_Na_i]	= 8.4;			//mM
				y[k + ((j - 1) * grid.neq_smc) + smc_K_i]	= 140.0;		//mM
				y[k + ((j - 1) * grid.neq_smc) + smc_Cl_i] 	= 59.4;			//mM
				y[k + ((j - 1) * grid.neq_smc) + smc_IP3]	= 0.0;
				y[k + ((j - 1) * grid.neq_smc) + smc_DAG]	= 0.0;
				y[k + ((j - 1) * grid.neq_smc) + smc_cGMP_i]= 0.0;
				y[k + ((j - 1) * grid.neq_smc) + smc_V_cGMP]= 0.0;
				y[k + ((j - 1) * grid.neq_smc) + smc_d_L]   = 1 / (1 + exp(-y[k + ((j - 1) * grid.neq_smc) + smc_Vm] / 8.3));
				y[k + ((j - 1) * grid.neq_smc) + smc_f_L]   = 1 / (1 + exp((y[k + ((j - 1) * grid.neq_smc) + smc_Vm] + 42.0) / 9.1));

				smc.flux(i, j, V_half_KCa) = -41.7 * log10(y[k + ((j-1) * grid.neq_smc) + smc_Ca]) - 128.2;
				y[k + ((j - 1) * grid.neq_smc) + smc_p_f]   = 1 / (1 + exp(-(y[k + ((j - 1) * grid.neq_smc) + smc_Vm] - smc.flux(i, j, V_half_KCa)) / 18.25));


				y[k + ((j - 1) * grid.neq_smc) + smc_p_K] 	= 1 / (1 + exp(-(y[k + ((j - 1) * grid.neq_smc) + smc_Vm] + 11.0) / 15.0));
				y[k + ((j - 1) * grid.neq_smc) + smc_q_1] 	= 1 / (1 + exp((y[k + ((j - 1) * grid.neq_smc) + smc_Vm] + 40.0) / 14.0));
				y[k + ((j - 1) * grid.neq_smc) + smc_q_2] 	= 1 / (1 + exp((y[k + ((j - 1) * grid.neq_smc) + smc_Vm] + 40.0) / 14.0));
				y[k + ((j - 1) * grid.neq_smc) + smc_R_01]  = 0.995;
				y[k + ((j - 1) * grid.neq_smc) + smc_R_10]  = 0.0033;
				y[k + ((j - 1) * grid.neq_smc) + smc_R_11]  = 4.0e-6;
				y[k + ((j - 1) * grid.neq_smc) + smc_h_IP3] = K_inhIP3 / (y[k + ((j - 1) * grid.neq_smc) + smc_Ca] + K_inhIP3);

				y[k + ((j - 1) * grid.neq_smc) + smc_PIP2] = PIP2_T - (1 + (k_deg_G / r_r_G)) * gama_G * y[k + ((j - 1) * grid.neq_smc) + smc_IP3];

				smc.flux(i, j, r_h_G) 	=   k_deg_G * gama_G * y[k + ((j - 1) * grid.neq_smc) + smc_IP3] / y[k + ((j - 1) * grid.neq_smc) + smc_PIP2];

				y[k + ((j - 1) * grid.neq_smc) + smc_R_S_G] = R_T_G * ksi_G;
				y[k + ((j - 1) * grid.neq_smc) + smc_R_S_P_G] = 0.0;
				y[k + (j * grid.neq_smc) + smc_G] = smc.flux(i, j, r_h_G) * (K_c_G *  y[k + ((j - 1) * grid.neq_smc) + smc_Ca]) / (alpha_G * y[k + ((j - 1) * grid.neq_smc) + smc_Ca]);

				smc.flux(i, j, delta_G) =   k_d_G * y[k + ((j - 1) * grid.neq_smc) + smc_G] / (k_a_G * (G_T_G - y[k + ((j - 1) * grid.neq_smc) + smc_G]));

				smc.NE(i, j) = 0.0;	//mM
				smc.NO(i, j) = 0.0;	//nM
				smc[i][j].I_stim = 0.0;	//pA
			}
		}

		for (int i = 0; i < (grid.num_smc_circumferentially + grid.num_ghost_cells);
				i++) {
			for (int j = 0; j < (grid.num_smc_axially + grid.num_ghost_cells);
					j++) {
				smc.var(i, j, smc_Vm)				= 0.0;
				smc.var(i, j, smc_d_L) 			= 0.0;
				smc.var(i, j, smc_f_L) 		 	= 0.0;
				smc.var(i, j, smc_p_f) 			= 0.0;
				smc.var(i, j, smc_p_s) 			= 0.0;
				smc.var(i, j, smc_q_1) 			= 0.0;
				smc.var(i, j, smc_q_2) 			= 0.0;
				smc.var(i, j, smc_p_K) 			= 0.0;
				smc.var(i, j, smc_Ca_u) 			= 0.0;
				smc.var(i, j, smc_Ca_r) 			= 0.0;
				smc.var(i, j, smc_R_10) 			= 0.0;
				smc.var(i, j, smc_R_11) 			= 0.0;
				smc.var(i, j, smc_R_01) 			= 0.0;
				smc.var(i, j, smc_h_IP3) 			= 0.0;
				smc.var(i, j, smc_R_S_G) 			= 0.0;
				smc.var(i, j, smc_R_S_P_G)		= 0.0;
				smc.var(i, j, smc_G) 				= 0.0;
				smc.var(i, j, smc_IP3) 			= 0.0;
				smc.var(i, j, smc_PIP2) 			= 0.0;
				smc.var(i, j, smc_V_cGMP)			= 0.0;
				smc.var(i, j, smc_cGMP_i)			= 0.0;
				smc.var(i, j, smc_Ca)				= 0.0;
				smc.var(i, j, smc_Na_i) 			= 0.0;
				smc.var(i, j, smc_K_i) 			= 0.0;
				smc.var(i, j, smc_Cl_i)			= 0.0;
				smc.var(i, j, smc_DAG)			= 0.0;

				for (int k = 1; k <= grid.num_fluxes_smc; k++) {
					smc[i][j].fluxes[k - 1] = 0.0;
				}
				for (int k = 1; k <= grid.num_coupling_species_smc; k++) {
					smc[i][j].homo_fluxes[k - 1] = 0.0;
					smc[i][j].hetero_fluxes[k - 1] = 0.0;
				}
			}
		}
}
/************************************************************************************************/
void tsoukias_smc(grid_parms grid, SMC_type&  smc)
/*																				   		      	*
 * This is the multicell version of the function tosukias_smc() used for single   		  	  	*
 * EC-SMC Unit simulation.														  		  		*
 * int i and j are the indicies of the SMC being simulated.						  		  		*
 * Int NO_Path and cGMP_Path are the two switches to enable or disable the		  		  		*
 * corresponding pathways.														  		  		*
 ************************************************************************************************/
{

// --------------------------------------------------------------------------------------------------------
//   SYSTEM OF EQNS OF A SINGLE VASCULAR SMOOTH MUSCLE CELL   E L E C T R O P H Y S I O L O G Y   M O D E L
// --------------------------------------------------------------------------------------------------------
	for (int i = 1; i <= grid.num_smc_circumferentially; i++) {
		for (int j = 1; j <= grid.num_smc_axially; j++) {

			smc.var(i, j, smc_DAG) = smc.var(i, j, smc_IP3);

//SMC single cell ionic currents

			smc.flux(i, j, E_K) = (R_const * T / (z_K * F)) * log(K_o / smc.var(i, j, smc_K_i));
			smc.flux(i, j, E_Ca) = (R_const * T / (z_Ca * F)) * log(Ca_o / smc.var(i, j, smc_Ca));
			smc.flux(i, j, E_Na) = (R_const * T / (z_Na * F)) * log(Na_o / smc.var(i, j, smc_Na_i));
			smc.flux(i, j, E_Cl) = (R_const * T / (z_Cl * F)) * log(Cl_o / smc.var(i, j, smc_Cl_i));

/**** Membrane Electrophysiology ****/

			//L-type voltage operated Ca channels

			smc.flux(i, j, dbar_L) = 1 / (1 + exp(-smc.var(i, j, smc_Vm) / 8.3));
			smc.flux(i, j, tau_dL) = 1.15 + 2.5 * exp(-P2((smc.var(i, j, smc_Vm) + 40.0) / 30.0));
			smc.flux(i, j, fbar_L) = 1 / (1 + exp((smc.var(i, j, smc_Vm) + 42.0) / 9.1));
			smc.flux(i, j, tau_fL) = 45.0 + 65.0 * exp(-P2((smc.var(i, j, smc_Vm) + 35.0) / 25.0));
			smc.flux(i, j, I_VOCC) =
			    1e6 * Am * P_VOCC * smc.var(i, j, smc_d_L) *
			    smc.var(i, j, smc_f_L) * smc.var(i, j, smc_Vm) *
			    (P2(z_Ca * F) / (R_const * T)) *
			    ((Ca_o -
			      (smc.var(i, j, smc_Ca) *
			       exp((smc.var(i, j, smc_Vm) * z_Ca * F) /
				   (R_const * T)))) / (1 - exp((smc.var(i, j, smc_Vm) * z_Ca * F) / (R_const * T))));

			//Large conductance Ca activated potassium channels

			smc.flux(i, j, R_NO) = grid.NO_path * smc.NO(i, j) / (smc.NO(i, j) + 200.00);
			smc.flux(i, j, R_cGMP) = grid.cGMP_path * P2(smc.var(i, j, smc_cGMP_i)) / (P2(smc.var(i, j, smc_cGMP_i)) + P2(1.5));
			smc.flux(i, j, V_half_KCa) = -41.7 * log10(smc.var(i, j, smc_Ca)) - 128.2;

			smc.flux(i, j, pbar_o) = 1 / (1 + exp(-(smc.var(i, j, smc_Vm) - smc.flux(i, j, V_half_KCa)) / 18.25));
			smc.flux(i, j, P_KCa) = 0.17 * smc.var(i, j, smc_p_f) + 0.83 * smc.var(i, j, smc_p_s);
			smc.flux(i, j, i_KCa) =
			    1e6 * P_BKCa * smc.var(i, j, smc_Vm) * (P2(F) /
								   (R_const * T)) *
			    ((K_o - (smc.var(i, j, smc_K_i) * exp(smc.var(i, j, smc_Vm) * F / (R_const * T)))) / (1 - exp(smc[i]
														 [j].vars[smc_Vm]
														 * F / (R_const * T))));
			smc.flux(i, j, I_BKCa) = Am * N_BKCa * smc.flux(i, j, P_KCa) * smc.flux(i, j, i_KCa);

			//Voltage dependent K channels

			smc.flux(i, j, qbar) = 1.0 / (1 + exp((smc.var(i, j, smc_Vm) + 40.0) / 14.0));
			smc.flux(i, j, tau_pK) = 61.5 * exp(-0.027 * smc.var(i, j, smc_Vm));
			smc.flux(i, j, pbar_K) = 1.0 / (1 + exp(-(smc.var(i, j, smc_Vm) + 11.0) / 15.0));
			smc.flux(i, j, I_Kv) =
			    g_Kv * smc.var(i, j, smc_p_K) * (0.45 *
							   smc.var(i, j, smc_q_1) +
							   0.55 * smc.var(i, j, smc_q_2)) * (smc.var(i, j, smc_Vm) - smc.flux(i, j, E_K));

			//Unspecified K Leak channels
			smc.flux(i, j, I_Kleak) = g_K_leak * (smc.var(i, j, smc_Vm) - smc.flux(i, j, E_K));

			//Non-selective cation channels
			smc.flux(i, j, Po_NSC) = 0.4344 + ((1 - 0.4344) / (1 + exp(-(smc.var(i, j, smc_Vm) - 47.12)
										/ 24.24)));
			smc.flux(i, j, INa_NSC) =
			    1e6 * Am *
			    ((smc.var(i, j, smc_DAG) /
			      (smc.var(i, j, smc_DAG) + K_NSC)) +
			     d_NSCmin) * smc.flux(i, j, Po_NSC) * P_NaNSC *
			    smc.var(i, j, smc_Vm) * (P2(F) / (R_const * T)) *
			    ((Na_o - (smc.var(i, j, smc_Na_i) * exp(smc.var(i, j, smc_Vm) * F / (R_const * T)))) / (1 - exp(smc[i]
														   [j].vars[smc_Vm]
														   * F / (R_const * T))));
			smc.flux(i, j, IK_NSC) =
			    1e6 * Am *
			    ((smc.var(i, j, smc_DAG) /
			      (smc.var(i, j, smc_DAG) + K_NSC)) +
			     d_NSCmin) * smc.flux(i, j, Po_NSC) * P_KNSC *
			    smc.var(i, j, smc_Vm) * (P2(F) / (R_const * T)) *
			    ((K_o - (smc.var(i, j, smc_K_i) * exp(smc.var(i, j, smc_Vm) * F / (R_const * T)))) / (1 - exp(smc[i]
														 [j].vars[smc_Vm]
														 * F / (R_const * T))));
			smc.flux(i, j, ICa_NSC) =
			    1e6 * Am * d_NSCmin * smc.flux(i, j, Po_NSC) *
			    P_CaNSC * smc.var(i, j, smc_Vm) * (P2(z_Ca * F) /
							      (R_const * T)) *
			    ((Ca_o -
			      (smc.var(i, j, smc_Ca) *
			       exp(smc.var(i, j, smc_Vm) * z_Ca * F /
				   (R_const * T)))) / (1 - exp(smc.var(i, j, smc_Vm) * z_Ca * F / (R_const * T))));
			smc.flux(i, j, I_NSC) = smc.flux(i, j, INa_NSC) + smc.flux(i, j, IK_NSC) + smc.flux(i, j, ICa_NSC);

			//Store operated non selective cation channels
			smc.flux(i, j, P_SOC) = 1.0 / (1 + (smc.var(i, j, smc_Ca_u) / K_SOC));
			smc.flux(i, j, ICa_SOC) = g_SOCCa * smc.flux(i, j, P_SOC) * (smc.var(i, j, smc_Vm) - smc.flux(i, j, E_Ca));
			smc.flux(i, j, INa_SOC) = g_SOCNa * smc.flux(i, j, P_SOC) * (smc.var(i, j, smc_Vm) - smc.flux(i, j, E_Na));
			smc.flux(i, j, I_SOC) = smc.flux(i, j, ICa_SOC) + smc.flux(i, j, INa_SOC);

			//Calcium activated Chloride channels
			smc.flux(i, j, alpha_Cl) =
			    pow(smc.var(i, j, smc_cGMP_i),
				n_ClcGMP) / (pow(smc.var(i, j, smc_cGMP_i), n_ClcGMP) + pow(K_ClcGMP, n_ClcGMP));
			smc.flux(i, j, K_ClCacGMP) = (1 - 0.9 * smc.flux(i, j, alpha_Cl)) * 400.00;
			smc.flux(i, j, P_Cl) =
			    (0.0132 * pow(smc.var(i, j, smc_Ca) * 1e6, n_ClCa) /
			     (pow(smc.var(i, j, smc_Ca) * 1e6, n_ClCa) +
			      pow(K_ClCa,
				  n_ClCa))) +
			    (smc.flux(i, j, alpha_Cl) *
			     pow(smc.var(i, j, smc_Ca) * 1e6,
				 n_ClCa) / (pow(smc.var(i, j, smc_Ca) * 1e6, n_ClCa) + pow(smc.flux(i, j, K_ClCacGMP), n_ClCa)));
			smc.flux(i, j, I_ClCa) = C_m * g_ClCa * smc.flux(i, j, P_Cl) * (smc.var(i, j, smc_Vm) - smc.flux(i, j, E_Cl));

			//Plasma membrane Ca pump
			smc.flux(i, j, I_PMCA) = I_bar_PMCA * (smc.var(i, j, smc_Ca)) / ((smc.var(i, j, smc_Ca)) + K_mPMCA);

			//Plasma membrane Na-Ca exchanger
			smc.flux(i, j, phi_F) = exp((F * gama_NCX * smc.var(i, j, smc_Vm)) / (R_const * T));
			smc.flux(i, j, phi_R) = exp((F * (gama_NCX - 1) * smc.var(i, j, smc_Vm)) / (R_const * T));
			smc.flux(i, j, R_NCX_cGMP) = 1 + 0.55 * (smc.var(i, j, smc_cGMP_i) / (smc.var(i, j, smc_cGMP_i) + 45.0));
			smc.flux(i, j, I_NCX) =
			    g_NCX * smc.flux(i, j, R_NCX_cGMP) *
			    (((P3(smc.var(i, j, smc_Na_i)) * Ca_o *
			       smc.flux(i, j, phi_F)) - (P3(Na_o) * smc.var(i, j, smc_Ca) * smc.flux(i, j, phi_R)))
			     / (1 + d_NCX * (P3(Na_o) * smc.var(i, j, smc_Ca) + P3(smc.var(i, j, smc_Na_i)) * Ca_o)));
			//Sodium potassium pump
			smc.flux(i, j, Q_10) = 1.87;
			smc.flux(i, j, Q) = pow(smc.flux(i, j, Q_10), ((T - 309.15) / 10.0));
			smc.flux(i, j, I_NaK) =
			    C_m * I_bar_NaK * smc.flux(i, j, Q) *
			    (pow(K_o, n_HKo) /
			     (pow(K_o, n_HKo) +
			      pow(K_dKo, n_HKo))) * (pow(smc.var(i, j, smc_Na_i),
							 n_HNai) /
						     (pow
						      (smc.var(i, j, smc_Na_i),
						       n_HNai) + pow(K_dNai,
								     n_HNai))) *
			    ((smc.var(i, j, smc_Vm) + 150.0) / (smc.var(i, j, smc_Vm) + 200.0));

			//Sodium potassium chloride cotransport
			smc.flux(i, j, R_NaKCl_cGMP) = 1.0 + 3.5 * (smc.var(i, j, smc_cGMP_i) / (smc.var(i, j, smc_cGMP_i) + 6.4));
			smc.flux(i, j, I_Cl_NaKCl) =
			    -(smc.flux(i, j, R_NaKCl_cGMP) * z_Cl * Am * L_NaKCl *
			      R_const * F * T) * log((Na_o / smc.var(i, j, smc_Na_i)) *
					       (K_o / smc.var(i, j, smc_K_i)) * P2(Cl_o / smc.var(i, j, smc_Cl_i)));
			smc.flux(i, j, I_Na_NaKCl) = -0.5 * smc.flux(i, j, I_Cl_NaKCl);
			smc.flux(i, j, I_K_NaKCl) = -0.5 * smc.flux(i, j, I_Cl_NaKCl);

	/**** Sarcoplasmic Reticulum ****/
			smc.flux(i, j, I_SERCA) = I_bar_SERCA * (smc.var(i, j, smc_Ca) * 1e3 / (smc.var(i, j, smc_Ca) * 1e3 + K_mUp));
			smc.flux(i, j, I_tr) = (smc.var(i, j, smc_Ca_u) - smc.var(i, j, smc_Ca_r)) * z_Ca * F * vol_u / tau_tr;
			smc.flux(i, j, I_rel) =
			    (P2(smc.var(i, j, smc_R_10)) +
			     R_leak) * (smc.var(i, j, smc_Ca_r) - smc.var(i, j, smc_Ca)) * z_Ca * F * vol_r / tau_rel;

			// Ryanodine receptor
			smc.flux(i, j, R_00) = 1.0 - smc.var(i, j, smc_R_01) - smc.var(i, j, smc_R_10) - smc.var(i, j, smc_R_11);

			// IP3 Receptor
			smc.flux(i, j, I_IP3) =
			    I_bar_IP3 * z_Ca * vol_Ca * F *
			    P3(smc.var(i, j, smc_h_IP3) *
			       (smc.var(i, j, smc_IP3) /
				(smc.var(i, j, smc_IP3) +
				 K_IP3)) * (smc.var(i, j, smc_Ca) /
					    (smc.var(i, j, smc_Ca) + K_actIP3))) * (smc.var(i, j, smc_Ca_u) - smc.var(i, j, smc_Ca));

	/**** Alpha Adrenoreceptor Activation and IP3 formation ****/
			smc.flux(i, j, r_h_G) =
			    alpha_G * (smc.var(i, j, smc_Ca) / (smc.var(i, j, smc_Ca) + K_c_G)) * smc.var(i, j, smc_G);
			smc.flux(i, j, rho_r_G) = smc.NE(i, j) * smc.var(i, j, smc_R_S_G) / (ksi_G * R_T_G * (K_1_G + smc.NE(i, j)));
			smc.flux(i, j, delta_G) = 0.0;	//deltaG_0;//0.0;

	/**** sGC activation and cGMP formation ****/
			smc.flux(i, j, A0_sGC) = (((k_1_sGC + k2_sGC) * kD_sGC) + (k_1_sGC * k_2_sGC)) / (k1_sGC * k3_sGC);
			smc.flux(i, j, A1_sGC) = ((k1_sGC + k3_sGC) * kD_sGC + (k2_sGC + k_2_sGC) * k1_sGC) / (k1_sGC * k3_sGC);
			smc.flux(i, j, V_bar_cGMP) =
			    V_cGMPmax * (B5_sGC * smc.NO(i, j) * 1e-6 +
					 P2(smc.NO(i, j) * 1e-6)) /
			    (smc.flux(i, j, A0_sGC) + smc.flux(i, j, A1_sGC) * smc.NO(i, j) * 1e-6 + P2(smc.NO(i, j) * 1e-6));
			if ((smc.flux(i, j, V_bar_cGMP) - smc.var(i, j, smc_V_cGMP)) >= 0) {
				smc.flux(i, j, tau_sGC) = tau_a_sGC;
			} else {
				smc.flux(i, j, tau_sGC) = tau_d_sGC;
			}

	/**** Ionic Balances ****/
			smc.flux(i, j, I_Catotm) =
			    smc.flux(i, j, ICa_SOC) + smc.flux(i, j, I_VOCC) -
			    2 * smc.flux(i, j, I_NCX) + smc.flux(i, j, I_PMCA) + smc.flux(i, j, ICa_NSC);
			smc.flux(i, j, I_Natotm) =
			    smc.flux(i, j, I_Na_NaKCl) + smc.flux(i, j, INa_SOC) +
			    3 * smc.flux(i, j, I_NaK) + 3 * smc.flux(i, j, I_NCX) + smc.flux(i, j, INa_NSC);
			smc.flux(i, j, I_Ktotm) =
			    smc.flux(i, j, I_K_NaKCl) + smc.flux(i, j, I_Kv) +
			    smc.flux(i, j, I_BKCa) + smc.flux(i, j, IK_NSC) + smc.flux(i, j, I_Kleak) - 2 * smc.flux(i, j, I_NaK);
			smc.flux(i, j, I_Cltotm) = smc.flux(i, j, I_Cl_NaKCl) + smc.flux(i, j, I_ClCa);

		}		//end j

	}			//end i

}
/**************************************************************************************/
/**/ void tsoukias_smc_derivatives(double* f, grid_parms grid,SMC_type&  smc) /**/
/**************************************************************************************/
{
	int k;
	for (int i = 1; i <= grid.num_smc_circumferentially; i++) {
		for (int j = 1; j <= grid.num_smc_axially; j++) {
			if (i > 1)
				k = ((i - 1) * grid.neq_smc_axially);
			else if (i == 1)
				k = 0;

			f[k + ((j - 1) * grid.neq_smc) + smc_Vm] = (-1e3
					* (smc.flux(i, j, I_VOCC) + smc.flux(i, j, I_Kv)
							+ smc.flux(i, j, I_BKCa) + smc.flux(i, j, I_Kleak)
							+ smc.flux(i, j, I_NSC) + smc.flux(i, j, I_SOC)
							+ smc.flux(i, j, I_ClCa) + smc.flux(i, j, I_PMCA)
							+ smc.flux(i, j, I_NaK) + smc.flux(i, j, I_NCX)) / C_m)
					+ smc.homo_flux(i, j, cpl_Vm) + smc.hetero_flux(i, j, cpl_Vm);

			f[k + ((j - 1) * grid.neq_smc) + smc_d_L] = 1e3
					* (smc.flux(i, j, dbar_L) - smc.var(i, j, smc_d_L))
					/ smc.flux(i, j, tau_dL);
			f[k + ((j - 1) * grid.neq_smc) + smc_f_L] = 1e3
					* (smc.flux(i, j, fbar_L) - smc.var(i, j, smc_f_L))
					/ smc.flux(i, j, tau_fL);

			f[k + ((j - 1) * grid.neq_smc) + smc_p_f] = 1e3
					* (smc.flux(i, j, pbar_o) - smc.var(i, j, smc_p_f)) / tau_pf;
			f[k + ((j - 1) * grid.neq_smc) + smc_p_s] = 1e3
					* (smc.flux(i, j, pbar_o) - smc.var(i, j, smc_p_s)) / tau_ps;

			f[k + ((j - 1) * grid.neq_smc) + smc_q_1] = 1e3
					* (smc.flux(i, j, qbar) - smc.var(i, j, smc_q_1)) / tau_q1;
			f[k + ((j - 1) * grid.neq_smc) + smc_q_2] = 1e3
					* (smc.flux(i, j, qbar) - smc.var(i, j, smc_q_2)) / tau_q2;
			f[k + ((j - 1) * grid.neq_smc) + smc_p_K] = 1e3
					* (smc.flux(i, j, pbar_K) - smc.var(i, j, smc_p_K))
					/ smc.flux(i, j, tau_pK);

			f[k + ((j - 1) * grid.neq_smc) + smc_Ca_u] = 1e3
					* (smc.flux(i, j, I_SERCA) - smc.flux(i, j, I_tr)
							- smc.flux(i, j, I_IP3)) / (z_Ca * F * vol_u);
			f[k + ((j - 1) * grid.neq_smc) + smc_Ca_r] = 1e3
					* ((smc.flux(i, j, I_tr) - smc.flux(i, j, I_rel))
							/ (z_Ca * F * vol_r))
					/ ((1
							+ (CSQN * K_CSQN
									/ P2(K_CSQN + smc.var(i, j, smc_Ca_r)))));

			f[k + ((j - 1) * grid.neq_smc) + smc_R_10] =
					1e3
							* ((Kr1 * P2(smc.var(i, j, smc_Ca))
									* smc.flux(i, j, R_00))
									- ((K_r1 + Kr2 * smc.var(i, j, smc_Ca))
											* smc.var(i, j, smc_R_10))
									+ (K_r2 * smc.var(i, j, smc_R_11)));
			f[k + ((j - 1) * grid.neq_smc) + smc_R_11] =
					1e3
							* ((Kr2 * smc.var(i, j, smc_Ca)
									* smc.var(i, j, smc_R_10))
									- ((K_r1 + K_r2) * smc.var(i, j, smc_R_11))
									+ (Kr1 * P2(smc.var(i, j, smc_Ca))
											* smc.var(i, j, smc_R_01)));
			f[k + ((j - 1) * grid.neq_smc) + smc_R_01] =
					1e3	* ((Kr2 * smc.var(i, j, smc_Ca) * smc.flux(i, j, R_00))
									+ (K_r1 * smc.var(i, j, smc_R_11))
									- ((K_r2 + Kr1 * P2(smc.var(i, j, smc_Ca)))* smc.var(i, j, smc_R_01)));
			f[k + ((j - 1) * grid.neq_smc) + smc_h_IP3] = 1e3 * k_on_IP3
					* (K_inhIP3
							- (smc.var(i, j, smc_Ca) + K_inhIP3)
									* smc.var(i, j, smc_h_IP3));

			f[k + ((j - 1) * grid.neq_smc) + smc_R_S_G] = 1e3
					* (k_r_G * ksi_G * R_T_G
							- (smc.var(i, j, smc_R_S_G)
									* (k_r_G
											+ (k_p_G * smc.NE(i, j)
													/ (K_1_G + smc.NE(i, j)))))
							- (k_r_G * smc.var(i, j, smc_R_S_P_G)));
			f[k + ((j - 1) * grid.neq_smc) + smc_R_S_P_G] = 1e3 * smc.NE(i, j)
					* ((k_p_G * smc.var(i, j, smc_R_S_G) / (K_1_G + smc.NE(i, j)))
							- (k_e_G * smc.var(i, j, smc_R_S_P_G)
									/ (K_2_G + smc.NE(i, j))));

			f[k + ((j - 1) * grid.neq_smc) + smc_G] = 1e3
					* (k_a_G
							* ((smc.flux(i, j, delta_G) + smc.flux(i, j, rho_r_G))
									* (G_T_G - smc.var(i, j, smc_G)))
							- (k_d_G * smc.var(i, j, smc_G)));
			f[k + ((j - 1) * grid.neq_smc) + smc_IP3] = 1e3
					* (((smc.flux(i, j, r_h_G) / gama_G) * smc.var(i, j, smc_PIP2))
							- (k_deg_G * smc.var(i, j, smc_IP3)))
					+ smc.homo_flux(i, j, cpl_IP3) + smc.hetero_flux(i, j, cpl_IP3);
			f[k + ((j - 1) * grid.neq_smc) + smc_PIP2] = 1e3
					* ((r_r_G * PIP2_T)
							- (smc.flux(i, j, r_h_G) + r_r_G)
									* smc.var(i, j, smc_PIP2)
							- (r_r_G * gama_G * smc.var(i, j, smc_IP3)));

			f[k + ((j - 1) * grid.neq_smc) + smc_V_cGMP] = 1e3
					* (smc.flux(i, j, V_bar_cGMP) - smc.var(i, j, smc_V_cGMP))
					/ smc.flux(i, j, tau_sGC);
			f[k + ((j - 1) * grid.neq_smc) + smc_cGMP_i] = 1e3
					* (smc.var(i, j, smc_V_cGMP)
							- k_pde_cGMP * P2(smc.var(i, j, smc_cGMP_i))
									/ (smc.var(i, j, smc_cGMP_i) + K_m_pde));

			f[k + ((j - 1) * grid.neq_smc) + smc_Ca] =
					1e3
							* ((-(smc.flux(i, j, I_Catotm) + smc.flux(i, j, I_SERCA)
									- smc.flux(i, j, I_rel) - smc.flux(i, j, I_IP3))
									/ (z_Ca * F * vol_Ca))
									/ (1
											+ ((S_CM * K_d)
													/ P2(K_d + smc.var(i, j, smc_Ca)))+
											((B_F * K_dB) / P2(K_dB + smc.var(i, j, smc_Ca))))) + smc.homo_flux(i, j, cpl_Ca) + smc.hetero_flux(i, j, cpl_Ca);

			f[k + ((j - 1) * grid.neq_smc) + smc_Na_i] = 1e3
					* (-smc.flux(i, j, I_Natotm) / (F * vol_i));
			f[k + ((j - 1) * grid.neq_smc) + smc_K_i] = 1e3
					* (-smc.flux(i, j, I_Ktotm) / (F * vol_i));
			f[k + ((j - 1) * grid.neq_smc) + smc_Cl_i] = 1e3
					* (-smc.flux(i, j, I_Cltotm) / (z_Cl * F * vol_i));
			f[k + ((j - 1) * grid.neq_smc) + smc_DAG] = f[k
					+ ((j - 1) * grid.neq_smc) + smc_IP3];
		}
	}
}
