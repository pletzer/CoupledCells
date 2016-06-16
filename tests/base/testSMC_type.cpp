#include "computelib.h"
#include <cassert>

int main() {
	SMC_type smc;
	int nc2 = 22;
	int na2 = 6;
	int nv = 5;
	int nf = 15;
	int ncs = 3;
	smc.build(nc2, na2, nv, nf, ncs);
	// Check setting.
	for (int i = 0; i < nc2; i++) {
		for (int j = 0; j < na2; j++) {
			for (int k = 0; k < nv; k++)
				smc.var(i, j, k) = 0.0;
			for (int k = 0; k < nf; k++)
				smc.flux(i, j, k) = 1.0;
			for (int k = 0; k < ncs; k++)
				smc.homo_flux(i, j, k) = 3.0;
			for (int k = 0; k < ncs; k++)
				smc.hetero_flux(i, j, k) = 4.0;
			smc.NO(i, j) = 5.0;
			smc.NE(i, j) = 6.0;
		}
	}
	// Check getting.
	for (int i = 0; i < nc2; i++) {
		for (int j = 0; j < na2; j++) {
			for (int k = 0; k < nv; k++)
				assert(smc.var(i, j, k) == 0.0);
			for (int k = 0; k < nf; k++)
				assert(smc.flux(i, j, k) == 1.0);
			for (int k = 0; k < ncs; k++)
				assert(smc.homo_flux(i, j, k) == 3.0);
			for (int k = 0; k < ncs; k++)
				assert(smc.hetero_flux(i, j, k) = 4.0);
			assert(smc.NO(i, j) == 5.0);
			assert(smc.NE(i, j) == 6.0);
		}
	}

	smc.destroy();

	std::cout << "SUCCESS\n;";
	return 0;
}